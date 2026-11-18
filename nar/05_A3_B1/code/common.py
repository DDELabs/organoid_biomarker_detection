"""Shared paths and statistics for A3 (statistical controls) and B1 (external validation).

Reuses: nar/02_organoid_only (result pickles, coefficient tables, ssGSEA caches, loaders),
nar/03_organoid_plus_patient (training-table loader incl. STORM GSE109211 rank scores),
robust/obd (rank_score, Reactome gene sets). Cox = statsmodels PHReg, log-rank = statsmodels survdiff.
"""
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402
from statsmodels.duration.survfunc import survdiff  # noqa: E402

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
REPO = FOLDER.parents[1]
DATA = FOLDER / "_data"
TAB = FOLDER / "tables"
FIG = FOLDER / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)
ORG = REPO / "nar" / "02_organoid_only"
ORG2 = ORG / "BLCA_OV_LUAD_OTHER"
SEED = 20261117

# Okabe-Ito colour-blind-safe palette
OK = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "green": "#009E73", "sky": "#56B4E9",
      "purple": "#CC79A7", "grey": "#777777", "black": "#000000", "yellow": "#F0E442"}


def plot_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                         "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
    return plt


def save_fig(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")


def letter(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def short(p, n=44):
    s = str(p).replace("REACTOME_", "").replace("_", " ").lower()
    return s if len(s) <= n else s[: n - 1] + "…"


# ------------------------------------------------------------------------------------------ statistics
def logrank(T, E, g):
    T, E, g = np.asarray(T, float), np.asarray(E, int), np.asarray(g, int)
    if g.sum() in (0, len(g)) or E.sum() == 0:
        return np.nan
    try:
        return float(survdiff(T, E, g)[1])
    except Exception:
        return np.nan


def cox(T, E, X, names=None, check=1):
    """Cox PH (Efron ties). Returns DataFrame coef/se/HR/lo/hi/p per column, or None if it fails.
    Only the first `check` terms (the terms of interest) must have a finite, sane standard error."""
    T = np.clip(np.asarray(T, float), 1e-6, None)
    E = np.asarray(E, int)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    names = names or [f"x{i}" for i in range(X.shape[1])]
    if E.sum() < 2:
        return None
    try:
        r = PHReg(T, X, status=E, ties="efron").fit(disp=False)
        b, se = np.asarray(r.params, float), np.asarray(r.bse, float)
        if not np.all(np.isfinite(se[:check])) or np.any(se[:check] > 20):
            return None
        z = b / se
        return pd.DataFrame({"coef": b, "se": se, "HR": np.exp(b), "lo": np.exp(b - 1.96 * se),
                             "hi": np.exp(b + 1.96 * se), "p": 2 * norm.sf(np.abs(z))}, index=names)
    except Exception:
        return None


def logit(y, X, names=None):
    import statsmodels.api as sm
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    names = names or [f"x{i}" for i in range(X.shape[1])]
    try:
        r = sm.Logit(y, sm.add_constant(X, has_constant="add")).fit(disp=0, maxiter=200)
        b, se = np.asarray(r.params[1:], float), np.asarray(r.bse[1:], float)
        if not np.all(np.isfinite(se)) or np.any(se > 50):
            return None
        return pd.DataFrame({"coef": b, "se": se, "OR": np.exp(b), "lo": np.exp(b - 1.96 * se),
                             "hi": np.exp(b + 1.96 * se), "p": 2 * norm.sf(np.abs(b / se))}, index=names)
    except Exception:
        return None


def power_logrank(D, hr, alpha=0.05, frac=0.5):
    """Schoenfeld: power of a two-sided log-rank test with D events, allocation frac, true HR."""
    if not D or not np.isfinite(hr):
        return np.nan
    zc = norm.ppf(1 - alpha / 2)
    m = np.sqrt(D * frac * (1 - frac)) * abs(np.log(hr))
    return float(norm.cdf(m - zc) + norm.cdf(-m - zc))


def mdhr(D, power=0.8, alpha=0.05, frac=0.5):
    """Minimal detectable HR (>1) at the given power (Schoenfeld)."""
    if not D:
        return np.nan
    return float(np.exp((norm.ppf(1 - alpha / 2) + norm.ppf(power)) / np.sqrt(D * frac * (1 - frac))))


def retrodesign(A, s, alpha=0.05, n_sim=200000, seed=SEED):
    """Gelman & Carlin 2014: power, type-S error and exaggeration ratio (type M) for true effect A, SE s."""
    z = norm.ppf(1 - alpha / 2)
    p_hi = 1 - norm.cdf(z - A / s)
    p_lo = norm.cdf(-z - A / s)
    pw = p_hi + p_lo
    typeS = p_lo / pw
    rng = np.random.default_rng(seed)
    est = A + s * rng.standard_normal(n_sim)
    sig = np.abs(est) > s * z
    exag = float(np.mean(np.abs(est[sig])) / A)
    return float(pw), float(typeS), exag


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    o = np.argsort(q)
    r = q[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.minimum(r, 1)
    out[ok] = res
    return out


def holm(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    o = np.argsort(q)
    r = np.maximum.accumulate(q[o] * (n - np.arange(n)))
    res = np.empty(n)
    res[o] = np.minimum(r, 1)
    out[ok] = res
    return out


def stage_num(s):
    t = str(s).upper().replace("STAGE", "").strip()
    for k, v in (("IV", 4), ("III", 3), ("II", 2), ("I", 1), ("4", 4), ("3", 3), ("2", 2), ("1", 1)):
        if t.startswith(k):
            return float(v)
    return np.nan


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(x, d=2):
    return "NA" if x is None or not np.isfinite(x) else f"{x:.{d}f}"


def fp(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2g}" if p < 0.001 else f"{p:.3f}"
