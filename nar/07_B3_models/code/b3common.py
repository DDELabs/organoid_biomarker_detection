"""Shared paths, data access and evaluation statistics for the B3 model comparison (nar/07_B3_models).

Reuses: robust/obd (gene sets, drug targets, rank score, survival log-rank, network RWR), the paper's STRING>700
network extracted by the organoid-only agent (nar/02_organoid_only/BLCA_OV_LUAD_OTHER/_data/moesm3/Figure_2a_f.tsv),
the bundles built by build_bundles.py from the earlier agents' caches, and the alignment of
nar/03_organoid_plus_patient (contrastive PCA + CORAL, code copied verbatim in `align_models`).
"""
import os
import pickle
import sys
import warnings
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr, rankdata  # noqa: E402

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
B3 = HERE.parent
REPO = B3.parents[1]
NAR = REPO / "nar"
DATA = B3 / "_data"
BUNDLES = DATA / "raw" / "bundles"
TABLES = B3 / "tables"
FIGS = B3 / "figures"
for d in (DATA, TABLES, FIGS):
    d.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "robust"))
from obd import survival as SV  # noqa: E402
from obd import bridge as BR  # noqa: E402
from obd.reference import gene_sets, canonical_symbol, drug_targets, drugbank_targets  # noqa: E402

CUTOFF = -1.2816
SEED = 20261117
MOESM3 = NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "_data" / "moesm3"
PROXCACHE = NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "_data" / "cache"
DEPMAP = NAR / "_cache" / "depmap" / "derived"
# drugs with no protein target of their own use the cisplatin targets / fingerprint (as the earlier agents did)
PROXY = {"CARBOPLATIN": "CISPLATIN", "OXALIPLATIN": "CISPLATIN"}

# ---------------------------------------------------------------------------------------------------- data
BUNDLE_ORDER = ["COAD_paper", "BLCA_paper", "COADREAD_vdW", "COAD_isella", "PAAD_tiriac", "PAAD_shi", "LIHC_ji", "LIHC_broutier",
                "BLCA_current", "OV_vias", "SARC_alshihabi", "LUAD_cells", "PDX_BRCA", "PDX_CRC", "PDX_NSCLC",
                "PDX_PDAC", "PDX_CM"]


def bundle(name):
    return pickle.load(open(BUNDLES / f"{name}.pkl", "rb"))


def combos():
    """All (bundle, drug) pairs with >= 15 treated patients, in report order."""
    out = []
    for b in BUNDLE_ORDER:
        f = BUNDLES / f"{b}.pkl"
        if not f.exists():
            continue
        d = pickle.load(open(f, "rb"))["drugs"]
        out += [(b, k) for k in d]
    return out


@lru_cache(None)
def reactome():
    return dict(gene_sets(("REACTOME",)))


@lru_cache(None)
def paper_network():
    """The paper's STRING>700 network (largest component), canonical symbols (as the organoid-only agent)."""
    import networkx as nx
    e = pd.read_csv(MOESM3 / "Figure_2a_f.tsv", sep="\t")
    a = e["Gene_Symbol_1"].astype(str).map(lambda s: canonical_symbol(s) or s)
    b = e["Gene_Symbol_2"].astype(str).map(lambda s: canonical_symbol(s) or s)
    g = nx.Graph()
    g.add_edges_from(zip(a, b))
    g.remove_edges_from(nx.selfloop_edges(g))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def targets(drug):
    d = PROXY.get(drug, drug)
    f = PROXCACHE / f"proximity_papernet_{d}.tsv.targets"
    if f.exists():
        return open(f).read().strip().split(",")
    t = drug_targets().get(d) or drugbank_targets().get(d) or []
    return list(t)


@lru_cache(None)
def fingerprints():
    """Drug 'fingerprint' = network-proximity z-score per Reactome pathway (paper table + recomputed drugs)."""
    from obd import network as N
    P = N.precomputed_proximity()
    cols = {c: P[c] for c in P.columns}
    for f in PROXCACHE.glob("proximity_papernet_*.tsv"):
        d = f.stem.replace("proximity_papernet_", "")
        if d.endswith("_CHECK"):
            continue
        cols[d] = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
    F = pd.DataFrame(cols)
    for d, p in PROXY.items():
        F[d] = F[p]
    return F


def fp_group(F):
    """Drugs with identical fingerprints (same targets) share a group id (for drug-blind CV)."""
    key = {c: tuple(np.round(F[c].fillna(0).values, 4)) for c in F.columns}
    ids = {}
    return pd.Series({c: ids.setdefault(k, len(ids)) for c, k in key.items()})


def rank_scores(expr):
    from obd.scoring import rank_score
    return rank_score(expr, reactome())


# ---------------------------------------------------------------------------------------------------- alignment
def coral(Xm, Xp, shrink=0.3):
    """Copied from nar/03_organoid_plus_patient/code/common.py (CORAL, shrunk covariances)."""
    def shrunk(X):
        C = np.cov(X, rowvar=False)
        d = np.diag(C).mean()
        return (1 - shrink) * C + shrink * d * np.eye(C.shape[0])

    def mpow(C, p):
        w, v = np.linalg.eigh(C)
        w = np.clip(w, 1e-8, None)
        return (v * w ** p) @ v.T

    Cm, Cp = shrunk(Xm), shrunk(Xp)
    mu = Xm.mean(0)
    return (Xm - mu) @ mpow(Cm, -0.5) @ mpow(Cp, 0.5)


def align_models(Xm, Xp, k=10, min_ratio=0.33, use_coral=True):
    """Copied from nar/03_organoid_plus_patient/code/common.py: BRIDGE contrastive PCA, then CORAL to tumours."""
    Xa, n_drop = BR.align(Xm, Xp, k=k, min_ratio=min_ratio)
    if use_coral and Xa.shape[0] >= 5:
        Xa = coral(Xa, Xp)
    Xa = (Xa - Xa.mean(0)) / np.where(Xa.std(0) > 0, Xa.std(0), 1)
    return Xa, n_drop


# ---------------------------------------------------------------------------------------------------- statistics
def zcols(X, mu=None, sd=None):
    X = np.asarray(X, float)
    if mu is None:
        mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    return np.nan_to_num((X - mu) / sd), mu, sd


def spearman_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 4 or np.std(x[ok]) == 0:
        return dict(rho=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=n)
    r, p = spearmanr(x[ok], y[ok])
    zr, se = np.arctanh(np.clip(r, -0.999, 0.999)), 1 / np.sqrt(max(n - 3, 1))
    return dict(rho=float(r), lo=float(np.tanh(zr - 1.96 * se)), hi=float(np.tanh(zr + 1.96 * se)), p=float(p), n=n)


def cindex(T, E, risk):
    """Harrell's C (higher risk = shorter survival), vectorised."""
    T, E, r = np.asarray(T, float), np.asarray(E, int), np.asarray(risk, float)
    ev = np.where(E == 1)[0]
    if len(ev) == 0:
        return np.nan
    later = T[None, :] > T[ev][:, None]
    num = ((r[ev][:, None] > r[None, :]) & later).sum() + 0.5 * ((r[ev][:, None] == r[None, :]) & later).sum()
    den = later.sum()
    return num / den if den else np.nan


def cindex_ci(T, E, risk, n_boot=200, seed=0):
    rng = np.random.default_rng(seed)
    T, E, r = np.asarray(T, float), np.asarray(E, int), np.asarray(risk, float)
    c = cindex(T, E, r)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(T), len(T))
        bs.append(cindex(T[i], E[i], r[i]))
    bs = np.array([b for b in bs if np.isfinite(b)])
    return c, (np.percentile(bs, 2.5) if len(bs) else np.nan), (np.percentile(bs, 97.5) if len(bs) else np.nan)


def cox_sd(T, E, score, strata=None):
    """HR per SD of the risk score (HR > 1: higher predicted resistance, shorter survival)."""
    from statsmodels.duration.hazard_regression import PHReg
    s = np.asarray(score, float)
    if np.std(s) == 0 or np.sum(E) < 2:
        return dict(HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=np.nan)
    s = (s - s.mean()) / s.std()
    try:
        r = PHReg(np.clip(np.asarray(T, float), 0, None), s[:, None], status=np.asarray(E, float), ties="efron").fit(disp=False)
        b, se = float(r.params[0]), float(r.bse[0])
        return dict(HR=np.exp(b), HR_lo=np.exp(b - 1.96 * se), HR_hi=np.exp(b + 1.96 * se), p=float(r.pvalues[0]))
    except Exception:
        return dict(HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=np.nan)


def median_split(T, E, score):
    s = np.asarray(score, float)
    hi = s > np.median(s)
    T, E = np.asarray(T, float), np.asarray(E, int)
    if hi.sum() in (0, len(hi)):
        return dict(lr_p=np.nan, split_HR=np.nan, split_HR_lo=np.nan, split_HR_hi=np.nan), hi
    p = SV.logrank(T[hi], E[hi], T[~hi], E[~hi])
    h = cox_sd(T, E, hi.astype(float))
    # cox_sd standardises the 0/1 indicator; recompute the plain HR for the split
    from statsmodels.duration.hazard_regression import PHReg
    try:
        r = PHReg(np.clip(T, 0, None), hi[:, None].astype(float), status=E.astype(float), ties="efron").fit(disp=False)
        b, se = float(r.params[0]), float(r.bse[0])
        h = dict(split_HR=np.exp(b), split_HR_lo=np.exp(b - 1.96 * se), split_HR_hi=np.exp(b + 1.96 * se))
    except Exception:
        h = dict(split_HR=np.nan, split_HR_lo=np.nan, split_HR_hi=np.nan)
    return dict(lr_p=p, **h), hi


def interaction(Tt, Et, st, Tu, Eu, su):
    """Treated-vs-untreated interaction: Cox on pooled patients, terms arm + score + arm x score (score scaled on pool)."""
    from statsmodels.duration.hazard_regression import PHReg
    if len(Tu) < 15 or np.sum(Eu) < 5:
        return dict(int_HR=np.nan, int_lo=np.nan, int_hi=np.nan, int_p=np.nan)
    s = np.r_[st, su].astype(float)
    s = (s - s.mean()) / (s.std() or 1)
    a = np.r_[np.ones(len(st)), np.zeros(len(su))]
    X = np.column_stack([a, s, a * s])
    try:
        r = PHReg(np.clip(np.r_[Tt, Tu].astype(float), 0, None), X, status=np.r_[Et, Eu].astype(float), ties="efron").fit(disp=False)
        b, se = r.params, r.bse
        return dict(int_HR=np.exp(b[2]), int_lo=np.exp(b[2] - 1.96 * se[2]), int_hi=np.exp(b[2] + 1.96 * se[2]),
                    int_p=float(r.pvalues[2]))
    except Exception:
        return dict(int_HR=np.nan, int_lo=np.nan, int_hi=np.nan, int_p=np.nan)


def short(p, n=46):
    s = p.replace("REACTOME_", "").replace("_", " ").lower()
    return s if len(s) <= n else s[: n - 1] + "…"


# colour-blind-safe (Okabe-Ito)
OK = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "green": "#009E73", "sky": "#56B4E9",
      "purple": "#CC79A7", "grey": "#777777", "black": "#000000", "yellow": "#F0E442"}


def plot_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def savefig(fig, stem, caption):
    fig.savefig(FIGS / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / f"{stem}.pdf", bbox_inches="tight")
    open(FIGS / f"{stem}.caption.md", "w").write(caption.strip() + "\n")
