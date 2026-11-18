"""Shared helpers for nar/13_gaps_bio (biology gap-filling agent).

Paths, figure style (Okabe-Ito colours, 300 dpi PNG + PDF + caption), Spearman with 95% CI, a fast numpy ssGSEA
(Barbie 2009, weight 0.25, the same formula gseapy uses for its raw ES) for random gene-set nulls, and survival
helpers reused from nar/02_organoid_only/code/engine.py.
"""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE.parent                         # nar/13_gaps_bio
REPO = OUT.parents[1]
NAR = REPO / "nar"
FIG, TAB, DATA = OUT / "figures", OUT / "tables", OUT / "_data"
for p in (FIG, TAB, DATA):
    p.mkdir(parents=True, exist_ok=True)
PRECLIN = REPO / "data" / "curated" / "preclinical"
VDW = PRECLIN / "colorectal_vandewetering2015"
O2 = NAR / "02_organoid_only"
DEPMAP = NAR / "_cache" / "depmap"
BH3 = "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS"
SEED = 20261118

OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "red": "#D55E00", "purple": "#CC79A7", "grey": "#999999", "black": "#000000"}


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def panel(ax, letter, x=-0.14):
    ax.text(x, 1.05, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def savefig(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")


def spearman_ci(x, y, n_boot=2000, seed=SEED):
    """Spearman rho, p, n, Fisher-z 95% CI (se = 1.06/sqrt(n-3), Fieller 1957) and percentile bootstrap CI."""
    x, y = pd.Series(np.asarray(x, float)), pd.Series(np.asarray(y, float))
    ok = (x.notna() & y.notna()).values
    x, y = x.values[ok], y.values[ok]
    n = len(x)
    if n < 5:
        return dict(rho=np.nan, p=np.nan, n=n, lo=np.nan, hi=np.nan, boot_lo=np.nan, boot_hi=np.nan)
    r, p = stats.spearmanr(x, y)
    z, se = np.arctanh(np.clip(r, -0.9999, 0.9999)), 1.06 / np.sqrt(max(n - 3, 1))
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if np.ptp(x[i]) == 0 or np.ptp(y[i]) == 0:
            continue
        bs.append(stats.spearmanr(x[i], y[i])[0])
    bs = np.array(bs)
    return dict(rho=float(r), p=float(p), n=n, lo=float(np.tanh(z - 1.96 * se)), hi=float(np.tanh(z + 1.96 * se)),
                boot_lo=float(np.nanpercentile(bs, 2.5)), boot_hi=float(np.nanpercentile(bs, 97.5)))


def ssgsea_slow(expr, sets, alpha=0.25):
    """Raw ssGSEA ES (sets x samples). expr genes x samples (DataFrame); sets: dict name -> list of genes.
    Within each sample genes are ranked (high to low); ES = sum over the walk of (P_hit - P_miss), as gseapy."""
    X = expr.values
    G, S = X.shape
    order = np.argsort(-X, axis=0, kind="stable")              # positions of genes, high expression first
    ranks = np.empty_like(order)
    for j in range(S):
        ranks[order[:, j], j] = np.arange(G)
    # gseapy ranks by rank-normalised value; weight = rank (G..1) ** alpha
    w_by_pos = (G - np.arange(G)).astype(float) ** alpha
    idx = {g: i for i, g in enumerate(expr.index)}
    out = np.zeros((len(sets), S))
    for k, (name, genes) in enumerate(sets.items()):
        gi = np.array([idx[g] for g in genes if g in idx])
        m = len(gi)
        for j in range(S):
            pos = np.sort(ranks[gi, j])
            hit = np.zeros(G)
            hit[pos] = w_by_pos[pos]
            ph = np.cumsum(hit) / hit.sum()
            miss = np.ones(G)
            miss[pos] = 0
            pm = np.cumsum(miss) / (G - m)
            out[k, j] = (ph - pm).sum()
    return pd.DataFrame(out, index=list(sets), columns=expr.columns)


def ssgsea_fast(expr, sets, alpha=0.25):
    """Vectorised closed form of ssgsea_slow (identical values): the sum of a running cumulative sum equals
    sum_j x_j * (G - j), so ES = sum_hit w*(G-pos)/sum_hit w - (G(G+1)/2 - sum_hit (G-pos)) / (G - m)."""
    X = expr.values
    G, S = X.shape
    order = np.argsort(-X, axis=0, kind="stable")
    ranks = np.empty_like(order)
    cols = np.arange(S)
    ranks[order, cols[None, :]] = np.arange(G)[:, None]
    idx = {g: i for i, g in enumerate(expr.index)}
    out = np.zeros((len(sets), S))
    tot = G * (G + 1) / 2.0
    for k, genes in enumerate(sets.values()):
        gi = np.array(sorted({idx[g] for g in genes if g in idx}))
        P = ranks[gi, :].astype(float)
        D = G - P
        W = D ** alpha
        out[k] = (W * D).sum(0) / W.sum(0) - (tot - D.sum(0)) / (G - len(gi))
    return pd.DataFrame(out, index=list(sets), columns=expr.columns)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def version_rows(paths):
    rows = []
    for p in paths:
        p = Path(p)
        if p.exists():
            rows.append(dict(file=str(p.relative_to(REPO)), bytes=p.stat().st_size,
                             sha256=sha(p) if p.stat().st_size < 3e9 else "skipped (large)"))
    return rows


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2g}" if p < 0.001 else f"{p:.3f}"


def o2_path():
    sys.path.insert(0, str(O2 / "code"))
