"""Shared paths and small helpers for nar/08_depmap (DepMap / cell-line mechanism agent).

Run every script with `python -I nar/08_depmap/code/<script>.py` from the repo root.
-I isolates the interpreter (no script dir / cwd on sys.path); we add the two code
folders we trust explicitly below.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                       # repo root
OUT = HERE.parent                            # nar/08_depmap
FIG = OUT / "figures"
TAB = OUT / "tables"
CACHE = ROOT / "nar" / "_cache" / "depmap"   # untrusted downloads (each in its own folder)
DERIVED = CACHE / "derived"                  # our own processed matrices
CUR = ROOT / "data" / "curated" / "preclinical"
for p in (FIG, TAB, DERIVED):
    p.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "robust"))
sys.path.insert(0, str(HERE))

BH3 = "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS"
SEED = 20261117

# Colour-blind-safe palette (Okabe-Ito)
OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "red": "#D55E00", "purple": "#CC79A7", "grey": "#999999", "black": "#000000"}


def savefig(fig, name, caption):
    """Save 300-dpi PNG + PDF and a plain-language caption file."""
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 100})
    return plt


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def load(name):
    return pd.read_pickle(DERIVED / f"{name}.pkl")


def residualize(y, covars):
    """OLS residuals of y on covariates (DataFrame, may contain categoricals already dummy-coded)."""
    X = np.column_stack([np.ones(len(y)), covars.values]) if covars is not None and covars.shape[1] else np.ones((len(y), 1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def lineage_dummies(lin, min_n=5):
    lin = pd.Series(lin).fillna("other")
    vc = lin.value_counts()
    lin = lin.where(lin.map(vc) >= min_n, "other")
    return pd.get_dummies(lin, drop_first=True, dtype=float)


def rank_z(x):
    """Rank-based inverse normal-like transform (ranks scaled to z) for robust partial correlations."""
    from scipy.stats import norm, rankdata
    r = rankdata(x)
    return norm.ppf((r - 0.5) / len(r))
