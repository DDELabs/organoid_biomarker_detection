"""Shared helpers for W4 (breast pathways). Reuses robust/run_train_combined.py (loader, L2 logistic with lambda=100, AUC)."""
import os, sys, warnings
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "robust"))
import run_train_combined as TC  # noqa: E402
HERE = Path(__file__).resolve().parents[1]
TAB, FIG, LOG = HERE / "tables", HERE / "figures", HERE / "logs"
SEED = 20261119
PROLIF_KW = ["CELL_CYCLE", "MITOTIC", "G1_S", "G2_M", "E2F", "DNA_REPLICATION", "PRE_REPLICATIVE", "CHECKPOINT"]


def breast_data():
    """Same selection as W3 / part1: treated breast patients with a response label, cohorts with >=8 of each class."""
    P, Z = TC.load()
    L = P[P["has_response"] & P["treated"].eq(1)]
    d = L[L["cancer_group"] == "breast"]
    ok = [c for c, x in d.groupby("cohort") if x["responder"].sum() >= 8 and (1 - x["responder"]).sum() >= 8]
    d = d[d["cohort"].isin(ok)]
    return d, Z.loc[d.index]


def auc(y, s):
    y = np.asarray(y)
    return roc_auc_score(y, s) if len(set(y)) > 1 else np.nan


def hm_se(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))


def prolif_cols(cols):
    return [c for c in cols if any(k in c for k in PROLIF_KW)]


def uni_score(X, y, coh):
    """Signed univariate strength per pathway from TRAINING data only: mean over training cohorts of (AUC - 0.5)."""
    out = []
    for c in np.unique(coh):
        m = coh == c
        r = np.apply_along_axis(rankdata, 0, X[m])
        n1 = y[m].sum(); n0 = m.sum() - n1
        out.append((r[y[m] == 1].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0) - 0.5)
    return np.mean(out, axis=0)
