"""Shared loading for the W3 breast benchmark. Reuses robust/run_train_combined.py (loader, L2 logistic, AUC)."""
import os, sys, warnings
os.environ.setdefault("OMP_NUM_THREADS", "2")
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "robust"))
import run_train_combined as TC  # noqa: E402
from obd import survlearn as SL  # noqa: E402
HERE = Path(__file__).resolve().parents[1]
TAB, FIG, LOG = HERE / "tables", HERE / "figures", HERE / "logs"
SEED = 20261118


def breast_data():
    """Same selection as part1 of run_train_combined: treated breast patients with a response label, cohorts with >=8 of each class."""
    P, Z = TC.load()
    L = P[P["has_response"] & P["treated"].eq(1)]
    d = L[L["cancer_group"] == "breast"]
    ok = [c for c, x in d.groupby("cohort") if x["responder"].sum() >= 8 and (1 - x["responder"]).sum() >= 8]
    d = d[d["cohort"].isin(ok)]
    return d, Z.loc[d.index], P
