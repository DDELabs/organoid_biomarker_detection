"""Shared helpers for nar/14_gaps_stats. Reuses nar/05_A3_B1/code/common.py (stats, plotting, palette)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
NAR = FOLDER.parent
REPO = NAR.parent
A3 = NAR / "05_A3_B1"
sys.path.insert(0, str(A3 / "code"))
import common as C  # noqa: E402

TAB = FOLDER / "tables"
FIG = FOLDER / "figures"
DATA = FOLDER / "_data"
for d in (TAB, FIG, DATA):
    d.mkdir(parents=True, exist_ok=True)
SEED = 20261118
OK = C.OK


def save_fig(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")


def load_patients():
    p = pd.concat([pd.read_csv(A3 / "_data" / f"a3_patients_{k}.tsv.gz", sep="\t") for k in ("gi", "other")],
                  ignore_index=True)
    p["proj"] = p.analysis_id.str.split("__").str[1]
    p["pool"] = p.proj.replace({"COAD": "COADREAD", "READ": "COADREAD"})
    p.loc[p.analysis_id.str.startswith("BLCA_lee2018_paper"), "pool"] = "BLCA_paper"
    return p


def master():
    return pd.read_csv(A3 / "tables" / "A3_master_results.tsv", sep="\t")


def logrank_fast(T, E, g):
    """Two-sided log-rank p (chi-square, 1 df), hypergeometric variance with ties; same as survdiff."""
    from scipy.stats import chi2
    T = np.asarray(T, float); E = np.asarray(E, int); g = np.asarray(g, int)
    if g.sum() in (0, len(g)) or E.sum() == 0:
        return np.nan
    o = np.argsort(T, kind="stable")
    T, E, g = T[o], E[o], g[o]
    ut, idx = np.unique(T, return_inverse=True)
    d = np.bincount(idx, weights=E)
    d1 = np.bincount(idx, weights=E * g)
    c = np.bincount(idx)
    c1 = np.bincount(idx, weights=g)
    n = c[::-1].cumsum()[::-1].astype(float)
    n1 = c1[::-1].cumsum()[::-1]
    m = d > 0
    d, d1, n, n1 = d[m], d1[m], n[m], n1[m]
    exp1 = d * n1 / n
    var = np.where(n > 1, d * (n1 / n) * (1 - n1 / n) * (n - d) / np.maximum(n - 1, 1), 0)
    v = var.sum()
    if v <= 0:
        return np.nan
    return float(chi2.sf((d1.sum() - exp1.sum()) ** 2 / v, 1))


DRUG_CLASS = {"CISPLATIN": "platinum", "CARBOPLATIN": "platinum", "OXALIPLATIN": "platinum",
              "FLUOROURACIL": "antimetabolite", "GEMCITABINE": "antimetabolite", "PEMETREXED": "antimetabolite",
              "METHOTREXATE": "antimetabolite", "CAPECITABINE": "antimetabolite",
              "PACLITAXEL": "microtubule", "DOCETAXEL": "microtubule", "VINBLASTINE": "microtubule", "VINORELBINE": "microtubule",
              "IRINOTECAN": "topoisomerase/anthracycline", "DOXORUBICIN": "topoisomerase/anthracycline",
              "ETOPOSIDE": "topoisomerase/anthracycline", "TOPOTECAN": "topoisomerase/anthracycline",
              "EPIRUBICIN": "topoisomerase/anthracycline",
              "DACARBAZINE": "alkylating", "TEMOZOLOMIDE": "alkylating", "CYCLOPHOSPHAMIDE": "alkylating",
              "IFOSFAMIDE": "alkylating", "TRABECTEDIN": "alkylating"}


def drug_class(d):
    return DRUG_CLASS.get(str(d).upper(), "targeted/other")


def short_aid(a):
    s, pr, d = a.split("__")
    s = (s.replace("COAD_vandewetering", "vdW").replace("BLCA_lee2018", "Lee").replace("PDX_gao2015_", "PDX ")
         .replace("LUAD_cellline_substitute", "LungCL").replace("OV_vias2023", "Vias").replace("SARC_alshihabi2024", "AlShihabi")
         .replace("COAD_isella_pdx", "IsellaPDX").replace("PAAD_", "").replace("LIHC_", ""))
    return f"{s} | {pr} | {d.lower()}"
