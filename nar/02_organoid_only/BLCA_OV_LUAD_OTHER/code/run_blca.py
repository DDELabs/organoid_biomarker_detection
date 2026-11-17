"""Bladder organoids (Lee 2018) -> TCGA-BLCA, paper's organoid-only strategy.

Version A ("paper source data"): organoid and TCGA ssGSEA NES exactly as published (MOESM3 sheets 'Figure 2e' and
'Figure 2h'), the paper's proximity z (sheet 'Figure 2a,d'), paper-era BCR biotab survival and drug records.
Version B ("current GDC data"): organoid TPM (GEO GSE103990, NCBI re-quantification) and current GDC STAR FPKM-UQ,
re-scored with ssGSEA; current GDC survival.
Organoid response in both: CoderData 2.1.0 refit AUC (the paper's own IC50 table for Lee 2018 is not public).
Usage: python run_blca.py [drug ...]
"""
import sys

import numpy as np
import pandas as pd

import common as K
import engine as E

OUT = K.PARENT / "BLCA_lee2018"
DRUGS = {  # TCGA drug: (organoid drug, proximity drug, note)
    "CISPLATIN": ("CISPLATIN", "CISPLATIN", ""),
    "GEMCITABINE": ("GEMCITABINE", "GEMCITABINE", ""),
    "METHOTREXATE": ("METHOTREXATE", "METHOTREXATE", "The user's own BLCA choice; only 12 TCGA patients got it."),
    "DOXORUBICIN": ("DOXORUBICIN", "DOXORUBICIN", ""),
    "VINBLASTINE": ("VINBLASTINE", "VINBLASTINE", "Only 12 TCGA patients got it (part of MVAC)."),
    "PACLITAXEL": ("PACLITAXEL", "PACLITAXEL", "Only 11 TCGA patients got it."),
    "CARBOPLATIN": ("CISPLATIN", "CISPLATIN", "No carboplatin screen in the organoids: cisplatin response and cisplatin proximity used as the platinum proxy."),
}


def isogenic():
    df = pd.read_csv(K.MOESM3 / "Figure_4c.tsv", sep="\t", skiprows=4, header=None, index_col=0)
    df = df[df.index.astype(str).str.startswith("REACTOME")]
    df.columns = ["R1", "R2", "R3", "S1", "S2", "S3"]
    return {"nes": df.astype(float), "resistant": ["R1", "R2", "R3"], "sensitive": ["S1", "S2", "S3"]}


def paper_version():
    o = pd.read_csv(K.MOESM3 / "Figure_2e.tsv", sep="\t", index_col=0)
    o.columns = [c.split("_org")[0].replace(".", "_") for c in o.columns]
    t = pd.read_csv(K.MOESM3 / "Figure_2h.tsv", sep="\t", index_col=0)
    surv = K.survival_biotab("BLCA")
    return o, t, surv


def current_version():
    oexpr, _ = K.load_preclinical("bladder_lee2018")
    texpr = K.tcga_expr("BLCA")
    osc, tsc = K.paired_scores("BLCA_lee2018", oexpr, texpr)
    return osc, tsc, K.survival_current("BLCA")


def response(drug):
    r = pd.read_csv(K.PRECLIN / "bladder_lee2018" / "response.tsv", sep="\t")
    return r[r["drug"] == drug].groupby("sample")["response"].median()


def run(drugs, versions=("paper", "current")):
    rows = []
    iso = isogenic()
    anyd = K.any_drug("BLCA")
    for v in versions:
        osc, tsc, surv = paper_version() if v == "paper" else current_version()
        for d in drugs:
            if v == "paper" and d != "CISPLATIN" and "--all-paper" not in sys.argv:
                continue
            od, pdrug, note = DRUGS[d]
            z, src = K.proximity(pdrug)
            lab = "Bladder organoids (Lee 2018), paper source-data NES" if v == "paper" else "Bladder organoids (Lee 2018), current GDC data"
            a = E.Analysis(set_id=f"BLCA_lee2018_{v}", set_label=lab, project="BLCA", drug=d, model_drug=od,
                           osc=osc, tsc=tsc, y=response(od), surv=surv, treated=K.treated("BLCA", d),
                           untreated=set(surv.index) - anyd, z=z, z_source=src,
                           outdir=OUT / v, response_label="AUC of fitted viability curve (CoderData)",
                           n_perm=1000, n_boot=int(K.os_env("NBOOT_CIS", 10000)) if d == "CISPLATIN" else 2000,
                           isogenic=iso if d == "CISPLATIN" else None, note=note, seed=0)
            res = a.run()
            res["version"] = v
            rows.append(res)
            print(v, d, {k: res.get(k) for k in ("n_models", "n_proximal", "biomarker", "treated_n", "treated_p",
                                                  "treated_hr", "untreated_p", "perm_empirical_p")}, flush=True)
    return rows


if __name__ == "__main__":
    drugs = [a for a in sys.argv[1:] if not a.startswith("--")] or list(DRUGS)
    rows = run(drugs)
    f = OUT / "tables" / "BLCA_summary.tsv"
    f.parent.mkdir(parents=True, exist_ok=True)
    new = pd.DataFrame(rows)
    if f.exists():
        old = pd.read_csv(f, sep="\t")
        old = old[~(old.set + old.drug).isin(new.set + new.drug)]
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(f, sep="\t", index=False)
