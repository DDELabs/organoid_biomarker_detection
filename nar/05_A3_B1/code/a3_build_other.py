"""A3 step 1b: per-patient table for the 34 bladder / ovarian / lung / sarcoma / PDX organoid-only analyses
(nar/02_organoid_only/BLCA_OV_LUAD_OTHER). Patient scores are rebuilt from the agent's cached ssGSEA matrices and
saved Ridge coefficient tables (nothing retrained), and the median-split log-rank p is checked against its summary.
"""
import importlib.util
import sys

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from common import DATA, ORG, ORG2, logrank

spec = importlib.util.spec_from_file_location("k_common", ORG2 / "code" / "common.py")
K = importlib.util.module_from_spec(spec)
sys.modules["k_common"] = K
spec.loader.exec_module(K)

TSC = {"BLCA_lee2018_current": "BLCA_lee2018", "OV_vias2023": "OV_vias2023", "SARC_alshihabi2024": "SARC_alshihabi2024",
       "LUAD_cellline_substitute": "LUAD_gdsc_ccle"}


def tsc_for(set_id, project):
    if set_id == "BLCA_lee2018_paper":
        return pd.read_csv(K.MOESM3 / "Figure_2h.tsv", sep="\t", index_col=0)
    if set_id.startswith("PDX_gao2015_"):
        tag = f"PDX_{set_id.split('_')[-1]}_{project}"
    else:
        tag = TSC[set_id]
    return pd.read_csv(K.CACHE / f"ssgsea_{tag}_tcga.tsv.gz", sep="\t", index_col=0)


def main():
    S = pd.read_csv(ORG2 / "SUMMARY" / "tables" / "cross_set_summary.tsv", sep="\t")
    rows, sigs, checks = [], [], []
    for _, r in S.iterrows():
        proj, drug, set_id = r["project"], r["drug"], r["set"]
        aid = f"{set_id}__{proj}__{drug}"
        tsc = tsc_for(set_id, proj)
        surv = K.survival_biotab(proj) if set_id == "BLCA_lee2018_paper" else K.survival_current(proj)
        co = pd.read_csv(ORG / r["folder"] / "tables" / f"{proj}_{drug}_coefficients.tsv", sep="\t")
        co = co.set_index("pathway")["coef_Ridge"]
        order = co.abs().sort_values(ascending=False, kind="stable").index
        top, top10 = order[0], list(order[:10])
        assert top == r["biomarker"], (aid, top, r["biomarker"])
        for k, p in enumerate(top10, 1):
            sigs.append({"analysis_id": aid, "rank": k, "pathway": p, "coef": float(co[p])})
        pats = [p for p in tsc.columns if p in surv.index]
        cols = [p for p in top10 if p in tsc.index]
        P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[cols, pats].T.values), index=pats, columns=cols)
        treated = K.treated(proj, drug)
        anyd = K.any_drug(proj)
        untreated = set(surv.index) - anyd
        tl = K.treatment_long(proj)
        tl = tl[tl["drug"] != "__ANY_RECORD__"]
        cod = tl.groupby("patient")["drug"].apply(lambda x: ";".join(sorted(set(x)))).to_dict()
        for who, grp in (("treated", treated), ("untreated", untreated)):
            pp = [p for p in pats if p in grp]
            if len(pp) < 6:
                continue
            T = np.clip(surv.loc[pp, "T"].values.astype(float), 0, None)
            E = surv.loc[pp, "E"].values.astype(int)
            s1 = co[top] * P.loc[pp, top].values
            sk = P.loc[pp, top10].values @ co[top10].values
            g = (s1 > np.median(s1)).astype(int)
            pr = logrank(T, E, g)
            checks.append({"analysis_id": aid, "group": who, "n": len(pp), "n_stored": r[f"{who}_n"],
                           "p_stored": r[f"{who}_p"], "p_recomputed": pr})
            for i, p in enumerate(pp):
                rows.append({"analysis_id": aid, "patient": p, "group": who, "T": T[i], "E": int(E[i]),
                             "score1": float(s1[i]), "split": int(g[i]), "score10": float(sk[i]),
                             "codrugs": cod.get(p, "")})
        print(aid, checks[-2:] if len(checks) > 1 else checks, flush=True)
    pd.DataFrame(rows).to_csv(DATA / "a3_patients_other.tsv.gz", sep="\t", index=False)
    pd.DataFrame(sigs).to_csv(DATA / "signatures_other.tsv", sep="\t", index=False)
    c = pd.DataFrame(checks)
    c.to_csv(DATA / "a3_check_other.tsv", sep="\t", index=False)
    print(c.to_string())


if __name__ == "__main__":
    main()
