"""A3 step 1a: per-patient table for the 23 GI organoid-only analyses (nar/02_organoid_only/_data/results/*.pkl).

For each analysis: treated and untreated TCGA patients, OS time/event, the frozen rank-1 Ridge score
(coef x z-scored pathway score, exactly as the organoid-only agent), its median split (re-derived and checked against
the stored split), the multi-pathway (Ridge top-10) score, and the co-drugs of each patient (for the regimen covariate).
Also writes the frozen signatures (pathways + Ridge coefficients) for B1.
"""
import pickle
import sys

import numpy as np
import pandas as pd

from common import DATA, ORG, logrank

import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location("gi_common", ORG / "code" / "common.py")
GI = importlib.util.module_from_spec(spec)
spec.loader.exec_module(GI)

RES = ORG / "_data" / "results"


def main():
    rows, sigs, checks = [], [], []
    treat_cache = {}
    for f in sorted(RES.glob("*.pkl")):
        d = pickle.load(open(f, "rb"))
        m = d["meta"]
        aid = f.stem
        coef = d["coefs"]["Ridge"]
        order = coef.abs().sort_values(ascending=False, kind="stable").index
        top = order[0]
        top10 = list(order[:10])
        for r, p in enumerate(order[:10], 1):
            sigs.append({"analysis_id": aid, "rank": r, "pathway": p, "coef": float(coef[p])})
        P = d["Pall"]
        proj = m["cohort"]
        parts = ["COAD"] if proj == "COAD_authors" else (["COAD", "READ"] if proj == "COADREAD" else [proj])
        cod = {}
        for pr in parts:
            if pr not in treat_cache:
                treat_cache[pr] = GI.treatment_table(pr)
            t = treat_cache[pr]
            for p, g in t.groupby("patient")["drug"]:
                cod.setdefault(p, set()).update(g)
        for who in ("treated", "untreated"):
            pats = d["tp"] if who == "treated" else d["up"]
            if (who, 1) not in d["km"]:
                continue
            T, E, g, used = d["km"][(who, 1)]
            s1 = coef[top] * P.loc[pats, top].values
            sk = P.loc[pats, top10].values @ coef[top10].values
            g2 = (s1 > np.median(s1)).astype(int)
            checks.append({"analysis_id": aid, "group": who, "n": len(pats), "split_identical": bool((g2 == g).all()),
                           "p_stored": logrank(T, E, g), "p_recomputed": logrank(T, E, g2)})
            for i, p in enumerate(pats):
                rows.append({"analysis_id": aid, "patient": p, "group": who, "T": float(T[i]), "E": int(E[i]),
                             "score1": float(s1[i]), "split": int(g2[i]), "score10": float(sk[i]),
                             "codrugs": ";".join(sorted(cod.get(p, set())))})
    pd.DataFrame(rows).to_csv(DATA / "a3_patients_gi.tsv.gz", sep="\t", index=False)
    pd.DataFrame(sigs).to_csv(DATA / "signatures_gi.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(DATA / "a3_check_gi.tsv", sep="\t", index=False)
    print(pd.DataFrame(checks).to_string())


if __name__ == "__main__":
    main()
