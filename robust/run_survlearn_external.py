"""Frozen test of the survlearn drug scores in independent trial cohorts (data/curated/trials).

For every drug with a fitted model, every patient who received that drug in an independent cohort
and has a recorded response gets the frozen score. AUROC is computed within cohort and pooled.
Control: the same score in patients of the same cohorts who did NOT receive the drug.
Usage: python robust/run_survlearn_external.py [--perm 2000]
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

from obd import RESULTS, scoring  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.curated import load_curated_trial, trial_catalog  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "survlearn"


def cohort_scores(cid):
    def compute():
        expr, _ = load_curated_trial(cid)
        gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
        return scoring.rank_score(expr, {k: v for k, v in gs.items() if len(v) >= 5})
    sc = scoring.cached(f"curated_{cid}_reactome_rank", compute).T
    return ((sc - sc.mean()) / sc.std().replace(0, 1)).fillna(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=2000)
    a = ap.parse_args()
    co = pd.read_csv(OUT / "coefficients.tsv.gz", sep="\t")
    cat = trial_catalog()
    cat = cat[cat["usable"].astype(str).str.lower().eq("true")]
    data = {}
    for cid in cat["cohort_id"]:
        try:
            _, clin = load_curated_trial(cid, canonical=False)
            if "responder" not in clin or "drugs" not in clin:
                continue
            data[cid] = (cohort_scores(cid), clin)
        except Exception as exc:  # noqa: BLE001
            print(cid, "skipped:", str(exc)[:80])
    rng = np.random.default_rng(3)
    rows, per = [], []
    for drug, c in co.groupby("drug"):
        c = c.set_index("pathway")
        vecs = {}
        if "drug_specific" in c and c["drug_specific"].notna().any():
            vecs["surv"] = -c["drug_specific"].dropna()
        if "notes" in c and c["notes"].notna().any():
            vecs["notes"] = c["notes"].dropna()
        for arm_kind in ("received", "not_received"):
            S, Y, G = {k: [] for k in vecs}, [], []
            for cid, (sc, clin) in data.items():
                got = clin["drugs"].fillna("").astype(str).str.upper().str.split(";").map(lambda x: drug in x)
                sel = clin.index[(got if arm_kind == "received" else ~got) & clin["responder"].notna()]
                sel = [s for s in sel if s in sc.index]
                if arm_kind == "not_received" and not got.any():
                    continue                                  # control only from cohorts that also used the drug
                y = clin.loc[sel, "responder"].astype(int)
                if len(sel) < 15 or y.value_counts().min() < 4 if len(sel) else True:
                    continue
                for k, v in vecs.items():
                    p = [x for x in v.index if x in sc.columns]
                    s = sc.loc[sel, p].values @ v[p].values
                    S[k].append(s)
                    if arm_kind == "received":
                        auc, _ = SL.stratified_auc(s, y.values, np.zeros(len(y)))
                        per.append({"drug": drug, "cohort": cid, "score": k, "n": len(y), "responders": int(y.sum()), "AUC": auc})
                Y.append(y.values)
                G.append(np.full(len(y), cid))
            if not Y:
                continue
            y, g = np.concatenate(Y), np.concatenate(G)
            sc_all = {k: np.concatenate(v) for k, v in S.items()}
            if len(vecs) == 2:                               # combined: mean of within-cohort z-scores
                z = {}
                for k, v in sc_all.items():
                    zz = np.zeros(len(v))
                    for cid in np.unique(g):
                        i = g == cid
                        zz[i] = (v[i] - v[i].mean()) / (v[i].std() or 1)
                    z[k] = zz
                sc_all["combined"] = (z["surv"] + z["notes"]) / 2
            for k, s in sc_all.items():
                auc, pairs = SL.stratified_auc(s, y, g)
                null = []
                for _ in range(a.perm):
                    yp = y.copy()
                    for cid in np.unique(g):
                        i = np.flatnonzero(g == cid)
                        yp[i] = y[rng.permutation(i)]
                    null.append(SL.stratified_auc(s, yp, g)[0])
                rows.append({"drug": drug, "patients": arm_kind, "score": k, "cohorts": len(np.unique(g)), "n": len(y),
                             "responders": int(y.sum()), "AUC": auc, "p_one_sided": (np.sum(np.array(null) >= auc) + 1) / (a.perm + 1)})
    E = pd.DataFrame(rows)
    E.to_csv(OUT / "external.tsv", sep="\t", index=False)
    pd.DataFrame(per).to_csv(OUT / "external_per_cohort.tsv", sep="\t", index=False)
    print(E.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
