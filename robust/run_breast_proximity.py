"""Breast response biomarkers defined the paper's way: network-proximity shortlist per drug first, model second.

For each drug given to breast patients in the training table:
  1. shortlist = pathways close to the drug's targets in the protein network (z <= -1.2816), from the paper's own
     table where the drug is in it, and from the STRING v12 recomputation;
  2. among patients who received the drug, train a ridge model on the shortlisted pathways only and test it on a
     cohort never seen (leave-one-cohort-out), against 200 random shortlists of the same size and all pathways;
  3. biomarker = the shortlisted pathway with the largest model weight (the paper's rule), reported with its own
     association with response across cohorts (hormone-receptor / HER2 adjusted, from run_breast_pathways.py).
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
import run_survlearn as RS  # noqa: E402
import run_train_combined as TC  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import network as N  # noqa: E402

OUT = RESULTS / "train_combined"
CUT = -1.2816
DRUGS = ["PACLITAXEL", "DOCETAXEL", "DOXORUBICIN", "EPIRUBICIN", "CYCLOPHOSPHAMIDE", "FLUOROURACIL", "CARBOPLATIN", "TRASTUZUMAB"]


def loco(d, Z, cols):
    """Pooled leave-one-cohort-out AUROC (mean over test cohorts, weighted by responder x non-responder pairs)."""
    num = den = 0.0
    y = d["responder"].astype(int)
    for c in d["cohort"].unique():
        te, tr = d.index[d["cohort"] == c], d.index[d["cohort"] != c]
        n1, n0 = int(y.loc[te].sum()), int((1 - y.loc[te]).sum())
        if n1 < 8 or n0 < 8 or y.loc[tr].nunique() < 2:
            continue
        a = TC.auc(y.loc[te].values, TC.fit(Z.loc[tr, cols].values, y.loc[tr]).decision_function(Z.loc[te, cols].values))[0]
        num += a * n1 * n0; den += n1 * n0
    return num / den if den else np.nan


def main():
    P, Z = TC.load()
    uni = pd.read_csv(OUT / "breast_model_pathways.tsv", sep="\t", index_col=0)
    uni.index = "REACTOME_" + uni.index
    paper = N.precomputed_proximity()
    B = P[(P["cancer_group"] == "breast") & P["has_response"] & P["treated"].eq(1)]
    rng = np.random.default_rng(0)
    rows, top = [], []
    for drug in DRUGS:
        d = B[B["drugs"].fillna("").str.split(";").apply(lambda s: drug in s)]
        d = d[d["cohort"].map(d.groupby("cohort")["responder"].agg(lambda s: s.sum() >= 8 and (1 - s).sum() >= 8))]
        if d["cohort"].nunique() < 2:
            continue
        lists = {}
        z = RS.proximity(drug)
        if z is not None:
            lists["STRING v12"] = z
        if drug in paper.columns:
            lists["paper table"] = paper[drug].dropna()
        a_all = loco(d, Z, list(Z.columns))
        for src, z in lists.items():
            cols = [c for c in Z.columns if c in set(z.index[z <= CUT])]
            a = loco(d, Z, cols)
            rnd = np.array([loco(d, Z, list(rng.choice(Z.columns, len(cols), replace=False))) for _ in range(200)])
            y = d["responder"].astype(int)
            w = pd.Series(TC.fit(Z.loc[d.index, cols].values, y).coef_[0], cols)
            s = uni.index.isin(cols)
            rows.append(dict(drug=drug, proximity=src, patients=len(d), cohorts=d["cohort"].nunique(), shortlist=len(cols), auc_shortlist=a,
                             auc_random_shortlists=rnd.mean(), share_random_as_good=float((rnd >= a).mean()), auc_all_pathways=a_all,
                             shortlist_pathways_linked_to_response=float((uni.loc[s, "adj_q"] < 0.05).mean()),
                             other_pathways_linked_to_response=float((uni.loc[~s, "adj_q"] < 0.05).mean())))
            print(rows[-1], flush=True)
            for rank, (p, v) in enumerate(w.reindex(w.abs().sort_values(ascending=False).index).head(5).items(), 1):
                top.append(dict(drug=drug, proximity=src, rank=rank, pathway=p.replace("REACTOME_", ""), proximity_z=float(z[p]), model_weight=v,
                                odds_ratio_per_sd=uni.loc[p, "adj_odds_ratio_per_sd"], same_direction=f"{int(uni.loc[p, 'adj_same_direction'])}/{int(uni.loc[p, 'adj_cohorts'])}",
                                q=uni.loc[p, "adj_q"]))
    R, T = pd.DataFrame(rows), pd.DataFrame(top)
    R.to_csv(OUT / "breast_proximity_shortlist_tests.tsv", sep="\t", index=False)
    T.to_csv(OUT / "breast_proximity_biomarkers.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260, "display.max_colwidth", 70, "display.max_rows", 200)
    print(R.round(3).to_string()); print(T.round(4).to_string())


if __name__ == "__main__":
    main()
