"""Which pathways drive the breast chemotherapy-response model?

For each pathway: (1) its weight in the model trained on all breast response cohorts, (2) how stable that weight is
when each cohort is left out in turn, (3) its own association with response in every cohort separately, with and
without adjusting for hormone-receptor and HER2 status, combined across cohorts.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
import run_train_combined as TC  # noqa: E402
from obd import RESULTS  # noqa: E402

OUT = RESULTS / "train_combined"


def main():
    P, Z = TC.load()
    d = P[(P["cancer_group"] == "breast") & P["has_response"] & P["treated"].eq(1)]
    ok = [c for c, x in d.groupby("cohort") if x["responder"].sum() >= 8 and (1 - x["responder"]).sum() >= 8]
    d = d[d["cohort"].isin(ok)]
    y = d["responder"].astype(int)
    X = Z.loc[d.index]
    w = pd.Series(TC.fit(X.values, y).coef_[0], X.columns)
    loo = pd.DataFrame({c: TC.fit(X[d["cohort"] != c].values, y[d["cohort"] != c]).coef_[0] for c in ok}, index=X.columns)
    rank = loo.abs().rank(ascending=False)
    rows = {}
    for p in X.columns:
        est = {"raw": [], "adj": []}
        for c in ok:
            i = d.index[d["cohort"] == c]
            x = X.loc[i, p]
            if x.std() < 1e-9:
                continue
            for kind in ("raw", "adj"):
                D = pd.DataFrame({"x": x})
                if kind == "adj":
                    for k in ("er", "her2"):
                        v = d.loc[i, k]
                        if v.notna().mean() > 0.9 and v.nunique() > 1:
                            D[k] = v.fillna(v.mode().iat[0])
                    if D.shape[1] == 1:
                        continue                      # no receptor data in this cohort
                try:
                    r = sm.Logit(y.loc[i].values.astype(float), sm.add_constant(D.astype(float))).fit(disp=0)
                    est[kind].append((r.params["x"], r.bse["x"]))
                except Exception:
                    pass
        o = {}
        for kind, e in est.items():
            b, s = np.array([x[0] for x in e]), np.array([x[1] for x in e])
            wt = 1 / s ** 2
            bm, se = (b * wt).sum() / wt.sum(), np.sqrt(1 / wt.sum())
            o[f"{kind}_cohorts"], o[f"{kind}_same_direction"] = len(b), int((np.sign(b) == np.sign(bm)).sum())
            o[f"{kind}_odds_ratio_per_sd"], o[f"{kind}_p"] = np.exp(bm), 2 * norm.sf(abs(bm / se))
        rows[p] = o
    A = pd.DataFrame(rows).T
    A.insert(0, "model_weight", w)
    A.insert(1, "weight_same_sign_in_leave_one_out", (np.sign(loo).T == np.sign(w)).sum())
    A.insert(2, "worst_rank_in_leave_one_out", rank.max(axis=1).astype(int))
    A["adj_q"] = multipletests(A["adj_p"], method="fdr_bh")[1]
    A["raw_q"] = multipletests(A["raw_p"], method="fdr_bh")[1]
    A.index = A.index.str.replace("REACTOME_", "")
    A.index.name = "pathway"
    A = A.reindex(A["model_weight"].abs().sort_values(ascending=False).index)
    A.to_csv(OUT / "breast_model_pathways.tsv", sep="\t")
    pd.set_option("display.width", 260, "display.max_colwidth", 62, "display.max_columns", 30)
    c = ["model_weight", "weight_same_sign_in_leave_one_out", "worst_rank_in_leave_one_out", "raw_odds_ratio_per_sd", "raw_same_direction", "raw_cohorts",
         "adj_odds_ratio_per_sd", "adj_same_direction", "adj_cohorts", "adj_q"]
    print(len(d), ok)
    print(A[c].head(30).round(3).to_string())
    print("weight share of top 20:", round(float(A["model_weight"].abs().head(20).sum() / A["model_weight"].abs().sum()), 3))
    print("pathways with raw q<0.05:", int((A["raw_q"] < 0.05).sum()), "| adj q<0.05:", int((A["adj_q"] < 0.05).sum()), "of", len(A))
    print(A.sort_values("adj_p")[c].head(20).round(4).to_string())


if __name__ == "__main__":
    main()
