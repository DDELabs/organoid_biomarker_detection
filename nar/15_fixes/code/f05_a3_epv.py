"""Audit fix I09: adjusted Cox models in nar/05 A3 with too few deaths per covariate.

For each of the 57 A3 analyses: deaths per estimated term (split + covariates). Rows with < 3 deaths per term are
flagged. For every analysis a REDUCED adjusted model is also fitted: split + at most floor(deaths/3) - 1 covariates,
taken in a fixed priority order (stage, age, purity, MSI, sex, regimen terms), so every reduced model has >= 3
deaths per term. Writes nar/05_A3_B1/tables/A3_adjusted_cox_epv.tsv and adds the flag + reduced-model columns to
nar/05_A3_B1/tables/A3_master_results.tsv (old copy in nar/15_fixes/backup/05/tables/).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

C5 = Path(__file__).resolve().parents[2] / "05_A3_B1" / "code"
sys.path.insert(0, str(C5))
import a3_analyze as A  # noqa: E402
from common import DATA, TAB, cox  # noqa: E402

PRIO = ["stage", "stage_missing", "age10", "age10_missing", "purity", "purity_missing", "msi", "msi_missing", "male"]
S = A.summaries()
pat = pd.concat([pd.read_csv(DATA / f"a3_patients_{k}.tsv.gz", sep="\t") for k in ("gi", "other")])
cl = A.clinical()
M = pd.read_csv(TAB / "A3_master_results.tsv", sep="\t")
rows = []
for _, s in S.iterrows():
    aid = s["analysis_id"]
    t = pat[(pat["analysis_id"] == aid) & (pat.group == "treated")].reset_index(drop=True)
    D = int(t.E.sum())
    X, used, reg = A.design(t, cl, s["project"], s["drug"])
    X = X[X.columns[X.std() > 0]]
    terms = 1 + X.shape[1]
    order = [c for c in PRIO if c in X.columns] + [c for c in X.columns if c not in PRIO]
    kmax = max(0, D // 3 - 1)
    keep = order[:kmax]
    r = dict(analysis_id=aid, deaths=D, n_terms_full=terms, deaths_per_term_full=D / terms,
             low_epv_full=D / terms < 3, reduced_covariates=";".join(keep), n_terms_reduced=1 + len(keep),
             deaths_per_term_reduced=D / (1 + len(keep)))
    ca = cox(t["T"], t["E"], np.column_stack([t["split"].values] + ([X[keep].values] if keep else [])), ["split"] + keep)
    if ca is not None:
        r.update(adjred_HR=ca.loc["split", "HR"], adjred_lo=ca.loc["split", "lo"], adjred_hi=ca.loc["split", "hi"],
                 adjred_p=ca.loc["split", "p"])
    rows.append(r)
E = pd.DataFrame(rows)
E = E.merge(M[["analysis_id", "hr", "adj_HR", "adj_lo", "adj_hi", "adj_p"]], on="analysis_id", how="left")
ok = lambda v: np.isfinite(v) & (v > 1e-3) & (v < 1e3)  # noqa: E731
E["full_adj_same_direction"] = np.where(ok(E.hr) & ok(E.adj_HR), np.sign(np.log(E.hr)) == np.sign(np.log(E.adj_HR)), np.nan)
E["reduced_adj_same_direction"] = np.where(ok(E.hr) & ok(E.adjred_HR), np.sign(np.log(E.hr)) == np.sign(np.log(E.adjred_HR)), np.nan)
E.to_csv(TAB / "A3_adjusted_cox_epv.tsv", sep="\t", index=False)
M = M.drop(columns=[c for c in M.columns if c in ("low_epv_full", "deaths_per_term_full", "adjred_HR", "adjred_lo",
                                                   "adjred_hi", "adjred_p", "reduced_covariates")])
M = M.merge(E[["analysis_id", "deaths_per_term_full", "low_epv_full", "reduced_covariates", "adjred_HR", "adjred_lo",
               "adjred_hi", "adjred_p"]], on="analysis_id", how="left")
M.to_csv(TAB / "A3_master_results.tsv", sep="\t", index=False)
f = E[ok(E.hr) & ok(E.adj_HR)]
g = f[~f.low_epv_full.astype(bool)]
h = E[ok(E.hr) & ok(E.adjred_HR)]
summ = pd.DataFrame([
    dict(set="full adjusted model, all estimable", n=len(f), same_direction=int(f.full_adj_same_direction.sum())),
    dict(set="full adjusted model, >= 3 deaths per term only", n=len(g), same_direction=int(g.full_adj_same_direction.sum())),
    dict(set="full adjusted model, < 3 deaths per term (flagged)", n=len(f) - len(g),
         same_direction=int(f.full_adj_same_direction.sum() - g.full_adj_same_direction.sum())),
    dict(set="reduced adjusted model (>= 3 deaths per term), all estimable", n=len(h), same_direction=int(h.reduced_adj_same_direction.sum())),
])
summ.to_csv(TAB / "A3_adjusted_cox_epv_summary.tsv", sep="\t", index=False)
print(summ.to_string())
print(E[E.low_epv_full][["analysis_id", "deaths", "n_terms_full", "adj_HR", "adjred_HR", "adjred_p", "reduced_covariates"]].round(3).to_string())
