"""Audit fixes for nar/09 (primary analyses stay frozen; these are added sensitivity analyses).

I01: Stage 2 GSE39582 5-FU test with the shared definition (nar/15_fixes/definitions/GSE39582_treatment.tsv):
     5FU/FUFOL/FOLFOX vs no chemotherapy, stage II-III only (09 primary: 5FU/FUFOL/FOLFOX/FOLFIRI, all stages).
I19: Stage 1 primary PFI interaction after dropping patients whose record of the drug has no start day
     (09 primary counts them as exposed from day 0).
Writes nar/09_registered_test/tables/stage2_sensitivity_GSE39582_harmonised.tsv and
       nar/09_registered_test/tables/stage1_sensitivity_no_missing_start.tsv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

C9 = Path(__file__).resolve().parents[2] / "09_registered_test" / "code"
sys.path.insert(0, str(C9))
import common as C  # noqa: E402
import s02_stage2 as S2  # noqa: E402
import s01_stage1 as S1  # noqa: E402

# ---------------- I01
H = C.frozen()
P, S = C.load_table()
P = P[P["cohort"] == "GSE39582"]
S = S.loc[P.index]
p0 = S2.cohort(P, S, "GSE39582", H)
rows = []
for label, arms, stages in (("09 primary definition (5FU|FUFOL|FOLFOX|FOLFIRI, all stages)", S2.FU_39582, None),
                            ("shared definition I01 (5FU|FUFOL|FOLFOX, stage II-III)", ["5FU", "FUFOL", "FOLFOX"], (2, 3)),
                            ("subset: 5FU|FUFOL (no oxaliplatin), stage II-III", ["5FU", "FUFOL"], (2, 3))):
    for ep in ("rfs", "os"):
        p = p0[p0["arm"].isin(arms + ["no chemotherapy"])]
        if stages:
            p = p[p["stage"].map(C.stage_num).isin(stages)]
        arm = p["arm"].isin(arms).astype(float)
        r = S2.fit_tte(p, p["score_FLUOROURACIL"], arm, ep)
        row = dict(definition=label, endpoint=ep.upper(), **{k: v for k, v in r.items() if k not in ("b", "se")})
        if r.get("fitted"):
            row.update(HR_int=np.exp(r["b"]), lo=np.exp(r["b"] - 1.96 * r["se"]), hi=np.exp(r["b"] + 1.96 * r["se"]),
                       p_one_sided=S2.one_sided(r["b"], r["se"], "tte"), p_two_sided=2 * stats.norm.sf(abs(r["b"] / r["se"])))
        rows.append(row)
T = pd.DataFrame(rows)
T.to_csv(C.TAB / "stage2_sensitivity_GSE39582_harmonised.tsv", sep="\t", index=False)
print(T[["definition", "endpoint", "n_drug_arm", "n_comparator", "events_drug_arm", "events_comparator", "HR_int", "lo", "hi", "p_one_sided"]].round(4).to_string())

# ---------------- I19
d_all, Zm, long = S1.build()
raw = S1.R.load_response_labels()
raw = raw[raw["therapy_type"].isin(S1.THERAPY) & ~raw["drug"].isin(S1.SUPPORTIVE) & raw["drug"].notna()]
incl = {u: S1.included_projects(d_all, u) for u in S1.DRUG_UNITS}
out = []
for drug in C.DRUGS:
    miss = set(raw.loc[raw["drug"].isin(S1.DRUG_UNITS[drug]) & raw["start_days"].isna(), "patient"])
    d = d_all[d_all["project"].isin(incl[drug])]
    for lab, dd in (("09 primary (missing start = day 0)", d), ("drop patients with a missing start day", d[~d["patient"].isin(miss)])):
        des = S1.Design(dd, drug, "PFI")
        r = des.fit(dd[f"score_{drug}"].values)
        s = S1.summarise(r)
        out.append(dict(drug=drug, analysis=lab, n=des.n, n_exposed=des.n_exposed, n_dropped=int(d["patient"].isin(miss).sum()) if "drop" in lab else 0,
                        events_exposed=des.D1, events_unexposed=des.D0, **s))
O = pd.DataFrame(out)
O.to_csv(C.TAB / "stage1_sensitivity_no_missing_start.tsv", sep="\t", index=False)
print(O.round(4).to_string())
