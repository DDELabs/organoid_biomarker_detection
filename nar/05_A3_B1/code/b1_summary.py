"""B1 summary counts (tables/B1_summary_counts.tsv) and the per-cancer key rows used in the report."""
import numpy as np
import pandas as pd
from common import TAB, bh
R = pd.read_csv(TAB / "B1_external_results.tsv", sep="\t")
R["primary_p"] = np.where(R.has_comparator, R.int_p, R.p)
R["primary_est"] = np.where(R.has_comparator, R.int_est, R.est)
R["BH_q_primary"] = bh(R["primary_p"])
R.to_csv(TAB / "B1_external_results.tsv", sep="\t", index=False)
rows = []
for g, d in [("all", R)] + list(R.groupby("group")):
    rows.append({"group": g, "tests": len(d), "treated_p<0.05": int((d.p < 0.05).sum()),
                 "treated_p<0.05_expected_dir": int(((d.p < 0.05) & (d.est > 1)).sum()),
                 "treated_expected_dir": int((d.est > 1).sum()),
                 "with_comparator": int(d.has_comparator.sum()), "interaction_p<0.05": int((d.int_p < 0.05).sum()),
                 "interaction_p<0.05_expected_dir": int(((d.int_p < 0.05) & (d.int_est > 1)).sum()),
                 "null_p<0.05": int((d.null_p_two_sided < 0.05).sum()),
                 "BH_q_primary<0.05": int((d.BH_q_primary < 0.05).sum()), "expected_by_chance_0.05": 0.05 * len(d)})
S = pd.DataFrame(rows)
S.to_csv(TAB / "B1_summary_counts.tsv", sep="\t", index=False)
print(S.to_string())
k = R[(R.p < 0.05) | (R.int_p < 0.05) | (R.null_p_two_sided < 0.05)]
print(k[["sig_id", "cohort", "endpoint", "comparison", "n", "events", "est", "lo", "hi", "p", "comp_est", "comp_p",
         "int_est", "int_lo", "int_hi", "int_p", "null_p_two_sided", "BH_q_primary"]].round(3).to_string())
