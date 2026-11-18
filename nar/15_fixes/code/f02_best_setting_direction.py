"""Audit fix I03: direction of the best-of-30 setting in nar/02 GI, plus a direction-aware version.

For each set x cohort x drug pickle in nar/02_organoid_only/_data/results:
- best setting (smallest treated log-rank p of the 30 = 3 learners x k 1-10) and its HR / direction;
- best setting among those whose HR > 1 (expected: predicted non-responders die more often);
- conservative shuffle share for that expected-direction best p: share of the 1,000 response shuffles whose best
  p of 30 (ANY direction; the stored shuffles keep no HR) is <= it. Because a direction-restricted best p in a shuffle
  can only be larger, this share is an upper bound on the direction-aware shuffle p.
Writes nar/02_organoid_only/tables/cross_set_best_direction.tsv (read by 02 report_sections.py).
"""
import pickle
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2] / "02_organoid_only"
rows = []
for f in sorted((ROOT / "_data" / "results").glob("*.pkl")):
    r = pickle.load(open(f, "rb"))
    ev = r["ev"].copy()
    folder, cohort, drug = f.stem.split("__")
    pb = np.asarray(r["perm_best"], float)
    # direction: HR > 1 = expected; when the HR is not estimable (0 deaths in one group) use the death proportions
    hr_ok = np.isfinite(ev["treated_HR"]) & (ev["treated_HR"] > 1e-3) & (ev["treated_HR"] < 1e3)
    rate_nr = ev["treated_events_nonresp"] / ev["treated_n_nonresp"].clip(lower=1)
    rate_r = ev["treated_events_resp"] / ev["treated_n_resp"].clip(lower=1)
    ev["expected"] = np.where(hr_ok, ev["treated_HR"] > 1, rate_nr > rate_r)
    bi = ev["treated_p"].idxmin()
    exp = ev[ev["expected"]]
    row = dict(folder=folder, cohort=cohort, model_drug=drug,
               best_setting=f"{ev.loc[bi, 'model']} k={ev.loc[bi, 'k']}", best_p=ev.loc[bi, "treated_p"],
               best_HR=ev.loc[bi, "treated_HR"], best_HR_lo=ev.loc[bi, "treated_HR_lo"], best_HR_hi=ev.loc[bi, "treated_HR_hi"],
               best_direction="expected" if ev.loc[bi, "expected"] else "opposite",
               shuffle_share_best_any_direction=(np.sum(pb <= ev.loc[bi, "treated_p"]) + 1) / (len(pb) + 1),
               settings_p_below_0_05=int((ev["treated_p"] < 0.05).sum()),
               settings_p_below_0_05_expected=int(((ev["treated_p"] < 0.05) & ev["expected"]).sum()),
               settings_p_below_0_05_opposite=int(((ev["treated_p"] < 0.05) & ~ev["expected"]).sum()))
    if len(exp):
        bj = exp["treated_p"].idxmin()
        row.update(best_expected_setting=f"{ev.loc[bj, 'model']} k={ev.loc[bj, 'k']}", best_expected_p=ev.loc[bj, "treated_p"],
                   best_expected_HR=ev.loc[bj, "treated_HR"],
                   shuffle_share_best_expected_upper_bound=(np.sum(pb <= ev.loc[bj, "treated_p"]) + 1) / (len(pb) + 1))
    rows.append(row)
T = pd.DataFrame(rows)
T.to_csv(ROOT / "tables" / "cross_set_best_direction.tsv", sep="\t", index=False)
print(T[["folder", "cohort", "model_drug", "best_setting", "best_p", "best_HR", "best_direction", "shuffle_share_best_any_direction",
         "best_expected_p", "shuffle_share_best_expected_upper_bound"]].round(4).to_string())
