"""Write the numeric sections of REPORT_GI.md (tables + rule-based 'What I think' lines) from files written by
plots.py / extras.py. The narrative sections of the report are written by hand around this output
(`tables/report_numbers.md`), so every number in the report comes from a file the code wrote."""
import numpy as np
import pandas as pd

from common import ROOT
from plots import fmt_p, hr_ok, short

S = pd.read_csv(ROOT / "tables" / "cross_set_summary.tsv", sep="\t")


def hr(h, lo, hi):
    return f"{h:.2f} ({lo:.2f}-{hi:.2f})" if hr_ok(h, lo, hi) else "not estimable"


def think(r):
    p, pu, sh = r.treated_p, r.untreated_p, r.shuffle_share_rank1
    right = r.treated_events_nonresp > r.treated_events_resp
    bits = []
    if p < 0.05:
        bits.append("the split separates treated patients" + (" in the expected direction" if right else
                                                               ", but predicted responders did worse (opposite direction)"))
        bits.append(f"{'about ' + str(int(round(sh * 100))) + '% of' if sh >= 0.01 else 'under 1% of'} shuffled-response runs did as well")
    elif p < 0.1:
        bits.append("a weak trend in treated patients" + ("" if right else " (opposite direction)"))
    else:
        bits.append("treated patients are not separated by this pathway")
    if np.isfinite(pu):
        bits.append("untreated patients are also separated, so part of it may be general prognosis" if pu < 0.05
                    else "untreated patients are not separated")
    if r.treated_events < 15:
        bits.append(f"only {int(r.treated_events)} deaths, so the estimate is fragile")
    return "; ".join(bits) + "."


rows = []
for r in S.itertuples():
    rows.append(f"| {r.folder} | {r.cohort} | {r.model_drug.lower()}{' *' if r.borderline_n else ''} | {r.n_models} | {r.n_proximal} | "
                f"{short(r.rank1_pathway, 50)} | {int(r.treated_n)} ({int(r.treated_events)}: {int(r.treated_events_resp)} vs {int(r.treated_events_nonresp)}) | "
                f"{hr(r.treated_HR, r.treated_HR_lo, r.treated_HR_hi)} | {fmt_p(r.treated_p)} | "
                f"{int(r.untreated_n) if np.isfinite(r.untreated_n) else 'NA'} ({int(r.untreated_events) if np.isfinite(r.untreated_events) else 'NA'}) | "
                f"{hr(r.untreated_HR, r.untreated_HR_lo, r.untreated_HR_hi)} | {fmt_p(r.untreated_p)} | {fmt_p(r.best_p_30_tests)} ({r.best_setting}) | "
                f"{fmt_p(r.bootstrap_share)} | {fmt_p(r.shuffle_share_rank1)} | {fmt_p(r.shuffle_share_best)} | {r.loo_rho_Ridge:.2f} | {think(r)} |")
hdr = ("| folder | cohort | drug | models | proximal | rank-1 pathway (Ridge) | treated n (deaths: resp vs non-resp) | HR treated (95% CI) | "
       "p treated | untreated n (deaths) | HR untreated | p untreated | best of 30 tests | bootstrap share | shuffle share rank-1 | "
       "shuffle share best-of-30 | LOO rho Ridge | What I think |\n" + "|---" * 18 + "|")
out = ["## Cross-set table (paper's cut-off z <= -1.2816; rank-1 Ridge pathway; median split; overall survival)", "",
       "* = 12-14 treated patients (below the ~15 target; kept and flagged).", "", hdr, *rows, ""]
cut = pd.read_csv(ROOT / "tables" / "cross_set_cutoff_sensitivity.tsv", sep="\t")
out += ["## Cut-off sensitivity (rank-1 Ridge, treated p; number of proximal pathways in brackets)", "",
        "| folder | cohort | drug | " + " | ".join(f"{c:g}" for c in sorted(cut.cutoff.unique(), reverse=True)) + " |",
        "|---" * (3 + cut.cutoff.nunique()) + "|"]
for (f, c, d), g in cut.groupby(["folder", "cohort", "drug"], sort=False):
    g = g.set_index("cutoff").sort_index(ascending=False)
    out.append(f"| {f} | {c} | {d.lower()} | " + " | ".join(
        (f"{fmt_p(x.treated_p)} ({int(x.n_pathways)})" if x.n_pathways > 0 else "0 pathways") for x in g.itertuples()) + " |")
(ROOT / "tables" / "report_numbers.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
