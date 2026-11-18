"""Per set x drug report blocks (numbers + rule-based 'What I think') for REPORT_GI.md, written to
tables/report_sections.md. All numbers come from tables/cross_set_summary.tsv (plots.py) and
_data/results/run_log_*.tsv (run_all.py)."""
import numpy as np
import pandas as pd

from common import DATA, ROOT
from plots import fmt_p, hr_ok, short

S = pd.read_csv(ROOT / "tables" / "cross_set_summary.tsv", sep="\t")
# audit fix I03: direction of the best-of-30 setting (nar/15_fixes/code/f02_best_setting_direction.py)
_B = ROOT / "tables" / "cross_set_best_direction.tsv"
if _B.exists():
    S = S.merge(pd.read_csv(_B, sep="\t").drop(columns=["best_setting"]), on=["folder", "cohort", "model_drug"], how="left")
COH = {"COAD_authors": "TCGA colon, authors' 2019 files (exact paper inputs)", "COAD": "TCGA colon (current GDC data)",
       "READ": "TCGA rectum (current GDC data)", "COADREAD": "TCGA colon + rectum pooled (current GDC data)",
       "PAAD": "TCGA pancreas", "LIHC": "TCGA liver"}


def hr(h, lo, hi):
    return f"HR {h:.2f} (95% CI {lo:.2f}-{hi:.2f})" if hr_ok(h, lo, hi) else "HR not estimable (no deaths in one group)"


def _h(x):
    return f"{x:.2f}" if np.isfinite(x) and 1e-3 < x < 1e3 else "not estimable"


def think(r):
    p, pu = r.treated_p, r.untreated_p
    right = r.treated_events_nonresp > r.treated_events_resp
    bits = []
    if p < 0.05:
        bits.append("The pathway splits drug-treated patients into groups with different survival"
                    + (", and the predicted non-responders die more often, as expected." if right else
                       ", but the predicted responders are the ones who die more often (the opposite direction)."))
        bits.append(f"In the response-shuffle test {r.shuffle_share_rank1 * 100:.1f}% of shuffled runs did as well, "
                    f"and {r.bootstrap_share * 100:.1f}% of random pathway sets did as well.")
    elif p < 0.1:
        bits.append("There is a weak trend in drug-treated patients" + ("." if right else " (opposite direction)."))
    else:
        bits.append("With this pathway, survival of the two predicted groups of treated patients is similar.")
    if np.isfinite(pu):
        if pu < 0.05:
            bits.append("The same split also separates untreated patients, so part of the signal may be general prognosis "
                        "rather than drug response.")
        elif p < 0.1:   # audit fix I18: only say 'drug-specific' when there is a treated signal to speak of
            bits.append("Untreated patients are not separated, which fits a drug-specific signal if one exists.")
        else:
            bits.append("Untreated patients are not separated either.")
    if r.best_p_30_tests < 0.05 and p >= 0.05:   # audit fix I03: direction-aware wording
        d = getattr(r, "best_direction", None)
        beat = (1 - r.shuffle_share_best) * 100
        if r.shuffle_share_best < 0.05 and d == "opposite":
            be = getattr(r, "best_expected_p", np.nan)
            bits.append(f"Some of the other 29 model settings reach p < 0.05 ({int(r.tests_p_below_0_05)} of 30). The best "
                        f"setting ({r.best_setting}, p = {fmt_p(r.best_p_30_tests)}) beats {beat:.1f}% of shuffled runs, "
                        f"but its HR is {_h(r.best_HR)} (opposite direction): predicted responders do worse, the opposite of the drug "
                        "hypothesis. " + (f"The best setting in the expected direction has p = {fmt_p(be)}; at most "
                        f"{r.shuffle_share_best_expected_upper_bound * 100:.1f}% of shuffles do as well (upper bound)."
                        if np.isfinite(be) else "No setting points in the expected direction."))
        elif r.shuffle_share_best < 0.05:
            bits.append(f"Some of the other 29 model settings reach p < 0.05 ({int(r.tests_p_below_0_05)} of 30). The best "
                        f"setting ({r.best_setting}, p = {fmt_p(r.best_p_30_tests)}) beats {beat:.1f}% of shuffled runs "
                        "and points in the expected direction (predicted non-responders do worse). The rank-1 setting the "
                        "paper's method would use does not reach p < 0.05.")
        else:
            bits.append(f"Some of the other 29 model settings reach p < 0.05 ({int(r.tests_p_below_0_05)} of 30), "
                        f"but {r.shuffle_share_best * 100:.1f}% of shuffled runs find a best p at least this small, "
                        "so picking the best setting is not enough on its own.")
    if r.treated_events < 15:
        bits.append(f"Only {int(r.treated_events)} deaths among treated patients, so any estimate here is fragile.")
    if r.loo_rho_Ridge < 0:        # audit fix I18
        bits.append(f"Inside the lab models the pathway model does not predict drug response (leave-one-out rho = {r.loo_rho_Ridge:.2f}, negative).")
    elif r.loo_rho_Ridge < 0.3:
        bits.append(f"Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = {r.loo_rho_Ridge:.2f}).")
    return " ".join(bits)


out = []
for (lab, folder), g in S.groupby(["set", "folder"], sort=False):
    out.append(f"### {lab}  (`{folder}/`)\n")
    for r in g.itertuples():
        out.append(f"#### {COH.get(r.cohort, r.cohort)} | {r.model_drug.lower()}"
                   + ("  (12-14 treated patients: borderline size)" if r.borderline_n else ""))
        out.append(f"- Lab models: {r.n_models} (response = {r.metric}); proximal pathways: {r.n_proximal}; "
                   f"rank-1 Ridge pathway: {short(r.rank1_pathway, 80)}.")
        out.append(f"- Drug-treated patients: n = {int(r.treated_n)}, deaths = {int(r.treated_events)} "
                   f"({int(r.treated_events_resp)} in predicted responders, {int(r.treated_events_nonresp)} in predicted "
                   f"non-responders); log-rank p = {fmt_p(r.treated_p)}; {hr(r.treated_HR, r.treated_HR_lo, r.treated_HR_hi)}.")
        if np.isfinite(r.untreated_n):
            out.append(f"- Untreated patients (same split): n = {int(r.untreated_n)}, deaths = {int(r.untreated_events)}; "
                       f"p = {fmt_p(r.untreated_p)}; {hr(r.untreated_HR, r.untreated_HR_lo, r.untreated_HR_hi)}.")
        out.append(f"- Other algorithms (rank-1): SVR p = {fmt_p(r.treated_p_SVR)}, LinearRegression p = "
                   f"{fmt_p(r.treated_p_LinearRegression)}. Best of all 30 settings: p = {fmt_p(r.best_p_30_tests)} "
                   f"({r.best_setting}); settings with p < 0.05: {int(r.tests_p_below_0_05)} of 30.")
        out.append(f"- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = "
                   f"{fmt_p(r.bootstrap_share)}. Response shuffle (1,000 runs): rank-1 share = {fmt_p(r.shuffle_share_rank1)}; "
                   f"best-of-30 share = {fmt_p(r.shuffle_share_best)}"
                   + (f" (best setting HR {_h(r.best_HR)}, {r.best_direction} direction)." if isinstance(getattr(r, "best_direction", None), str) else "."))
        out.append(f"- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge {r.loo_rho_Ridge:.2f} "
                   f"(p = {fmt_p(r.loo_p_Ridge)}), SVR {r.loo_rho_SVR:.2f}, LinearRegression {r.loo_rho_LinearRegression:.2f}.")
        out.append(f"- What I think: {think(r)}\n")

logs = [pd.read_csv(f, sep="\t") for f in sorted((DATA / "results").glob("run_log_*.tsv"))]
if logs:
    L = pd.concat(logs)
    sk = L[L.get("skip").notna()] if "skip" in L else L.iloc[0:0]
    out.append("### Combinations that could not be analysed\n")
    out.append("| set | cohort | drug | reason |\n|---|---|---|---|")
    for r in sk.itertuples():
        out.append(f"| {r.set_label} | {r.cohort} | {r.model_drug.lower()} | {r.skip} |")
(ROOT / "tables" / "report_sections.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
