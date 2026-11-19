"""Figures for W2 (forest plots). Reads only tables written by run_v2.py and the v1 results; writes PNG + PDF + caption into nar/18_other_cohorts/figures."""
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[3]; T = ROOT/"nar/18_other_cohorts/tables"; F = ROOT/"nar/18_other_cohorts/figures"; F.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "font.family": "DejaVu Sans"})
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#7f7f7f"
def hm_ci(a, n1, n0):
    q1, q2 = a/(2-a), 2*a*a/(1+a)
    se = np.sqrt((a*(1-a)+(n1-1)*(q1-a*a)+(n0-1)*(q2-a*a))/(n1*n0)); return max(0, a-1.96*se), min(1, a+1.96*se)
old = pd.read_csv(ROOT/"robust/results/train_combined/part1_response_leave_one_cohort_out.tsv", sep="\t")
ci = [hm_ci(r.auc_pathways, r.responders, r.n-r.responders) for r in old.itertuples()]
old["lo"], old["hi"] = [c[0] for c in ci], [c[1] for c in ci]
p1 = pd.read_csv(T/"part1_v2.tsv", sep="\t"); new = p1[p1.table == "v2_rectal_split"].copy()
new = new.rename(columns={"auc_ci_low": "lo", "auc_ci_high": "hi"})
old_names = set(old[old.cancer == "breast"].test_cohort)

def save(fig, name, caption):
    fig.savefig(F/f"{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(F/f"{name}.pdf", bbox_inches="tight"); plt.close(fig)
    (F/f"{name}.caption.md").write_text(caption.strip()+"\n")

def forest_panel(ax, rows, title):
    """rows: list of (label, auc, lo, hi, kind, n, resp); kind in old_v1, old_v2, new"""
    y = 0; ticks = []; labs = []
    for lab, grp in rows:
        for kind, a, lo, hi, n, r in grp:
            col = {"old_v1": GREY, "old_v2": BLUE, "new": ORANGE}[kind]
            off = {"old_v1": 0.16, "old_v2": -0.16, "new": 0}[kind]
            ax.errorbar(a, y+off, xerr=[[a-lo], [hi-a]], fmt="o" if kind != "old_v1" else "s", ms=4 if kind != "old_v1" else 3, color=col, mfc=col if kind != "old_v1" else "white", lw=1)
        ticks.append(y); labs.append(f"{lab}  (n={grp[-1][4]}, resp={grp[-1][5]})"); y += 1
    ax.axvline(0.5, color="k", lw=0.8, ls=":"); ax.set_yticks(ticks); ax.set_yticklabels(labs); ax.invert_yaxis(); ax.set_xlim(0.1, 1.02)
    ax.set_xlabel("AUC (area under the ROC curve; 0.5 = coin toss)"); ax.set_title(title, loc="left")

def build_rows(cancer):
    rows = []
    n_ = new[new.cancer == cancer]; o_ = old[old.cancer == cancer]
    names = list(n_.sort_values("auc_pathways").test_cohort)
    for c in names:
        g = []; r2 = n_[n_.test_cohort == c].iloc[0]
        r1 = o_[o_.test_cohort == c]
        if len(r1):
            r1 = r1.iloc[0]; g.append(("old_v1", r1.auc_pathways, r1.lo, r1.hi, int(r1.n), int(r1.responders)))
            g.append(("old_v2", r2.auc_pathways, r2.lo, r2.hi, int(r2.n), int(r2.responders)))
        else:
            g.append(("new", r2.auc_pathways, r2.lo, r2.hi, int(r2.n), int(r2.responders)))
        rows.append((c, g))
    return rows

# ---- Figure 1: breast
rows = build_rows("breast")
fig, ax = plt.subplots(figsize=(8.5, 0.28*len(rows)+1.6))
forest_panel(ax, rows, "a  Breast cancer: leave one cohort out")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="s", color=GREY, mfc="white", ls="", label="original 13 cohorts, v1 training set"), Line2D([], [], marker="o", color=BLUE, ls="", label="original 13 cohorts, v2 training set (more cohorts)"), Line2D([], [], marker="o", color=ORANGE, ls="", label="new cohort (added in W2)")], loc="lower left", fontsize=8, frameon=False)
save(fig, "fig1_breast_loco_forest", """
**Figure 1. Breast cancer, one cohort held out at a time.** Each row is one patient cohort. For each row, the model was trained on all the other breast cohorts and then scored the held-out cohort. The dot is the AUC (1.0 = perfect, 0.5 = coin toss); the line is the 95% confidence interval. Grey open squares: the original 13 cohorts when the model was trained on the v1 table. Blue dots: the same cohorts after more breast cohorts were added to the training pool. Orange dots: new cohorts that were tested for the first time. The text next to each name gives the number of patients and the number of responders. What I think: the original cohorts keep their AUC pattern when the pool grows, and the new cohorts show a similar range.
""")
# ---- Figure 2: other cancers with at least one new cohort
others = ["ovarian", "rectal chemoradiotherapy", "bladder", "oesophagus", "colorectal", "melanoma", "lung"]
fig, axs = plt.subplots(len(others), 1, figsize=(8.5, 14), gridspec_kw={"height_ratios": [max(len(build_rows(c)), 1) for c in others]})
for i, (c, ax) in enumerate(zip(others, axs)):
    forest_panel(ax, build_rows(c), f"{'abcdefg'[i]}  {c.capitalize()}")
    if i < len(others)-1: ax.set_xlabel("")
axs[0].legend(handles=[Line2D([], [], marker="s", color=GREY, mfc="white", ls="", label="existing cohort, v1 training set"), Line2D([], [], marker="o", color=BLUE, ls="", label="existing cohort, v2 training set"), Line2D([], [], marker="o", color=ORANGE, ls="", label="new cohort (W2/W1)")], loc="lower left", fontsize=8, frameon=False)
fig.tight_layout()
save(fig, "fig2_other_cancers_loco_forest", """
**Figure 2. Other cancers, one cohort held out at a time.** Same layout as Figure 1, one panel per cancer type. Rectal chemoradiotherapy (radiation plus chemotherapy before surgery) is shown as its own group. The dot is the AUC with its 95% confidence interval. Grey open squares: the cohort scored with the v1 training table; blue: the same cohort with the v2 table; orange: a cohort that is new. What I think: cohorts with few responders have wide intervals, so the width of each line matters as much as the dot.
""")
# ---- Figure 3: randomised drug-specific response tests
p2 = pd.read_csv(T/"part2_v2.tsv", sep="\t"); p2 = p2[p2.table == "new_drug_specific_response"].copy()
orig = pd.read_csv(T/"part2_v2.tsv", sep="\t"); orig = orig[orig.table == "original_part2_rerun"]
allp2 = pd.concat([p2, orig])
fig, ax = plt.subplots(figsize=(9, 0.55*len(allp2)+1.8))
for i, r in enumerate(allp2.itertuples()):
    for arm, a, n, k, col, off in (("drug", r.auc_in_drug_arm, r.n_drug, r.resp_drug, ORANGE, -0.15), ("control", r.auc_in_control_arm, r.n_control, r.resp_control, GREY, 0.15)):
        if pd.isna(a) or k in (0, n): continue
        lo, hi = hm_ci(a, k, n-k); ax.errorbar(a, i+off, xerr=[[a-lo], [hi-a]], fmt="o", color=col, ms=4, lw=1)
ax.set_yticks(range(len(allp2))); ax.set_yticklabels([f"{t[:70]}\np(arm difference)={d:.3f}" for t, d in zip(allp2.test, allp2.p_auc_difference)], fontsize=7.5)
ax.invert_yaxis(); ax.axvline(0.5, color="k", ls=":", lw=0.8); ax.set_xlim(0, 1.02); ax.set_xlabel("AUC inside each arm (0.5 = coin toss)")
ax.legend(handles=[Line2D([], [], marker="o", color=ORANGE, ls="", label="drug arm"), Line2D([], [], marker="o", color=GREY, ls="", label="control arm")], loc="lower right", fontsize=8, frameon=False)
ax.set_title("Randomised comparisons: does the score separate responders in the drug arm and the control arm?", loc="left")
save(fig, "fig3_randomised_drug_specific_forest", """
**Figure 3. Randomised trials: AUC in the drug arm and in the control arm.** For each randomised comparison (NOAH trastuzumab, CHER-LOB lapatinib arms, APEX bortezomib, and the earlier I-SPY2 arms), the score is tested separately in the patients who got the drug and in those who got the control. Orange: drug arm. Grey: control arm. Lines are 95% confidence intervals. The p-value under each label compares the two AUCs. A score that is specific to the drug would show a higher AUC in the drug arm. What I think: the per-arm patient numbers are small, so intervals are wide; read the dots together with the n and p in REPORT.md.
""")
# ---- Figure 4: survival interaction forest
p3 = pd.read_csv(T/"part3_v2.tsv", sep="\t"); r = pd.read_csv(T/"new_randomised_survival_tests.tsv", sep="\t")
r = r.rename(columns={"test": "setting"}); r["test_cohort"] = r["endpoint"]
d = pd.concat([p3, r[["setting", "test_cohort", "n", "treated", "events", "interaction_hazard_ratio", "hr_ci_low", "hr_ci_high", "p_interaction"]]], ignore_index=True)
fig, ax = plt.subplots(figsize=(9, 0.28*len(d)+1.8)); y = 0; ticks = []; labs = []; last = None
for r_ in d.itertuples():
    if r_.setting != last and last is not None: y += 0.6
    last = r_.setting; col = ORANGE if "randomised" in r_.setting else BLUE
    ax.errorbar(r_.interaction_hazard_ratio, y, xerr=[[r_.interaction_hazard_ratio-r_.hr_ci_low], [r_.hr_ci_high-r_.interaction_hazard_ratio]], fmt="o", color=col, ms=4, lw=1)
    ticks.append(y); labs.append(f"{r_.setting[:48]} | {r_.test_cohort}  (n={r_.n}, events={r_.events})"); y += 1
ax.set_xscale("log"); ax.axvline(1, color="k", ls=":", lw=0.8); ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("Hazard ratio of treatment x score interaction (log scale; 1 = no interaction)"); ax.set_title("Survival: does the score change how much the drug helps?", loc="left")
save(fig, "fig4_survival_interaction_forest", """
**Figure 4. Survival interaction forest plot.** Each row is one cohort (or one endpoint of a randomised trial). The dot is the hazard ratio (HR) of the treatment-by-score interaction from a Cox model, with its 95% confidence interval. HR below 1 means the drug looks more helpful when the score is higher; above 1 means the opposite. The model for each row was trained without the cohort being scored. Blue: observational cohorts (treated compared with untreated). Orange: randomised trials (cross-fitted inside the trial). What I think: most intervals include 1; the number of events behind each row is printed next to its name and sets how wide the line is.
""")
print("figures done")
