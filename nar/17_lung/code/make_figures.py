"""Figures for W1 (lung): reads only the tables written by run_lung.py. PNG 300 dpi + PDF; captions are written next to them."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.duration.survfunc import SurvfuncRight

ROOT = Path(__file__).resolve().parents[3]
TAB, FIG = ROOT / "nar/17_lung/tables", ROOT / "nar/17_lung/figures"
FIG.mkdir(exist_ok=True, parents=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                     "legend.fontsize": 8.5, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
C = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9", "grey": "#666666"}


def save(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")
    plt.close(fig)


def letter(ax, s):
    ax.text(-0.02, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom", ha="right")


def forest(ax, rows, ycol, hr, lo, hi, color, off=0.0, label=None, marker="o"):
    ax.errorbar(rows[hr], rows[ycol] + off, xerr=[rows[hr] - rows[lo], rows[hi] - rows[hr]], fmt=marker, color=color, ms=5, lw=1.4,
                capsize=2, label=label)


# ------------------------------------------------------------------ figure 1: adjuvant interaction
B = pd.read_csv(TAB / "lung_adjuvant_interaction.tsv", sep="\t")
B = B[B["interaction_hr"].notna()].copy()
cohorts = [c for c in B["test_cohort"].unique() if not c.startswith("ALL")] + [c for c in B["test_cohort"].unique() if c.startswith("ALL")]
pos = {c: len(cohorts) - 1 - i for i, c in enumerate(cohorts)}
B["y"] = B["test_cohort"].map(pos)
sets = list(B["pathways"].unique())
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]
for s, col, off, mk in zip(sets, [C["blue"], C["orange"]], [0.13, -0.13], ["o", "s"]):
    forest(ax, B[B["pathways"] == s], "y", "interaction_hr", "ci_lo", "ci_hi", col, off, s, mk)
ax.axvline(1, color="k", lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"])
lab = []
for c in cohorts:
    r = B[B["test_cohort"] == c].iloc[0]
    nm = "All tested (pooled)" if c.startswith("ALL") else c + (" (randomised)" if r.get("randomised") is True or str(r.get("randomised")) == "True" else "")
    lab.append(f"{nm}\nn={int(r['n'])}, events={int(r['events'])}")
ax.set_yticks([pos[c] for c in cohorts]); ax.set_yticklabels(lab)
ax.set_xlabel("Chemo x score interaction, hazard ratio per 1 SD of predicted benefit (log scale)\nbelow 1 = chemo works better when the score is high")
ax.legend(loc="lower right", frameon=False); letter(ax, "a")
ax = axs[1]
H = B[B["pathways"] == sets[0]]
H = H[~H["test_cohort"].str.startswith("ALL")]
yy = H["test_cohort"].map(pos)
ax.errorbar(H["hr_high_benefit_half"], yy + 0.13, xerr=[H["hr_high_benefit_half"] - H["high_lo"], H["high_hi"] - H["hr_high_benefit_half"]], fmt="o", color=C["green"], label="predicted high benefit (top half)", ms=5, lw=1.4, capsize=2)
ax.errorbar(H["hr_low_benefit_half"], yy - 0.13, xerr=[H["hr_low_benefit_half"] - H["low_lo"], H["low_hi"] - H["hr_low_benefit_half"]], fmt="s", color=C["vermillion"], label="predicted low benefit (bottom half)", ms=5, lw=1.4, capsize=2)
ax.axvline(1, color="k", lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"])
ax.set_yticks([pos[c] for c in cohorts]); ax.set_yticklabels([""] * len(cohorts)); ax.set_ylim(*axs[0].get_ylim())
ax.set_xlabel("Hazard ratio, chemo vs no chemo (log scale)\nbelow 1 = chemo lengthens survival")
ax.legend(loc="lower right", frameon=False); letter(ax, "b")
fig.tight_layout()
save(fig, "lung_forest_interaction", """
**Does a pathway score tell us who benefits from chemotherapy after lung surgery?**
Each row is one lung cohort that was kept out of training (the model was learned on the other cohorts).
**a**, the "interaction" hazard ratio (HR, how much the chemo effect changes per 1 standard deviation of the score; the dashed line at 1 means no change)
with 95% confidence interval. Blue: the score uses all 673 Reactome pathways. Orange: the score uses only the 47 pathways close to cisplatin's targets in the STRING network (the paper's shortlist).
The bottom row pools the tested cohorts. GSE14814 is the randomised trial, so it is the cleanest test.
**b**, the plain chemo-vs-no-chemo hazard ratio inside each half of the patients, split at the median predicted benefit (all-pathway score). If the score works, green should sit left of red.
n = patients, events = deaths in the held-out cohort.
""")

# ------------------------------------------------------------------ figure 2: AUC by cohort
A = pd.read_csv(TAB / "lung_loco_response.tsv", sep="\t")
A["ci_hi"], A["ci_lo"] = A["ci_hi"].clip(upper=1), A["ci_lo"].clip(lower=0)
cs = list(dict.fromkeys(A.sort_values(["group", "test_cohort"])["test_cohort"]))
models = list(dict.fromkeys(A["model"]))
cols = dict(zip(models, [C["blue"], C["orange"], C["green"], C["purple"]]))
fig, ax = plt.subplots(figsize=(8.6, 5.6))
for j, m in enumerate(models):
    d = A[A["model"] == m]
    y = d["test_cohort"].map({c: len(cs) - 1 - i for i, c in enumerate(cs)}) + (1.5 - j) * 0.16
    ax.errorbar(d["auc"], y, xerr=[d["auc"] - d["ci_lo"], d["ci_hi"] - d["auc"]], fmt="o", color=cols[m], ms=4.5, lw=1.3, capsize=2, label=m)
ax.axvline(0.5, color="k", lw=0.8, ls="--")
g = A.drop_duplicates("test_cohort").set_index("test_cohort")
ax.set_yticks(range(len(cs))); ax.set_yticklabels([f"{c}\n{g.loc[c, 'group'].split(' (')[0]}, n={int(g.loc[c, 'n'])}, responders={int(g.loc[c, 'responders'])}" for c in cs][::-1])
ax.set_xlabel("AUC for predicting response in the held-out cohort (0.5 = chance; whiskers = 95% CI, cut at 0 and 1)")
ax.set_xlim(0, 1)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), frameon=False, ncol=1)
save(fig, "lung_auc_by_cohort", """
**Can pathway scores predict drug response in lung patients?** One row per lung cohort with a response label. Each dot is the area under the ROC curve
(AUC: the chance that a random responder gets a higher score than a random non-responder; 0.5 is a coin flip, 1 is perfect) with its 95% confidence interval.
Four ridge-regression models were each tested on the cohort without ever seeing it: one trained on the other cohorts of the same treatment group
(immunotherapy cohorts together, TCGA-LUAD and TCGA-LUSC together), one trained on all breast response cohorts, one trained on every other cohort of any cancer,
and one trained on every cancer except lung. GSE33072 has no same-group partner. GSE126044 (16 patients) and GSE207422 (24 patients, pathological response after neoadjuvant treatment) are very small.
""")

# ------------------------------------------------------------------ figure 3: untreated prognosis
U = pd.read_csv(TAB / "lung_untreated_prognosis.tsv", sep="\t")
U = U[U["hr_per_sd"].notna()].copy()
names = [c for c in U["test_cohort"].unique() if not c.startswith("ALL")] + ["ALL TESTED COHORTS (inverse-variance pool)"]
pos = {c: len(names) - 1 - i for i, c in enumerate(names)}
U["y"] = U["test_cohort"].map(pos)
fig, axs = plt.subplots(1, 2, figsize=(11, 6.2), sharey=True)
for ax, sc, L in zip(axs, ["prognostic score (main-effect part)", "benefit score (drug-specific part)"], "ab"):
    for s, col, off, mk in zip(sets, [C["blue"], C["orange"]], [0.13, -0.13], ["o", "s"]):
        forest(ax, U[(U["score"] == sc) & (U["pathways"] == s)], "y", "hr_per_sd", "ci_lo", "ci_hi", col, off, s, mk)
    ax.axvline(1, color="k", lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_xticks([0.5, 0.7, 1, 1.5, 2]); ax.set_xticklabels(["0.5", "0.7", "1", "1.5", "2"])
    ax.set_title(sc.split(" (")[0].capitalize()); ax.set_xlabel("Hazard ratio of death (or relapse) per 1 SD of the score\n(log scale; below 1 = higher score, longer survival)")
    letter(ax, L)
d0 = U.drop_duplicates("test_cohort").set_index("test_cohort")
axs[0].set_yticks([pos[c] for c in names])
axs[0].set_yticklabels([("All tested (pooled)" if c.startswith("ALL") else f"{c}\n{'no chemo recorded' if d0.loc[c, 'group'].startswith('no') else 'chemo not recorded'}, n={int(d0.loc[c, 'n'])}, events={int(d0.loc[c, 'events'])}") for c in names])
axs[0].legend(loc="lower right", frameon=False)
fig.tight_layout()
save(fig, "lung_untreated_prognosis", """
**Is the score just a general prognosis marker?** Each row is a lung cohort, looking only at patients who did not receive chemotherapy
(either recorded as "no chemo", or surgery-only cohorts in which treatment was not written down). The score was learned on the other cohorts with the chemo/no-chemo data
and then applied here. The dot is the hazard ratio (HR) per 1 standard deviation of the score with its 95% confidence interval; below 1 means a higher score goes with longer survival.
**a**, the score's general (main-effect) part: it asks "does the score predict outcome whatever the treatment?". **b**, the drug-specific part, the one used for chemo benefit in the other figure.
Blue: all 673 pathways. Orange: the 47 cisplatin-near pathways. The bottom row pools the cohorts. GSE8894 gives relapse-free survival; the others give overall survival.
""")

# ------------------------------------------------------------------ figure 4: GSE14814 KM
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_lung as RL  # noqa: E402
P, Z = RL.TC.load()
fr = RL.frame(P, ["GSE14814"])
sc = pd.read_pickle(TAB / "_scores_gse14814.pkl")
key = [k for k in sc if k[1] == "GSE14814" and k[0].startswith("all")][0]
fr["score"] = sc[key]
fr["half"] = np.where(fr["score"] >= fr["score"].median(), "predicted high benefit", "predicted low benefit")
fig, axs = plt.subplots(1, 2, figsize=(9.4, 4.0), sharey=True)
for ax, h, L in zip(axs, ["predicted high benefit", "predicted low benefit"], "ab"):
    for t, lab, col in ((1, "chemo (cisplatin + vinorelbine)", C["vermillion"]), (0, "observation", C["blue"])):
        d = fr[(fr["half"] == h) & (fr["treated"] == t)]
        s = SurvfuncRight(d["time"].values / 12, d["event"].values.astype(int))
        ax.step(np.r_[0, s.surv_times], np.r_[1, s.surv_prob], where="post", color=col, lw=1.8, label=f"{lab}, n={len(d)}, deaths={int(d['event'].sum())}")
    ax.set_title(h.capitalize()); ax.set_xlabel("Years since surgery"); ax.set_ylim(0, 1.02); ax.legend(loc="lower left", frameon=False); letter(ax, L)
axs[0].set_ylabel("Share of patients alive")
fig.tight_layout()
save(fig, "lung_gse14814_km", """
**The randomised trial (GSE14814, JBR.10): survival by treatment and predicted benefit.** Patients were randomly given chemotherapy or observation after lung surgery.
The score was learned only on the other lung cohorts (never on these patients), using all 673 pathways, then patients were split at the median score.
**a**, the half predicted to benefit more from chemo; **b**, the other half. If the score picked who benefits, the red line would sit clearly above the blue one in **a** and not in **b**.
The numbers of patients and deaths are in the legends. The test statistics are in the table lung_adjuvant_interaction.tsv.
""")
print("figures done")
