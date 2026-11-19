"""Figures (PNG 300 dpi + PDF + caption) from tables/.  Run after the chain finishes."""
import sys, pickle
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = Path(__file__).resolve().parents[1]; TAB = OUT / "tables"; FIG = OUT / "figures"; DATA = OUT / "_data"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#888888"

def save(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption + "\n"); plt.close(fig)

def forest(ax, df, label, hr, lo, hi, color=BLUE, p=None):
    y = np.arange(len(df))[::-1]
    ax.errorbar(df[hr], y, xerr=[df[hr] - df[lo], df[hi] - df[hr]], fmt="o", color=color, ecolor=color, capsize=2, ms=4)
    ax.axvline(1, color=GREY, lw=.8, ls="--"); ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(df[label], fontsize=7)

# Fig 1: lung cluster x drug values + cluster classifier
v = pd.read_csv(TAB / "luad_cluster_by_drug_value_listed.tsv", sep="\t", index_col=0)[["C1", "C2", "C3", "C4"]]
cv = pd.read_csv(TAB / "lung_cluster_classifier_cv.tsv", sep="\t") if (TAB / "lung_cluster_classifier_cv.tsv").exists() else None
fig, ax = plt.subplots(1, 2, figsize=(10, 6.5), gridspec_kw={"width_ratios": [1.3, 1]})
im = ax[0].imshow(v.values.astype(float), cmap="RdBu_r", vmin=-.8, vmax=.8, aspect="auto")
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(v.columns); ax[0].set_yticks(range(len(v))); ax[0].set_yticklabels(v.index, fontsize=7)
for i in range(v.shape[0]):
    for j in range(4):
        if pd.notna(v.iat[i, j]): ax[0].text(j, i, f"{v.iat[i,j]:.2f}", ha="center", va="center", fontsize=6)
ax[0].set_title("a  Value listed per cluster and drug\n(blank = no entry)"); fig.colorbar(im, ax=ax[0], shrink=.5, label="value (log2FC)")
if cv is not None:
    r = cv.iloc[0]; ax[1].bar(["chance", "shuffled labels", "pathway scores\n(5-fold CV)"], [.25, r.shuffled_balanced_accuracy_mean, r.balanced_accuracy], color=[GREY, GREY, BLUE])
    ax[1].set_ylabel("balanced accuracy"); ax[1].set_title(f"b  Can pathway scores tell the 4 clusters apart?\nshuffle p = {r.permutation_p:.3f}")
save(fig, "fig1_lung_cluster_drug_table", "a) The drug table gives one number per cluster and drug (102 patient-derived cell lines; 13 in C1, 45 in C2, 19 in C3, 25 in C4), not one per line. Blank cells have no entry. b) Cross-validated accuracy of a classifier that predicts the cluster from Reactome pathway scores, compared with chance and with shuffled cluster labels.")

# forest of apply tables
def apply_forest(plan, name, title):
    f = TAB / f"apply_{plan}.tsv"
    if not f.exists(): return
    d = pd.read_csv(f, sep="\t"); d = d[(d.score == "rank1") & d.HR.notna()]
    if d.empty: return
    d["lab"] = d.organoid_drug + " -> " + d.regimen + " | " + d.cohort + f" (n={0})"
    d["lab"] = [f"{a} > {b} | {c} n={int(n)}/{int(e)} ev" for a, b, c, n, e in zip(d.organoid_drug, d.regimen, d.cohort, d.n_treated, d.events_treated)]
    fig, ax = plt.subplots(figsize=(8, max(3, .22 * len(d) + 1))); forest(ax, d, "lab", "HR", "HR_lo", "HR_hi")
    ax.set_xlabel("hazard ratio, predicted less sensitive vs more sensitive (95% CI)"); ax.set_title(title)
    save(fig, name, f"{title}. Each row is one organoid drug used as the model, one patient group and one cohort. HR above 1 = patients predicted less sensitive did worse (the direction the method expects). Table: tables/apply_{plan}.tsv.")
apply_forest("lung", "fig2_lung_apply_forest", "Lung: organoid-only model (paper method, top pathway) in treated patients")
apply_forest("ovary", "fig5_ovary_apply_forest", "Ovary: PARP-drug model in platinum / paclitaxel treated patients")

# joint
def joint_fig(plan, name):
    import glob
    fs = glob.glob(str(DATA / "joint" / plan / "tables" / "*__eval.tsv"))
    if not fs: return
    d = pd.concat([pd.read_csv(f, sep="\t") for f in fs]); d = d[(d.cohort == "POOLED") & (d.outcome == "survival") & d.HR.notna()] if "HR" in d else pd.DataFrame()
    if d.empty: return
    d.to_csv(TAB / f"joint_{plan}_pooled.tsv", sep="\t", index=False)
    d["lab"] = d.unit.str.replace("_vs_platinum_proxy", " > platinum") + " | " + d.model
    fig, ax = plt.subplots(figsize=(8, max(3, .2 * len(d) + 1))); forest(ax, d, "lab", "HR", "HR_low" if "HR_low" in d else "lo", "HR_high" if "HR_high" in d else "hi")
    ax.set_xlabel("pooled HR per 1 SD of predicted benefit (95% CI); below 1 = expected direction")
    save(fig, name, f"{plan}: organoid + patient (joint) models, leave-one-cohort-out, pooled over held-out cohorts. Table: tables/joint_{plan}_pooled.tsv.")
try: joint_fig("lung", "fig3_lung_joint_pooled")
except Exception as e: print("joint lung fig", e)
try: joint_fig("ovary", "fig6_ovary_joint_pooled")
except Exception as e: print("joint ovary fig", e)

# cluster expected value forest
f = TAB / "lung_cluster_expected_value_survival.tsv"
if f.exists():
    d = pd.read_csv(f, sep="\t"); d = d[d.HR_treated.notna()]
    d["lab"] = [f"{a} {b} {c} n={int(n)}" for a, b, c, n in zip(d.organoid_drug, d.regimen, d.cohort, d.n_treated)]
    fig, ax = plt.subplots(figsize=(8, max(3, .22 * len(d) + 1))); forest(ax, d, "lab", "HR_treated", "HR_lo", "HR_hi")
    ax.set_xlabel("HR, predicted less sensitive vs more sensitive (95% CI)"); ax.set_title("Lung: expected drug value from predicted cluster")
    save(fig, "fig4_lung_cluster_expected_value", "Patient cluster probabilities (classifier trained on the 102 cell lines) turned into an expected drug value (sum of probability times the cluster value) and split at the median in treated patients. Table: tables/lung_cluster_expected_value_survival.tsv.")

# ovary: ic50 + signature
r = pd.read_csv(TAB / "ovary_ic50_raw_7lines.tsv", sep="\t", index_col=0) if (TAB / "ovary_ic50_raw_7lines.tsv").exists() else None
s = TAB / "ovary_signature_survival.tsv"
if r is not None:
    fig, ax = plt.subplots(1, 2, figsize=(10, 4)); r.plot.bar(ax=ax[0], color=[BLUE, ORANGE]); ax[0].set_ylabel("IC50 (as given)"); ax[0].set_title("a  IC50 per line"); ax[0].tick_params(axis="x", rotation=60)
    if s.exists():
        d = pd.read_csv(s, sep="\t"); d = d[d.set.str.startswith("6")]
        d["lab"] = [f"{a} {c} {e} n={int(n)}" for a, c, e, n in zip(d.drug, d.cohort, d.endpoint, d.n_treated)]
        forest(ax[1], d, "lab", "HR", "HR_lo", "HR_hi"); ax[1].set_title("b  Signature in platinum-treated patients"); ax[1].set_xlabel("HR (95% CI)")
    save(fig, "fig7_ovary_ic50_signature", "a) Olaparib and niraparib IC50 for the 7 lines (MCF7 is a breast line). b) Signature built from the 10 pathways most positively and 10 most negatively correlated with ln IC50 across 6 ovarian lines, tested in platinum-treated patients; HR above 1 = high score did worse.")

# liver heatmap
f = TAB / "liver_log10_ic50_centred_per_drug.tsv"
if f.exists():
    l = pd.read_csv(f, sep="\t", index_col=0); fig, ax = plt.subplots(figsize=(4.5, 4.5)); im = ax.imshow(l.values, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(l.columns); ax.set_yticks(range(len(l))); ax.set_yticklabels(l.index, fontsize=7); fig.colorbar(im, label="log10 IC50 minus drug mean")
    save(fig, "fig8_liver_ic50", "Median IC50 of 12 drugs in 3 liver organoids (72.3, 2E, CC3), centred on each drug's mean. Red = less sensitive than the other two. No expression data exist for these organoids in the uploaded files.")
