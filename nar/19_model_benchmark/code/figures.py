import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 100})
PRED = TAB / "preds"
R = pd.read_csv(TAB / "breast_loco_all_models.tsv", sep="\t"); C = pd.read_csv(TAB / "model_comparison.tsv", sep="\t")
RS = pd.read_csv(TAB / "random_pathway_sets_auc.tsv", sep="\t")
BASE = ["ER_HER2_only", "Proliferation_score"]
LABEL = {"L2_fixed": "L2 (existing, lambda=100)", "L2": "L2 (tuned)", "L1": "L1", "ElasticNet": "Elastic net", "SVM_linear": "Linear SVM",
         "SVM_RBF": "RBF SVM", "RandomForest": "Random forest", "XGBoost": "XGBoost", "MLP": "MLP (1 layer)", "NetSmoothL2": "Network-smoothed L2",
         "GCN_prior": "Graph network (GCN, proximity prior)", "ER_HER2_only": "ER + HER2 only", "Proliferation_score": "Proliferation score",
         "L2_shortlist": "L2, proximity shortlist", "L1_shortlist": "L1, proximity shortlist", "ElasticNet_shortlist": "Elastic net, proximity shortlist",
         "SVM_linear_shortlist": "Linear SVM, proximity shortlist", "L2_fixed_shortlist": "L2 lambda=100, proximity shortlist", "RandomSets_mean": "Random pathway sets (mean of 200)"}
lab = lambda m: LABEL.get(m, m)
# add random-set mean as a pseudo model for heatmap / bars
rm = RS.groupby("test_cohort").auc.agg(["mean", lambda x: np.percentile(x, 2.5), lambda x: np.percentile(x, 97.5)])
rm.columns = ["auc", "lo", "hi"]
A = R.pivot(index="model", columns="test_cohort", values="auc")
A.loc["RandomSets_mean"] = rm["auc"]
A["Mean"] = A.mean(axis=1)
order = A["Mean"].sort_values(ascending=False).index
A = A.loc[order]
# (1) heatmap
fig, ax = plt.subplots(figsize=(12, 0.36 * len(A) + 2))
im = ax.imshow(A.values, cmap="viridis", vmin=0.45, vmax=0.9, aspect="auto")
ax.set_xticks(range(A.shape[1])); ax.set_xticklabels(A.columns, rotation=60, ha="right")
ax.set_yticks(range(A.shape[0])); ax.set_yticklabels([lab(m) for m in A.index])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        v = A.values[i, j]
        if np.isfinite(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if v < 0.68 else "black")
ax.axvline(A.shape[1] - 1.5, color="white", lw=3)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("AUC (0.5 = coin flip)")
ax.text(-0.09, 1.02, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")
ax.set_title("AUC for each model on each held-out breast cohort", loc="left")
fig.tight_layout()
def save(fig, stem):
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight"); fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight"); plt.close(fig)
save(fig, "fig_model_heatmap")
# (2) mean AUC + CI
C2 = C.copy(); C2["lab"] = C2.model.map(lab)
rrow = dict(model="RandomSets_mean", mean_auc=rm.auc.mean(), mean_ci_low=np.nan, mean_ci_high=np.nan, lab=lab("RandomSets_mean"))
C2 = pd.concat([C2, pd.DataFrame([rrow])], ignore_index=True).sort_values("mean_auc")
ref = C.set_index("model").loc["L2_fixed"]
fig, ax = plt.subplots(figsize=(8, 0.34 * len(C2) + 1.8))
colr = ["#999999" if m in BASE or m == "RandomSets_mean" else ("#0072B2" if m.endswith("shortlist") else ("#D55E00" if m == "L2_fixed" else "#009E73")) for m in C2.model]
ax.barh(range(len(C2)), C2.mean_auc, color=colr, alpha=0.9)
ok = C2.mean_ci_low.notna()
ax.errorbar(C2.mean_auc[ok], np.arange(len(C2))[ok], xerr=[(C2.mean_auc - C2.mean_ci_low)[ok], (C2.mean_ci_high - C2.mean_auc)[ok]], fmt="none", ecolor="black", capsize=2, lw=1)
ax.axvline(ref.mean_auc, color="#D55E00", ls="--", lw=1); ax.axvline(0.5, color="k", lw=0.8)
for i, (v, m) in enumerate(zip(C2.mean_auc, C2.model)): ax.text(0.452, i, f"{v:.3f}", va="center", fontsize=8, color="white")
ax.set_yticks(range(len(C2))); ax.set_yticklabels(C2.lab); ax.set_xlim(0.45, 0.8); ax.set_xlabel("Mean AUC over the 13 held-out cohorts (bars: 95% CI)")
ax.text(-0.5, 1.02, "b", transform=ax.transAxes, fontsize=14, fontweight="bold")
fig.tight_layout(); save(fig, "fig_model_mean_auc")
# (3) ROC of the 3 best (non-baseline) models + L2 existing
cand = C[~C.model.isin(BASE + ["L2_fixed"])].sort_values("mean_auc", ascending=False).model.tolist()[:3]
grid = np.linspace(0, 1, 101)
def curves(m):
    p = pd.read_csv(PRED / f"{m}.tsv.gz", sep="\t", index_col=0); out = {}
    for c, g in p.groupby("cohort"):
        f, t, _ = roc_curve(g.y, g.score); out[c] = np.interp(grid, f, t)
    return pd.DataFrame(out)
ref_c = curves("L2_fixed").mean(axis=1)
fig, axs = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True)
for k, (ax, m) in enumerate(zip(axs, cand)):
    cc = curves(m)
    for c in cc: ax.plot(grid, cc[c], color="#56B4E9", lw=0.7, alpha=0.6)
    ax.plot(grid, cc.mean(axis=1), color="#0072B2", lw=2.5, label=f"{lab(m)} (mean of cohorts)")
    ax.plot(grid, ref_c, color="#D55E00", lw=1.8, ls="--", label="L2 existing (mean)")
    ax.plot([0, 1], [0, 1], color="k", lw=0.7, ls=":")
    ax.set_xlabel("False positive rate"); ax.text(-0.02, 1.04, "abc"[k], transform=ax.transAxes, fontsize=14, fontweight="bold")
    ax.set_title(lab(m), fontsize=9, loc="left"); ax.legend(loc="lower right", fontsize=7, frameon=False)
axs[0].set_ylabel("True positive rate")
fig.tight_layout(); save(fig, "fig_roc_top3")
print("top3:", cand)
