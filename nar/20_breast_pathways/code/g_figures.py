import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, textwrap
from common import *
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "figure.dpi": 100, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})
COL = {"proliferation": "#D55E00", "immune": "#0072B2", "signalling": "#009E73", "metabolism": "#E69F00", "other": "#8c8c8c"}
def save(fig, name):
    for e in ("png", "pdf"): fig.savefig(FIG / f"{name}.{e}", bbox_inches="tight")
    plt.close(fig)
def lab(x, n=58): return x if len(x) <= n else x[:n - 1] + "…"
def panel(ax, l): ax.text(-0.12, 1.06, l, transform=ax.transAxes, fontsize=12, fontweight="bold")
A = pd.read_csv(TAB / "auc_vs_k.tsv", sep="\t"); G = pd.read_csv(TAB / "minimal_subset.tsv", sep="\t"); full = A[(A.ranking == "weight") & (A.k == 673)].mean_auc.iloc[0]
# 1 AUC vs k
fig, ax = plt.subplots(1, 2, figsize=(10, 4)); a = ax[0]
for rk, c, nm, mk in [("weight", "#0072B2", "ranked by model weight", "o"), ("univariate", "#009E73", "ranked by one-at-a-time strength", "s"), ("random", "#8c8c8c", "random pathways", "^")]:
    g = A[A.ranking == rk].sort_values("k"); a.plot(g.k, g.mean_auc, marker=mk, color=c, label=nm, ms=4); a.fill_between(g.k, g.ci_low, g.ci_high, color=c, alpha=.12)
a.axhline(full, color="k", ls="--", lw=.8); a.axhline(full - .01, color="k", ls=":", lw=.8); a.text(1.1, full + .004, f"all 673 pathways: {full:.3f}", fontsize=8); a.text(1.1, full - .017, "within 0.01 of all", fontsize=8)
a.set_xscale("log"); a.set_xticks([1, 2, 5, 10, 20, 50, 100, 200, 673]); a.set_xticklabels([1, 2, 5, 10, 20, 50, 100, 200, 673]); a.set_xlabel("number of pathways kept (log scale)"); a.set_ylabel("mean AUC over 13 held-out cohorts"); a.legend(frameon=False, fontsize=8, loc="lower right"); a.set_ylim(.5, .79); panel(a, "a")
b = ax[1]; b.plot(G.step, G.mean_test_auc, "o-", color="#D55E00", ms=4); b.fill_between(G.step, G.ci_low, G.ci_high, color="#D55E00", alpha=.15)
b.axhline(full, color="k", ls="--", lw=.8); b.axhline(full - .01, color="k", ls=":", lw=.8); b.set_xlabel("pathways added one by one (greedy, chosen in training cohorts only)"); b.set_ylabel("mean AUC over 13 held-out cohorts"); b.set_ylim(.5, .79); panel(b, "b")
save(fig, "fig_auc_vs_k")
# 2 top30 weights
S = pd.read_csv(TAB / "top30_pathways.tsv", sep="\t", index_col=0).iloc[::-1]
fig, ax = plt.subplots(figsize=(8.5, 8)); ax.barh(range(30), S.L2_weight, color=[COL[t] for t in S.theme]); ax.set_yticks(range(30)); ax.set_yticklabels([lab(i.replace("_", " ").lower()) for i in S.index], fontsize=8)
ax.axvline(0, color="k", lw=.6); ax.set_xlabel("model weight (positive = higher score goes with response)")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in COL.values()], labels=list(COL), frameon=False, loc="lower right", title="theme"); save(fig, "fig_top30_weights")
# 3 themes
F = pd.read_csv(TAB / "family_weights.tsv", sep="\t", index_col=0)
fig, ax = plt.subplots(1, 2, figsize=(10, 4)); ax[0].pie(F.total_abs_weight, labels=[f"{i}\n{v:.1%}" for i, v in zip(F.index, F.share_of_total_abs_weight)], colors=[COL[t] for t in F.index], startangle=90, counterclock=False, textprops={"fontsize": 8}); panel(ax[0], "a")
x = np.arange(len(F)); ax[1].bar(x - .2, F.share_of_total_abs_weight * 100, .4, color=[COL[t] for t in F.index], label="share of total |weight|"); ax[1].bar(x + .2, F.share_if_equal_per_pathway * 100, .4, color="white", edgecolor="k", hatch="//", label="share of pathways (what equal weights would give)")
ax[1].set_xticks(x); ax[1].set_xticklabels(F.index); ax[1].set_ylabel("percent"); ax[1].legend(frameon=False, fontsize=8); panel(ax[1], "b"); save(fig, "fig_themes")
# 4 stability heatmap
FF = pd.read_csv(TAB / "stability_by_fold.tsv", sep="\t", index_col=0); T30 = pd.read_csv(TAB / "top30_pathways.tsv", sep="\t", index_col=0)
H = FF.loc[T30.index].copy(); H.insert(0, "all data\n(100 boots)", T30.stability.values)
fig, ax = plt.subplots(figsize=(9, 8)); im = ax.imshow(H.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns, rotation=75, ha="right", fontsize=8); ax.set_yticks(range(30)); ax.set_yticklabels([lab(i.replace("_", " ").lower(), 55) for i in H.index], fontsize=8)
for t, tl in zip(ax.get_yticklabels(), T30.theme): t.set_color(COL[tl] if tl != "other" else "#555555")
ax.axvline(.5, color="w", lw=2); cb = plt.colorbar(im, ax=ax, shrink=.6); cb.set_label("fraction of bootstrap fits that kept the pathway (L1 model)"); ax.set_xlabel("held-out cohort (training = the other 12; 20 bootstraps each)"); save(fig, "fig_stability")
# 5 forest
O = pd.read_csv(TAB / "subgroup_auc_pooled.tsv", sep="\t"); order = ["All patients", "ER+/HER2-", "HER2+", "TNBC"]; ms = [("full", "#0072B2", "All 673 pathways"), ("top100", "#009E73", "Top 100 pathways"), ("resid", "#E69F00", "673, proliferation regressed out"), ("prolif", "#D55E00", "Proliferation score only")]
fig, ax = plt.subplots(figsize=(8, 5)); y0 = 0; yt = []; yl = []
for gname in order:
    for j, (m, c, nm) in enumerate(ms):
        r = O[(O.group == gname) & (O.model == m)].iloc[0]; yy = y0 + j * .7
        ax.errorbar(r.mean_auc, yy, xerr=[[r.mean_auc - r.ci_low], [r.ci_high - r.mean_auc]], fmt="o", color=c, capsize=3, label=nm if gname == order[0] else None)
    yt.append(y0 + 1.05); yl.append(f"{gname}\n({int(r.n_patients)} pts, {int(r.n_cohorts)} cohorts)"); y0 += 3.6
ax.axvline(.5, color="k", ls=":", lw=.8); ax.set_yticks(yt); ax.set_yticklabels(yl); ax.invert_yaxis(); ax.set_xlabel("mean AUC over held-out cohorts (95% CI)"); ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.14), ncol=2); save(fig, "fig_subgroup_forest")
# 6 proliferation
P = pd.read_csv(TAB / "proliferation_per_cohort.tsv", sep="\t"); fig, ax = plt.subplots(1, 2, figsize=(10, 4)); x = np.arange(len(P))
for j, (k, c, nm) in enumerate([("auc_full", "#0072B2", "all 673"), ("auc_resid", "#E69F00", "673, proliferation regressed out"), ("auc_prolif", "#D55E00", "proliferation score only")]): ax[0].bar(x + (j - 1) * .27, P[k], .27, color=c, label=nm)
ax[0].set_xticks(x); ax[0].set_xticklabels(P.test_cohort, rotation=70, ha="right", fontsize=7); ax[0].axhline(.5, color="k", ls=":", lw=.8); ax[0].set_ylim(.45, .9); ax[0].set_ylabel("AUC in held-out cohort"); ax[0].legend(frameon=False, fontsize=7, loc="upper left", ncol=1); panel(ax[0], "a")
ax[1].bar(x - .2, P.spearman_model_vs_prolif, .4, color="#0072B2", label="full model vs proliferation"); ax[1].bar(x + .2, P.spearman_resid_vs_prolif, .4, color="#E69F00", label="regressed-out model vs proliferation")
ax[1].set_xticks(x); ax[1].set_xticklabels(P.test_cohort, rotation=70, ha="right", fontsize=7); ax[1].axhline(0, color="k", lw=.6); ax[1].set_ylabel("Spearman correlation"); ax[1].legend(frameon=False, fontsize=7); panel(ax[1], "b"); save(fig, "fig_proliferation")
