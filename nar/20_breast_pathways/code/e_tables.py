"""Build the summary tables for Q1-Q5 from the files written by a-d."""
import statsmodels.api as sm
from scipy.stats import norm, binomtest
from themes import *
d, Z = breast_data(); y = d["responder"].astype(int).values; coh = d["cohort"].values; cols = list(Z.columns); X = Z.values
cohorts = sorted(set(coh)); short = [c.replace("REACTOME_", "") for c in cols]
def mean_ci(per):   # per: DataFrame with columns auc, n1, n0 -> plain mean over cohorts, CI from summed variances (as W3)
    v = np.array([hm_se(a, n1, n0) ** 2 for a, n1, n0 in zip(per.auc, per.n1, per.n0)]); m = per.auc.mean(); se = np.sqrt(v.sum()) / len(v)
    return m, m - 1.96 * se, m + 1.96 * se
cnt = pd.DataFrame({"n1": pd.Series(y).groupby(coh).sum(), "n0": (1 - pd.Series(y)).groupby(coh).sum()})
# ---------------- themes
W = pd.read_csv(TAB / "l2_weights_per_fold.tsv", sep="\t", index_col=0)
old = pd.read_csv(ROOT / "robust/results/train_combined/breast_model_pathways.tsv", sep="\t", index_col=0)
w = W["full"].copy(); w.index = short
print("max |difference| vs robust/ weights:", float((w - old["model_weight"].reindex(w.index)).abs().max()))
T = pd.DataFrame({"weight": w.values, "theme": [theme(c) for c in cols], "broad_cellcycle": [broad_prolif(c) for c in cols]}, index=short)
T.index.name = "pathway"; T.to_csv(TAB / "pathway_themes.tsv", sep="\t")
tot = T.weight.abs().sum()
F = T.assign(a=T.weight.abs()).groupby("theme").agg(n_pathways=("a", "size"), total_abs_weight=("a", "sum"), net_weight=("weight", "sum"))
F["share_of_total_abs_weight"] = F.total_abs_weight / tot
F["share_if_equal_per_pathway"] = F.n_pathways / len(T)
F["share_in_top30"] = T.assign(a=T.weight.abs()).sort_values("a", ascending=False).head(30).groupby("theme").size().reindex(F.index).fillna(0) / 30
# per-fold range of the proliferation / immune shares
sh = {}
for f in cohorts + ["full"]:
    a = W[f].abs().values; th = np.array([theme(c) for c in cols]); bp = np.array([broad_prolif(c) for c in cols])
    sh[f] = {t: a[th == t].sum() / a.sum() for t in F.index} | {"broad_cellcycle": a[bp].sum() / a.sum()}
SH = pd.DataFrame(sh).T
F["share_min_over_LOO_folds"] = SH.loc[cohorts, F.index].min(); F["share_max_over_LOO_folds"] = SH.loc[cohorts, F.index].max()
F = F.sort_values("total_abs_weight", ascending=False); F.to_csv(TAB / "family_weights.tsv", sep="\t")
SH.to_csv(TAB / "theme_share_per_fold.tsv", sep="\t")
print(F.round(3).to_string()); print("broad cell-cycle share:", SH.loc["full", "broad_cellcycle"], SH.loc[cohorts, "broad_cellcycle"].min(), SH.loc[cohorts, "broad_cellcycle"].max())
# ---------------- Q1 auc vs k
R = pd.read_csv(TAB / "auc_vs_k_per_cohort.tsv", sep="\t").join(cnt, on="test_cohort")
rows = []
for (rk, k), g in R.groupby(["ranking", "k"]):
    m, lo, hi = mean_ci(g); rows.append(dict(ranking=rk, k=k, mean_auc=m, ci_low=lo, ci_high=hi, min_cohort_auc=g.auc.min(), max_cohort_auc=g.auc.max()))
A = pd.DataFrame(rows)
full = A[(A.ranking == "weight") & (A.k == 673)].mean_auc.iloc[0]
A["diff_vs_all_673"] = A.mean_auc - full; A["within_0.01_of_all"] = A.diff_vs_all_673 >= -0.01
A.sort_values(["ranking", "k"]).to_csv(TAB / "auc_vs_k.tsv", sep="\t", index=False); print(A.sort_values(["ranking", "k"]).round(3).to_string())
# ---------------- Q2 stability
B = np.load(TAB / "_l1_boot_coefs.npy"); lab = pd.read_csv(TAB / "_l1_boot_labels.tsv", sep="\t").iloc[:, 0].values
sel = B != 0
fa = sel[lab == "all"].mean(0)
fold_freq = pd.DataFrame({c: sel[lab == c].mean(0) for c in cohorts}, index=short); fold_freq.to_csv(TAB / "stability_by_fold.tsv", sep="\t")
sgn_all = np.sign(B[lab == "all"]); same = np.array([(sgn_all[:, j] == np.sign(w.values[j])).sum() / max((sgn_all[:, j] != 0).sum(), 1) for j in range(len(cols))])
L2s = np.sign(W[cohorts].values); sign_ok = (L2s == np.sign(W["full"].values)[:, None]).sum(1)
S = pd.DataFrame({"theme": T.theme.values, "L2_weight": w.values, "sign": np.where(w.values > 0, "+ (higher score, more response)", "- (higher score, less response)"),
                  "raw_odds_ratio_per_sd": old["raw_odds_ratio_per_sd"].reindex(w.index).values, "raw_q": old["raw_q"].reindex(w.index).values,
                  "adj_odds_ratio_per_sd_er_her2": old["adj_odds_ratio_per_sd"].reindex(w.index).values,
                  "L1_selected_in_100_bootstraps": fa, "L1_selected_in_fold_bootstraps_mean": fold_freq.mean(1).values,
                  "L1_same_sign_when_selected": same, "L2_same_sign_in_13_LOO_folds": sign_ok}, index=short)
S.index.name = "pathway"; S["stability"] = S.L1_selected_in_100_bootstraps
S = S.reindex(S.L2_weight.abs().sort_values(ascending=False).index); S["rank_by_abs_weight"] = np.arange(1, len(S) + 1)
S.to_csv(TAB / "stability_selection.tsv", sep="\t")
S.head(30).to_csv(TAB / "top30_pathways.tsv", sep="\t")
print("selected in >=80% of the 100 bootstraps:", int((fa >= .8).sum()), "| >=50%:", int((fa >= .5).sum()), "| mean selected per fit:", sel[lab == "all"].sum(1).mean())
print("of top 30 by weight: median stability", S.head(30).stability.median(), "| >=80%:", int((S.head(30).stability >= .8).sum()))
print(S.head(30)[["theme", "L2_weight", "raw_odds_ratio_per_sd", "stability", "L2_same_sign_in_13_LOO_folds"]].round(3).to_string())
# ---------------- Q3 greedy
G = pd.read_csv(TAB / "greedy_path_per_fold.tsv", sep="\t").join(cnt, on="test_cohort"); rows = []
for st, g in G.groupby("step"):
    g = g.assign(auc=g.test_auc); m, lo, hi = mean_ci(g); rows.append(dict(step=st, mean_test_auc=m, ci_low=lo, ci_high=hi, mean_inner_auc=g.inner_auc.mean()))
GS = pd.DataFrame(rows); GS["diff_vs_all_673"] = GS.mean_test_auc - full; GS["within_0.01_of_all"] = GS.diff_vs_all_673 >= -0.01
# honest stopping rule per fold: smallest step whose INNER auc is within 0.01 of that fold's best inner auc
hon = []
for c, g in G.groupby("test_cohort"):
    best = g.inner_auc.max(); st = int(g[g.inner_auc >= best - 0.01].step.min()); hon.append(dict(test_cohort=c, auc=float(g[g.step == st].test_auc.iloc[0]), step_chosen=st))
H = pd.DataFrame(hon).join(cnt, on="test_cohort"); m, lo, hi = mean_ci(H)
GS.to_csv(TAB / "minimal_subset.tsv", sep="\t", index=False); H.to_csv(TAB / "minimal_subset_honest_stop.tsv", sep="\t", index=False)
print(GS.round(3).to_string()); print("honest stop rule: mean AUC %.3f (%.3f-%.3f); steps chosen %s" % (m, lo, hi, H.step_chosen.tolist()))
freq = G[G.step <= 5].groupby("added").test_cohort.nunique().sort_values(ascending=False); freq.rename("n_folds_among_first_5_picks").to_csv(TAB / "greedy_first5_pick_counts.tsv", sep="\t")
print(freq.head(10))
# ---------------- Q4 summary
P4 = pd.read_csv(TAB / "proliferation_per_cohort.tsv", sep="\t").join(cnt, on="test_cohort"); rows = []
for k in ["full", "prolif", "resid", "noprolif", "top20", "top50", "top100", "top200"]:
    g = P4.assign(auc=P4[f"auc_{k}"]); m, lo, hi = mean_ci(g)
    rows.append(dict(model=k, mean_auc=m, ci_low=lo, ci_high=hi, cohorts_better_than_full=int((P4[f"auc_{k}"] > P4.auc_full).sum()), mean_diff_vs_full=(P4[f"auc_{k}"] - P4.auc_full).mean()))
Q = pd.DataFrame(rows); Q.to_csv(TAB / "proliferation_tests.tsv", sep="\t", index=False); print(Q.round(3).to_string())
for nm in ("full", "resid"):
    z = P4[f"lr_z_{nm}"].dropna(); zc = z.sum() / np.sqrt(len(z)); print(nm, "LR combined z", round(zc, 2), "p", 2 * norm.sf(abs(zc)), "n cohorts", len(z), "cohorts p<.05:", int((P4.filter(like=f"lr_p_{nm}").bfill(axis=1).iloc[:, 0] < .05).sum()))
print("spearman model vs prolif: mean", P4.spearman_model_vs_prolif.mean(), P4.spearman_model_vs_prolif.min(), P4.spearman_model_vs_prolif.max(), "| resid:", P4.spearman_resid_vs_prolif.mean())
