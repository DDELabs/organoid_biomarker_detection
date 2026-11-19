"""Q5: does it hold in ER+/HER2-, HER2+ and triple-negative patients? Uses out-of-fold scores from d_proliferation.py (test cohort never in training).
Group from the er/her2 flags (TNBC = ER- and HER2-; HER2+ = HER2+ with any ER). Cohorts with no receptor data (GSE22093, GSE23988, TCGA-BRCA) are left out here.
A cohort-subgroup is used when it has >=5 responders and >=5 non-responders. BRIGHTNESS has only TNBC patients (flag ER-/HER2- for all)."""
import statsmodels.api as sm
from scipy.stats import norm
from common import *
S = pd.read_csv(TAB / "loco_scores.tsv.gz", sep="\t", index_col=0)
d, Z = breast_data(); assert (S.index == d.index).all()
def grp(r):
    if r.er == 1 and r.her2 == 0: return "ER+/HER2-"
    if r.her2 == 1: return "HER2+"
    if r.er == 0 and r.her2 == 0: return "TNBC"
    return "unknown"
S["grp"] = S.apply(grp, axis=1)
models = {"full": "All 673 pathways", "top100": "Top 100 pathways", "resid": "673, proliferation regressed out", "prolif": "Proliferation score only"}
rows = []
for gname, g in list(S.groupby("grp")) + [("All patients", S)]:
    if gname == "unknown": continue
    for c, x in g.groupby("cohort"):
        n1, n0 = int(x.y.sum()), int((1 - x.y).sum())
        if n1 < 5 or n0 < 5: continue
        for m in models:
            a = auc(x.y, x[m]); rows.append(dict(group=gname, cohort=c, model=m, n=len(x), responders=n1, auc=a, se=hm_se(a, n1, n0)))
R = pd.DataFrame(rows); R.to_csv(TAB / "subgroup_auc_per_cohort.tsv", sep="\t", index=False)
out = []
for (gname, m), g in R.groupby(["group", "model"]):
    mean = g.auc.mean(); se = np.sqrt((g.se ** 2).sum()) / len(g)
    out.append(dict(group=gname, model=m, n_cohorts=len(g), n_patients=int(g.n.sum()), mean_auc=mean, ci_low=mean - 1.96 * se, ci_high=mean + 1.96 * se,
                    cohorts_auc_above_0_5=int((g.auc > .5).sum()), min_cohort_auc=g.auc.min(), max_cohort_auc=g.auc.max()))
O = pd.DataFrame(out); O.to_csv(TAB / "subgroup_auc_pooled.tsv", sep="\t", index=False); print(O.round(3).to_string())
# top-30 pathways: odds ratio per SD within each subgroup (meta over cohorts)
T = pd.read_csv(TAB / "top30_pathways.tsv", sep="\t", index_col=0)
y = d["responder"].astype(int); res = []
for gname in ["ER+/HER2-", "HER2+", "TNBC"]:
    for p in T.index:
        est = []
        for c, x in S[S.grp == gname].groupby("cohort"):
            if x.y.sum() < 5 or (1 - x.y).sum() < 5: continue
            v = Z.loc[x.index, "REACTOME_" + p].values
            try:
                r = sm.Logit(x.y.values.astype(float), sm.add_constant(v)).fit(disp=0); est.append((r.params[1], r.bse[1]))
            except Exception: pass
        b = np.array([e[0] for e in est]); s = np.array([e[1] for e in est]); wt = 1 / s ** 2; bm = (b * wt).sum() / wt.sum(); se = np.sqrt(1 / wt.sum())
        res.append(dict(group=gname, pathway=p, n_cohorts=len(b), odds_ratio_per_sd=np.exp(bm), ci_low=np.exp(bm - 1.96 * se), ci_high=np.exp(bm + 1.96 * se), p=2 * norm.sf(abs(bm / se))))
Rr = pd.DataFrame(res); Rr["overall_raw_or"] = Rr.pathway.map(T.raw_odds_ratio_per_sd)
Rr["same_direction_as_overall"] = np.sign(np.log(Rr.odds_ratio_per_sd)) == np.sign(np.log(Rr.overall_raw_or))
Rr.to_csv(TAB / "subgroup_top30_odds_ratios.tsv", sep="\t", index=False)
print(Rr.groupby("group").agg(n=("pathway", "size"), same_direction=("same_direction_as_overall", "sum"), p_lt_05=("p", lambda v: int((v < .05).sum()))))
print("top30 flip between weight sign and overall raw OR direction:", int((np.sign(T.L2_weight) != np.sign(np.log(T.raw_odds_ratio_per_sd))).sum()))
