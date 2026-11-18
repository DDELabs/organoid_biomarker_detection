import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import *
from scipy.stats import norm, wilcoxon, binomtest
PRED = TAB / "preds"
REF = "L2_fixed"
d, Z, P = breast_data()

def midrank(x):
    o = np.argsort(x); z = x[o]; n = len(x); T = np.zeros(n); i = 0
    while i < n:
        j = i
        while j < n and z[j] == z[i]: j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1; i = j
    out = np.empty(n); out[o] = T; return out

def delong(y, S):
    """y: 0/1; S: (k, n) scores. Returns aucs (k,), covariance (k,k)  [DeLong 1988, fast version Sun & Xu 2014]."""
    pos, neg = S[:, y == 1], S[:, y == 0]; m, n = pos.shape[1], neg.shape[1]; k = S.shape[0]
    tx = np.array([midrank(r) for r in pos]); ty = np.array([midrank(r) for r in neg]); tz = np.array([midrank(r) for r in S])
    aucs = tz[:, :0].sum(1) if False else (tz[:, y == 1].sum(1) / m / n) - (m + 1.0) / 2.0 / n
    v01 = (tz[:, y == 1] - tx) / n; v10 = 1.0 - (tz[:, y == 0] - ty) / m
    return aucs, np.cov(v01) / m + np.cov(v10) / n if k > 1 else np.atleast_2d(np.cov(v01) / m + np.cov(v10) / n)

names = sorted(f.name[:-7] for f in PRED.glob("*.tsv.gz"))
preds = {n: pd.read_csv(PRED / f"{n}.tsv.gz", sep="\t", index_col=0) for n in names}
er = d[["er", "her2"]]
cohorts = sorted(d["cohort"].unique()); K = len(cohorts)
rows = []
for n in names:
    for c in cohorts:
        pc = preds[n][preds[n].cohort == c]
        if len(pc) == 0: continue
        r = preds[REF][preds[REF].cohort == c].loc[pc.index]       # reference on exactly the same patients
        yy = pc.y.values
        a, cv = delong(yy, np.vstack([pc.score.values, r.score.values]))
        se = np.sqrt(cv[0, 0]); dd = a[0] - a[1]; sed = np.sqrt(max(cv[0, 0] + cv[1, 1] - 2 * cv[0, 1], 0))
        # AUC inside ER/HER2 groups (patients with both known; pooled within-group pairs)
        k = pc.index.intersection(er.dropna().index)
        wg = np.nan
        if len(k) > 20:
            g = (er.loc[k, "er"].astype(int).astype(str) + er.loc[k, "her2"].astype(int).astype(str)).values
            wg = SL.stratified_auc(pc.loc[k, "score"].values, pc.loc[k, "y"].values, g, 5)[0]
        rows.append(dict(model=n, test_cohort=c, n=len(pc), responders=int(yy.sum()), auc=a[0], se=se,
                         ci_low=max(a[0] - 1.96 * se, 0), ci_high=min(a[0] + 1.96 * se, 1), p_vs_chance=2 * norm.sf(abs(a[0] - .5) / se),
                         auc_within_er_her2=wg, ref_auc_same_patients=a[1], diff_vs_L2=dd, diff_se=sed,
                         diff_ci_low=dd - 1.96 * sed, diff_ci_high=dd + 1.96 * sed,
                         p_paired_vs_L2=(2 * norm.sf(abs(dd) / sed) if sed > 0 else np.nan)))
R = pd.DataFrame(rows)
R.to_csv(TAB / "breast_loco_all_models.tsv", sep="\t", index=False)

# --- model comparison: mean across cohorts
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.minimum(q, 1)
cmp_ = []
for n, g in R.groupby("model"):
    nc = len(g); mean = g.auc.mean(); se = np.sqrt((g.se ** 2).sum()) / nc
    md = g.diff_vs_L2.mean(); sed = np.sqrt((g.diff_se ** 2).sum()) / nc
    dif = g.diff_vs_L2.values
    try: wp = wilcoxon(dif).pvalue if n != REF and np.any(dif != 0) else np.nan
    except Exception: wp = np.nan
    cmp_.append(dict(model=n, cohorts=nc, mean_auc=mean, mean_ci_low=mean - 1.96 * se, mean_ci_high=mean + 1.96 * se,
                     mean_diff_vs_L2=md, diff_ci_low=md - 1.96 * sed, diff_ci_high=md + 1.96 * sed,
                     p_mean_diff=(2 * norm.sf(abs(md) / sed) if sed > 0 else np.nan),
                     cohorts_better_than_L2=int((dif > 0).sum()), cohorts_worse_than_L2=int((dif < 0).sum()),
                     cohorts_sig_better_p05=int(((dif > 0) & (g.p_paired_vs_L2 < .05)).sum()),
                     cohorts_sig_worse_p05=int(((dif < 0) & (g.p_paired_vs_L2 < .05)).sum()),
                     p_wilcoxon_across_cohorts=wp, mean_auc_within_er_her2=g.auc_within_er_her2.mean(), cohorts_with_within=int(g.auc_within_er_her2.notna().sum())))
C = pd.DataFrame(cmp_).sort_values("mean_auc", ascending=False)
m = C.model != REF
C.loc[m, "q_mean_diff_BH"] = bh(C.loc[m, "p_mean_diff"].values)
C.to_csv(TAB / "model_comparison.tsv", sep="\t", index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(C.round(4).to_string(index=False))
