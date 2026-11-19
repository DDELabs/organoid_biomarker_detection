"""Q4: proliferation. All steps inside LOO folds (12 training cohorts). Saves per-patient out-of-fold scores for later use."""
import statsmodels.api as sm
from scipy.stats import spearmanr, norm, chi2
from common import *
d, Z = breast_data(); y = d["responder"].astype(int).values; coh = d["cohort"].values; cols = list(Z.columns); X = Z.values
cohorts = sorted(set(coh)); pi = np.array([c in set(prolif_cols(cols)) for c in cols]); print("proliferation pathways:", pi.sum())
pd.DataFrame({"pathway": [c for c, f in zip(cols, pi) if f]}).to_csv(TAB / "proliferation_pathways.tsv", sep="\t", index=False)
KS = [20, 50, 100, 200]
S = pd.DataFrame(index=d.index, columns=["full", "prolif", "resid", "noprolif"] + [f"top{k}" for k in KS], dtype=float)
def zs(v, cc):   # standardise within cohort
    v = pd.Series(v); return v.groupby(cc).transform(lambda x: (x - x.mean()) / (x.std() + 1e-12)).values
for c in cohorts:
    tr, te = coh != c, coh == c
    pr = zs(X[:, pi].mean(1), coh)                                  # proliferation score (mean of the pathways, z within cohort)
    sign = np.sign(np.mean([np.mean(uni_score(pr[tr][:, None], y[tr], coh[tr]))]))   # direction from training only
    S.loc[d.index[te], "prolif"] = sign * pr[te]
    w = TC.fit(X[tr], y[tr]); S.loc[d.index[te], "full"] = w.decision_function(X[te])
    for k in KS:
        idx = np.argsort(-np.abs(w.coef_[0]))[:k]
        S.loc[d.index[te], f"top{k}"] = TC.fit(X[tr][:, idx], y[tr]).decision_function(X[te][:, idx])
    b = (X[tr] * pr[tr][:, None]).sum(0) / (pr[tr] ** 2).sum()      # slope of each pathway on proliferation (training only)
    Xr = X - np.outer(pr, b)
    S.loc[d.index[te], "resid"] = TC.fit(Xr[tr], y[tr]).decision_function(Xr[te])
    S.loc[d.index[te], "noprolif"] = TC.fit(X[tr][:, ~pi], y[tr]).decision_function(X[te][:, ~pi])
    print(c, flush=True)
S.assign(cohort=coh, y=y, er=d["er"].values, her2=d["her2"].values).to_csv(TAB / "loco_scores.tsv.gz", sep="\t")
rows = []
for c in cohorts:
    m = coh == c; r = dict(test_cohort=c, n=int(m.sum()), responders=int(y[m].sum()))
    for k in ["full", "prolif", "resid", "noprolif"] + [f"top{k}" for k in KS]:
        a = auc(y[m], S.loc[d.index[m], k]); r[f"auc_{k}"] = a
    r["spearman_model_vs_prolif"] = spearmanr(S.loc[d.index[m], "full"], S.loc[d.index[m], "prolif"])[0]
    r["spearman_resid_vs_prolif"] = spearmanr(S.loc[d.index[m], "resid"], S.loc[d.index[m], "prolif"])[0]
    # likelihood-ratio tests: does the score add to proliferation (+ ER + HER2 where recorded)?
    for nm in ("full", "resid"):
        D0 = pd.DataFrame({"p": S.loc[d.index[m], "prolif"].values}); has = []
        for k in ("er", "her2"):
            v = d.loc[d.index[m], k]
            if v.notna().mean() > 0.9 and v.nunique() > 1:
                D0[k] = v.fillna(v.mode().iat[0]).values; has.append(k)
        D1 = D0.assign(s=S.loc[d.index[m], nm].values)
        try:
            f0 = sm.Logit(y[m], sm.add_constant(D0.astype(float))).fit(disp=0); f1 = sm.Logit(y[m], sm.add_constant(D1.astype(float))).fit(disp=0)
            lr = 2 * (f1.llf - f0.llf); r[f"lr_p_{nm}_adds_to_prolif" + ("_er_her2" if has else "_only")] = chi2.sf(lr, 1)
            r[f"lr_z_{nm}"] = np.sign(f1.params["s"]) * np.sqrt(max(lr, 0)); r["adjusted_for"] = "+".join(["prolif"] + has)
        except Exception as e:
            r[f"lr_err_{nm}"] = str(e)[:40]
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(TAB / "proliferation_per_cohort.tsv", sep="\t", index=False)
print(R.round(3).to_string())
