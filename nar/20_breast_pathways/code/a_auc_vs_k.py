"""Q1: AUC vs number of pathways kept. Ranking and refit use TRAINING cohorts only (12 of 13); test cohort untouched."""
from common import *
d, Z = breast_data(); y = d["responder"].astype(int).values; coh = d["cohort"].values; cols = list(Z.columns); X = Z.values
cohorts = sorted(set(coh)); KS = [1, 2, 3, 5, 10, 20, 50, 100, 200, 673]
rng = np.random.default_rng(SEED); rows = []; chosen = {}
for c in cohorts:
    tr, te = coh != c, coh == c
    w = TC.fit(X[tr], y[tr]).coef_[0]
    u = uni_score(X[tr], y[tr], coh[tr])
    rk = {"weight": np.argsort(-np.abs(w)), "univariate": np.argsort(-np.abs(u))}
    for k in KS:
        for name, order in rk.items():
            idx = order[:k]
            m = TC.fit(X[tr][:, idx], y[tr]); s = m.decision_function(X[te][:, idx])
            rows.append(dict(test_cohort=c, ranking=name, k=k, n=int(te.sum()), auc=auc(y[te], s)))
            if k == 1 and name == "weight": chosen[c] = cols[idx[0]]
            if k == 1 and name == "univariate": chosen[c + "|uni"] = cols[idx[0]]
        if k < 673:   # random sets of the same size (20 draws) as a reference
            a = []
            for _ in range(20):
                idx = rng.choice(len(cols), k, replace=False)
                a.append(auc(y[te], TC.fit(X[tr][:, idx], y[tr]).decision_function(X[te][:, idx])))
            rows.append(dict(test_cohort=c, ranking="random", k=k, n=int(te.sum()), auc=float(np.mean(a))))
    print(c, flush=True)
R = pd.DataFrame(rows); R.to_csv(TAB / "auc_vs_k_per_cohort.tsv", sep="\t", index=False)
pd.Series(chosen).rename("best_single_pathway").to_csv(TAB / "best_single_pathway_per_fold.tsv", sep="\t")
