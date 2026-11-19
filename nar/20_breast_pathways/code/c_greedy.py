"""Q3: nested greedy forward selection. Inside each outer fold (12 training cohorts only):
 candidate pool = top 60 pathways by training-only univariate strength; at each step add the pathway that gives the best
 inner AUC (4 inner folds that keep whole cohorts together; mean of per-cohort AUC). Up to 15 steps. The outer test cohort is
 scored with the model refit on the chosen set. Also an honest stopping rule: smallest step within 0.01 of the best inner AUC."""
from joblib import Parallel, delayed
from common import *
import warnings; warnings.filterwarnings("ignore")
d, Z = breast_data(); y = d["responder"].astype(int).values; coh = d["cohort"].values; cols = list(Z.columns); X = Z.values
cohorts = sorted(set(coh)); POOL, STEPS = 60, 15

def inner_auc(Xs, ys, cs, folds):
    a = []
    for f in folds:
        te = np.isin(cs, f)
        s = TC.fit(Xs[~te], ys[~te]).decision_function(Xs[te])
        for c in f:
            m = cs[te] == c
            a.append(auc(ys[te][m], s[m]))
    return np.nanmean(a)

def outer(c):
    tr, te = coh != c, coh == c
    Xt, yt, ct = X[tr], y[tr], coh[tr]
    u = uni_score(Xt, yt, ct); pool = list(np.argsort(-np.abs(u))[:POOL])
    tc = sorted(set(ct)); r = np.random.default_rng(SEED); r.shuffle(tc)
    folds = [tc[i::4] for i in range(4)]
    sel, path = [], []
    for step in range(STEPS):
        best, ba = None, -1
        for j in pool:
            if j in sel: continue
            a = inner_auc(Xt[:, sel + [j]], yt, ct, folds)
            if a > ba: best, ba = j, a
        sel.append(best)
        s = TC.fit(Xt[:, sel], yt).decision_function(X[te][:, sel])
        path.append(dict(test_cohort=c, step=step + 1, added=cols[best], inner_auc=ba, test_auc=auc(y[te], s)))
    return path

res = Parallel(n_jobs=2)(delayed(outer)(c) for c in cohorts)
R = pd.DataFrame([r for p in res for r in p]); R.to_csv(TAB / "greedy_path_per_fold.tsv", sep="\t", index=False)
print("done")
