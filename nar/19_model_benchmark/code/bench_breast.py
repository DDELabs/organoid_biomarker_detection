"""W3 breast benchmark. Usage: python bench_breast.py MODEL [MODEL ...]   (one process each; max 2 at a time)
Every model is trained on the 12 other cohorts and tested on the held-out one. Hyper-parameters are chosen by
cohort-grouped 4-fold cross-validation INSIDE the 12 training cohorts only (score = mean within-cohort AUC).
Writes tables/preds/<model>.tsv.gz (one row per test patient) and tables/preds/<model>_hp.tsv (chosen settings)."""
import sys, time, itertools; sys.path.insert(0, __file__.rsplit("/", 1)[0])
import os
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
from common import *
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC, SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
import xgboost as xgb
PRED = TAB / "preds"; PRED.mkdir(exist_ok=True, parents=True)
d, Z, P = breast_data()
y = d["responder"].astype(int); coh = d["cohort"]
flags = pd.read_csv(TAB / "pathway_flags.tsv", sep="\t", index_col=0)
cols = list(Z.columns)
SL_COLS = [c for c in cols if flags.loc[c, "shortlist_paper"]]
PROLIF = [c for c in cols if flags.loc[c, "proliferation"]]
J = np.load(TAB / "pathway_jaccard.npy")


# ---------------------------------------------------------------- model factories: name -> (grid list of dict, fit(Xtr,ytr,hp)->score_fn)
def lr(**kw):
    return lambda X, y_, hp: LogisticRegression(max_iter=3000, **{**kw, **hp}).fit(X, y_).decision_function
def lsvm(X, y_, hp): return LinearSVC(C=hp["C"], max_iter=5000, dual=True).fit(X, y_).decision_function
def rf(X, y_, hp):
    m = RandomForestClassifier(n_estimators=300, min_samples_leaf=hp["leaf"], max_features="sqrt", n_jobs=1, random_state=SEED).fit(X, y_)
    return lambda Xt: m.predict_proba(Xt)[:, 1]
def xg(X, y_, hp):
    m = xgb.XGBClassifier(n_estimators=200, learning_rate=0.05, max_depth=hp["depth"], subsample=0.8, colsample_bytree=0.3,
                          reg_lambda=hp["lam"], n_jobs=1, random_state=SEED, verbosity=0).fit(X, y_)
    return lambda Xt: m.predict_proba(Xt)[:, 1]
def svr(X, y_, hp): return SVC(kernel="rbf", C=hp["C"], gamma=hp["g"] / X.shape[1], cache_size=500).fit(X, y_).decision_function
def mlp(X, y_, hp):
    m = MLPClassifier((32,), alpha=hp["alpha"], max_iter=300, early_stopping=True, validation_fraction=0.15, n_iter_no_change=15, random_state=SEED).fit(X, y_)
    return lambda Xt: m.predict_proba(Xt)[:, 1]

# network smoothing on a pathway-overlap graph: X_s = (1-s) X + s * A_norm X   (A = 10 nearest neighbours by gene overlap)
def _smooth_matrix(idx, s):
    Jj = J[np.ix_(idx, idx)].copy(); np.fill_diagonal(Jj, 0)
    A = np.zeros_like(Jj)
    for i in range(len(idx)):
        nn_ = np.argsort(-Jj[i])[:10]; A[i, nn_] = Jj[i, nn_]
    A = np.maximum(A, A.T); A = A / np.maximum(A.sum(1, keepdims=True), 1e-12)
    return (1 - s) * np.eye(len(idx)) + s * A
ALLIDX = np.arange(len(cols))
def smoothlr(X, y_, hp):
    S = _smooth_matrix(ALLIDX, hp["s"]).T
    m = LogisticRegression(C=hp["C"], max_iter=3000).fit(X @ S, y_)
    return lambda Xt: m.decision_function(Xt @ S)

# GCN from nets.py (graph = pathway overlap graph, node prior = paper-table network proximity)
_gcn = {}
def gcn(X, y_, hp):
    sys.path.insert(0, str(ROOT / "nar/07_B3_models/code"))
    import torch, scipy.sparse as sp
    import nets
    torch.set_num_threads(1)
    if "g" not in _gcn:
        Jj = J.copy(); np.fill_diagonal(Jj, 0)
        A = np.zeros_like(Jj)
        for i in range(len(Jj)):
            nn_ = np.argsort(-Jj[i])[:10]; A[i, nn_] = 1
        A = np.maximum(A, A.T) + np.eye(len(Jj))
        dg = A.sum(1); S = sp.coo_matrix(A / np.sqrt(np.outer(dg, dg)))
        At = torch.sparse_coo_tensor(np.vstack([S.row, S.col]), S.data, S.shape, dtype=torch.float32).coalesce()
        pr = (-flags["min_z_paper"].fillna(0).values); pr = (pr - pr.mean()) / pr.std()   # closer to targets = larger prior
        _gcn["g"] = (At, torch.tensor(pr, dtype=torch.float32))
    At, prior = _gcn["g"]
    m = nets.gcn_fit(At, prior, X, y_.values.astype(float) if hasattr(y_, "values") else y_.astype(float), steps=hp["steps"], wd=hp["wd"], seed=SEED)
    return lambda Xt: nets.gcn_predict(m, Xt)

G = lambda **k: [dict(zip(k, v)) for v in itertools.product(*k.values())]
MODELS = {
 # name: (columns, grid, fit)
 "L2_fixed":        (cols, [dict(C=0.01)], lr()),                    # the existing model, lambda = 100, no tuning
 "L2":              (cols, G(C=[0.0003, 0.001, 0.003, 0.01, 0.03, 0.1]), lr()),
 "L1":              (cols, G(C=[0.003, 0.01, 0.03, 0.1, 0.3, 1.0]), lr(penalty="l1", solver="liblinear")),
 "ElasticNet":      (cols, G(C=[0.01, 0.03, 0.1], l1_ratio=[0.2, 0.5, 0.8]), lr(penalty="elasticnet", solver="saga", tol=1e-3)),
 "SVM_linear":      (cols, G(C=[1e-4, 1e-3, 1e-2]), lsvm),
 "SVM_RBF":         (cols, G(C=[0.1, 1, 10], g=[0.3, 1.0]), svr),
 "RandomForest":    (cols, G(leaf=[5, 25]), rf),
 "XGBoost":         (cols, G(depth=[2, 4], lam=[1, 30]), xg),
 "MLP":             (cols, G(alpha=[0.1, 1, 10]), mlp),
 "NetSmoothL2":     (cols, G(s=[0.3, 0.6], C=[0.003, 0.01, 0.03]), smoothlr),
 "GCN_prior":       (cols, G(wd=[1e-4, 1e-2], steps=[150]), gcn),
 # network-proximity shortlist (183 pathways near paclitaxel/docetaxel/doxorubicin/5-FU targets, paper table)
 "L2_shortlist":         (SL_COLS, G(C=[0.0003, 0.001, 0.003, 0.01, 0.03, 0.1]), lr()),
 "L1_shortlist":         (SL_COLS, G(C=[0.003, 0.01, 0.03, 0.1, 0.3, 1.0]), lr(penalty="l1", solver="liblinear")),
 "ElasticNet_shortlist": (SL_COLS, G(C=[0.01, 0.03, 0.1], l1_ratio=[0.2, 0.5, 0.8]), lr(penalty="elasticnet", solver="saga", tol=1e-3)),
 "SVM_linear_shortlist": (SL_COLS, G(C=[1e-4, 1e-3, 1e-2]), lsvm),
 "L2_fixed_shortlist":   (SL_COLS, [dict(C=0.01)], lr()),
}


def inner_score(fitfn, Xtr, ytr, gtr, hp):
    """mean within-cohort AUC over cohort-grouped 4 folds, inside the training cohorts"""
    aucs = []
    for a, b in GroupKFold(4).split(Xtr, ytr, gtr):
        f = fitfn(Xtr[a], ytr.iloc[a], hp)
        s = np.asarray(f(Xtr[b])).ravel() if True else None
        for c in np.unique(gtr[b]):
            m = gtr[b] == c
            yy = ytr.iloc[b].values[m]
            if len(set(yy)) == 2: aucs.append(roc_auc_score(yy, s[m]))
    return float(np.mean(aucs))


def run(name):
    ccols, grid, fitfn = MODELS[name]
    X = Z[ccols].values
    out, hps = [], []
    t0 = time.time()
    for c in sorted(coh.unique()):
        te = (coh == c).values; tr = ~te
        Xtr, ytr, gtr = X[tr], y[tr], coh[tr].values
        if len(grid) > 1:
            sc = [inner_score(fitfn, Xtr, ytr, gtr, hp) for hp in grid]
            hp = grid[int(np.argmax(sc))]
        else:
            hp, sc = grid[0], [np.nan]
        f = fitfn(Xtr, ytr, hp)
        s = np.asarray(f(X[te])).ravel()
        out.append(pd.DataFrame(dict(cohort=c, score=s, y=y[te].values), index=d.index[te]))
        hps.append(dict(test_cohort=c, **hp, inner_cv_auc=float(np.nanmax(sc))))
        print(name, c, hp, f"{roc_auc_score(y[te], s):.3f}", f"{time.time()-t0:.0f}s", flush=True)
    pd.concat(out).to_csv(PRED / f"{name}.tsv.gz", sep="\t")
    pd.DataFrame(hps).to_csv(PRED / f"{name}_hp.tsv", sep="\t", index=False)


if __name__ == "__main__":
    for n in (sys.argv[1:] or list(MODELS)):
        run(n)
