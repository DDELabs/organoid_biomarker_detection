"""Model zoo for the W3 'positives' benchmark (regression version of the breast zoo).

Every model learns a continuous target y (organoid drug resistance, or a patient 'badness' score) from pathway
scores X and returns a score for new samples (higher = predicted more resistant / worse outcome).
Hyper-parameters are chosen ONLY by cross-validation inside the training data given to `fit_model`.

Models: L2 (ridge), L1 (lasso), ElasticNet, SVM_linear, SVM_RBF, RandomForest, XGBoost, MLP, NetSmoothL2, GCN_prior,
plus the 'shortlist' variants (network-proximity pathways only) of the four linear models.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.svm import SVR, LinearSVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold, GroupKFold
import xgboost as xgb

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
SEED = 20261119
J_ALL = np.load(HERE / "tables" / "pathway_jaccard.npy")
J_NAMES = list(pd.read_csv(HERE / "tables" / "pathway_flags.tsv", sep="\t", index_col=0).index)
J_IDX = {n: i for i, n in enumerate(J_NAMES)}

G = lambda **k: [dict(zip(k, v)) for v in itertools.product(*k.values())]


def jaccard_sub(feats):
    """Jaccard overlap among the given pathways (zero for pathways missing from the 673-pathway table)."""
    ix = [J_IDX.get(f, -1) for f in feats]
    J = np.zeros((len(feats), len(feats)))
    ok = [i for i, k in enumerate(ix) if k >= 0]
    if ok:
        J[np.ix_(ok, ok)] = J_ALL[np.ix_([ix[i] for i in ok], [ix[i] for i in ok])]
    np.fill_diagonal(J, 0)
    return J


def knn_adj(J, k=10):
    A = np.zeros_like(J)
    for i in range(len(J)):
        nn_ = np.argsort(-J[i])[:k]; A[i, nn_] = J[i, nn_]
    return np.maximum(A, A.T)


def smooth_matrix(J, s):
    A = knn_adj(J); A = A / np.maximum(A.sum(1, keepdims=True), 1e-12)
    return (1 - s) * np.eye(len(J)) + s * A


# ------------------------------------------------------------------------------------------------ fitters: (X, y, hp, ctx) -> predict fn
def f_ridge(X, y, hp, ctx):
    m = Ridge(alpha=hp["alpha"]).fit(X, y); return m.predict
def f_lasso(X, y, hp, ctx):
    m = Lasso(alpha=hp["alpha"], max_iter=20000).fit(X, y); return m.predict
def f_enet(X, y, hp, ctx):
    m = ElasticNet(alpha=hp["alpha"], l1_ratio=hp["l1"], max_iter=20000).fit(X, y); return m.predict
def f_lsvr(X, y, hp, ctx):
    m = LinearSVR(C=hp["C"], epsilon=0.0, loss="epsilon_insensitive", dual=True, max_iter=20000, random_state=SEED).fit(X, y); return m.predict
def f_rbf(X, y, hp, ctx):
    m = SVR(kernel="rbf", C=hp["C"], gamma=hp["g"] / X.shape[1]).fit(X, y); return m.predict
def f_rf(X, y, hp, ctx):
    m = RandomForestRegressor(n_estimators=300, min_samples_leaf=hp["leaf"], max_features="sqrt", n_jobs=1, random_state=SEED).fit(X, y); return m.predict
def f_xgb(X, y, hp, ctx):
    m = xgb.XGBRegressor(n_estimators=200, learning_rate=0.05, max_depth=hp["depth"], subsample=0.8, colsample_bytree=0.3,
                         reg_lambda=hp["lam"], n_jobs=1, random_state=SEED, verbosity=0).fit(X, y); return m.predict
def f_mlp(X, y, hp, ctx):
    es = len(y) >= 60
    m = MLPRegressor(hidden_layer_sizes=(32,), alpha=hp["alpha"], max_iter=300, early_stopping=es, validation_fraction=0.15 if es else 0.1,
                     n_iter_no_change=15, random_state=SEED).fit(X, y); return m.predict
def f_smooth(X, y, hp, ctx):
    S = smooth_matrix(ctx["J"], hp["s"]).T
    m = Ridge(alpha=hp["alpha"]).fit(X @ S, y); return lambda Xt: m.predict(Xt @ S)
_gcn = {}
def f_gcn(X, y, hp, ctx):
    import torch, scipy.sparse as sp
    sys.path.insert(0, str(ROOT / "nar/07_B3_models/code"))
    import nets
    torch.set_num_threads(1)
    key = id(ctx)
    if key not in _gcn:
        A = (knn_adj(ctx["J"]) > 0).astype(float) + np.eye(len(ctx["J"]))
        A = (A > 0).astype(float)
        dg = A.sum(1); S = sp.coo_matrix(A / np.sqrt(np.outer(dg, dg)))
        At = torch.sparse_coo_tensor(np.vstack([S.row, S.col]), S.data, S.shape, dtype=torch.float32).coalesce()
        pr = np.asarray(ctx["prior"], float); pr = (pr - pr.mean()) / (pr.std() or 1.0)
        _gcn.clear(); _gcn[key] = (At, torch.tensor(pr, dtype=torch.float32))
    At, prior = _gcn[key]
    m = nets.gcn_fit(At, prior, X, np.asarray(y, float), steps=hp["steps"], wd=hp["wd"], seed=SEED)
    return lambda Xt: nets.gcn_predict(m, Xt)


def zoo(small):
    """name -> (feature set 'all'|'short', grid, fitter).  `small` = few training samples (< 60): smaller leaves."""
    leaf = [1, 3] if small else [5, 25]
    A = dict(alpha=[1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0])
    return {
        "L2": ("all", G(**A), f_ridge),
        "L1": ("all", G(alpha=[0.003, 0.01, 0.03, 0.1, 0.3]), f_lasso),
        "ElasticNet": ("all", G(alpha=[0.01, 0.03, 0.1], l1=[0.2, 0.5, 0.8]), f_enet),
        "SVM_linear": ("all", G(C=[1e-5, 1e-4, 1e-3, 1e-2]), f_lsvr),
        "SVM_RBF": ("all", G(C=[0.1, 1.0, 10.0], g=[0.3, 1.0]), f_rbf),
        "RandomForest": ("all", G(leaf=leaf), f_rf),
        "XGBoost": ("all", G(depth=[2, 4], lam=[1, 30]), f_xgb),
        "MLP": ("all", G(alpha=[0.1, 1.0, 10.0]), f_mlp),
        "NetSmoothL2": ("all", G(s=[0.3, 0.6], alpha=[10.0, 100.0, 1000.0, 10000.0]), f_smooth),
        "GCN_prior": ("all", G(wd=[1e-4, 1e-2], steps=[150]), f_gcn),
        "L2_shortlist": ("short", G(**A), f_ridge),
        "L1_shortlist": ("short", G(alpha=[0.003, 0.01, 0.03, 0.1, 0.3]), f_lasso),
        "ElasticNet_shortlist": ("short", G(alpha=[0.01, 0.03, 0.1], l1=[0.2, 0.5, 0.8]), f_enet),
        "SVM_linear_shortlist": ("short", G(C=[1e-5, 1e-4, 1e-3, 1e-2]), f_lsvr),
    }


def _spear(a, b):
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return -1.0 if np.std(a) < 1e-12 else 0.0
    r = spearmanr(a, b)[0]
    return -1.0 if not np.isfinite(r) else float(r)


def inner_score(fit, X, y, hp, ctx, folds, groups=None):
    """Cross-validated Spearman between out-of-fold predictions and y, inside the training data only.
    Leave-one-out / K-fold on samples: Spearman of the pooled out-of-fold predictions.
    Grouped folds (whole cohorts held out): mean Spearman within each held-out cohort."""
    oof = np.full(len(y), np.nan); per = []
    for a, b in folds:
        f = fit(X[a], y[a], hp, ctx)
        p = np.asarray(f(X[b])).ravel(); oof[b] = p
        if groups is not None:
            for g in np.unique(groups[b]):
                m = groups[b] == g
                if m.sum() >= 8: per.append(_spear(p[m], y[b][m]))
    return float(np.mean(per)) if groups is not None and per else _spear(oof, y)


def make_folds(n, groups=None, seed=0):
    if groups is not None:
        k = min(4, len(np.unique(groups)))
        return list(GroupKFold(k).split(np.zeros(n), groups=groups))
    if n <= 20:
        return [(np.setdiff1d(np.arange(n), [i]), np.array([i])) for i in range(n)]
    return list(KFold(5, shuffle=True, random_state=seed).split(np.arange(n)))


def fit_model(name, Xtr, ytr, feats, short_idx, groups=None, prior=None, folds=None):
    """Select hyper-parameters by inner CV in the training data, refit on all of it.
    Returns (predict_fn_on_full_feature_matrix, chosen hp, inner score)."""
    small = len(ytr) < 60
    fs, grid, fit = zoo(small)[name]
    cols = np.arange(Xtr.shape[1]) if fs == "all" else np.asarray(short_idx)
    X = Xtr[:, cols]
    J = jaccard_sub([feats[i] for i in cols])
    ctx = dict(J=J, prior=(np.zeros(len(cols)) if prior is None else np.asarray(prior)[cols]))
    if folds is None:
        folds = make_folds(len(ytr), groups)
    if len(grid) > 1:
        sc = [inner_score(fit, X, ytr, hp, ctx, folds, groups) for hp in grid]
        hp = grid[int(np.argmax(sc))]; best = float(np.max(sc))
    else:
        hp, best = grid[0], np.nan
    f = fit(X, ytr, hp, ctx)
    return (lambda Xt: np.asarray(f(Xt[:, cols])).ravel()), hp, best


MODEL_NAMES = list(zoo(False))
