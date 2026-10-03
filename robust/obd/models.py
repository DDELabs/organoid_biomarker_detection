"""Organoid response models: the original baseline and the robust signature learner."""
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet, LinearRegression, Ridge, RidgeCV
from sklearn.svm import SVR

warnings.filterwarnings("ignore", category=ConvergenceWarning)


def zscore(x):
    """Column-wise z-score of a samples x features array (population SD, as StandardScaler)."""
    x = np.asarray(x, float)
    sd = x.std(axis=0)
    sd[sd == 0] = 1.0
    return (x - x.mean(axis=0)) / sd


# ----------------------------------------------------------------- baseline (paper)
BASELINE_MODELS = {
    "Ridge": lambda: RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)),
    "SVR": lambda: SVR(kernel="linear"),
    "LinearRegression": lambda: LinearRegression(),
}


def baseline_coefficients(X, y, model):
    """Kong et al.: fit on z-scored pathway scores, return coefficients."""
    return np.ravel(BASELINE_MODELS[model]().fit(zscore(X), y).coef_)


# ----------------------------------------------------------------- robust learner
def _ensemble_coefs(X, y, rng=None):
    """Coefficients of four penalised/robust linear learners on standardised data."""
    Xz = zscore(X)
    yz = (y - y.mean()) / (y.std() or 1.0)
    return np.vstack([
        Ridge(alpha=1.0).fit(Xz, yz).coef_,
        np.ravel(SVR(kernel="linear", C=1.0).fit(Xz, yz).coef_),
        ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=5000).fit(Xz, yz).coef_,
        np.array([spearmanr(Xz[:, j], yz)[0] if Xz[:, j].std() > 0 else 0.0 for j in range(Xz.shape[1])]),
    ])


def stability_selection(X, y, top_k=7, n_boot=200, frac=0.8, seed=0):
    """Selection frequency of each feature among the top_k (by |coef|) of every learner,
    over random subsamples of organoids (Meinshausen & Buhlmann 2010 style).

    Returns (frequency, sign agreement, mean coefficient) arrays over features.
    """
    rng = np.random.default_rng(seed)
    n, p = X.shape
    m = max(4, int(round(frac * n)))
    freq = np.zeros(p)
    signs = np.zeros(p)
    coef_sum = np.zeros(p)
    total = 0
    for _ in range(n_boot):
        idx = rng.choice(n, m, replace=False)
        C = _ensemble_coefs(X[idx], y[idx])
        for c in C:
            top = np.argsort(-np.abs(c))[:top_k]
            freq[top] += 1
            signs[top] += np.sign(c[top])
            coef_sum += c
            total += 1
    sel = np.maximum(freq, 1)
    return freq / total, np.abs(signs) / sel, coef_sum / total


def robust_signature(X, y, features, top_k=7, n_boot=200, min_freq=0.5, min_sign=0.9, seed=0):
    """Pathways selected stably and with a consistent sign; weights = mean ensemble coef."""
    freq, sign_agree, coef = stability_selection(X, y, top_k, n_boot, seed=seed)
    table = pd.DataFrame({"selection_freq": freq, "sign_agreement": sign_agree, "mean_coef": coef}, index=features)
    chosen = table[(table.selection_freq >= min_freq) & (table.sign_agreement >= min_sign)]
    if chosen.empty:  # fall back to the most stable pathways
        chosen = table.sort_values("selection_freq", ascending=False).head(max(2, top_k // 2))
    weights = chosen["mean_coef"]
    return weights, table


# ----------------------------------------------------------------- domain alignment
def precise_directions(Xs, Xt, n_components=10, n_shared=5):
    """PRECISE (Mourragui 2019): principal vectors shared by source (organoid) and
    target (patient) pathway-score spaces. Returns a features x n_shared basis."""
    def pcs(X):
        Xz = zscore(X)
        _, _, vt = np.linalg.svd(Xz, full_matrices=False)
        return vt[:n_components].T

    Ps, Pt = pcs(Xs), pcs(Xt)
    u, s, vt = np.linalg.svd(Ps.T @ Pt)
    shared = (Ps @ u + Pt @ vt.T) / 2.0
    shared /= np.linalg.norm(shared, axis=0)
    return shared[:, :n_shared], s[:n_shared]


def precise_weights(Xs, y, basis, alpha=1.0):
    """Ridge on organoid data projected onto the shared basis, back-projected to pathways."""
    Z = zscore(Xs) @ basis
    w = Ridge(alpha=alpha).fit(Z, (y - y.mean()) / (y.std() or 1.0)).coef_
    return basis @ w


# ----------------------------------------------------------------- organoid CV
def loocv(X, y, fit_predict):
    """Leave-one-organoid-out predictions; fit_predict(Xtrain, ytrain, Xtest) -> yhat.
    Any feature selection must happen inside fit_predict (avoids selection bias)."""
    pred = np.empty(len(y))
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        pred[i] = fit_predict(X[tr], y[tr], X[~tr])[0]
    rho, p = spearmanr(pred, y)
    return pred, rho, p


def standardise_with(train, test):
    mu, sd = train.mean(axis=0), train.std(axis=0)
    sd[sd == 0] = 1.0
    return (test - mu) / sd
