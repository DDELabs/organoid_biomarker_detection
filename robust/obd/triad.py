"""TRIAD: patient-anchored, prior-guided, pan-cancer drug-response biomarker discovery.

Why: every organoid-only analysis failed on the patient side (few treated deaths per
cancer, non-randomised overall survival). TRIAD makes the PATIENT the unit of learning
and uses the organoid / cell-line / network evidence as priors:

  label     RECIST best response (CR/PR vs SD/PD) recorded per drug in TCGA, pooled
            across all cancer types treated with the drug (thousands of labels)
  features  per-sample Reactome rank scores (platform- and cohort-independent),
            centred within cancer type, so cancer identity cannot drive predictions
  priors    (i) network proximity of each pathway to the drug targets (per-feature
                penalty: closer pathways are shrunk less)
            (ii) a pre-clinical coefficient vector (organoid NIT and/or GDSC cell lines)
                 that the patient model is shrunk TOWARD instead of toward zero
  model     prior-centred ridge logistic regression:
              min  -loglik(beta) + lambda * sum_j (beta_j - m_j)^2 / w_j^2
  validate  leave-one-cancer-out (LOCO) AUROC, a within-cancer label-permutation null,
            and frozen application to external trial cohorts
Each prior is ablated (none / network / pre-clinical / both), so the value of the
organoid evidence is measured directly on held-out patients.
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.metrics import roc_auc_score


def center_within(scores, groups):
    """Centre and scale features within each group (cancer type). scores: samples x features."""
    out = scores.copy()
    for g, idx in scores.groupby(groups).groups.items():
        block = scores.loc[idx]
        sd = block.std().replace(0, 1).fillna(1)
        out.loc[idx] = (block - block.mean()) / sd
    return out.fillna(0.0)


def fit_prior_logistic(X, y, lam=1.0, w=None, m=None, sample_weight=None):
    """Ridge logistic regression centred on prior mean m with per-feature scale w."""
    n, p = X.shape
    w = np.ones(p) if w is None else np.asarray(w, float)
    m = np.zeros(p) if m is None else np.asarray(m, float)
    sw = np.ones(n) if sample_weight is None else np.asarray(sample_weight, float)
    inv = 1.0 / np.maximum(w, 1e-3) ** 2

    def f(theta):
        b0, b = theta[0], theta[1:]
        z = b0 + X @ b
        ll = sw * (y * z - np.logaddexp(0, z))
        pen = lam * np.sum(inv * (b - m) ** 2)
        pr = 1 / (1 + np.exp(-z))
        g0 = -np.sum(sw * (y - pr))
        gb = -X.T @ (sw * (y - pr)) + 2 * lam * inv * (b - m)
        return -ll.sum() + pen, np.r_[g0, gb]

    theta0 = np.r_[np.log((y.mean() + 1e-3) / (1 - y.mean() + 1e-3)), m]
    res = minimize(f, theta0, jac=True, method="L-BFGS-B")
    return res.x[0], res.x[1:]


def loco(X, y, groups, lam=1.0, w=None, m=None, min_test=10):
    """Leave-one-cancer-out predictions. Returns per-sample scores and per-cancer AUROC."""
    X, y, groups = np.asarray(X, float), np.asarray(y, float), np.asarray(groups)
    pred = np.full(len(y), np.nan)
    per = {}
    for g in np.unique(groups):
        te = groups == g
        if te.sum() < min_test or len(np.unique(y[te])) < 2 or len(np.unique(y[~te])) < 2:
            continue
        # balance cancers in training so large cohorts do not dominate
        counts = pd.Series(groups[~te]).value_counts()
        sw = 1.0 / pd.Series(groups[~te]).map(counts).values
        sw = sw * len(sw) / sw.sum()
        b0, b = fit_prior_logistic(X[~te], y[~te], lam, w, m, sw)
        pred[te] = b0 + X[te] @ b
        per[g] = {"n": int(te.sum()), "responders": int(y[te].sum()), "auc": roc_auc_score(y[te], pred[te])}
    ok = ~np.isnan(pred)
    pooled = roc_auc_score(y[ok], pred[ok]) if ok.sum() and len(np.unique(y[ok])) == 2 else np.nan
    # mean within-cancer AUC is the honest summary (pooled AUC can reflect base-rate shifts)
    mean_auc = float(np.mean([v["auc"] for v in per.values()])) if per else np.nan
    return pred, per, pooled, mean_auc


def prior_only_auc(X, y, groups, m):
    """No patient fitting at all: score = X @ prior (e.g. frozen organoid signature)."""
    s = np.asarray(X, float) @ np.asarray(m, float)
    per = {}
    for g in np.unique(groups):
        te = np.asarray(groups) == g
        if te.sum() >= 10 and len(np.unique(np.asarray(y)[te])) == 2:
            per[g] = roc_auc_score(np.asarray(y)[te], s[te])
    return float(np.mean(list(per.values()))) if per else np.nan, per


def permutation_null(X, y, groups, lam, w, m, n=200, seed=0):
    """Mean within-cancer LOCO AUROC with labels shuffled within each cancer."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    groups = np.asarray(groups)
    vals = []
    for _ in range(n):
        yp = y.copy()
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            yp[idx] = rng.permutation(yp[idx])
        vals.append(loco(X, yp, groups, lam, w, m)[3])
    return np.array(vals)
