"""Patient-trained, drug-specific pathway biomarkers from survival + treatment (+ response notes).

For one drug the model is a Cox regression stratified by cancer type:

    hazard = baseline_cancer(t) * exp( covariates + general_risk + treated
                                       + sum_k a_k * pathway_k                 (prognostic part)
                                       + sum_k d_k * treated * pathway_k )     (drug-specific part)

The pathways k are the drug's network-proximal shortlist (Kong et al. 2020 step 1). The
coefficients d_k are the drug's biomarker: `benefit = -sum_k d_k * pathway_k` is high when the
drug is expected to help more. Everything is cross-fitted by patient, so every patient is scored
by a model that never saw that patient.
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize


# --------------------------------------------------------------------------- penalised Cox
class StratCox:
    """Cox partial likelihood (Breslow ties), stratified, with a per-coefficient ridge penalty."""

    def __init__(self, time, event, strata):
        time, event = np.asarray(time, float), np.asarray(event, float)
        codes = pd.factorize(np.asarray(strata))[0]
        self.order = np.lexsort((-time, codes))                 # within stratum: time descending
        t, s = time[self.order], codes[self.order]
        self.e = event[self.order]
        n = len(t)
        new_s = np.r_[True, s[1:] != s[:-1]]
        self.start = np.flatnonzero(new_s)                      # first row of each stratum
        self.sid = np.cumsum(new_s) - 1
        grp_end = np.r_[(t[1:] != t[:-1]) | new_s[1:], True]    # last row of each tie group
        nxt = np.flatnonzero(grp_end)
        self.last = nxt[np.searchsorted(nxt, np.arange(n))]     # risk set of row i = rows start..last[i]
        self.n_events = float(self.e.sum())

    def _cum(self, v):                                          # cumulative sum restarted per stratum
        c = np.cumsum(v)
        base = np.r_[0.0, c][self.start]
        return c - base[self.sid]

    def _rcum(self, v):                                         # reverse cumulative sum per stratum
        c = np.cumsum(v[::-1])[::-1]
        end = np.r_[self.start[1:], len(v)]
        tail = np.r_[c, 0.0][end]
        return c - tail[self.sid]

    def loss_grad(self, b, X, pen):
        eta = X @ b
        m = np.maximum.reduceat(eta, self.start)[self.sid]       # per-stratum shift: no underflow in any risk set
        w = np.exp(eta - m)
        S0 = self._cum(w)[self.last]
        ll = np.sum(self.e * (eta - m - np.log(S0)))
        a = np.zeros(len(w))
        np.add.at(a, self.last, self.e / S0)                    # each event adds 1/S0 to its risk set
        c = self._rcum(a)
        g = X.T @ (self.e - w * c)
        f, gr = (-ll + 0.5 * np.sum(pen * b * b)) / self.n_events, (-g + pen * b) / self.n_events
        if not np.isfinite(f):                                   # wild line-search step: tell the optimiser to back off
            return 1e10, np.zeros_like(b)
        return f, gr

    def fit(self, X, pen, b0=None):
        X = np.asarray(X, float)[self.order]
        pen = np.broadcast_to(np.asarray(pen, float), (X.shape[1],))
        b0 = np.zeros(X.shape[1]) if b0 is None else b0
        with np.errstate(all="ignore"):
            r = minimize(self.loss_grad, b0, args=(X, pen), jac=True, method="L-BFGS-B",
                         options={"maxiter": 500, "ftol": 1e-10, "gtol": 1e-7})
        return r.x

    def loglik(self, b, X):
        return -self.loss_grad(b, np.asarray(X, float)[self.order], np.zeros(X.shape[1]))[0] * self.n_events


def cox_test(df, cols, time="time", event="event", strata="project"):
    """Unpenalised stratified Cox for a few columns; returns DataFrame(coef, se, z, p)."""
    from statsmodels.duration.hazard_regression import PHReg
    d = df[cols + [time, event, strata]].dropna()
    X = d[cols].astype(float)
    X = X.loc[:, X.std() > 1e-9]
    r = PHReg(d[time].values, X.values, status=d[event].values, strata=d[strata].values, ties="breslow").fit(disp=0)
    out = pd.DataFrame({"coef": r.params, "se": r.bse, "z": r.params / r.bse, "p": r.pvalues}, index=X.columns)
    out.attrs["n"], out.attrs["events"] = len(d), int(d[event].sum())
    return out


# --------------------------------------------------------------------------- helpers
def zscore_within(M, proj):
    g = M.groupby(proj)
    return ((M - g.transform("mean")) / g.transform("std").replace(0, 1)).fillna(0)


def folds_by_patient(index, k=5, seed=0):
    rng = np.random.default_rng(seed)
    return pd.Series(rng.permutation(len(index)) % k, index=index)


def shortlist(z, all_pathways, cutoff=-1.2816, min_k=15, fallback_k=30):
    """Network-proximal pathways for a drug. Returns (pathways, how)."""
    if z is None:
        return list(all_pathways), "no targets: all pathways"
    z = z.reindex(all_pathways).dropna()
    s = list(z.index[z <= cutoff])
    if len(s) >= min_k:
        return s, "proximal"
    return list(z.sort_values().index[:fallback_k]), f"fewer than {min_k} proximal: {fallback_k} closest"


def general_risk(Z, cdr, time, event, fold, lam=200.0):
    """Cross-fitted pan-cancer prognostic score from all pathways (out-of-fold)."""
    out = pd.Series(np.nan, index=Z.index)
    for f in sorted(fold.unique()):
        tr, te = fold.index[fold != f], fold.index[fold == f]
        m = StratCox(cdr.loc[tr, time], cdr.loc[tr, event], cdr.loc[tr, "project"])
        b = m.fit(Z.loc[tr].values, lam)
        out.loc[te] = Z.loc[te].values @ b
    return out


def fit_interaction(D, Zs, base_cols, lam_main, lam_int):
    """Fit the per-drug model. D: frame with time,event,project,treated + base_cols; Zs: shortlist scores."""
    T = D["treated"].values[:, None]
    X = np.hstack([D[base_cols].values, D[["treated"]].values, Zs.values, Zs.values * T])
    k, nb = Zs.shape[1], len(base_cols) + 1
    pen = np.r_[np.full(nb, 1e-3), np.full(k, lam_main), np.full(k, lam_int)]
    m = StratCox(D["time"], D["event"], D["project"])
    b = m.fit(X, pen)
    return pd.Series(b[nb:nb + k], Zs.columns), pd.Series(b[nb + k:], Zs.columns)


def response_model(Zs, y, lam):
    """Ridge logistic of response (1 = CR/PR) on shortlist pathways; returns coefficient Series."""
    from sklearn.linear_model import LogisticRegression
    if y.nunique() < 2 or y.value_counts().min() < 5:
        return None
    m = LogisticRegression(C=1.0 / lam, max_iter=2000).fit(Zs.values, y.values)
    return pd.Series(m.coef_[0], Zs.columns)


def stratified_auc(score, y, groups, min_each=3):
    """AUROC pooled over within-group comparisons only (never compares patients of different cancers)."""
    from scipy.stats import rankdata
    num = den = 0.0
    for g in pd.unique(groups):
        i = groups == g
        s, yy = score[i], y[i]
        n1, n0 = int(yy.sum()), int((1 - yy).sum())
        if n1 < min_each or n0 < min_each:
            continue
        r = rankdata(s)
        num += r[yy == 1].sum() - n1 * (n1 + 1) / 2
        den += n1 * n0
    return (num / den, int(den)) if den else (np.nan, 0)
