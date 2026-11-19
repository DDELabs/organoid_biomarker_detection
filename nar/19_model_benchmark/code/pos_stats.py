"""Statistics for the W3 positives benchmark: metrics, paired bootstrap against the original method.

Score convention everywhere: higher score = predicted more resistant / worse outcome.
Survival tasks : C-index (higher = better ranking of shorter survival), HR per SD (Cox), log-rank of the median split.
Binary tasks   : AUC of the (oriented) score for the labelled class, one-sided Mann-Whitney p (as in the original test).
Paired bootstrap: the SAME resampled patients are used for every model, so differences are paired.
"""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, rankdata, norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import survdiff


def cparts(T, E, r):
    """Harrell C numerator / denominator (higher risk should die earlier)."""
    T, E, r = np.asarray(T, float), np.asarray(E, int), np.asarray(r, float)
    ev = np.where(E == 1)[0]
    if len(ev) == 0:
        return 0.0, 0.0
    later = T[None, :] > T[ev][:, None]
    num = ((r[ev][:, None] > r[None, :]) & later).sum() + 0.5 * ((r[ev][:, None] == r[None, :]) & later).sum()
    return float(num), float(later.sum())


def cindex(T, E, r):
    n, d = cparts(T, E, r)
    return n / d if d else np.nan


def cox(T, E, s, strata=None, cov=None):
    """Cox HR per SD of s (optionally stratified and adjusted). Returns dict(coef, se, hr, lo, hi, p)."""
    s = np.asarray(s, float)
    if np.std(s) < 1e-12 or np.sum(E) < 2:
        return dict(coef=np.nan, se=np.nan, hr=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
    z = (s - s.mean()) / s.std()
    X = z[:, None] if cov is None else np.column_stack([z, cov])
    try:
        r = PHReg(np.asarray(T, float), X, status=np.asarray(E, float), strata=None if strata is None else np.asarray(strata), ties="efron").fit(disp=False)
        b, se = float(r.params[0]), float(r.bse[0])
        return dict(coef=b, se=se, hr=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=float(2 * norm.sf(abs(b / se))))
    except Exception:
        return dict(coef=np.nan, se=np.nan, hr=np.nan, lo=np.nan, hi=np.nan, p=np.nan)


def split_stats(T, E, s):
    """Median split (above median = high predicted resistance): log-rank p and HR (95% CI)."""
    s = np.asarray(s, float); T = np.asarray(T, float); E = np.asarray(E, int)
    hi = (s > np.median(s)).astype(int)
    if hi.sum() in (0, len(hi)) or E.sum() == 0:
        return dict(lr_p=np.nan, split_hr=np.nan, split_lo=np.nan, split_hi=np.nan)
    try:
        lp = float(survdiff(T, E, hi)[1])
    except Exception:
        lp = np.nan
    try:
        r = PHReg(T, hi[:, None].astype(float), status=E.astype(float), ties="efron").fit(disp=False)
        b, se = float(r.params[0]), float(r.bse[0])
        return dict(lr_p=lp, split_hr=np.exp(b), split_lo=np.exp(b - 1.96 * se), split_hi=np.exp(b + 1.96 * se))
    except Exception:
        return dict(lr_p=lp, split_hr=np.nan, split_lo=np.nan, split_hi=np.nan)


def auc_stats(y, s, alternative="greater"):
    """AUC = P(score of a labelled-1 sample > score of a labelled-0 sample); one-sided MWU p (greater)."""
    y, s = np.asarray(y, int), np.asarray(s, float)
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < 2 or n0 < 2:
        return dict(auc=np.nan, p=np.nan)
    r = rankdata(s)
    auc = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    p = float(mannwhitneyu(s[y == 1], s[y == 0], alternative=alternative).pvalue)
    return dict(auc=float(auc), p=p)


def auc_fast(y, s):
    y = np.asarray(y, int); n1 = y.sum(); n0 = len(y) - n1
    if n1 < 2 or n0 < 2:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot_idx(n, strata, B, seed):
    """B resamples of patient indices, drawn within strata (cohorts) so that every cohort keeps its size."""
    rng = np.random.default_rng(seed)
    strata = np.zeros(n, int) if strata is None else pd.factorize(np.asarray(strata))[0]
    groups = [np.where(strata == g)[0] for g in np.unique(strata)]
    out = []
    for _ in range(B):
        out.append(np.concatenate([rng.choice(g, len(g), replace=True) for g in groups]))
    return out


def paired_boot(scores, fns, n, strata, B=500, seed=1):
    """scores: DataFrame (samples x models); fns: dict name -> f(full_score_array, resample_idx)->float.
    Returns dict name -> DataFrame (B x models) of the metric on each resample (same resamples for all models)."""
    idxs = boot_idx(n, strata, B, seed)
    out = {k: np.full((B, scores.shape[1]), np.nan) for k in fns}
    M = scores.values
    for b, ix in enumerate(idxs):
        for j in range(M.shape[1]):
            for k, f in fns.items():
                try:
                    out[k][b, j] = f(M[:, j], ix)
                except Exception:
                    pass
    return {k: pd.DataFrame(v, columns=scores.columns) for k, v in out.items()}


def ci(x, lo=2.5, hi=97.5):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return (np.percentile(x, lo), np.percentile(x, hi)) if len(x) > 10 else (np.nan, np.nan)


def boot_p(d):
    """Two-sided bootstrap p-value of a paired difference (share of resamples on each side of 0)."""
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    if len(d) < 10:
        return np.nan
    return float(min(1.0, 2 * (min((d <= 0).sum(), (d >= 0).sum()) + 1) / (len(d) + 1)))


def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan)
    ok = np.isfinite(p); q = p[ok]
    if not len(q):
        return out
    o = np.argsort(q); r = q[o] * len(q) / (np.arange(len(q)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(len(q)); qq[o] = np.minimum(r, 1.0); out[ok] = qq
    return out
