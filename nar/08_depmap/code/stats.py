"""Small statistics helpers: Spearman, lineage-adjusted partial Spearman, BH FDR, random-signature nulls."""
import numpy as np
import pandas as pd
from scipy import stats

from common import lineage_dummies, rank_z, residualize


def spearman(x, y):
    x, y = pd.Series(x), pd.Series(y)
    ok = x.notna() & y.notna()
    if ok.sum() < 5:
        return np.nan, np.nan, int(ok.sum())
    r, p = stats.spearmanr(x[ok], y[ok])
    return float(r), float(p), int(ok.sum())


def partial_spearman(x, y, lineage, extra=None, min_lineage=5):
    """Spearman correlation of x and y after removing lineage (and optional extra covariates).

    Both variables are rank-normalised, residualised on lineage dummies (+extra) and correlated;
    p from the t distribution with n - 2 - k degrees of freedom.
    """
    df = pd.DataFrame({"x": x, "y": y, "lin": lineage})
    if extra is not None:
        df = df.join(extra)
    df = df.dropna()
    n = len(df)
    if n < 10:
        return np.nan, np.nan, n
    cov = lineage_dummies(df["lin"], min_lineage)
    if extra is not None:
        cov = pd.concat([cov, df[extra.columns].astype(float)], axis=1)
    rx = residualize(rank_z(df["x"].values), cov)
    ry = residualize(rank_z(df["y"].values), cov)
    r = float(np.corrcoef(rx, ry)[0, 1])
    k = cov.shape[1]
    dof = n - 2 - k
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), dof)), n


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def sig_score(scores, weights):
    """Frozen signature applied unchanged: z-score each pathway across the samples, weighted sum, z-score."""
    w = weights[[p for p in weights.index if p in scores.index]]
    if w.empty:
        return None
    Z = scores.loc[w.index].T
    Z = (Z - Z.mean()) / Z.std(ddof=0).replace(0, 1)
    s = Z.values @ w.values
    return pd.Series((s - s.mean()) / (s.std() or 1.0), index=scores.columns)


def random_signature_null(scores, k, y, n=1000, seed=11):
    """Spearman rho between y and n random k-pathway random-sign signatures."""
    rng = np.random.default_rng(seed)
    pool = np.array(scores.index)
    out = []
    for _ in range(n):
        pick = rng.choice(pool, k, replace=False)
        w = pd.Series(rng.choice([-1.0, 1.0], k), index=pick)
        s = sig_score(scores, w)
        out.append(stats.spearmanr(s.values, y.loc[s.index].values)[0])
    return np.array(out)
