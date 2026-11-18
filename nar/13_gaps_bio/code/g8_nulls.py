"""Shared pieces for G8 (mutation-concordance nulls and confounder adjustment)."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, norm, rankdata
from sklearn.linear_model import RidgeCV

SILENT = {"synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "upstream_gene_variant",
          "downstream_gene_variant", "intergenic_region", "non_coding_transcript_exon_variant", "splice_region_variant"}


def rz(x):
    r = rankdata(x)
    return norm.ppf((r - 0.5) / len(r))


def mwu(score, mut, alt):
    a, b = score[score.index.isin(mut)], score[~score.index.isin(mut)]
    return mannwhitneyu(a, b, alternative=alt).pvalue if len(a) >= 3 else np.nan


def tmb(mutfile, patients=None):
    m = pd.read_csv(mutfile, sep="\t", usecols=["sample", "gene", "effect"])
    m = m[m["sample"].str[13:15] == "01"]
    m["patient"] = m["sample"].str[:12]
    prof = set(m.patient)
    t = m[~m.effect.isin(SILENT)].groupby("patient").size().reindex(sorted(prof)).fillna(0)
    return t, m


def all_pathway_null(Pz, prof, mut, alt, obs_p):
    """Every pathway in both directions (random pathway with random sign): share with p <= observed."""
    S = Pz.loc[Pz.index.isin(prof)]
    ps = []
    for c in S.columns:
        for s in (1, -1):
            ps.append(mwu(s * S[c], mut, alt))
    ps = np.array(ps)
    return float(np.nanmean(ps <= obs_p + 1e-15)), len(ps), ps


def pipeline_null(X, y, Pz, pool, m, prof, mut, alt, obs_p, n=2000, seed=0):
    """Random sets of m pathways (the proximal-set size) -> Ridge on the lab models -> rank-1 pathway with its sign
    -> one-sided Mann-Whitney in TCGA. Mimics the whole procedure with non-proximal pathways."""
    rng = np.random.default_rng(seed)
    ys = (y - y.mean()) / y.std()
    S = Pz.loc[Pz.index.isin(prof)]
    ps = np.empty(n)
    cache = {}
    for b in range(n):
        idx = rng.choice(len(pool), m, replace=False)
        cols = [pool[i] for i in idx]
        c = RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(X[cols].values, ys).coef_
        i = int(np.argmax(np.abs(c)))
        key = (cols[i], np.sign(c[i]))
        if key not in cache:
            cache[key] = mwu(np.sign(c[i]) * S[cols[i]], mut, alt)
        ps[b] = cache[key]
    return float(np.nanmean(ps <= obs_p + 1e-15)), ps


def adjust(score, mut, cov, expect_sign):
    """OLS of rank-normalised score on mutation + covariates, and logistic of mutation on score + covariates.
    expect_sign: -1 if mutants are expected to have a LOWER predicted-resistance score."""
    df = pd.DataFrame({"score": score}).join(cov, how="inner").dropna()
    df["mut"] = df.index.isin(mut).astype(float)
    df["score_rz"] = rz(df["score"].values)
    covs = [c for c in cov.columns]
    X = sm.add_constant(df[["mut"] + covs].astype(float))
    o = sm.OLS(df["score_rz"], X).fit()
    b, p2 = o.params["mut"], o.pvalues["mut"]
    p1 = p2 / 2 if np.sign(b) == expect_sign else 1 - p2 / 2
    out = dict(n=len(df), n_mut=int(df.mut.sum()), ols_beta_mut=b, ols_lo=o.conf_int().loc["mut", 0],
               ols_hi=o.conf_int().loc["mut", 1], ols_p_one_sided=p1)
    try:
        L = sm.Logit(df["mut"], sm.add_constant(df[["score_rz"] + covs].astype(float))).fit(disp=0)
        bl, pl = L.params["score_rz"], L.pvalues["score_rz"]
        out.update(logit_OR_per_SD=np.exp(bl), logit_lo=np.exp(L.conf_int().loc["score_rz", 0]),
                   logit_hi=np.exp(L.conf_int().loc["score_rz", 1]),
                   logit_p_one_sided=pl / 2 if np.sign(bl) == expect_sign else 1 - pl / 2)
    except Exception as e:  # perfect separation etc.
        out["logit_error"] = str(e)[:60]
    return out
