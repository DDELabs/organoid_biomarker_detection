"""Apply an organoid-only (paper-style) model to the patient cohorts of robust/results/training_table_v2.
Shared by the lung, ovary and liver scripts.  One CPU."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u22
from u22 import C, N, pd, np, REPO, DATA
import engine as EN
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import survdiff
from scipy.stats import chi2, mannwhitneyu, norm

TT2 = REPO / "robust" / "results" / "training_table_v2"
CUT = -1.2816
_cache = {}


def v2():
    if "P" in _cache:
        return _cache["P"], _cache["S"]
    P = pd.read_csv(TT2 / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    S = pd.concat([pd.read_csv(TT2 / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")])
    S = S.loc[P.index]
    P["duplicate_flag"] = P["possible_duplicate_of"].notna()
    _cache["P"], _cache["S"] = P, S
    return P, S


def zwithin(S):
    sd = S.std(0).replace(0, np.nan)
    return ((S - S.mean(0)) / sd).fillna(0.0)


def endpoint(d, prefer=("os",)):
    for a in prefer:
        ok = d[f"{a}_months"].notna() & d[f"{a}_event"].notna() & (d[f"{a}_months"] > 0)
        if ok.sum() >= 10:
            return d.loc[ok, f"{a}_months"].astype(float).values, d.loc[ok, f"{a}_event"].astype(int).values, ok, a.upper()
    return None, None, None, None


def logrank(T, E, g):
    g = np.asarray(g)
    if len(set(g)) < 2 or np.sum(E) == 0:
        return np.nan
    try:
        return float(survdiff(T, E, g)[1])
    except Exception:
        return np.nan


def cox_binary(T, E, g):
    """HR of group 1 vs 0 with 95% CI."""
    g = np.asarray(g, float)
    if len(set(g)) < 2 or E[g == 0].sum() == 0 or E[g == 1].sum() == 0:
        return np.nan, np.nan, np.nan
    try:
        r = PHReg(T, g.reshape(-1, 1), status=E, ties="efron").fit(disp=False)
        b, se = float(r.params[0]), float(r.bse[0])
        return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))
    except Exception:
        return np.nan, np.nan, np.nan


def interaction(T, E, g, arm):
    """Cox model group + arm + group*arm: HR (95% CI) and p of the product term."""
    X = np.column_stack([g, arm, np.asarray(g) * np.asarray(arm)]).astype(float)
    if len(set(arm)) < 2 or len(set(g)) < 2:
        return np.nan, np.nan, np.nan, np.nan
    try:
        r = PHReg(T, X, status=E, ties="efron").fit(disp=False)
        b, se = float(r.params[2]), float(r.bse[2])
        return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(2 * norm.sf(abs(b / se)))
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def shortlist(drug, columns):
    z, col = u22.proximity(drug)
    if z is None:
        return [], col
    z = z.reindex(columns).dropna()
    return list(z.index[z <= CUT]), col


def organoid_fit(Sorg, y, feats):
    """Ridge on z-scored organoid pathway scores (paper's RidgeCV); returns coef (Series) over feats."""
    X = Sorg[feats]
    X = (X - X.mean()) / X.std().replace(0, np.nan)
    X = X.fillna(0.0)
    c, b = EN.fit("Ridge", X.loc[y.index].values, y.values)
    return pd.Series(c, index=feats)


def rank1_score(Zp, coef):
    i = coef.abs().idxmax()
    return Zp[i] * coef[i], i


def multi_score(Zp, coef):
    return Zp[coef.index] @ coef


def split(score):
    return (np.asarray(score) > np.median(score)).astype(int)


def evaluate_cohort(Zp, T_idx, treated_idx, comp_idx, endpoint_pref, P, coef, how="rank1"):
    """Returns dict with treated and comparison statistics for one score definition."""
    score = rank1_score(Zp, coef)[0] if how == "rank1" else multi_score(Zp, coef)
    out = {}
    d = P
    res = {}
    for label, idx in (("treated", treated_idx), ("comparison", comp_idx)):
        if idx is None or len(idx) == 0:
            continue
        x = d.loc[idx]
        T, E, ok, ep = endpoint(x, endpoint_pref)
        if T is None:
            continue
        sc = score.loc[x.index[ok.values]]
        g = split(score.loc[idx].values)          # split at the median of the arm itself (paper: median of treated patients)
        g = pd.Series(g, index=idx).loc[x.index[ok.values]].values
        hr = cox_binary(T, E, g)
        res[label] = dict(n=len(T), events=int(E.sum()), logrank_p=logrank(T, E, g), HR=hr[0], HR_lo=hr[1], HR_hi=hr[2], endpoint=ep, T=T, E=E, g=g)
    return res, score
