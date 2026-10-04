"""Frozen-signature validation in independent patient cohorts.

A signature (pathway -> weight, learned on organoids only) is applied without any
refitting. Pathway rank scores are computed per sample (platform-independent) and
z-scored within the cohort. Tests:
  * treated patients: Cox HR per SD of predicted resistance (adjusted where possible)
  * predictive test: score x treated interaction among comparable patients
  * response cohorts: AUC for separating responders (low score) from non-responders
  * random-signature null: the same tests for random pathway sets of equal size
    (Venet et al. 2011: random signatures are often prognostic in cancer)
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from . import survival as S
from .models import zscore


def signature_score(path_scores, weights):
    """path_scores: pathways x samples; weights: Series over pathways."""
    w = weights[[p for p in weights.index if p in path_scores.index]]
    if w.empty:
        return None
    Z = pd.DataFrame(zscore(path_scores.loc[w.index].T.values), index=path_scores.columns, columns=w.index)
    s = Z.values @ w.values
    return pd.Series((s - s.mean()) / (s.std() or 1.0), index=path_scores.columns)


def _cov(clin):
    df = pd.DataFrame(index=clin.index)
    if "stage" in clin and clin["stage"].notna().mean() > 0.7:
        df["stage_III"] = (clin["stage"] == "III").astype(float)
        df["stage_IV"] = (clin["stage"] == "IV").astype(float)
        df.loc[clin["stage"].isna(), ["stage_III", "stage_IV"]] = np.nan
    for c in ("age", "sex"):
        if c in clin and clin[c].notna().mean() > 0.7:
            df[c] = clin[c] / (10.0 if c == "age" else 1.0)
    return df


def survival_tests(score, clin, treated, untreated, time_col, event_col):
    out = {}
    d = _cov(clin)
    d["score"] = score.reindex(clin.index)
    d["months"], d["event"] = clin[time_col], clin[event_col]
    d = d.dropna(subset=["score", "months", "event"])
    covs = [c for c in d.columns if c not in ("score", "months", "event")]
    tr = d.loc[d.index.intersection(treated)]
    if len(tr) >= 20 and tr["event"].sum() >= 5:
        r = S.cox(tr, ["score"] + covs)
        out.update({"n_treated": len(tr), "events_treated": int(tr["event"].sum()),
                    "HR_treated": r.loc["score", "HR"], "HR_low": r.loc["score", "HR_low"],
                    "HR_high": r.loc["score", "HR_high"], "p_treated": r.loc["score", "p"],
                    "adjusted_for": ",".join(c for c in r.attrs.get("covariates", []) if c != "score"),
                    "c_index_treated": S.c_index(tr["months"], tr["event"], tr["score"])})
    if untreated is not None:
        un = d.loc[d.index.intersection(untreated)]
        both = pd.concat([tr.assign(treated=1.0), un.assign(treated=0.0)])
        if len(un) >= 20:
            both["score_x_treated"] = both["score"] * both["treated"]
            r = S.cox(both, ["score", "treated", "score_x_treated"] + covs)
            out.update({"n_untreated": len(un), "HR_untreated": r.loc["score", "HR"], "p_untreated": r.loc["score", "p"],
                        "interaction_HR": r.loc["score_x_treated", "HR"], "interaction_p": r.loc["score_x_treated", "p"]})
    return out


def response_test(score, responder):
    """AUC that a LOW predicted-resistance score identifies responders."""
    s = score.reindex(responder.index).dropna()
    r = responder.loc[s.index].astype(int)
    if r.nunique() < 2:
        return {}
    u, p = mannwhitneyu(-s[r == 1], -s[r == 0], alternative="greater")
    return {"n": len(s), "responders": int(r.sum()), "AUC": u / (r.sum() * (len(r) - r.sum())), "p_one_sided": p}


def random_null(path_scores, size, stat_fn, n=1000, pool=None, seed=11):
    """Distribution of stat_fn over random equal-size, random-sign signatures."""
    rng = np.random.default_rng(seed)
    pool = list(path_scores.index) if pool is None else [p for p in pool if p in path_scores.index]
    vals = []
    for _ in range(n):
        pick = rng.choice(pool, size, replace=False)
        w = pd.Series(rng.choice([-1.0, 1.0], size), index=pick)
        v = stat_fn(signature_score(path_scores, w))
        if v is not None and np.isfinite(v):
            vals.append(v)
    return np.array(vals)
