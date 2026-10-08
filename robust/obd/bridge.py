"""BRIDGE: patient-derived models -> patient tumours, pathway-level drug-response transfer.

Domains (pre-clinical): PDX in vivo (Novartis PDXE), patient-derived organoids, cell lines.
For each domain and drug:
  1. Reactome rank scores of the models, standardised within tumour type.
  2. Alignment: principal axes of the model data that carry little variance in patient tumours
     (variance ratio patient/model < `min_ratio`, e.g. loss of stroma/immune, culture effects)
     are projected out (contrastive-PCA style; linear, so pathways stay interpretable).
  3. Sensitivity effect vector: ridge (RidgeCV) of within-tumour-type rank-normalised sensitivity
     on pathways, rescaled to unit sd.
Domain vectors are averaged within domain class (pdx / organoid / cell_line; sqrt(n) weights) and then
across classes, giving the pre-clinical prior beta_pre.
The patient model is TRIAD's prior-centred ridge logistic with prior mean s * TAU * beta_pre; the
transfer strength s is learned per drug by leave-one-cancer-out (nested for the reported AUROC).
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from sklearn.linear_model import RidgeCV

from . import scoring
from .curated import PRECLIN, load_curated_preclinical
from .reference import gene_sets
from .triad import fit_prior_logistic, loco

DOMAINS = {"pdx_novartis_gao2015": "pdx",
           "pancreas_shi2022": "organoid", "pancreas_tiriac2018": "organoid", "liver_ji2023": "organoid",
           "bladder_lee2018": "organoid", "sarcoma_alshihabi2024": "organoid", "liver_broutier2017": "organoid",
           "ctrpv2_ccle": "cell_line", "prism_repurposing": "cell_line"}
TAU = 0.05
GRID = (0.0, 1.0, 2.0, 4.0, 8.0)


def _direction(name):
    r = pd.read_csv(PRECLIN / name / "response.tsv", sep="\t", usecols=["direction"], nrows=1)
    return -1.0 if str(r["direction"].iloc[0]).startswith("lower") else 1.0


def _groups(name, samples):
    f = PRECLIN / name / "samples.tsv"
    if f.exists():
        s = pd.read_csv(f, sep="\t", index_col=0)
        col = next((c for c in ("tumour_type", "lineage", "tissue", "cancer", "primary_disease") if c in s), None)
        if col is not None:
            return s[col].reindex(samples).fillna("NA").astype(str).values
    return np.array(["all"] * len(samples))


def domain_scores(name, pathways):
    expr, resp = load_curated_preclinical(name)
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    sc = scoring.cached(f"BRIDGE_{name}_reactome_rank", scoring.rank_score, expr,
                        {k: v for k, v in gs.items() if len(v) >= 5})
    return sc.reindex(pathways).T, resp


def _within(X, groups):
    X = X.copy()
    for g in np.unique(groups):
        i = groups == g
        mu, sd = X[i].mean(0), X[i].std(0)
        X[i] = (X[i] - mu) / np.where(sd > 0, sd, 1)
    return np.nan_to_num(X)


def align(Xm, Xp, k=10, min_ratio=0.33):
    """Project out model-specific axes: top-k PCs of the model data whose variance in patients is low."""
    u, s, vt = np.linalg.svd(Xm - Xm.mean(0), full_matrices=False)
    V = vt[:k]
    vm = (Xm @ V.T).var(0)
    vp = (Xp @ V.T).var(0)
    drop = V[(vp / np.maximum(vm, 1e-9)) < min_ratio]
    if len(drop):
        Xm = Xm - (Xm @ drop.T) @ drop
    return Xm, len(drop)


def domain_vector(name, drug, pathways, Xp_ref, min_n=12):
    """Unit-sd sensitivity effect vector of one domain for one drug (None if not screened)."""
    try:
        sc, resp = domain_scores(name, pathways)
    except Exception:
        return None
    if drug not in resp:
        return None
    r = resp[drug].dropna()
    r = r[r.index.isin(sc.index)]
    if len(r) < min_n:
        return None
    grp = _groups(name, r.index)
    X = _within(sc.loc[r.index].values.astype(float), grp)
    y = _direction(name) * r.values.astype(float)  # higher = more sensitive
    yz = np.empty(len(y))
    for g in np.unique(grp):
        i = grp == g
        yz[i] = norm.ppf((rankdata(y[i]) - 0.5) / i.sum()) if i.sum() > 2 else 0.0
    X, n_drop = align(X, Xp_ref)
    b = RidgeCV(alphas=np.logspace(0, 4, 9)).fit(X, yz).coef_
    return {"domain": name, "class": DOMAINS[name], "n": len(r), "dropped_axes": n_drop,
            "beta": pd.Series(b / (b.std() or 1), index=pathways)}


def combine(vectors):
    """Average within domain class (sqrt-n weights), then across classes."""
    by = {}
    for v in vectors:
        by.setdefault(v["class"], []).append(v)
    cls = {}
    for c, vs in by.items():
        w = np.array([np.sqrt(v["n"]) for v in vs])
        b = sum(wi * v["beta"] for wi, v in zip(w, vs)) / w.sum()
        cls[c] = b / (b.std() or 1)
    if not cls:
        return None, {}
    b = sum(cls.values()) / len(cls)
    return b / (b.std() or 1), cls


def prior_mean(beta, s, n_cov):
    return np.r_[s * TAU * beta.values, np.zeros(n_cov)]


def choose_s(X, y, groups, lam, beta, n_cov, k):
    """Transfer strength by LOCO mean within-cancer AUROC over GRID."""
    best, res = None, {}
    for s in GRID:
        auc = loco(X, y, groups, lam, None, prior_mean(beta, s, n_cov) if beta is not None else None, n_score=k)[3]
        res[s] = auc
        if np.isfinite(auc) and (best is None or auc > res[best] + 1e-9):
            best = s
    return (0.0 if best is None else best), res


def nested_loco(X, y, groups, lam, beta, n_cov, k, min_test=10):
    """Outer leave-one-cancer-out; s chosen by inner LOCO on the training cancers only."""
    from sklearn.metrics import roc_auc_score
    pred = np.full(len(y), np.nan)
    per, chosen = {}, {}
    for g in np.unique(groups):
        te = groups == g
        if te.sum() < min_test or len(np.unique(y[te])) < 2:
            continue
        tr = ~te
        if beta is None or len(np.unique(groups[tr])) < 3:
            s = 0.0
        else:
            s, _ = choose_s(X[tr], y[tr], groups[tr], lam, beta, n_cov, k)
        counts = pd.Series(groups[tr]).value_counts()
        sw = 1.0 / pd.Series(groups[tr]).map(counts).values
        sw = sw * len(sw) / sw.sum()
        m = prior_mean(beta, s, n_cov) if beta is not None else None
        _, b = fit_prior_logistic(X[tr], y[tr], lam, None, m, sw)
        pred[te] = X[te, :k] @ b[:k]
        per[g] = roc_auc_score(y[te], pred[te])
        chosen[g] = s
    return float(np.mean(list(per.values()))) if per else np.nan, per, chosen
