"""Shared helpers for the organoid + patient joint model (nar/03_organoid_plus_patient).

Reuses robust/obd: rank scoring (scoring.rank_score), Reactome gene sets, curated loaders,
StratCox (penalised stratified Cox from survlearn), BRIDGE's contrastive alignment (bridge.align),
and the survival helpers (log-rank, Kaplan-Meier, C-index).
"""
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import mannwhitneyu, norm, rankdata

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
REPO = FOLDER.parents[1]
sys.path.insert(0, str(REPO / "robust"))
warnings.filterwarnings("ignore")

from obd import bridge as BR  # noqa: E402
from obd import survival as SV  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.cohorts import to_canonical  # noqa: E402
from obd.curated import PRECLIN, load_curated_preclinical, load_curated_trial  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from obd.scoring import rank_score  # noqa: E402

TT = REPO / "robust" / "results" / "training_table"
CACHE = FOLDER / "cache"
TABLES = FOLDER / "tables"
FIGS = FOLDER / "figures"
for d in (CACHE, TABLES, FIGS, FIGS / "supp"):
    d.mkdir(parents=True, exist_ok=True)

SEED = 20240611
CUT = -1.2816


# ----------------------------------------------------------------------------------------------- patients
_P = {}


def patients():
    """Combined training table (20,766 samples) + pathway scores; duplicates dropped; GSE109211 added."""
    if "P" in _P:
        return _P["P"], _P["S"]
    P = pd.read_csv(TT / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    S = pd.concat([pd.read_csv(TT / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")])
    S = S.loc[P.index]
    P = P[P["possible_duplicate_of"].isna()].copy()
    S = S.loc[P.index]
    # STORM (GSE109211, sorafenib vs placebo) is not in the table: score it the same way
    expr, clin = load_curated_trial("GSE109211")
    sc = cached_scores("GSE109211", expr, canonical_done=True)
    clin.index = [f"GSE109211:{i}" for i in clin.index]
    sc.index = [f"GSE109211:{i}" for i in sc.index]
    add = pd.DataFrame(index=clin.index, columns=P.columns)
    add["cohort"], add["source"], add["cancer"], add["cancer_group"] = "GSE109211", "GEO (added here)", "hepatocellular", "liver"
    add["arm"] = clin["arm"].values
    add["drugs"] = np.where(clin["arm"].eq("sorafenib"), "SORAFENIB", "PLACEBO")
    add["treated"] = clin["arm"].eq("sorafenib").astype(float).values
    add["responder"] = clin["responder"].astype(float).values
    add["randomised"] = True
    add["has_response"] = True
    P = pd.concat([P, add])
    S = pd.concat([S, sc.reindex(index=add.index, columns=S.columns)])
    ok = S.notna().any(axis=1)                     # one STORM sample in the clinical file has no expression
    P, S = P[ok.values], S[ok.values]
    _P["P"], _P["S"] = P, S
    return P, S


def zscore_within(S, groups):
    return SL.zscore_within(S, pd.Series(np.asarray(groups), index=S.index))


def endpoint(d):
    """Per cohort: first of PFS, RFS, DFS, OS that is recorded. Returns time, event, name."""
    for a in ("pfs", "rfs", "dfs", "os"):
        ok = d[f"{a}_months"].notna() & d[f"{a}_event"].notna()
        if ok.sum() >= 0.6 * len(d) and d.loc[ok, f"{a}_event"].sum() >= 5:
            return d[f"{a}_months"].astype(float), d[f"{a}_event"].astype(float), a.upper()
    return None, None, None


def stage_num(s):
    t = str(s).upper().replace("STAGE", "").strip()
    for k, v in (("IV", 4), ("III", 3), ("II", 2), ("I", 1), ("4", 4), ("3", 3), ("2", 2), ("1", 1)):
        if t.startswith(k):
            return v
    return np.nan


# ----------------------------------------------------------------------------------------------- lab models
def _gs(index):
    gs = {k: [g for g in v if g in index] for k, v in gene_sets(("REACTOME",)).items()}
    return {k: v for k, v in gs.items() if len(v) >= 5}


def cached_scores(name, expr, canonical_done=False):
    """Reactome rank scores (samples x pathways), cached under nar/03.../cache."""
    f = CACHE / f"scores_{name}.tsv.gz"
    if f.exists():
        out = pd.read_csv(f, sep="\t", index_col=0)
        out.index = out.index.astype(str)
        return out
    e = expr if canonical_done else to_canonical(expr)
    sc = rank_score(e, _gs(set(e.index))).T
    sc.index = sc.index.astype(str)
    sc.to_csv(f, sep="\t")
    return sc


def vias2023():
    """Ovarian HGSOC organoids (Vias et al. 2023 eLife), fetched from the authors' repo (see REPORT)."""
    raw = FOLDER / "data" / "vias2023_raw"
    lib = {"19902": "119148", "19904": "54327", "19905": "118976", "19906": "119178", "19907": "119127",
           "19917": "23868", "19921": "119058", "19925": "54288", "19937": "54059", "19938": "151723",
           "19939": "151761", "19940": "54276", "19941": "32077", "19942": "151773"}
    t = pd.read_csv(raw / "20191218_ViasM_BJ_orgaBrs_tpm.csv")
    cols = {f"JBLAB{k}": v for k, v in lib.items()}
    t = t[["symbol"] + list(cols)].dropna(subset=["symbol"]).rename(columns=cols)
    expr = np.log2(t.groupby("symbol").sum() + 1)
    expr = expr.loc[expr.sum(1) > 0]
    a = pd.read_csv(raw / "20200419-AUC-organoids.csv")
    a["id"] = a["id"].astype(str)
    resp = a.set_index("id")[[c for c in a.columns if c.startswith("auc_ll5.")]]
    resp.columns = [c.replace("auc_ll5.", "").upper() for c in resp.columns]
    resp = resp.groupby(level=0).median()
    keep = [s for s in resp.index if s in expr.columns]
    return expr[keep], resp.loc[keep]


def lab_set(name, drug, lineage=None):
    """Return (scores samples x pathways, response Series higher = MORE sensitive, group labels)."""
    if name == "ovarian_vias2023":
        expr, resp = vias2023()
        sc = cached_scores(name, expr)
        direction = -1.0                                   # AUC of viability: lower = more sensitive
    else:
        expr, resp = load_curated_preclinical(name)
        sc = cached_scores(name, expr, canonical_done=True)
        direction = BR._direction(name)
    if drug not in resp:
        return None
    r = resp[drug].dropna()
    r = r[r.index.isin(sc.index)]
    if lineage is not None:
        s = pd.read_csv(PRECLIN / name / "samples.tsv", sep="\t", index_col=0)
        r = r[s["lineage"].reindex(r.index).eq(lineage).values]
    grp = pd.Series(lineage or "all", index=r.index)
    y = direction * r.astype(float)
    return sc.loc[r.index], y, grp


def lab_expression_scores(name, lineage=None):
    """All models of a set with expression (for the alignment figures), regardless of drug."""
    if name == "ovarian_vias2023":
        expr, _ = vias2023()
        return cached_scores(name, expr)
    expr, _ = load_curated_preclinical(name)
    sc = cached_scores(name, expr, canonical_done=True)
    if lineage is not None:
        s = pd.read_csv(PRECLIN / name / "samples.tsv", sep="\t", index_col=0)
        sc = sc[s["lineage"].reindex(sc.index).eq(lineage).values]
    return sc


def rank_normal(y):
    y = np.asarray(y, float)
    return norm.ppf((rankdata(y) - 0.5) / len(y))


# ----------------------------------------------------------------------------------------------- alignment
def coral(Xm, Xp, shrink=0.3):
    """CORAL (Sun et al. 2016): re-colour model features with the tumour covariance (shrunk to identity)."""
    def shrunk(X):
        C = np.cov(X, rowvar=False)
        d = np.diag(C).mean()
        return (1 - shrink) * C + shrink * d * np.eye(C.shape[0])

    def mpow(C, p):
        w, v = np.linalg.eigh(C)
        w = np.clip(w, 1e-8, None)
        return (v * w ** p) @ v.T

    Cm, Cp = shrunk(Xm), shrunk(Xp)
    mu = Xm.mean(0)
    return (Xm - mu) @ mpow(Cm, -0.5) @ mpow(Cp, 0.5)


def align_models(Xm, Xp, k=10, min_ratio=0.33, use_coral=True):
    """Celligner-like step: remove model-specific axes (BRIDGE contrastive PCA), then CORAL to tumours."""
    Xa, n_drop = BR.align(Xm, Xp, k=k, min_ratio=min_ratio)
    if use_coral and Xa.shape[0] >= 5:
        Xa = coral(Xa, Xp)
    Xa = (Xa - Xa.mean(0)) / np.where(Xa.std(0) > 0, Xa.std(0), 1)
    return Xa, n_drop


# ----------------------------------------------------------------------------------------------- joint model
class PatientBlock:
    """One patient cohort used for training: 'cox' (treated survival, stratum = cohort) or 'logit' (response)."""

    def __init__(self, kind, X, time=None, event=None, y=None, strata=None):
        self.kind, self.X = kind, np.asarray(X, float)
        if kind == "cox":
            self.time, self.event = np.asarray(time, float), np.asarray(event, float)
            self.strata = np.asarray(strata)
            self.m = SL.StratCox(self.time, self.event, self.strata)
            self.Xo = self.X[self.m.order]
            self.weight = float(self.event.sum())
        else:
            self.y = np.asarray(y, float)
            self.weight = float(len(self.y))

    def subset(self, idx):
        if self.kind == "cox":
            return PatientBlock("cox", self.X[idx], self.time[idx], self.event[idx], strata=self.strata[idx])
        return PatientBlock("logit", self.X[idx], y=self.y[idx])


def _block_loss(bl, beta, a):
    if bl.kind == "cox":
        f, g = bl.m.loss_grad(-beta, bl.Xo, np.zeros(len(beta)))   # benefit score: higher = lower hazard
        return f, -g, 0.0
    eta = a + bl.X @ beta
    p = 1 / (1 + np.exp(-eta))
    eps = 1e-12
    f = -np.mean(bl.y * np.log(p + eps) + (1 - bl.y) * np.log(1 - p + eps))
    r = (p - bl.y) / len(bl.y)
    return f, bl.X.T @ r, r.sum()


def fit_joint(Xo, yo, blocks, w, lam):
    """Shared pathway weights beta for both tasks.

    loss = w * 0.5 * mean((yo - c - Xo beta)^2) + (1 - w) * patient loss + 0.5 * lam * |beta|^2
    patient loss = event/patient-weighted mean of the blocks' per-event Cox loss or per-patient log loss.
    w = 1: organoid-only; w = 0: patient-only.
    """
    F = (Xo.shape[1] if Xo is not None else blocks[0].X.shape[1])
    nlog = sum(1 for b in blocks if b.kind == "logit")
    W = sum(b.weight for b in blocks) if blocks else 1.0

    def fun(theta):
        beta, c, a = theta[:F], theta[F], theta[F + 1:]
        f = 0.5 * lam * beta @ beta
        g = lam * beta.copy()
        gc, ga = 0.0, np.zeros(nlog)
        if w > 0 and Xo is not None:
            r = (c + Xo @ beta) - yo
            f += w * 0.5 * np.mean(r * r)
            g += w * Xo.T @ r / len(r)
            gc = w * r.mean()
        j = 0
        if w < 1:
            for bl in blocks:
                aj = a[j] if bl.kind == "logit" else 0.0
                fb, gb, gab = _block_loss(bl, beta, aj)
                s = (1 - w) * bl.weight / W
                f += s * fb
                g += s * gb
                if bl.kind == "logit":
                    ga[j] = s * gab
                    j += 1
        return f, np.r_[g, gc, ga]

    theta0 = np.zeros(F + 1 + nlog)
    with np.errstate(all="ignore"):
        r = minimize(fun, theta0, jac=True, method="L-BFGS-B", options={"maxiter": 800, "gtol": 1e-7})
    return r.x[:F], r.x[F + 1:]


def heldout_ll(train_blocks, full_blocks, test_masks, beta, a):
    """Held-out log-likelihood: Cox by Verweij-van Houwelingen (full minus train), logistic directly."""
    ll, j = 0.0, 0
    for tb, fb, te in zip(train_blocks, full_blocks, test_masks):
        if fb.kind == "cox":
            ll += fb.m.loglik(-beta, fb.X) - tb.m.loglik(-beta, tb.X)
        else:
            eta = a[j] + fb.X[te] @ beta
            p = np.clip(1 / (1 + np.exp(-eta)), 1e-9, 1 - 1e-9)
            ll += np.sum(fb.y[te] * np.log(p) + (1 - fb.y[te]) * np.log(1 - p))
            j += 1
    return ll


def inner_folds(blocks, k=3, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for bl in blocks:
        n = bl.X.shape[0]
        out.append(rng.permutation(n) % k)
    return out


LAMS = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
WS = (0.25, 0.5, 0.75)


def tune(Xo, yo, blocks, ws, lams=LAMS, k=3, seed=0):
    """Choose (w, lam) by inner patient cross-validation (held-out partial likelihood)."""
    folds = inner_folds(blocks, k, seed)
    best, res = None, {}
    for w in ws:
        for lam in lams:
            ll = 0.0
            for f in range(k):
                tr = [fo != f for fo in folds]
                if any(b.kind == "cox" and b.subset(t).event.sum() < 2 for b, t in zip(blocks, tr)):
                    continue
                tb = [b.subset(t) for b, t in zip(blocks, tr)]
                beta, a = fit_joint(Xo, yo, tb, w, lam)
                ll += heldout_ll(tb, blocks, [~t for t in tr], beta, a)
            res[(w, lam)] = ll
            if best is None or ll > res[best]:
                best = (w, lam)
    return best, res


def organoid_only(Xo, yo):
    from sklearn.linear_model import RidgeCV
    m = RidgeCV(alphas=np.logspace(-1, 4, 11)).fit(Xo, yo)
    return m.coef_, float(m.alpha_)


# ----------------------------------------------------------------------------------------------- evaluation
def cox_hr(score, time, event, strata=None):
    """HR per SD of the benefit score (statsmodels PHReg). HR < 1: higher score, longer survival."""
    from statsmodels.duration.hazard_regression import PHReg
    s = np.asarray(score, float)
    s = (s - s.mean()) / (s.std() or 1)
    r = PHReg(np.asarray(time, float), s[:, None], status=np.asarray(event, float),
              strata=None if strata is None else np.asarray(strata), ties="efron").fit(disp=False)
    b, se = float(r.params[0]), float(r.bse[0])
    return dict(HR=np.exp(b), HR_low=np.exp(b - 1.96 * se), HR_high=np.exp(b + 1.96 * se), p=float(r.pvalues[0]))


def concord(score, time, event):
    """Harrell's C for benefit score (higher score should mean LONGER survival). Returns (C, num, den)."""
    t, e, r = np.asarray(time, float), np.asarray(event, int), -np.asarray(score, float)
    num = den = 0.0
    for i in np.where(e == 1)[0]:
        later = t > t[i]
        den += later.sum()
        num += (r[i] > r[later]).sum() + 0.5 * (r[i] == r[later]).sum()
    return (num / den if den else np.nan), num, den


def auc_test(score, y):
    s, y = np.asarray(score, float), np.asarray(y, int)
    n1, n0 = y.sum(), (1 - y).sum()
    if n1 < 2 or n0 < 2:
        return dict(AUROC=np.nan, p=np.nan, num=0.0, den=0.0)
    r = rankdata(s)
    num = r[y == 1].sum() - n1 * (n1 + 1) / 2
    return dict(AUROC=num / (n1 * n0), p=float(mannwhitneyu(s[y == 1], s[y == 0]).pvalue), num=num, den=float(n1 * n0))


def logrank_median(score, time, event):
    s = np.asarray(score, float)
    hi = s > np.median(s)
    t, e = np.asarray(time, float), np.asarray(event, int)
    return SV.logrank(t[hi], e[hi], t[~hi], e[~hi]), hi


def interaction_cox(score, arm, time, event):
    from statsmodels.duration.hazard_regression import PHReg
    s = np.asarray(score, float)
    s = (s - s.mean()) / (s.std() or 1)
    a = np.asarray(arm, float)
    X = np.column_stack([a, s, a * s])
    r = PHReg(np.asarray(time, float), X, status=np.asarray(event, float), ties="efron").fit(disp=False)
    b, se = r.params, r.bse
    return dict(int_HR=np.exp(b[2]), int_low=np.exp(b[2] - 1.96 * se[2]), int_high=np.exp(b[2] + 1.96 * se[2]),
                int_p=float(r.pvalues[2]), arm_HR=np.exp(b[0]), score_HR_in_comparison=np.exp(b[1]))


def interaction_logit(score, arm, y):
    import statsmodels.api as sm
    s = np.asarray(score, float)
    s = (s - s.mean()) / (s.std() or 1)
    a = np.asarray(arm, float)
    X = sm.add_constant(np.column_stack([a, s, a * s]))
    r = sm.Logit(np.asarray(y, float), X).fit(disp=0)
    b, se = r.params, r.bse
    return dict(int_OR=np.exp(b[3]), int_low=np.exp(b[3] - 1.96 * se[3]), int_high=np.exp(b[3] + 1.96 * se[3]),
                int_p=float(r.pvalues[3]), arm_OR=np.exp(b[1]), score_OR_in_comparison=np.exp(b[2]))


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
