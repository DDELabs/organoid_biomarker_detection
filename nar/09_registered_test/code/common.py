"""Shared helpers: frozen hypotheses (hash-checked), training table, within-cohort scores, Cox with
entry/strata/weights, small stats utilities."""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
ROBUST = REPO / "robust"
TT = ROBUST / "results" / "training_table"
FIG, TAB, DATA = HERE / "figures", HERE / "tables", HERE / "_data"
for d in (FIG, TAB, DATA):
    d.mkdir(exist_ok=True)
sys.path.insert(0, str(ROBUST))
sys.path.insert(0, str(HERE / "code"))

FROZEN = REPO / "nar" / "08_depmap" / "tables" / "frozen_hypotheses.tsv"
FROZEN_SHA = "a0bdd620293e2559667cf7f4c218894ae3e70cf1eb404558a084d28aee878ab7"
PURITY = ROBUST / "revamp_audit" / "data" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
SEED = 20261117
DRUGS = ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "GEMCITABINE", "PACLITAXEL", "PEMETREXED", "SORAFENIB"]
LABEL = {"FLUOROURACIL": "5-FU", "OXALIPLATIN": "Oxaliplatin", "IRINOTECAN": "Irinotecan", "GEMCITABINE": "Gemcitabine",
         "PACLITAXEL": "Paclitaxel", "PEMETREXED": "Pemetrexed", "SORAFENIB": "Sorafenib"}
PROLIF = ["REACTOME_CELL_CYCLE_MITOTIC", "REACTOME_DNA_REPLICATION", "REACTOME_MITOTIC_M_M_G1_PHASES",
          "REACTOME_G2_M_CHECKPOINTS", "REACTOME_S_PHASE"]
# Okabe-Ito colour-blind-safe palette
CB = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00", "purple": "#CC79A7",
      "sky": "#56B4E9", "yellow": "#F0E442", "grey": "#7F7F7F", "black": "#000000"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def frozen():
    got = sha256(FROZEN)
    if got != FROZEN_SHA:
        raise SystemExit(f"frozen_hypotheses.tsv hash changed: {got}")
    f = pd.read_csv(FROZEN, sep="\t")
    return {d: dict(zip(g["pathway"], np.where(g["sign"].eq("+"), 1.0, -1.0))) for d, g in f.groupby("drug")}


def load_table():
    P = pd.read_csv(TT / "patients.tsv.gz", sep="\t", low_memory=False, index_col="sample_id")
    S = pd.concat([pd.read_csv(TT / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")])
    P = P[P["possible_duplicate_of"].isna()]
    S = S.reindex(P.index)
    return P, S


def zwithin(S, groups):
    g = S.groupby(groups.values)
    return (S - g.transform("mean")) / g.transform("std")


def signed_score(Z, weights, groups):
    """Mean of sign x standardised pathway score, re-standardised within group."""
    cols = list(weights)
    s = (Z[cols] * pd.Series(weights)).mean(axis=1)
    g = s.groupby(groups.values)
    return (s - g.transform("mean")) / g.transform("std")


def stage_num(s):
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return np.nan
    t = str(s).upper().replace("STAGE", "").strip()
    for k, v in (("IV", 4), ("III", 3), ("II", 2), ("I", 1)):
        if t.startswith(k):
            return v
    if t[:1].isdigit():
        v = int(t[0])
        return v if 1 <= v <= 4 else np.nan
    return np.nan


# ------------------------------------------------------------------ Cox (Breslow) with entry/strata/weights
def _prep(time, event, X, entry, strata, w):
    n = len(time)
    entry = np.zeros(n) if entry is None else np.asarray(entry, float)
    strata = np.zeros(n, int) if strata is None else pd.factorize(np.asarray(strata))[0]
    w = np.ones(n) if w is None else np.asarray(w, float)
    out = []
    for s in np.unique(strata):
        m = strata == s
        t, e, en, ww, XX = np.asarray(time, float)[m], np.asarray(event, float)[m], entry[m], w[m], X[m]
        if e.sum() == 0:
            continue
        T = np.unique(t[e > 0])
        oe, oen = np.argsort(t, kind="stable"), np.argsort(en, kind="stable")
        ie = np.searchsorted(t[oe], T, "left")
        ien = np.searchsorted(en[oen], T, "left")   # rows with entry >= T are not at risk
        # weighted events per event time
        ev = e > 0
        D = np.bincount(np.searchsorted(T, t[ev]), weights=ww[ev], minlength=len(T))
        sx = (ww[ev, None] * XX[ev]).sum(0)
        out.append((XX, ww, oe, oen, ie, ien, D, sx))
    return out


def _suffix(a):
    return np.concatenate([np.cumsum(a[::-1], axis=0)[::-1], np.zeros((1,) + a.shape[1:])])


def cox_fit(time, event, X, entry=None, strata=None, weights=None, max_iter=50, tol=1e-9):
    X = np.asarray(X, float)
    p = X.shape[1]
    parts = _prep(time, event, X, entry, strata, weights)
    b = np.zeros(p)
    ll_old, step = -np.inf, np.zeros(p)
    for it in range(max_iter):
        ll, g, H = 0.0, np.zeros(p), np.zeros((p, p))
        for XX, ww, oe, oen, ie, ien, D, sx in parts:
            eta = XX @ b
            r = ww * np.exp(eta - eta.max())
            rx = r[:, None] * XX
            rxx = rx[:, :, None] * XX[:, None, :]
            S0 = _suffix(r[oe])[ie] - _suffix(r[oen])[ien]
            S1 = _suffix(rx[oe])[ie] - _suffix(rx[oen])[ien]
            S2 = _suffix(rxx[oe])[ie] - _suffix(rxx[oen])[ien]
            S0 = np.maximum(S0, 1e-300)
            ll += (sx @ b) - (D * (np.log(S0) + eta.max())).sum()
            m1 = S1 / S0[:, None]
            g += sx - (D[:, None] * m1).sum(0)
            H -= (D[:, None, None] * (S2 / S0[:, None, None] - m1[:, :, None] * m1[:, None, :])).sum(0)
        if (not np.isfinite(ll)) or (ll < ll_old - 1e-8 and it > 0):
            # step-halving: go back half of the previous step
            step = step / 2
            b = b - step
            if np.abs(step).max() < 1e-10:
                break
            continue
        try:
            step = np.linalg.solve(-H, g)
        except np.linalg.LinAlgError:
            try:
                step = np.linalg.lstsq(-H, g, rcond=None)[0]
            except np.linalg.LinAlgError:
                nan = np.full(p, np.nan)
                return nan, nan, np.full((p, p), np.nan)
        if not np.all(np.isfinite(step)):
            nan = np.full(p, np.nan)
            return nan, nan, np.full((p, p), np.nan)
        b = b + step
        if abs(ll - ll_old) < tol and np.abs(step).max() < 1e-6:
            break
        ll_old = ll
    try:
        cov = np.linalg.inv(-H)
    except Exception:
        nan = np.full(p, np.nan)
        return b, nan, np.full((p, p), np.nan)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return b, se, cov


def wald(b, se):
    z = b / se
    return z, 2 * stats.norm.sf(abs(z))


def holm(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    out = np.full_like(p, np.nan)
    pv = p[ok]
    o = np.argsort(pv)
    m = len(pv)
    adj = np.maximum.accumulate((m - np.arange(m)) * pv[o])
    r = np.empty(m)
    r[o] = np.minimum(adj, 1)
    out[ok] = r
    return out


def dersimonian_laird(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float)
    w = 1 / se ** 2
    fe = (w * y).sum() / w.sum()
    Q = (w * (y - fe) ** 2).sum()
    k = len(y)
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / c) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    re = (wr * y).sum() / wr.sum()
    se_re = np.sqrt(1 / wr.sum())
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 and k > 1 else 0.0
    return {"fe": fe, "se_fe": np.sqrt(1 / w.sum()), "re": re, "se_re": se_re, "tau2": tau2, "I2": I2, "Q": Q, "k": k,
            "p_re": 2 * stats.norm.sf(abs(re / se_re)), "p_fe": 2 * stats.norm.sf(abs(fe * np.sqrt(w.sum())))}


Z = stats.norm.ppf


def power_two_sided(log_eff, se, alpha=0.05):
    d = abs(log_eff) / se
    z = Z(1 - alpha / 2)
    return float(stats.norm.cdf(d - z) + stats.norm.cdf(-d - z))


def power_one_sided(log_eff, se, alpha=0.025):
    return float(stats.norm.cdf(abs(log_eff) / se - Z(1 - alpha)))


def mdhr(se, alpha=0.05, power=0.8, two_sided=True):
    za = Z(1 - alpha / 2) if two_sided else Z(1 - alpha)
    return float(np.exp((za + Z(power)) * se))


def save_fig(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")
