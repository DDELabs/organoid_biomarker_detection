"""Shared helpers for the revamp audit (read-only use of the existing pipeline).

Nothing here modifies existing code or results; outputs go to robust/revamp_audit/out/.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROBUST = HERE.parent
REPO = ROBUST.parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROBUST))

from obd import tcga_response as R  # noqa: E402
from run_atlas_relaxed import CR_LIKE, NR_LIKE, SUPPORTIVE, THERAPY  # noqa: E402

Z = stats.norm.ppf
PHI = stats.norm.cdf
LOG_HR = abs(np.log(0.7))  # 0.3567 per SD


def systemic_long():
    """Drug records restricted exactly as run_final_biomarker_survival.py does."""
    long = R.load_response_labels()
    return long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()].copy()


def expr_patients(project):
    """Patients with a cached Reactome rank-score profile (the set every TCGA analysis can use)."""
    import gzip
    f = ROBUST / "results" / "cache" / f"TCGA_{project}_reactome_rank_all.tsv.gz"
    if not f.exists():
        return set()
    with gzip.open(f, "rt") as fh:
        return set(fh.readline().rstrip("\n").split("\t")[1:])


# ---------------------------------------------------------------- power formulas
def events_main(hr=0.7, alpha=0.05, power=0.8, sd=1.0, r2=0.0):
    """Schoenfeld/Hsieh-Lavori: events for a continuous covariate (per-SD HR), two-sided alpha.
    r2 = squared multiple correlation of the marker with adjustment covariates (variance inflation)."""
    return (Z(1 - alpha / 2) + Z(power)) ** 2 / ((np.log(hr) * sd) ** 2 * (1 - r2))


def events_interaction(hr_int=0.7, alpha=0.05, power=0.8, frac_events_treated=0.5):
    """Events for a treatment x continuous-marker interaction (marker SD 1 in each arm).
    var(beta_int) ~= 1/D_t + 1/D_u = 1/(D * f * (1-f))."""
    f = frac_events_treated
    return (Z(1 - alpha / 2) + Z(power)) ** 2 / (np.log(hr_int) ** 2 * f * (1 - f))


def power_from_var(log_effect, var, alpha=0.05):
    if not np.isfinite(var) or var <= 0:
        return np.nan
    se = np.sqrt(var)
    z = Z(1 - alpha / 2)
    d = abs(log_effect) / se
    return float(PHI(d - z) + PHI(-d - z))


def power_main(D, hr=0.7, alpha=0.05):
    return power_from_var(np.log(hr), 1.0 / D, alpha) if D and D > 0 else np.nan


def power_interaction(Dt, Du, hr=0.7, alpha=0.05):
    if not Dt or not Du or Dt <= 0 or Du <= 0:
        return np.nan
    return power_from_var(np.log(hr), 1.0 / Dt + 1.0 / Du, alpha)


def mdhr(D, alpha=0.05, power=0.8):
    """Minimum detectable per-SD HR (< 1) with D events."""
    if not D or D <= 0:
        return np.nan
    return float(np.exp(-(Z(1 - alpha / 2) + Z(power)) / np.sqrt(D)))


def auc_to_delta(auc):
    """Binormal equal-variance: AUC = Phi(delta / sqrt 2); delta ~ log OR per within-class SD."""
    return np.sqrt(2) * Z(auc)


def power_auc(n1, n0, auc=0.65, alpha=0.05):
    """Power to detect AUROC `auc` vs 0.5 (two-sample z on the standardised mean difference)."""
    if not n1 or not n0 or n1 < 2 or n0 < 2:
        return np.nan
    return power_from_var(auc_to_delta(auc), 1.0 / n1 + 1.0 / n0, alpha)


def n_for_auc(auc, prevalence=0.5, alpha=0.05, power=0.8):
    d = auc_to_delta(auc)
    k = (Z(1 - alpha / 2) + Z(power)) ** 2 / d ** 2  # = 1/n1 + 1/n0 must be <= 1/k
    return k * (1 / prevalence + 1 / (1 - prevalence))


def jung_alpha(m, m1, q=0.1, power=0.8):
    """Per-test alpha giving expected FDR q with m tests, m1 true effects, each at `power` (Jung 2005)."""
    m0 = m - m1
    return m1 * power * q / (m0 * (1 - q))


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    if ok.sum():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty_like(pv)
        out[o] = np.minimum(r, 1)
        q[ok] = out
    return q


def reverse_km_median(time, event):
    """Median potential follow-up (reverse Kaplan-Meier: censoring treated as the event)."""
    t = np.asarray(time, float)
    c = 1 - np.asarray(event, float)
    ok = np.isfinite(t) & np.isfinite(c) & (t >= 0)
    t, c = t[ok], c[ok]
    if len(t) == 0:
        return np.nan
    o = np.argsort(t)
    t, c = t[o], c[o]
    s, n = 1.0, len(t)
    for i, (ti, ci) in enumerate(zip(t, c)):
        if ci:
            s *= 1 - 1 / (n - i)
            if s <= 0.5:
                return float(ti)
    return np.nan
