"""Survival statistics without lifelines (log-rank, Kaplan-Meier, Cox, C-index)."""
import numpy as np
import pandas as pd
from scipy.stats import chi2


def logrank(t1, e1, t2, e2):
    """Two-group log-rank test p-value (matches lifelines.logrank_test)."""
    t = np.r_[t1, t2].astype(float)
    e = np.r_[e1, e2].astype(int)
    g = np.r_[np.zeros(len(t1)), np.ones(len(t2))]
    o = ex = v = 0.0
    for u in np.unique(t[e == 1]):
        at_risk = t >= u
        n, n1 = at_risk.sum(), (at_risk & (g == 0)).sum()
        dead = (t == u) & (e == 1)
        d, d1 = dead.sum(), (dead & (g == 0)).sum()
        o += d1
        ex += d * n1 / n
        if n > 1:
            v += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    return float(chi2.sf((o - ex) ** 2 / v, 1)) if v > 0 else 1.0


def km_at(t, e, at=60.0):
    """Kaplan-Meier survival probability at `at` months."""
    t, e = np.asarray(t, float), np.asarray(e, int)
    s = 1.0
    for u in np.unique(t[e == 1]):
        if u > at:
            break
        s *= 1 - ((t == u) & (e == 1)).sum() / (t >= u).sum()
    return s


def km_curve(t, e):
    t, e = np.asarray(t, float), np.asarray(e, int)
    xs, ys, s = [0.0], [1.0], 1.0
    for u in np.unique(t[e == 1]):
        s *= 1 - ((t == u) & (e == 1)).sum() / (t >= u).sum()
        xs.append(u)
        ys.append(s)
    xs.append(t.max())
    ys.append(s)
    return np.array(xs), np.array(ys)


def c_index(t, e, risk):
    """Harrell's C: higher risk should mean shorter survival."""
    t, e, r = np.asarray(t, float), np.asarray(e, int), np.asarray(risk, float)
    num = den = 0.0
    for i in np.where(e == 1)[0]:
        later = t > t[i]
        den += later.sum()
        num += (r[i] > r[later]).sum() + 0.5 * (r[i] == r[later]).sum()
    return num / den if den else np.nan


def cox(df, covariates, duration="months", event="event"):
    """Cox PH fit (statsmodels PHReg). Returns DataFrame of HR, CI, p per covariate."""
    from statsmodels.duration.hazard_regression import PHReg

    # drop adjustment covariates that are mostly missing (e.g. no AJCC stage in GBM)
    covariates = [covariates[0]] + [c for c in covariates[1:] if df[c].notna().mean() >= 0.7]
    d = df[[duration, event] + covariates].dropna()
    x = d[covariates].astype(float)
    keep = x.std() > 0
    # binary adjustment covariates need events in both levels, else the HR is infinite
    for c in covariates[1:]:
        v = x[c]
        if keep[c] and set(v.unique()) <= {0.0, 1.0}:
            ev = d[event]
            cells = [ev[v == 1].sum(), (1 - ev[v == 1]).sum(), ev[v == 0].sum(), (1 - ev[v == 0]).sum()]
            if min(cells) == 0 or min((v == 1).sum(), (v == 0).sum()) < 5:
                keep[c] = False
    x = x.loc[:, keep]
    res = PHReg(d[duration].values, x.values, status=d[event].values, ties="efron").fit(disp=False)
    if not np.all(np.isfinite(res.params)) and x.shape[1] > 1:
        return cox(df, covariates[:1], duration, event)  # fall back to unadjusted
    ci = np.exp(np.clip(res.conf_int(), -50, 50))
    out = pd.DataFrame({"HR": np.exp(res.params), "HR_low": ci[:, 0], "HR_high": ci[:, 1],
                        "p": res.pvalues, "n": len(d), "events": int(d[event].sum())}, index=x.columns)
    out.attrs["covariates"] = list(x.columns)
    return out


def covariate_frame(clinical, patients):
    """Numeric adjustment covariates: stage dummies (ref = I/II), age, sex."""
    c = clinical.loc[patients].copy()
    out = pd.DataFrame(index=c.index)
    out["stage_III"] = (c["stage"] == "III").astype(float)
    out["stage_IV"] = (c["stage"] == "IV").astype(float)
    out.loc[c["stage"].isna(), ["stage_III", "stage_IV"]] = np.nan
    out["age"] = c["age"] / 10.0
    out["sex"] = c["sex"]
    return out
