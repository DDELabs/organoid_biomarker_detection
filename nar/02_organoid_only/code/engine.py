"""The organoid-only strategy of Kong et al. 2020, run end to end for one lab-model set x drug x patient cohort.

Steps (python/single_pathway_prediction.py and python/multiple_pathway_prediction.py, unchanged):
  1. proximal pathways = Reactome pathways with network-proximity z <= cutoff (-1.2816) to the drug targets;
  2. pathway scores (ssGSEA) z-scored across lab models (StandardScaler) and across all patients of the cohort;
  3. Ridge (RidgeCV, cv=3, alpha 0.1-0.9) / linear SVR / LinearRegression fitted on lab-model drug response;
  4. biomarker = top-|coefficient| pathway(s); predicted response = sum(coef * patient score);
  5. drug-treated patients split at the median prediction (<= median = predicted responder); log-rank test (OS).
Added for reading the result (none change steps 1-5): hazard ratio (Cox) with 95% CI and events per group,
the same biomarker in untreated patients, leave-one-out fit in the models, the paper's random-pathway bootstrap
(Fig. 4e/f), a shuffle (permutation) of the lab-model response, and alternative proximity cut-offs.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import survdiff

from fastridge import ridgecv_coef

MLS = ("Ridge", "SVR", "LinearRegression")
KMAX = 10


def fit(ml, X, y):
    if ml == "Ridge":
        c, b, _ = ridgecv_coef(X, y)
        return c, b
    if ml == "SVR":
        m = SVR(kernel="linear").fit(X, y)
        return m.coef_[0], m.intercept_[0]
    m = LinearRegression().fit(X, y)
    return m.coef_, m.intercept_


def logrank(T, E, g):
    if g.sum() in (0, len(g)) or E.sum() == 0:
        return 1.0
    try:
        return float(survdiff(T, E, g)[1])
    except np.linalg.LinAlgError:          # zero variance (e.g. all deaths when only one group is at risk)
        return np.nan


def surv_stats(T, E, g, hr=True):
    """g = 1 predicted non-responder (prediction above the median), 0 predicted responder."""
    out = {"n": int(len(T)), "events": int(E.sum()), "n_resp": int((g == 0).sum()), "n_nonresp": int((g == 1).sum()),
           "events_resp": int(E[g == 0].sum()), "events_nonresp": int(E[g == 1].sum()), "p": logrank(T, E, g),
           "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan}
    if hr and 0 < g.sum() < len(g) and E[g == 0].sum() > 0 and E[g == 1].sum() > 0:
        try:
            r = PHReg(T, g.reshape(-1, 1).astype(float), status=E, ties="efron").fit(disp=False)
            b, se = float(r.params[0]), float(r.bse[0])
            out.update(HR=np.exp(b), HR_lo=np.exp(b - 1.96 * se), HR_hi=np.exp(b + 1.96 * se))
        except Exception:
            pass
    return out


def split(pred):
    return (pred > np.median(pred)).astype(int)


class Combo:
    """One lab-model set x drug x patient cohort."""

    def __init__(self, osc, y, tsc, sur, treated, untreated, z, cutoff=-1.2816, meta=None):
        self.meta = dict(meta or {})
        self.y_all = y.dropna()
        self.samples = sorted(s for s in self.y_all.index if s in osc.columns)
        self.y = self.y_all.loc[self.samples].values.astype(float)
        common = [p for p in osc.index if p in tsc.index]
        self.osc, self.z_all = osc.loc[common, self.samples], z
        pats = [p for p in tsc.columns if p in sur.index]
        self.tsc = tsc.loc[common, pats]
        # patient scaling over every patient of the cohort with survival (as the paper)
        self.P = pd.DataFrame(StandardScaler().fit_transform(self.tsc.T.values), index=pats, columns=common)
        # organoid scaling over the models with a response (as the paper)
        self.X = pd.DataFrame(StandardScaler().fit_transform(self.osc.T.values), index=self.samples, columns=common)
        self.tp = [p for p in pats if p in treated]
        self.up = [p for p in pats if p in untreated]
        self.T = sur.loc[self.tp, "months"].values.astype(float)
        self.E = sur.loc[self.tp, "event"].values.astype(int)
        self.Tu = sur.loc[self.up, "months"].values.astype(float)
        self.Eu = sur.loc[self.up, "event"].values.astype(int)
        self.common = common
        self.set_cutoff(cutoff)

    def set_cutoff(self, cutoff):
        self.cutoff = cutoff
        z = self.z_all.reindex(self.common)
        self.feats = [p for p in self.common if pd.notna(z[p]) and z[p] <= cutoff]
        return self.feats

    # ---------------------------------------------------------------- core
    def coefs(self, feats=None, y=None):
        feats = self.feats if feats is None else feats
        y = self.y if y is None else y
        X = self.X[feats].values
        return {ml: fit(ml, X, y)[0] for ml in MLS}

    def predict_split(self, feats, coef, k, Pmat):
        order = np.argsort(-np.abs(coef), kind="stable")[:k]
        pred = Pmat[:, order] @ coef[order]
        return split(pred), [feats[i] for i in order], coef[order]

    def evaluate(self, feats=None, y=None, untreated=True, hr=True):
        """Every (model, k) test: k = 1 is the paper's single-pathway biomarker, k = 2..10 multi-pathway."""
        feats = self.feats if feats is None else feats
        cf = self.coefs(feats, y)
        Pt = self.P.loc[self.tp, feats].values
        Pu = self.P.loc[self.up, feats].values
        rows = []
        for ml in MLS:
            for k in range(1, min(KMAX, len(feats)) + 1):
                g, used, cu = self.predict_split(feats, cf[ml], k, Pt)
                r = {"model": ml, "k": k, "pathways": ";".join(used), "coefs": ";".join(f"{c:.4g}" for c in cu)}
                r.update({f"treated_{a}": b for a, b in surv_stats(self.T, self.E, g, hr).items()})
                if untreated and len(self.up):
                    gu = self.predict_split(feats, cf[ml], k, Pu)[0]
                    r.update({f"untreated_{a}": b for a, b in surv_stats(self.Tu, self.Eu, gu).items()})
                rows.append(r)
        return pd.DataFrame(rows), cf

    def groups(self, ml="Ridge", k=1, feats=None, who="treated"):
        feats = self.feats if feats is None else feats
        cf = self.coefs(feats)[ml]
        Pm = self.P.loc[self.tp if who == "treated" else self.up, feats].values
        return self.predict_split(feats, cf, k, Pm)

    def loo(self, feats=None):
        """Leave-one-model-out prediction with all proximal pathways (scaler and model refitted per fold)."""
        feats = self.feats if feats is None else feats
        raw = self.osc.loc[feats].T.values
        out = {}
        for ml in MLS:
            pred = np.zeros(len(self.y))
            for i in range(len(self.y)):
                tr = np.arange(len(self.y)) != i
                sc = StandardScaler().fit(raw[tr])
                c, b = fit(ml, sc.transform(raw[tr]), self.y[tr])
                pred[i] = (sc.transform(raw[i:i + 1]) @ c)[0] + b
            rho, p = spearmanr(self.y, pred) if len(self.y) > 2 else (np.nan, np.nan)
            out[ml] = (pred, rho, p)
        return out

    # ---------------------------------------------------------------- nulls
    def rank1_p(self, feats, y):
        c = ridgecv_coef(self.X[feats].values, y)[0]
        i = int(np.argmax(np.abs(c)))
        pred = c[i] * self.P.loc[self.tp, feats[i]].values
        return logrank(self.T, self.E, split(pred))

    def single_pathway_p(self, sign=1):
        """Log-rank p of every pathway used alone (median split of sign * score). The coefficient's sign only
        decides on which side the median patient falls, so p is precomputed for both signs; reused by the
        bootstrap and by Supp. Fig. 10."""
        Pa = sign * self.P.loc[self.tp].values
        return pd.Series([logrank(self.T, self.E, split(Pa[:, j])) for j in range(len(self.common))], index=self.common)

    def top_rank(self, coef, feats, pvec, iso_p=None, alpha=0.05):
        """Paper's Fig. 4e statistic: first rank (by |coefficient|) whose pathway gives survival p < alpha in treated
        patients (and, when isogenic cell lines exist, isogenic t-test p < alpha). len(feats) + 1 = none found."""
        order = np.argsort(-np.abs(coef), kind="stable")
        for r, i in enumerate(order, 1):
            f = feats[i]
            if pvec[f] < alpha and (iso_p is None or (f in iso_p.index and iso_p[f] < alpha)):
                return r
        return len(feats) + 1

    def bootstrap(self, n=10000, seed=0, iso_p=None):
        """Paper Fig. 4e/f: random pathways matched in number to the proximal set -> Ridge. Returns, per random set,
        the rank-1 pathway's log-rank p and the paper's 'top rank' statistic."""
        rng = np.random.default_rng(seed)
        m = len(self.feats)
        Xa = self.X.values
        pp, pn = self.single_pathway_p(1).values, self.single_pathway_p(-1).values
        iso_ok = np.ones(len(self.common), bool) if iso_p is None else iso_p.reindex(self.common).fillna(1).values < 0.05
        ps, tr = np.empty(n), np.empty(n, dtype=int)
        for b in range(n):
            idx = rng.choice(len(self.common), m, replace=False)
            c = ridgecv_coef(Xa[:, idx], self.y)[0]
            o = np.argsort(-np.abs(c), kind="stable")
            order = idx[o]
            pv = np.where(c[o] > 0, pp[order], pn[order])
            ps[b] = pv[0]
            hit = np.nonzero((pv < 0.05) & iso_ok[order])[0]
            tr[b] = hit[0] + 1 if len(hit) else m + 1
        return ps, tr

    def permutation(self, n=1000, seed=1):
        """Shuffle the lab-model response n times; rerun rank-1 (Ridge) and the full 3 models x k=1..10 grid."""
        rng = np.random.default_rng(seed)
        r1, best = np.empty(n), np.empty(n)
        for b in range(n):
            yp = rng.permutation(self.y)
            ev, _ = self.evaluate(y=yp, untreated=False, hr=False)
            r1[b] = ev.loc[(ev.model == "Ridge") & (ev.k == 1), "treated_p"].iloc[0]
            best[b] = ev["treated_p"].min()
        return r1, best


def km_curve(T, E):
    from statsmodels.duration.survfunc import SurvfuncRight
    if len(T) == 0:
        return np.array([0.0]), np.array([1.0]), np.array([1.0]), np.array([1.0])
    sf = SurvfuncRight(T, E)
    t = np.r_[0, sf.surv_times]
    s = np.r_[1, sf.surv_prob]
    se = np.r_[0, sf.surv_prob_se]
    with np.errstate(divide="ignore", invalid="ignore"):
        # log-log CI as lifelines
        v = se ** 2 / (s ** 2 * np.log(s) ** 2)
        lo = s ** np.exp(1.96 * np.sqrt(v))
        hi = s ** np.exp(-1.96 * np.sqrt(v))
    lo[0] = hi[0] = 1
    lo, hi = np.nan_to_num(lo, nan=0.0), np.nan_to_num(hi, nan=1.0)
    return t, s, lo, hi
