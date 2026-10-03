"""One drug x cancer study: baseline reproduction, robust signature, validation."""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from . import models as M
from . import survival as S


@dataclass
class Study:
    cancer: str
    drug: str
    org_scores: pd.DataFrame      # pathways x organoids
    response: pd.Series           # organoid -> IC50 (higher = resistant)
    pat_scores: pd.DataFrame      # pathways x patients
    clinical: pd.DataFrame        # patient -> months, event, stage, age, sex
    treated: set                  # patients treated with the drug
    features: list                # candidate (network-proximal) pathways
    label: str = ""
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        samples = [s for s in self.response.dropna().index if s in self.org_scores.columns]
        feats = [f for f in self.features if f in self.org_scores.index and f in self.pat_scores.index]
        pats = [p for p in self.pat_scores.columns if p in self.clinical.index]
        self.samples, self.features, self.patients = samples, feats, pats
        self.X = self.org_scores.loc[feats, samples].T.values.astype(float)
        self.y = self.response[samples].values.astype(float)
        # patients are standardised within the patient cohort, as in the original pipeline
        self.P = pd.DataFrame(M.zscore(self.pat_scores.loc[feats, pats].T.values), index=pats, columns=feats)
        self.treated_pats = [p for p in pats if p in self.treated]

    def summary(self):
        return {"cancer": self.cancer, "drug": self.drug, "data": self.label,
                "organoids": len(self.samples), "candidate_pathways": len(self.features),
                "patients_with_survival": len(self.patients), "treated_patients": len(self.treated_pats)}

    # ------------------------------------------------------------ patient side
    def patient_score(self, weights, pats=None):
        pats = self.treated_pats if pats is None else pats
        return self.P.loc[pats, weights.index].values @ weights.values

    def evaluate(self, weights, interaction=True):
        """Validate a weighted pathway signature (score = predicted IC50) in patients."""
        tp = self.treated_pats
        score = self.patient_score(weights, tp)
        cl = self.clinical.loc[tp]
        t, e = cl["months"].values, cl["event"].values
        resp = score <= np.median(score)
        out = {"n_treated": len(tp), "events": int(e.sum()),
               "logrank_p": S.logrank(t[resp], e[resp], t[~resp], e[~resp]),
               "os5_responder": S.km_at(t[resp], e[resp]), "os5_nonresponder": S.km_at(t[~resp], e[~resp]),
               "c_index": S.c_index(t, e, score)}
        df = S.covariate_frame(self.clinical, tp)
        df["score"] = (score - score.mean()) / (score.std() or 1.0)
        df["months"], df["event"] = t, e
        try:
            u = S.cox(df, ["score"]).loc["score"]
            a = S.cox(df, ["score", "stage_III", "stage_IV", "age", "sex"]).loc["score"]
            out.update({"cox_HR": u.HR, "cox_p": u.p, "adj_HR": a.HR, "adj_HR_low": a.HR_low,
                        "adj_HR_high": a.HR_high, "adj_p": a.p, "adj_n": int(a.n)})
        except Exception as exc:  # singular fits on tiny cohorts
            out.update({"cox_HR": np.nan, "cox_p": np.nan, "adj_HR": np.nan, "adj_p": np.nan, "error": str(exc)})
        if interaction:
            out.update(self.interaction(weights))
        return out

    def interaction(self, weights):
        """Predictive vs prognostic: Cox on all patients with score x treated interaction."""
        pats = self.patients
        score = self.patient_score(weights, pats)
        df = S.covariate_frame(self.clinical, pats)
        df["score"] = (score - score.mean()) / (score.std() or 1.0)
        df["treated"] = [float(p in self.treated) for p in pats]
        df["score_x_treated"] = df["score"] * df["treated"]
        df["months"] = self.clinical.loc[pats, "months"].values
        df["event"] = self.clinical.loc[pats, "event"].values
        try:
            r = S.cox(df, ["score", "treated", "score_x_treated", "stage_III", "stage_IV", "age", "sex"])
            return {"prognostic_HR_untreated": r.loc["score", "HR"], "prognostic_p": r.loc["score", "p"],
                    "interaction_HR": r.loc["score_x_treated", "HR"], "interaction_p": r.loc["score_x_treated", "p"]}
        except Exception:
            return {"interaction_HR": np.nan, "interaction_p": np.nan}

    # ------------------------------------------------------------ baseline (paper)
    def baseline(self, y=None, ks=range(2, 11)):
        """Kong et al. top-k pathways for Ridge / SVR / OLS, median split log-rank."""
        y = self.y if y is None else y
        tp = self.treated_pats
        t = self.clinical.loc[tp, "months"].values
        e = self.clinical.loc[tp, "event"].values
        Pt = self.P.loc[tp].values
        rows = []
        for model in M.BASELINE_MODELS:
            c = M.baseline_coefficients(self.X, y, model)
            order = np.argsort(-np.abs(c), kind="stable")
            for k in ks:
                idx = order[:k]
                score = Pt[:, idx] @ c[idx]
                r = score <= np.median(score)
                rows.append({"ML": model, "k": k, "pathways": ";".join(self.features[i] for i in idx),
                             "logrank_p": S.logrank(t[r], e[r], t[~r], e[~r]),
                             "os5_responder": S.km_at(t[r], e[r]), "os5_nonresponder": S.km_at(t[~r], e[~r])})
        return pd.DataFrame(rows)

    # ------------------------------------------------------------ robust signature
    def signature(self, y=None, X=None, top_k=7, n_boot=200, seed=0, align=False):
        X = self.X if X is None else X
        y = self.y if y is None else y
        weights, table = M.robust_signature(X, y, self.features, top_k=top_k, n_boot=n_boot, seed=seed)
        if align:  # PRECISE: re-estimate weights on organoid/patient shared directions
            basis, _ = M.precise_directions(X, self.P.values, n_components=min(10, X.shape[0] - 2), n_shared=5)
            w = pd.Series(M.precise_weights(X, y, basis), index=self.features)
            weights = w[weights.index]
            table["precise_coef"] = w
        return weights, table

    def organoid_cv(self, top_k=7, n_boot=50):
        """LOOCV of the full robust learner (selection inside each fold) vs the baseline."""
        def robust_fp(Xtr, ytr, Xte):
            w, _ = M.robust_signature(Xtr, ytr, list(range(Xtr.shape[1])), top_k=top_k, n_boot=n_boot)
            idx = w.index.values
            return M.standardise_with(Xtr, Xte)[:, idx] @ w.values

        def baseline_fp(Xtr, ytr, Xte):
            c = M.baseline_coefficients(Xtr, ytr, "Ridge")
            idx = np.argsort(-np.abs(c))[:top_k]
            return M.standardise_with(Xtr, Xte)[:, idx] @ c[idx]

        _, r_rho, r_p = M.loocv(self.X, self.y, robust_fp)
        _, b_rho, b_p = M.loocv(self.X, self.y, baseline_fp)
        return {"robust_loocv_spearman": r_rho, "robust_loocv_p": r_p,
                "baseline_loocv_spearman": b_rho, "baseline_loocv_p": b_p}

    # ------------------------------------------------------------ permutation nulls
    def permutation_null(self, n_perm=500, top_k=7, n_boot=60, n_jobs=4, seed=7):
        """Shuffle organoid IC50 and rerun everything downstream.

        robust: statistic = adjusted Cox z of the signature score (one-sided, HR>1).
        baseline: statistic = min log-rank p over 3 models x k=2..10, i.e. the
        multiplicity of the original 'pick the best k' readout.
        """
        rng = np.random.default_rng(seed)
        perms = [rng.permutation(self.y) for _ in range(n_perm)]

        def one(yp, i):
            w, _ = self.signature(y=yp, top_k=top_k, n_boot=n_boot, seed=i)
            ev = self.evaluate(w, interaction=False)
            z = np.log(ev["adj_HR"]) if np.isfinite(ev.get("adj_HR", np.nan)) else 0.0
            return z, ev["logrank_p"], self.baseline(y=yp)["logrank_p"].min()

        res = Parallel(n_jobs=n_jobs)(delayed(one)(yp, i) for i, yp in enumerate(perms))
        return pd.DataFrame(res, columns=["log_adj_HR", "robust_logrank_p", "baseline_min_p"])
