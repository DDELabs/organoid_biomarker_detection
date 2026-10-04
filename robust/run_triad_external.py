"""External validation of TRIAD drug models in clinical-trial cohorts (frozen models).

Each model variant (patient_only, network, cell_line_prior, organoid_prior, TRIAD_all,
and the frozen pre-clinical priors themselves) is fitted on ALL TCGA patients labelled for
the drug, frozen, and applied to independent trial cohorts:

  response cohorts     AUROC of the response score for pCR / response in the treated arm
  controlled cohorts   predictive test: score x arm interaction
                         binary response -> logistic  responder ~ score + arm + score:arm
                         survival        -> Cox        OS ~ score + arm + score:arm (+ age, stage)
Null: the model's coefficient vector permuted across pathways (1,000x), which keeps the
signature's distribution but destroys its biology.

Usage: python robust/run_triad_external.py [--null 1000]
"""
import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd import survival as S  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import trial_cohorts as TC  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "triad" / "external"

# drug -> list of (cohort, test kind, treated-arm selector, control-arm selector)
TESTS = {
    "PACLITAXEL": [("GSE25066_hatzis", "auc", None, None),
                   ("GSE194040_ispy2", "auc", lambda c: c["arm"] == "Ctr", None),
                   ("GSE164458_brightness", "auc", lambda c: c["arm"].str.startswith("paclitaxel"), None)],
    "DOXORUBICIN": [("GSE25066_hatzis", "auc", None, None),
                    ("GSE194040_ispy2", "auc", lambda c: c["arm"] == "Ctr", None)],
    "CARBOPLATIN": [("GSE164458_brightness", "logit_interaction", lambda c: c["arm"].str.contains("carboplatin"),
                     lambda c: c["arm"].str.startswith("paclitaxel")),
                    ("GSE194040_ispy2", "logit_interaction", lambda c: c["arm"] == "VC", lambda c: c["arm"] == "Ctr")],
    "CISPLATIN": [("GSE14814_jbr10", "cox_interaction", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "observation"),
                  ("GSE68465_lung_adj", "cox_interaction", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "observation"),
                  ("GSE37745_lung_adj", "cox_interaction", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "observation")],
}


def fitted_models(drug, graph):
    """Fit every TRIAD variant on all TCGA patients; return {name: Series over pathways}."""
    gs = gene_sets(("REACTOME",))
    lab = R.patient_labels(drug)
    lab = lab[lab["responder"].notna()]
    M, proj = RT.patient_matrix(sorted(lab["project"].unique()))
    pats = [p for p in lab.index if p in M.index]
    lab = lab.loc[pats]
    feats = list(M.columns)
    Xp = T.center_within(M.loc[pats, feats], proj.loc[pats])
    setting = pd.get_dummies(lab["setting_proxy"]).drop(columns=["early"], errors="ignore").astype(float)
    X = np.hstack([Xp.values, setting.values])
    y = lab["responder"].astype(float).values
    k = len(feats)
    net, _ = RT.network_w(drug, gs, graph)
    w_net = np.r_[net.reindex(feats).fillna(0.25).values if net is not None else np.ones(k), 10 * np.ones(setting.shape[1])]
    w_flat = np.r_[np.ones(k), 10 * np.ones(setting.shape[1])]
    g_prior, _ = RT.gdsc_prior(drug, gs)
    o_prior, _ = RT.organoid_priors(drug, gs)

    def m_from(*priors):
        ps = [p.reindex(feats).fillna(0) for p in priors if p is not None]
        if not ps:
            return None
        v = sum(p / p.std() for p in ps) / len(ps)
        return np.r_[-RT.TAU * v.values / v.std(), np.zeros(setting.shape[1])]

    spec = {"patient_only": (w_flat, None), "network": (w_net, None), "TRIAD_all": (w_net, m_from(g_prior, o_prior))}
    if g_prior is not None:
        spec["cell_line_prior"] = (w_flat, m_from(g_prior))
    if o_prior is not None:
        spec["organoid_prior"] = (w_flat, m_from(o_prior))
    out = {}
    for name, (w, m) in spec.items():
        _, b = T.fit_prior_logistic(X, y, RT.LAM, w, m)
        out[name] = pd.Series(b[:k], index=feats)
    if g_prior is not None:
        out["frozen_cell_line"] = -g_prior.reindex(feats).fillna(0)
    if o_prior is not None:
        out["frozen_organoid"] = -o_prior.reindex(feats).fillna(0)
    return out


def cohort_scores(cid):
    expr, clin = TC.load_trial(cid, canonical=True)
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    sc = scoring.cached(f"TRIAL_{cid}_reactome_rank", scoring.rank_score, expr, {k: v for k, v in gs.items() if len(v) >= 5})
    return sc, clin


def apply(coef, sc, samples):
    feats = [f for f in coef.index if f in sc.index]
    Z = sc.loc[feats, samples].T
    Z = (Z - Z.mean()) / Z.std().replace(0, 1)
    return pd.Series(Z.fillna(0).values @ coef[feats].values, index=samples)


def statistic(kind, score, clin, treated, control):
    """Return (effect, p, extra) for one test; effect oriented so > 0.5 / > 1 = model works."""
    if kind == "auc":
        c = clin.loc[treated].dropna(subset=["responder"])
        s = score.loc[c.index]
        return roc_auc_score(c["responder"], s), None, {"n": len(c), "responders": int(c["responder"].sum())}
    idx = treated | control
    c = clin.loc[idx].copy()
    c["score"] = (score.loc[c.index] - score.loc[c.index].mean()) / score.loc[c.index].std()
    c["arm_t"] = treated.loc[c.index].astype(float)
    c["sxa"] = c["score"] * c["arm_t"]
    if kind == "logit_interaction":
        import statsmodels.api as sm
        c = c.dropna(subset=["responder"])
        fit = sm.Logit(c["responder"].astype(float), sm.add_constant(c[["score", "arm_t", "sxa"]])).fit(disp=0)
        auc_t = roc_auc_score(c.loc[c.arm_t == 1, "responder"], c.loc[c.arm_t == 1, "score"])
        return float(np.exp(fit.params["sxa"])), float(fit.pvalues["sxa"]), {"n": len(c), "AUC_treated_arm": auc_t}
    # survival: higher response score should reduce hazard more in the treated arm (HR_interaction < 1)
    c = c.rename(columns={"os_months": "months", "os_event": "event"}).dropna(subset=["months", "event"])
    covs = [x for x in ("age",) if x in c and c[x].notna().mean() > 0.7]
    r = S.cox(c, ["sxa", "score", "arm_t"] + covs)
    return float(r.loc["sxa", "HR"]), float(r.loc["sxa", "p"]), {"n": len(c), "events": int(c["event"].sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--drugs", nargs="*", default=list(TESTS))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    graph = {}

    def g():
        if "g" not in graph:
            graph["g"] = N.load_string()
        return graph["g"]

    rows = []
    rng = np.random.default_rng(5)
    for drug in a.drugs:
        models = fitted_models(drug, g)
        for cid, kind, tsel, csel in TESTS[drug]:
            sc, clin = cohort_scores(cid)
            treated = tsel(clin) if tsel else pd.Series(True, index=clin.index)
            control = csel(clin) if csel else pd.Series(False, index=clin.index)
            samples = list(clin.index[treated | control])
            for name, coef in models.items():
                score = apply(coef, sc, samples)
                eff, p, extra = statistic(kind, score, clin, treated.loc[samples], control.loc[samples])
                null = []
                for _ in range(a.null):
                    perm = pd.Series(rng.permutation(coef.values), index=coef.index)
                    null.append(statistic(kind, apply(perm, sc, samples), clin, treated.loc[samples], control.loc[samples])[0])
                null = np.array(null)
                good = (null >= eff) if kind == "auc" or kind == "logit_interaction" else (null <= eff)
                rows.append({"drug": drug, "cohort": cid, "test": kind, "model": name, "effect": eff, "p": p,
                             "null_p": float((1 + good.sum()) / (1 + len(null))), **extra})
                print(rows[-1], flush=True)
            pd.DataFrame(rows).to_csv(OUT / "triad_external_validation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
