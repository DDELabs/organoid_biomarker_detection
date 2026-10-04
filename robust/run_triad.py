"""TRIAD pan-cancer analysis: patient-anchored drug-response biomarkers with priors.

Pre-registered choices (fixed before looking at results):
  label     first recorded RECIST response per patient for the drug (CR/PR = 1, SD/PD = 0)
  adjust    treatment setting (early / post-progression / late) as covariates that are
            fitted but NOT part of the biomarker score
  features  673 Reactome rank scores, centred and scaled within cancer type
  lambda    10; prior strength tau = 0.05 log-odds per SD of the prior coefficient vector
  primary   model with all available priors (network + cell line + organoid)
  metric    leave-one-cancer-out (LOCO) AUROC averaged over held-out cancers
  null      100 label permutations within cancer type (same LOCO pipeline)

Usage: python robust/run_triad.py [--drugs FLUOROURACIL CISPLATIN ...] [--perm 100]
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from obd import RESULTS  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import nextgen as G  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "triad"
DRUGS = ["CISPLATIN", "FLUOROURACIL", "CARBOPLATIN", "PACLITAXEL", "GEMCITABINE", "DOXORUBICIN",
         "DOCETAXEL", "TEMOZOLOMIDE", "OXALIPLATIN", "ETOPOSIDE", "CAPECITABINE"]
LAM, TAU = 10.0, 0.05


def patient_matrix(projects):
    blocks = []
    for p in projects:
        f = RESULTS / "cache" / f"TCGA_{p}_reactome_rank_all.tsv.gz"
        if f.exists():
            blocks.append(pd.read_csv(f, sep="\t", index_col=0).T.assign(_project=p))
    M = pd.concat(blocks)
    proj = M.pop("_project")
    M = M.loc[:, M.notna().all()]
    return M, proj


def ridge_prior(scores, response, min_n=12, alpha=10.0):
    """Standardised ridge coefficients (pathways) predicting resistance in a pre-clinical set."""
    from sklearn.linear_model import Ridge
    y = response.dropna()
    s = [x for x in y.index if x in scores.columns]
    if len(s) < min_n:
        return None
    X = scores[s].T
    X = (X - X.mean()) / X.std().replace(0, 1)
    yz = (y[s] - y[s].mean()) / y[s].std()
    return pd.Series(Ridge(alpha=alpha).fit(X.fillna(0).values, yz.values).coef_, index=scores.index)


def gdsc_prior(drug, gs):
    from obd.preclinical_sets import _gdsc_expression, load_gdsc
    try:
        _, resp = load_gdsc(drugs=[drug])
    except Exception:
        return None, 0
    if drug not in resp:
        return None, 0
    sc = scoring.cached("GDSC_all_lines_reactome_rank",
                        lambda: scoring.rank_score(_gdsc_expression(), {k: v for k, v in gs.items()}))
    return ridge_prior(sc, resp[drug], min_n=50), int(resp[drug].notna().sum())


def organoid_priors(drug, gs):
    """Ridge priors from every organoid set that screened the drug (mean of available)."""
    from obd.organoid_sets import load_organoid_set
    from obd.preclinical_sets import load_organoid_set_extra
    sets = {"vdw_coad": lambda: C.coad_organoids(),
            "tiriac_paad": lambda: load_organoid_set("pancreas_tiriac2018", canonical=True),
            "licob_lihc": lambda: C.licob_organoids(),
            "vias_ov": lambda: load_organoid_set_extra("ovarian_vias2023")}
    priors, used = [], []
    for name, fn in sets.items():
        try:
            expr, resp = fn()
        except Exception:
            continue
        if drug not in resp:
            continue
        sc = scoring.cached(f"ORG_{name}_reactome_rank",
                            lambda: scoring.rank_score(expr, {k: [g for g in v if g in expr.index] for k, v in gs.items()}))
        pr = ridge_prior(sc, resp[drug])
        if pr is not None:
            priors.append(pr / pr.std())
            used.append(f"{name}(n={int(resp[drug].notna().sum())})")
    if not priors:
        return None, []
    return pd.concat(priors, axis=1).mean(axis=1), used


def network_w(drug, gs, graph):
    from run_multicancer import targets_for
    t = targets_for(drug)
    if not t:
        return None, []
    z = scoring.cached(f"proximity_STRINGv12_{drug}", lambda: N.proximity_z(graph(), t, gs).to_frame("z"))["z"]
    w = 1.0 / (1.0 + np.exp((z - (-1.2816)) / 0.5))
    return 0.25 + 0.75 * w, t   # never fully exclude a pathway


def run_drug(drug, n_perm, graph, primary_only=False):
    gs = gene_sets(("REACTOME",))
    lab = R.patient_labels(drug)
    lab = lab[lab["responder"].notna()]
    M, proj = patient_matrix(sorted(lab["project"].unique()))
    pats = [p for p in lab.index if p in M.index]
    lab = lab.loc[pats]
    feats = list(M.columns)
    Xp = T.center_within(M.loc[pats, feats], proj.loc[pats])
    setting = pd.get_dummies(lab["setting_proxy"]).drop(columns=["early"], errors="ignore").astype(float)
    X = np.hstack([Xp.values, setting.values])
    k = len(feats)
    y = lab["responder"].astype(float).values
    groups = lab["project"].values

    net, targets = network_w(drug, gs, graph)
    w_net = np.r_[net.reindex(feats).fillna(0.25).values if net is not None else np.ones(k), 10 * np.ones(setting.shape[1])]
    w_flat = np.r_[np.ones(k), 10 * np.ones(setting.shape[1])]
    g_prior, n_lines = gdsc_prior(drug, gs)
    o_prior, org_used = organoid_priors(drug, gs)

    def m_from(*priors):
        ps = [p.reindex(feats).fillna(0) for p in priors if p is not None]
        if not ps:
            return None
        v = sum(p / p.std() for p in ps) / len(ps)
        return np.r_[-TAU * v.values / v.std(), np.zeros(setting.shape[1])]  # resistance -> lower response

    models = {"patient_only": (w_flat, None), "network": (w_net, None)}
    if g_prior is not None:
        models["cell_line_prior"] = (w_flat, m_from(g_prior))
    if o_prior is not None:
        models["organoid_prior"] = (w_flat, m_from(o_prior))
    models["TRIAD_all"] = (w_net, m_from(g_prior, o_prior))

    rows, per_cancer = [], {}
    for name, (w, m) in models.items():
        pred, per, pooled, mean_auc = T.loco(X, y, groups, LAM, w, m, n_score=k)
        _, per_k, kfold_auc = T.kfold(X, y, groups, LAM, w, m, n_score=k)
        rows.append({"model": name, "mean_within_cancer_AUC": mean_auc, "kfold_within_cancer_AUC": kfold_auc,
                     "pooled_AUC": pooled, "cancers": len(per)})
        per_cancer[name] = per
    # frozen pre-clinical models applied directly to patients (no patient fitting)
    for name, pr in (("frozen_cell_line", g_prior), ("frozen_organoid", o_prior)):
        if pr is not None:
            auc, _ = T.prior_only_auc(Xp.values, y, groups, -pr.reindex(feats).fillna(0).values)
            rows.append({"model": name, "mean_within_cancer_AUC": auc, "kfold_within_cancer_AUC": auc})
    res = pd.DataFrame(rows)
    # permutation null for the pre-registered primary model and the patient-only model
    for name in (("TRIAD_all",) if primary_only else ("TRIAD_all", "patient_only")):
        w, m = models[name]
        for scheme, col in ((("loco", "mean_within_cancer_AUC"),) if primary_only else
                            (("loco", "mean_within_cancer_AUC"), ("kfold", "kfold_within_cancer_AUC"))):
            null = T.permutation_null(X, y, groups, LAM, w, m, n=n_perm, n_score=k, scheme=scheme)
            obs = res.loc[res.model == name, col].iloc[0]
            res.loc[res.model == name, f"perm_p_{scheme}"] = (1 + (null >= obs).sum()) / (1 + len(null))
            res.loc[res.model == name, f"null_mean_{scheme}"] = null.mean()
    # final primary model on all patients: top pathways
    w, m = models["TRIAD_all"]
    b0, b = T.fit_prior_logistic(X, y, LAM, w, m)
    coef = pd.Series(b[:k], index=feats).sort_values()
    out = OUT / drug if not primary_only else OUT / "confirmatory" / drug
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "models.tsv", sep="\t", index=False)
    pd.DataFrame({m: {c: v["auc"] for c, v in pc.items()} for m, pc in per_cancer.items()}).to_csv(out / "per_cancer_auc.tsv", sep="\t")
    coef.to_csv(out / "triad_coefficients.tsv", sep="\t", header=["coef_log_odds_response"])
    meta = {"drug": drug, "patients": len(y), "responders": int(y.sum()), "cancers": sorted(set(groups)),
            "targets": targets, "gdsc_lines": n_lines, "organoid_sets": org_used,
            "top_response_pathways": list(coef.index[-10:][::-1]), "top_resistance_pathways": list(coef.index[:10])}
    json.dump(meta, open(out / "meta.json", "w"), indent=1)
    print(drug, meta["patients"], "pts |", res.round(3).to_string(index=False), flush=True)
    return res.assign(drug=drug, patients=len(y), responders=int(y.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", nargs="*", default=DRUGS)
    ap.add_argument("--perm", type=int, default=100)
    ap.add_argument("--primary-only", action="store_true", help="confirmatory run: primary model, LOCO null only")
    a = ap.parse_args()
    graph = {}

    def g():
        if "g" not in graph:
            graph["g"] = N.load_string()
        return graph["g"]

    OUT.mkdir(parents=True, exist_ok=True)
    allres = []
    for d in a.drugs:
        try:
            allres.append(run_drug(d, a.perm, g, a.primary_only))
        except Exception as exc:
            print(d, "FAILED", repr(exc), flush=True)
        if allres and not a.primary_only:
            path = OUT / "triad_summary.tsv"
            new = pd.concat(allres)
            old = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame(columns=["drug"])
            pd.concat([old[~old.drug.isin(new.drug)], new]).to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main()
