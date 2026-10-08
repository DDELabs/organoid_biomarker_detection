"""BRIDGE: train on patient-derived models (PDX, organoids, cell lines), align to patient tumours,
anchor on patient response, and report per-drug pathway biomarkers with honest validation.

Per drug (TCGA RECIST labels, >= 2 cancers):
  zero-shot            pre-clinical prior alone scores patients (no patient fitting); mean within-cancer
                       AUROC, null = prior permuted across pathways (random signature, 500x)
  zero-shot per class  the same for the pdx / organoid / cell_line vectors separately
  patient_only         LOCO AUROC of the patient-only ridge logistic (s = 0)
  BRIDGE (nested)      LOCO AUROC with transfer strength s chosen inside each training fold
  s_full               transfer strength chosen on all cancers (used for the final fit)
Pathways: final model, 100 patient bootstraps -> z; 'consensus' = |z| >= 2 and same sign as the
pre-clinical prior. External: frozen final score in every curated/harvested trial cohort whose
regimen contains the drug (AUROC + random-signature null).

Outputs: robust/results/bridge/{summary.tsv, pathways.tsv.gz, domains.tsv, external.tsv, BRIDGE_RESULTS.md}
Usage: python robust/run_bridge.py [--boot 100] [--null 500] [--jobs 4]
"""
import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import bridge as B  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.atlas import DRUG_CLASS  # noqa: E402

OUT = RESULTS / "bridge"
LAM = 10.0
DRUGS = sorted(set(DRUG_CLASS) | {"CETUXIMAB", "DACARBAZINE", "TAMOXIFEN"})


def patient_data(drug):
    lab = R.patient_labels(drug)
    lab = lab[lab["responder"].notna()]
    M, proj = RT.patient_matrix(sorted(lab["project"].unique()))
    lab = lab[lab.index.isin(M.index)]
    Xp = T.center_within(M.loc[lab.index], proj.loc[lab.index])
    setting = pd.get_dummies(lab["setting_proxy"]).drop(columns=["early"], errors="ignore").astype(float)
    return lab, Xp, setting


def zero_shot(Xp, y, groups, beta, n_null, rng):
    auc, _ = T.prior_only_auc(Xp, y, groups, beta.values)
    null = [T.prior_only_auc(Xp, y, groups, rng.permutation(beta.values))[0] for _ in range(n_null)]
    null = np.array([v for v in null if np.isfinite(v)])
    return auc, float((1 + (null >= auc).sum()) / (1 + len(null))) if np.isfinite(auc) else np.nan


def run_drug(drug, n_boot, n_null):
    rng = np.random.default_rng(sum(map(ord, drug)))
    lab, Xpd, setting = patient_data(drug)
    y = lab["responder"].astype(float).values
    groups = lab["project"].values
    ok_c = [g for g in np.unique(groups) if (groups == g).sum() >= 10 and 0 < y[groups == g].sum() < (groups == g).sum()]
    if len(ok_c) < 2:
        return None
    pathways = list(Xpd.columns)
    Xp = Xpd.values
    vecs = [v for v in (B.domain_vector(n, drug, pathways, Xp) for n in B.DOMAINS) if v is not None]
    beta, cls = B.combine(vecs)
    X = np.hstack([Xp, setting.values])
    k, n_cov = len(pathways), setting.shape[1]
    row = {"drug": drug, "class": DRUG_CLASS.get(drug, "targeted/hormonal"), "patients": len(y),
           "responders": int(y.sum()), "cancers_evaluable": len(ok_c),
           "domains": ";".join(f"{v['domain']}({v['n']})" for v in vecs)}
    # patient-only
    row["patient_only_loco"] = T.loco(X, y, groups, LAM, None, None, n_score=k)[3]
    if beta is not None:
        row["zero_shot_auc"], row["zero_shot_null_p"] = zero_shot(Xp, y, groups, beta, n_null, rng)
        for c, b in cls.items():
            row[f"zero_shot_{c}"], row[f"zero_shot_{c}_p"] = zero_shot(Xp, y, groups, b, n_null // 2, rng)
        row["bridge_nested_loco"], per, chosen = B.nested_loco(X, y, groups, LAM, beta, n_cov, k)
        s_full, grid = B.choose_s(X, y, groups, LAM, beta, n_cov, k)
        row["s_chosen_per_fold"] = ";".join(f"{g}:{s:g}" for g, s in chosen.items())
        row["s_full"] = s_full
        row.update({f"loco_s{s:g}": a for s, a in grid.items()})
    else:
        s_full = 0.0
    m = B.prior_mean(beta, s_full, n_cov) if beta is not None else None
    _, bfull = T.fit_prior_logistic(X, y, LAM, None, m)
    pats = lab.index.values
    pw = pathway_inference(drug, X, y, k, pathways, beta, cls, n_boot, rng)
    pw["bridge_coef"] = bfull[:k]
    dom = pd.DataFrame([{"drug": drug, "domain": v["domain"], "class": v["class"], "n": v["n"],
                         "dropped_axes": v["dropped_axes"]} for v in vecs])
    print(drug, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in row.items()
                 if kk in ("patients", "patient_only_loco", "zero_shot_auc", "zero_shot_null_p", "bridge_nested_loco", "s_full")},
          flush=True)
    return row, pw, dom, pd.Series(bfull[:k], index=pathways), len(pats)


def pathway_inference(drug, X, y, k, pathways, beta, cls, n_boot, rng):
    """Prior-free patient bootstrap z (a prior-centred fit would shrink every bootstrap towards the same
    prior and inflate z). consensus = |patient z| >= 2 and same sign as the pre-clinical prior."""
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max():
            continue
        boots.append(T.fit_prior_logistic(X[i], y[i], LAM, None, None)[1][:k])
    boots = np.array(boots)
    z = boots.mean(0) / (boots.std(0) + 1e-12)
    pw = pd.DataFrame({"drug": drug, "pathway": pathways, "patient_boot_z": z,
                       "q": T_bh(2 * norm.sf(np.abs(z))),
                       "preclinical_beta": beta.values if beta is not None else np.nan})
    for c, b in cls.items():
        pw[f"beta_{c}"] = b.values
    pw["consensus"] = (np.abs(z) >= 2) & (np.sign(z) == np.sign(pw["preclinical_beta"]))
    return pw


def pathways_only(drugs, n_boot):
    rows = []
    for drug in drugs:
        rng = np.random.default_rng(sum(map(ord, drug)))
        lab, Xpd, setting = patient_data(drug)
        y = lab["responder"].astype(float).values
        pathways = list(Xpd.columns)
        vecs = [v for v in (B.domain_vector(n, drug, pathways, Xpd.values) for n in B.DOMAINS) if v is not None]
        beta, cls = B.combine(vecs)
        X = np.hstack([Xpd.values, setting.values])
        rows.append(pathway_inference(drug, X, y, len(pathways), pathways, beta, cls, n_boot, rng))
        print(drug, "pathways done", flush=True)
    return pd.concat(rows, ignore_index=True)


def T_bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    r = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def external(coefs, n_null, rng):
    import run_atlas_external as AE
    rows = []
    for tag, loader in AE.cohorts():
        try:
            expr, clin = loader()
        except Exception:
            continue
        if "responder" not in clin or clin["responder"].notna().sum() < 20 or clin["responder"].nunique() < 2:
            continue
        sc = AE.pathway_scores(tag.replace(":", "_"), expr)
        clin = clin.loc[[s for s in clin.index if s in sc.columns]]
        clin = clin[clin["responder"].notna()]
        drugs = AE.drugs_of(clin)
        for drug, vec in coefs.items():
            sel = drugs.str.split(";").apply(lambda x: drug in x)
            c = clin[sel]
            if len(c) < 20 or c["responder"].nunique() < 2:
                continue
            yv = c["responder"].astype(int)
            s = AE.score(vec, sc, list(c.index))
            auc = roc_auc_score(yv, s)
            p = AE.null_p(lambda zz: roc_auc_score(yv, zz), vec, sc, list(c.index), auc, n_null, rng)
            rows.append({"cohort": tag, "drug": drug, "n": len(c), "responders": int(yv.sum()), "AUC": auc, "null_p": p})
            print(rows[-1], flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=100)
    ap.add_argument("--null", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--drugs", nargs="*", default=DRUGS)
    ap.add_argument("--pathways-only", action="store_true", help="recompute pathways.tsv.gz for drugs in summary.tsv")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.pathways_only:
        old = pd.read_csv(OUT / "pathways.tsv.gz", sep="\t")
        drugs = list(pd.read_csv(OUT / "summary.tsv", sep="\t")["drug"])
        pw = pathways_only(drugs, a.boot)
        if "coef" in old:
            pw = pw.merge(old[["drug", "pathway", "coef"]].rename(columns={"coef": "bridge_coef"}), on=["drug", "pathway"], how="left")
        pw.to_csv(OUT / "pathways.tsv.gz", sep="\t", index=False)
        return
    # warm the score caches serially (avoids parallel writers)
    _, Xp0, _ = patient_data("PACLITAXEL")
    for n in B.DOMAINS:
        try:
            B.domain_scores(n, list(Xp0.columns))
        except Exception as exc:
            print(n, "unavailable:", exc)
    res = Parallel(n_jobs=a.jobs)(delayed(run_drug)(d, a.boot, a.null) for d in a.drugs)
    res = [r for r in res if r is not None]
    summ = pd.DataFrame([r[0] for r in res])
    summ.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    pw = pd.concat([r[1] for r in res], ignore_index=True)
    pw.to_csv(OUT / "pathways.tsv.gz", sep="\t", index=False)
    pd.concat([r[2] for r in res], ignore_index=True).to_csv(OUT / "domains.tsv", sep="\t", index=False)
    coefs = {r[0]["drug"]: r[3] for r in res}
    ext = external(coefs, a.null, np.random.default_rng(5))
    ext.to_csv(OUT / "external.tsv", sep="\t", index=False)
    print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
