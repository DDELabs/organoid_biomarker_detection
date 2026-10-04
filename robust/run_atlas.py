"""ATLAS: pathway determinants of drug response across all drugs and cancers (TCGA).

Pre-registered choices
  records    first RECIST response per patient x drug, 21 drugs with >= 30 labelled patients
             (leucovorin excluded as a modulator); weight 1/(drugs in the regimen)
  features   673 Reactome rank scores, centred and scaled within cancer type
  covariates treatment setting (early / post-progression / late) and proliferation
             (REACTOME_CELL_CYCLE score); fitted, never part of pathway effects
  model      multi-task ridge logistic: shared + class + drug layers, lambda = 10
  variants   separate  : drug layer only (one model per drug, no borrowing)
             atlas     : shared + class + drug layers
             atlas_prior: atlas + GDSC prior means + network penalty weights on drug layers
  held-out   leave-one-cancer-out; metric = mean AUROC over (drug, cancer) pairs with
             >= 10 records and both outcomes
  null       labels permuted within (drug, cancer), same LOCO pipeline
  pathways   effects from the full-data fit; stability = 100 patient bootstraps;
             z = mean / sd, BH-FDR across all (layer, pathway) tests

Usage: python robust/run_atlas.py [--perm 20] [--boot 100]
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
from scipy.stats import norm  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import atlas as A  # noqa: E402
from obd import network as N  # noqa: E402
from obd import survival as S  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "atlas"
DRUGS = [d for d in A.DRUG_CLASS]
LAM = 10.0
PROLIF = "REACTOME_CELL_CYCLE"


def build_records():
    recs = []
    for d in DRUGS:
        lab = R.patient_labels(d)
        lab = lab[lab["responder"].notna()].copy()
        lab["drug"] = d
        recs.append(lab.reset_index().rename(columns={"index": "patient"}))
    rec = pd.concat(recs, ignore_index=True)
    if "patient" not in rec:
        rec = rec.rename(columns={rec.columns[0]: "patient"})
    M, proj = RT.patient_matrix(sorted(rec["project"].unique()))
    rec = rec[rec["patient"].isin(M.index)].reset_index(drop=True)
    pats = rec["patient"].unique()
    Xc = T.center_within(M.loc[pats], proj.loc[pats])
    X = Xc.loc[rec["patient"]].values
    setting = pd.get_dummies(rec["setting_proxy"]).drop(columns=["early"], errors="ignore").astype(float)
    cov = np.column_stack([setting.values, Xc.loc[rec["patient"], PROLIF].values])
    w = 1.0 / rec["n_drugs_in_regimen"].clip(lower=1).fillna(1).values
    return rec, X, cov, w, list(Xc.columns), Xc


def priors(pathways, graph):
    gs = gene_sets(("REACTOME",))
    pm, pw = {}, {}
    for d in DRUGS:
        g, _ = RT.gdsc_prior(d, gs)
        if g is not None:
            v = g.reindex(pathways).fillna(0)
            pm[d] = -RT.TAU * v / v.std()
        try:
            net, _ = RT.network_w(d, gs, graph)
        except Exception:
            net = None
        pw[d] = net
    return pm, pw


def variant(name, pathways):
    if name == "separate":
        return A.AtlasDesign(pathways, DRUGS, scale_shared=0.0, scale_class=0.0, scale_drug=1.0)
    return A.AtlasDesign(pathways, DRUGS)


def loco(rec, X, cov, w, pathways, vname, y=None, pm=None, pw=None):
    y = rec["responder"].astype(float).values if y is None else y
    score = np.full(len(rec), np.nan)
    from joblib import Parallel, delayed

    def fold(c):
        te = (rec["project"] == c).values
        des = variant(vname, pathways)
        _, b = A.fit(des, X[~te], rec["drug"].values[~te], cov[~te], y[~te], LAM, w[~te], pm, pw)
        return te, A.predict(des, b, X[te], rec["drug"].values[te], cov[te])

    for te, s in Parallel(n_jobs=4)(delayed(fold)(c) for c in rec["project"].unique()):
        score[te] = s
    rows = []
    for (d, c), g in rec.assign(s=score, y=y).groupby(["drug", "project"]):
        if len(g) >= 10 and g["y"].nunique() == 2:
            rows.append({"drug": d, "project": c, "n": len(g), "auc": roc_auc_score(g["y"], g["s"])})
    pairs = pd.DataFrame(rows)
    return pairs, (pairs["auc"].mean() if len(pairs) else np.nan), score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=20)
    ap.add_argument("--boot", type=int, default=100)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    graph = {}

    def g():
        if "g" not in graph:
            graph["g"] = N.load_string()
        return graph["g"]

    rec, X, cov, w, pathways, Xc = build_records()
    print("records", len(rec), "patients", rec.patient.nunique(), "cancers", rec.project.nunique(), flush=True)
    pm, pw = priors(pathways, g)

    # 1. held-out performance of the three variants
    summary, per_pair = [], []
    for vname, use_prior in (("separate", False), ("atlas", False), ("atlas_prior", True)):
        pairs, mean_auc, _ = loco(rec, X, cov, w, pathways, "separate" if vname == "separate" else "atlas",
                                  pm=pm if use_prior else None, pw=pw if use_prior else None)
        per_pair.append(pairs.assign(model=vname))
        per_drug = pairs.groupby("drug")["auc"].mean()
        summary.append({"model": vname, "mean_pair_auc": mean_auc, "pairs": len(pairs),
                        **{f"auc_{d}": per_drug.get(d, np.nan) for d in DRUGS}})
        print(vname, "mean AUC over drug x cancer pairs: %.3f (%d pairs)" % (mean_auc, len(pairs)), flush=True)
    pd.concat(per_pair).to_csv(OUT / "loco_pairs.tsv", sep="\t", index=False)
    summ = pd.DataFrame(summary)

    # 2. permutation null for the primary (atlas_prior) model
    rng = np.random.default_rng(1)
    null = []
    for i in range(a.perm):
        yp = rec["responder"].astype(float).values.copy()
        for _, idx in rec.groupby(["drug", "project"]).groups.items():
            idx = np.asarray(list(idx))
            yp[idx] = rng.permutation(yp[idx])
        null.append(loco(rec, X, cov, w, pathways, "atlas", y=yp, pm=pm, pw=pw)[1])
        print("perm", i, round(null[-1], 3), flush=True)
    obs = summ.loc[summ.model == "atlas_prior", "mean_pair_auc"].iloc[0]
    summ["perm_p_primary"] = (1 + sum(v >= obs for v in null)) / (1 + len(null))
    summ["null_mean"] = np.mean(null) if null else np.nan
    summ.to_csv(OUT / "model_comparison.tsv", sep="\t", index=False)

    # 3. full fit + bootstrap stability of pathway effects (primary model)
    y = rec["responder"].astype(float).values
    des = variant("atlas", pathways)
    _, b = A.fit(des, X, rec["drug"].values, cov, y, LAM, w, pm, pw)
    eff = des.effects(b)
    boots = []
    pats = rec["patient"].unique()
    for i in range(a.boot):
        pick = rng.choice(pats, len(pats), replace=True)
        cnt = pd.Series(pick).value_counts()
        ww = w * rec["patient"].map(cnt).fillna(0).values
        keep = ww > 0
        _, bb = A.fit(des, X[keep], rec["drug"].values[keep], cov[keep], y[keep], LAM, ww[keep], pm, pw)
        boots.append(des.effects(bb))
        if i % 10 == 0:
            print("boot", i, flush=True)
    stack = np.stack([e.values for e in boots])
    mean, sd = stack.mean(0), stack.std(0) + 1e-12
    z = pd.DataFrame(mean / sd, index=eff.index, columns=eff.columns)
    sign_cons = pd.DataFrame((np.sign(stack) == np.sign(mean)).mean(0), index=eff.index, columns=eff.columns)
    p = pd.DataFrame(2 * norm.sf(np.abs(z.values)), index=z.index, columns=z.columns)
    flat = p.stack()
    order = flat.sort_values()
    bh = (order * len(order) / np.arange(1, len(order) + 1))[::-1].cummin()[::-1].clip(upper=1)
    q = bh.reindex(flat.index).unstack()
    long = pd.concat({"effect": eff, "boot_mean": pd.DataFrame(mean, index=eff.index, columns=eff.columns),
                      "z": z, "p": p, "q": q, "sign_consistency": sign_cons}, axis=1).stack(level=1)
    long.index.names = ["pathway", "layer"]
    long.reset_index().to_csv(OUT / "pathway_effects.tsv.gz", sep="\t", index=False)

    # per-drug total effects (shared + class + drug), with bootstrap z
    tot = {}
    for d in DRUGS:
        cols = [list(eff.columns).index(c) for c in ("shared", f"class:{A.DRUG_CLASS[d]}", f"drug:{d}")]
        t = stack[:, :, cols].sum(2)
        tot[d] = pd.DataFrame({"effect": t.mean(0), "z": t.mean(0) / (t.std(0) + 1e-12)}, index=eff.index)
    pd.concat(tot, names=["drug", "pathway"]).to_csv(OUT / "drug_total_effects.tsv.gz", sep="\t")

    # 4. overall response score (shared layer) vs survival in all treated patients (PFI from CDR)
    cdr = R.load_cdr()
    shared = pd.Series(Xc.values @ eff["shared"].reindex(Xc.columns).fillna(0).values, index=Xc.index)
    surv = cdr.reindex(shared.index)
    tcol = next(c for c in ("PFI.time", "PFI_time", "pfi_time") if c in surv.columns)
    ecol = next(c for c in ("PFI", "pfi") if c in surv.columns)
    df = pd.DataFrame({"score": (shared - shared.mean()) / shared.std(), "months": surv[tcol] / 30.44,
                       "event": surv[ecol]}).dropna()
    df = df.join(rec.drop_duplicates("patient").set_index("patient")["project"]).dropna()
    from statsmodels.duration.hazard_regression import PHReg
    fit = PHReg(df["months"].values, df[["score"]].values, status=df["event"].values,
                strata=df["project"].values).fit(disp=False)
    overall = {"n": len(df), "events": int(df.event.sum()), "HR_per_SD_higher_response_score": float(np.exp(fit.params[0])),
               "p": float(fit.pvalues[0])}
    json.dump({"records": len(rec), "patients": int(rec.patient.nunique()), "cancers": int(rec.project.nunique()),
               "drugs": DRUGS, "overall_response_score_PFI": overall,
               "null": null}, open(OUT / "meta.json", "w"), indent=1)
    print(summ[["model", "mean_pair_auc", "pairs", "perm_p_primary", "null_mean"]].round(3).to_string(index=False))
    print("overall response score vs PFI:", overall)


if __name__ == "__main__":
    main()
