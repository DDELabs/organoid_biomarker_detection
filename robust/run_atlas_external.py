"""ATLAS external validation: frozen pathway scores in independent trial cohorts.

The primary ATLAS model (shared + class + drug layers, GDSC + network priors) is fitted
once on all TCGA records and frozen. For every external cohort:

  regimen score   mean total effect (shared + class + drug) over the regimen's drugs
                  that ATLAS models (PLATINUM -> platinum class; others ignored)
  shared score    'overall response' layer only
  specific score  class + drug layers only (drug-specific part, shared removed)

Tests
  response cohorts       AUROC (higher score -> response) per cohort / arm
  randomised designs     drug-specificity: logit(response) ~ specific_score * arm, or
                         Cox OS ~ specific_score * arm (chemo vs observation)
  null                   effect vectors permuted across pathways (500x)

Usage: python robust/run_atlas_external.py [--null 500]
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
import statsmodels.api as sm  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import run_atlas as RA  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import atlas as A  # noqa: E402
from obd import curated as CU  # noqa: E402
from obd import network as N  # noqa: E402
from obd import survival as S  # noqa: E402
from obd import trial_cohorts as TC  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "atlas" / "external"
ALIAS = {"PLATINUM": None, "EPIRUBICIN": "EPIRUBICIN", "XELODA": "CAPECITABINE"}


def frozen_effects():
    cache = RESULTS / "atlas" / "full_fit_effects.tsv.gz"
    if cache.exists():
        return pd.read_csv(cache, sep="\t", index_col=0)
    graph = {}
    rec, X, cov, w, pathways, _ = RA.build_records()
    pm, pw = RA.priors(pathways, lambda: graph.setdefault("g", N.load_string()))
    des = RA.variant("atlas", pathways)
    _, b = A.fit(des, X, rec["drug"].values, cov, rec["responder"].astype(float).values, RA.LAM, w, pm, pw)
    eff = des.effects(b)
    eff.to_csv(cache, sep="\t")
    return eff


def regimen_vectors(eff, drugs):
    """(total, shared, specific) effect vectors for a regimen (list of drug names)."""
    comps, spec = [], []
    for d in drugs:
        d = ALIAS.get(d, d)
        if d is None:   # generic platinum -> class layer only
            c = eff["class:platinum"]
            comps.append(eff["shared"] + c)
            spec.append(c)
        elif d in A.DRUG_CLASS and f"drug:{d}" in eff:
            c = eff[f"class:{A.DRUG_CLASS[d]}"] + eff[f"drug:{d}"]
            comps.append(eff["shared"] + c)
            spec.append(c)
    if not comps:
        return None, eff["shared"], None
    return sum(comps) / len(comps), eff["shared"], sum(spec) / len(spec)


def score(vec, sc, samples):
    feats = [f for f in vec.index if f in sc.index]
    Z = sc.loc[feats, samples].T
    Z = (Z - Z.mean()) / Z.std().replace(0, 1)
    s = pd.Series(Z.fillna(0).values @ vec[feats].values, index=samples)
    return (s - s.mean()) / (s.std() or 1)


def pathway_scores(tag, expr):
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    return scoring.cached(f"ATLASEXT_{tag}_reactome_rank", scoring.rank_score, expr, {k: v for k, v in gs.items() if len(v) >= 5})


def null_p(stat_fn, vec, sc, samples, obs, n, rng, greater=True):
    vals = []
    for _ in range(n):
        v = pd.Series(rng.permutation(vec.values), index=vec.index)
        r = stat_fn(score(v, sc, samples))
        if r is not None and np.isfinite(r):
            vals.append(r)
    vals = np.array(vals)
    return float((1 + ((vals >= obs) if greater else (vals <= obs)).sum()) / (1 + len(vals)))


def cohorts():
    """(tag, loader) for curated GEO cohorts and previously harvested trial cohorts."""
    out = [(f"geo:{c}", (lambda c=c: CU.load_curated_trial(c))) for c in CU.trial_catalog()["cohort_id"]]
    try:
        lst = TC.list_trials(usable_only=True)
        ids = list(lst["cohort_id"]) if hasattr(lst, "columns") else list(lst)
        for c in ids:
            out.append((f"trial:{c}", (lambda c=c: TC.load_trial(c, canonical=True))))
    except Exception as exc:
        print("harvested trials unavailable:", exc)
    return out


def drugs_of(clin):
    col = "drugs" if "drugs" in clin else "drug"
    return clin[col].fillna("").astype(str).str.upper().str.replace(" ", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=500)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    eff = frozen_effects()
    rng = np.random.default_rng(3)
    rows = []
    for tag, loader in cohorts():
        try:
            expr, clin = loader()
        except Exception as exc:
            print(tag, "load failed", exc)
            continue
        if "responder" not in clin or clin["responder"].notna().sum() < 20 or clin["responder"].nunique() < 2:
            continue
        sc = pathway_scores(tag.replace(":", "_"), expr)
        clin = clin.loc[[s for s in clin.index if s in sc.columns]]
        arms = clin["arm"].fillna("all") if "arm" in clin else pd.Series("all", index=clin.index)
        for arm in ["ALL"] + sorted(arms.unique()):
            c = clin if arm == "ALL" else clin[arms == arm]
            c = c[c["responder"].notna()]
            if len(c) < 20 or c["responder"].nunique() < 2:
                continue
            regimen = sorted({x for s in drugs_of(c) for x in s.replace(",", ";").split(";") if x})
            total, shared, spec = regimen_vectors(eff, regimen)
            y = c["responder"].astype(int)
            for kind, vec in (("regimen_total", total), ("shared", shared)):
                if vec is None:
                    continue
                s = score(vec, sc, list(c.index))
                auc = roc_auc_score(y, s)
                p = null_p(lambda z: roc_auc_score(y, z), vec, sc, list(c.index), auc, a.null, rng)
                rows.append({"cohort": tag, "arm": arm, "score": kind, "drugs_modelled": ";".join(
                    d for d in regimen if ALIAS.get(d, d) is None or ALIAS.get(d, d) in A.DRUG_CLASS),
                    "n": len(c), "responders": int(y.sum()), "AUC": auc, "null_p": p})
                print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(OUT / "atlas_external_auc.tsv", sep="\t", index=False)

    # ---------------- drug-specificity tests in randomised / controlled designs
    spec_rows = []
    designs = [
        ("geo:GSE20271", "paclitaxel added (T/FAC vs FAC)", lambda c: c["arm"] == "T/FAC", lambda c: c["arm"] == "FAC", ["PACLITAXEL"], "logit"),
        ("geo:GSE41998", "paclitaxel vs ixabepilone after AC", lambda c: c["arm"] == "AC->paclitaxel", lambda c: c["arm"] == "AC->ixabepilone", ["PACLITAXEL"], "logit"),
        ("trial:GSE164458_brightness", "carboplatin added", lambda c: c["arm"].str.contains("carboplatin"), lambda c: c["arm"].str.startswith("paclitaxel"), ["CARBOPLATIN"], "logit"),
        ("trial:GSE194040_ispy2", "veliparib+carboplatin added", lambda c: c["arm"] == "VC", lambda c: c["arm"] == "Ctr", ["CARBOPLATIN"], "logit"),
        ("geo:GSE42127", "adjuvant platinum vs observation (OS)", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "observation", ["PLATINUM"], "cox"),
        ("trial:GSE14814_jbr10", "cisplatin/vinorelbine vs observation (OS)", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "observation", ["CISPLATIN", "VINORELBINE"], "cox"),
        ("geo:GSE103479", "adjuvant 5-FU vs surgery (OS)", lambda c: c["arm"] == "adjuvant chemo", lambda c: c["arm"] == "surgery only", ["FLUOROURACIL"], "cox"),
    ]
    loaders = dict(cohorts())
    for tag, label, tsel, csel, regimen, kind in designs:
        if tag not in loaders:
            continue
        expr, clin = loaders[tag]()
        sc = pathway_scores(tag.replace(":", "_"), expr)
        clin = clin.loc[[s for s in clin.index if s in sc.columns]]
        t, ctl = tsel(clin).fillna(False), csel(clin).fillna(False)
        c = clin[t | ctl].copy()
        c["arm_t"] = t[t | ctl].astype(float)
        total, shared, spec = regimen_vectors(eff, regimen)
        for sname, vec in (("specific", spec), ("total", total), ("shared", shared)):
            s = score(vec, sc, list(c.index))

            def stat(z, c=c):
                d = c.assign(score=z.loc[c.index].values, sxa=z.loc[c.index].values * c["arm_t"].values)
                if kind == "logit":
                    d = d.dropna(subset=["responder"])
                    f = sm.Logit(d["responder"].astype(float), sm.add_constant(d[["score", "arm_t", "sxa"]])).fit(disp=0)
                    return float(f.params["sxa"]), float(f.pvalues["sxa"])
                d = d.rename(columns={"os_months": "months", "os_event": "event"}).dropna(subset=["months", "event"])
                r = S.cox(d, ["sxa", "score", "arm_t"])
                return -float(np.log(r.loc["sxa", "HR"])), float(r.loc["sxa", "p"])  # >0 = benefit grows with score

            eff_obs, p_obs = stat(s)
            p_null = null_p(lambda z: stat(z)[0], vec, sc, list(c.index), eff_obs, a.null, rng)
            spec_rows.append({"cohort": tag, "design": label, "score": sname, "n": len(c),
                              "interaction_log_effect": eff_obs, "wald_p": p_obs, "null_p": p_null,
                              "direction": "expected" if eff_obs > 0 else "opposite"})
            print(spec_rows[-1], flush=True)
    pd.DataFrame(spec_rows).to_csv(OUT / "atlas_drug_specificity.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
