"""External validation of frozen 5-FU signatures in independent colorectal cohorts.

Signatures are learned on the 19 van de Wetering organoids only, then applied unchanged:
  paper_top7     Kong et al. Ridge top-7 (paper's ssGSEA data)
  robust_v1      robust stability-selected signature (paper data, scorecard 5/5)
  v1_rank        robust signature re-learned on rank scores (current pipeline)
  nit_soft       next-gen: soft network prior, raw response
  nit_deconf     next-gen: soft network prior + deconfounded response
Cohorts (robust/obd/geo_cohorts.py): GSE39582 (adjuvant, RFS/OS), GSE14333 (RFS),
TCGA-READ (OS), GSE28702 (FOLFOX response).

Usage: python robust/run_external_validation.py [--null 1000]
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
from obd import geo_cohorts as GEO  # noqa: E402
from obd import models as M  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd import validation as V  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from obd.study import Study  # noqa: E402

OUT = RESULTS / "external_validation"
DRUG = "FLUOROURACIL"


def signatures():
    sigs = {}
    _, resp = C.coad_organoids()
    feats = N.proximal([N.precomputed_proximity()[DRUG]])
    paper = Study("COAD", DRUG, C.paper_scores("COAD", "organoid"), resp[DRUG], C.paper_scores("COAD", "TCGA"),
                  C.tcga_biotab_clinical("COAD"), C.tcga_biotab_drugs("COAD")[DRUG], feats)
    c = M.baseline_coefficients(paper.X, paper.y, "Ridge")
    top = np.argsort(-np.abs(c))[:7]
    sigs["paper_top7"] = pd.Series(c[top], index=[paper.features[i] for i in top])
    rep = json.load(open(RESULTS / "COAD_FLUOROURACIL" / "report.json"))
    sigs["robust_v1"] = pd.Series(rep["versions"]["paper"]["signature"])

    import run_nextgen as NG
    ctx = NG.load("COAD", DRUG, "vdw", None, lambda: None)
    for name, kw in (("v1_rank", {}), ("nit_soft", dict(soft=True)), ("nit_deconf", dict(soft=True, deconf=True))):
        w, _ = NG.make_study(ctx, **kw).signature()
        sigs[name] = w
    return sigs, ctx


def cohort_scores(gse):
    expr, clin = GEO.load_geo_cohort(gse, canonical=True)
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    sc = scoring.cached(f"external_{gse}_reactome_rank", scoring.rank_score, expr, gs)
    return sc, clin


def designs(gse, clin):
    """Patient groups and endpoints per cohort."""
    if gse == "GSE39582":
        base = clin[clin["stage"].isin(["II", "III"])]
        return base, set(base.index[base["fu_based"] == 1]), set(base.index[base["chemo"] == 0]), [("rfs_months", "rfs_event"), ("os_months", "os_event")]
    if gse == "GSE14333":
        tr = set(GEO.fu_treated(clin, presumed=True))
        return clin, tr, set(clin.index[clin["chemo"] == 0]), [("rfs_months", "rfs_event")]
    if gse == "TCGA-READ":
        un = clin.index[clin["chemo"] == 0]
        return clin, set(GEO.fu_treated(clin)), set(un) if len(un) >= 20 else None, [("os_months", "os_event")]
    return clin, None, None, []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=1000)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sigs, _ = signatures()
    pd.concat({k: v for k, v in sigs.items()}, names=["signature", "pathway"]).rename("weight").to_csv(OUT / "frozen_signatures.tsv", sep="\t")
    rows, nulls = [], []
    for gse in GEO.AVAILABLE:
        sc, clin = cohort_scores(gse)
        base, tr, un, endpoints = designs(gse, clin)
        for name, w in sigs.items():
            score = V.signature_score(sc, w)
            if score is None:
                continue
            for tcol, ecol in endpoints:
                res = V.survival_tests(score.reindex(base.index), base, tr, un, tcol, ecol)
                rows.append({"cohort": gse, "endpoint": tcol.split("_")[0].upper(), "signature": name,
                             "n_pathways_used": int(sum(p in sc.index for p in w.index)), **res})
            if gse == "GSE28702":
                r = V.response_test(score, clin["responder"])
                rows.append({"cohort": gse, "endpoint": "FOLFOX response", "signature": name, **r})
        # random-signature null on the primary design
        if gse == "GSE39582" and args.null:
            for tcol, ecol in endpoints:
                ep = tcol.split("_")[0].upper()

                def stat_tr(s, tcol=tcol, ecol=ecol):
                    r = V.survival_tests(s.reindex(base.index), base, tr, None, tcol, ecol)
                    return abs(np.log(r["HR_treated"])) if "HR_treated" in r else None

                def stat_int(s, tcol=tcol, ecol=ecol):
                    r = V.survival_tests(s.reindex(base.index), base, tr, un, tcol, ecol)
                    return abs(np.log(r["interaction_HR"])) if "interaction_HR" in r else None

                for name, w in sigs.items():
                    k = int(sum(p in sc.index for p in w.index))
                    obs = [r for r in rows if r["cohort"] == gse and r["signature"] == name and r["endpoint"] == ep][0]
                    n_tr = V.random_null(sc, k, stat_tr, n=args.null)
                    n_in = V.random_null(sc, k, stat_int, n=args.null, seed=12)
                    nulls.append({"endpoint": ep, "signature": name, "k": k,
                                  "random_null_p_treated": float((1 + (n_tr >= abs(np.log(obs["HR_treated"]))).sum()) / (1 + len(n_tr))),
                                  "random_null_p_interaction": float((1 + (n_in >= abs(np.log(obs["interaction_HR"]))).sum()) / (1 + len(n_in)))})
                    print("null", nulls[-1], flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "external_validation.tsv", sep="\t", index=False)
    pd.DataFrame(nulls).to_csv(OUT / "random_signature_null_GSE39582.tsv", sep="\t", index=False)
    cols = ["cohort", "endpoint", "signature", "n_treated", "events_treated", "HR_treated", "HR_low", "HR_high", "p_treated",
            "HR_untreated", "interaction_HR", "interaction_p", "AUC", "p_one_sided"]
    print(res[[c for c in cols if c in res]].round(3).to_string())


if __name__ == "__main__":
    main()
