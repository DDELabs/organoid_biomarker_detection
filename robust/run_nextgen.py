"""NIT (next-generation) analysis: step-wise ablation and full robust runs.

For each study the fixes are added one at a time, so each effect is visible:
  v1            robust pipeline as before (hard proximity cut-off, raw response, OS from diagnosis)
  +landmark     survival re-based at treatment start (immortal-time bias removed)
  +prolif_cov   patient proliferation added to the adjusted Cox model
  +deconfound   organoid response residualised on proliferation + general sensitivity
  +soft_prior   all pathways, weighted by network proximity instead of a hard cut-off
  +transfer     shrink toward a pan-cancer cell-line model of the same drug (if available)
Then the full model is run through the robust runner (permutation null, jackknife,
bootstrap, scorecard) and written to robust/results/nextgen/<STUDY>/.

Usage: python robust/run_nextgen.py [--only COAD_FLUOROURACIL LIHC_SORAFENIB] [--perm 500] [--ablation-only]
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
from obd import network as N  # noqa: E402
from obd import nextgen as G  # noqa: E402
from obd import scoring  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from obd.runner import run  # noqa: E402
from obd.study import Study  # noqa: E402

STUDIES = [
    ("COAD", "FLUOROURACIL", "vdw", None),
    ("LIHC", "SORAFENIB", "licob", None),
    ("LUAD", "CISPLATIN", "gdsc", ["lung: NSCLC"]),
    ("BRCA", "DOXORUBICIN", "gdsc", ["breast"]),
    ("BLCA", "GEMCITABINE", "gdsc", ["bladder"]),
    ("PAAD", "GEMCITABINE", "gdsc", ["pancreas"]),
    ("STAD", "CISPLATIN", "gdsc", ["stomach"]),
    ("PAAD", "GEMCITABINE", "tiriac", None),      # real PDOs (Tiriac 2018)
    ("PAAD", "FLUOROURACIL", "tiriac", None),
    ("OV", "PACLITAXEL", "vias", None),            # real PDOs (Vias 2023)
]
ORGANOID_SOURCES = ("vdw", "licob", "tiriac", "vias")


def study_tag(project, drug, source):
    return f"{project}_{drug}" if source in ("vdw", "licob", "gdsc") else f"{project}_{drug}_{source}"
OUT = RESULTS / "nextgen"


def load(project, drug, source, tissues, g):
    """Everything a study needs, with pathway rank scores on a shared gene universe."""
    from run_multicancer import targets_for
    if source == "vdw":
        pre_expr, pre_resp = C.coad_organoids()
        prox = N.precomputed_proximity()[drug]          # the paper's STRING proximity
    else:
        if source == "licob":
            pre_expr, pre_resp = C.licob_organoids()
        elif source == "tiriac":
            from obd.organoid_sets import load_organoid_set
            pre_expr, pre_resp = load_organoid_set("pancreas_tiriac2018", canonical=True)
        elif source == "vias":
            from obd.preclinical_sets import load_organoid_set_extra
            pre_expr, pre_resp = load_organoid_set_extra("ovarian_vias2023")
        else:
            pre_expr, pre_resp = C.gdsc_cell_lines(tissues)
        prox = scoring.cached(f"proximity_STRINGv12_{drug}",
                              lambda: N.proximity_z(g(), targets_for(drug), gene_sets(("REACTOME",))).to_frame("z"))["z"]
    tcga = C.tcga_star_fpkm_uq(project)
    genes = pre_expr.index.intersection(tcga.index)
    gs = {k: [x for x in v if x in genes] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    tag = study_tag(project, drug, source)
    pre_sc = scoring.cached(f"{tag}_{source}_reactome_rank", scoring.rank_score, pre_expr.loc[genes], gs)
    pat_sc = scoring.cached(f"{tag}_TCGA_reactome_rank", scoring.rank_score, tcga.loc[genes], gs)
    ctx = dict(project=project, drug=drug, source=source, pre_expr=pre_expr, pre_resp=pre_resp, prox=prox,
               pre_sc=pre_sc, pat_sc=pat_sc, tcga=tcga, clinical=C.tcga_xena_clinical(project),
               treated=C.tcga_biotab_drugs(project).get(drug, set()),
               start=C.tcga_biotab_drug_start(project, drug))
    ctx["pre_prolif"] = G.proliferation_score(pre_expr)
    ctx["pat_prolif"] = G.proliferation_score(tcga).to_frame("proliferation")
    ctx["general"] = G.general_sensitivity(pre_resp, drug)
    ctx["prior"] = None
    if source in ORGANOID_SOURCES:  # transfer prior from GDSC2 (all tissues); not for cell-line studies
        try:
            from obd.preclinical_sets import load_gdsc
            cl_expr, cl_resp = load_gdsc(drugs=[drug], datasets=("GDSC2",), metric="ln_ic50")
            cg = cl_expr.index.intersection(genes)
            cl_sc = scoring.cached(f"GDSC2_all_{tag}_rank", scoring.rank_score, cl_expr.loc[cg],
                                   {k: [x for x in v if x in cg] for k, v in gs.items()})
            ctx["prior"] = G.transfer_prior(cl_sc, cl_resp[drug], list(pre_sc.index))
            ctx["prior_n_lines"] = int(cl_resp[drug].notna().sum())
        except Exception as exc:
            print("no transfer prior:", exc)
    return ctx


def make_study(ctx, landmark=False, prolif_cov=False, deconf=False, soft=False, transfer=False, label=""):
    y = ctx["pre_resp"][ctx["drug"]]
    info = {}
    if deconf:
        y, info = G.deconfound(ctx["pre_resp"], ctx["drug"], prolif=ctx["pre_prolif"], general=ctx["general"])
    if soft:
        w = G.network_weights(ctx["prox"])
        feats = list(w.index)
    else:
        w, feats = None, N.proximal([ctx["prox"]])
    st = Study(ctx["project"], ctx["drug"], ctx["pre_sc"], y, ctx["pat_sc"], ctx["clinical"], ctx["treated"], feats,
               label=label, notes={"model": ctx["source"], "deconfound": info},
               feature_weights=w if soft else None,
               prior_coef=ctx["prior"] if transfer else None,
               patient_covars=ctx["pat_prolif"] if prolif_cov else None,
               landmark_days=ctx["start"] if landmark else None)
    return st


ABLATION = [
    ("v1", {}),
    ("+landmark", dict(landmark=True)),
    ("+prolif_cov", dict(landmark=True, prolif_cov=True)),
    ("+deconfound", dict(landmark=True, prolif_cov=True, deconf=True)),
    ("+soft_prior", dict(landmark=True, prolif_cov=True, deconf=True, soft=True)),
    ("+transfer", dict(landmark=True, prolif_cov=True, deconf=True, soft=True, transfer=True)),
]


def ablation(ctx):
    rows = []
    for name, kw in ABLATION:
        if kw.get("transfer") and ctx["prior"] is None:
            continue
        st = make_study(ctx, label=name, **kw)
        w, _ = st.signature()
        ev = st.evaluate(w)
        cv = st.organoid_cv(n_boot=30) if name in ("v1", "+deconfound", "+soft_prior") else {}
        rows.append({"step": name, "features": len(st.features), "treated": ev["n_treated"], "events": ev["events"],
                     "adj_HR": ev.get("adj_HR"), "adj_p": ev.get("adj_p"), "c_index": ev["c_index"],
                     "interaction_p": ev.get("interaction_p"), "loocv_rho": cv.get("robust_loocv_spearman"),
                     "signature": ";".join(x.replace("REACTOME_", "") for x in w.index)})
        print(f"  {name:12s} feats={len(st.features):4d} HR={ev.get('adj_HR', np.nan):.2f} p={ev.get('adj_p', np.nan):.3f}")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--perm", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--ablation-only", action="store_true")
    args = ap.parse_args()
    graph = {}

    def g():
        if "g" not in graph:
            graph["g"] = N.load_string()
        return graph["g"]

    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for project, drug, source, tissues in STUDIES:
        tag = study_tag(project, drug, source)
        if args.only and tag not in args.only:
            continue
        print(f"\n===== {tag} =====")
        ctx = load(project, drug, source, tissues, g)
        ab = ablation(ctx)
        ab.to_csv(OUT / f"ablation_{tag}.tsv", sep="\t", index=False)
        if args.ablation_only:
            continue
        full = dict(landmark=True, prolif_cov=True, deconf=True, soft=True, transfer=ctx["prior"] is not None)
        st = make_study(ctx, label="nit", **full)
        alt = make_study(ctx, label="nit_nodeconf", **{**full, "deconf": False})
        rep = run({"nit": st, "nit_nodeconf": alt}, ctx["prox"], OUT / tag, n_perm=args.perm, n_jobs=args.jobs)
        e = rep["versions"]["nit"]["robust"]
        summary.append({"study": tag, "model": source, "deconfound_R2": st.notes["deconfound"].get("confounder_R2"),
                        "transfer_prior": ctx["prior"] is not None, "treated": e["n_treated"], "events": e["events"],
                        "adj_HR": e.get("adj_HR"), "adj_p": e.get("adj_p"), "c_index": e["c_index"],
                        "interaction_p": e.get("interaction_p"), "perm_p": rep["permutation"]["robust_adjHR_empirical_p"],
                        "direction": rep["verdict"]["direction"], "status": rep["verdict"]["signature_status"],
                        "signature": ";".join(rep["versions"]["nit"]["signature"])})
        path = OUT / "nextgen_summary.tsv"
        old = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame(columns=["study"])
        new = pd.DataFrame(summary)
        pd.concat([old[~old.study.isin(new.study)], new]).to_csv(path, sep="\t", index=False)
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
