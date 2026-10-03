"""Run the robust pipeline on other cancer types.

Pre-clinical models:
  * LICOB liver cancer organoids (Ji et al. 2023, iLICOB GitHub release): real PDOs
  * GDSC cell lines (Garnett et al. 2012) as a stand-in where no open organoid
    pharmacogenomic set is reachable. Results are labelled 'cell-line proxy'.
Patients: current GDC STAR FPKM-UQ + survival (Xena mirror), BCR biotab drug tables.
Network: drug-target to Reactome proximity recomputed on STRING v12 (score >= 700),
with random-walk propagation reported alongside.

Usage: python robust/run_multicancer.py [--only LIHC_SORAFENIB] [--perm 500]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402

from obd import RESULTS  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd.reference import drug_targets, drugbank_targets, gene_sets  # noqa: E402
from obd.runner import run  # noqa: E402
from obd.study import Study  # noqa: E402

# (TCGA project, drug, pre-clinical source, GDSC tissue filter)
STUDIES = [
    ("LIHC", "SORAFENIB", "licob", None),
    ("LUAD", "CISPLATIN", "gdsc", ["lung: NSCLC"]),
    ("BRCA", "DOXORUBICIN", "gdsc", ["breast"]),
    ("BLCA", "GEMCITABINE", "gdsc", ["bladder"]),
    ("PAAD", "GEMCITABINE", "gdsc", ["pancreas"]),
    ("STAD", "CISPLATIN", "gdsc", ["stomach"]),
]


def targets_for(drug):
    """Curated organoid-drug targets (as in the paper) united with DrugBank targets."""
    return sorted(set(drug_targets().get(drug, [])) | set(drugbank_targets().get(drug, [])))


def run_one(project, drug, source, tissues, g, n_perm, jobs):
    tag = f"{project}_{drug}"
    print(f"\n===== {tag} ({source}) =====")
    if source == "licob":
        pre_expr, pre_resp = C.licob_organoids()
        model = "patient-derived organoids (LICOB)"
    else:
        pre_expr, pre_resp = C.gdsc_cell_lines(tissues)
        model = "cell-line proxy (GDSC 2012)"
    tcga = C.tcga_star_fpkm_uq(project)
    clinical = C.tcga_xena_clinical(project)
    treated = C.tcga_biotab_drugs(project).get(drug, set())

    # gene universe shared by both cohorts, as in the original run_ssGSEA.py
    genes = pre_expr.index.intersection(tcga.index)
    gs = {k: [x for x in v if x in genes] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}

    targets = targets_for(drug)
    prox = scoring.cached(f"proximity_STRINGv12_{drug}", lambda: N.proximity_z(g, targets, gene_sets(("REACTOME",))).to_frame("z"))["z"]
    rwr = scoring.cached(f"rwr_STRINGv12_{drug}", lambda: N.rwr_z(g, targets, gene_sets(("REACTOME",)), n_random=100).to_frame("z"))["z"]
    feats = N.proximal([prox])

    pre_ss = scoring.cached(f"{tag}_{source}_reactome_ssgsea", scoring.ssgsea, pre_expr.loc[genes], gs)
    pat_ss = scoring.cached(f"{tag}_TCGA_reactome_ssgsea", scoring.ssgsea, tcga.loc[genes], gs)
    pre_rk = scoring.cached(f"{tag}_{source}_reactome_rank", scoring.rank_score, pre_expr.loc[genes], gs)
    pat_rk = scoring.cached(f"{tag}_TCGA_reactome_rank", scoring.rank_score, tcga.loc[genes], gs)

    notes = {"model": model, "targets": targets, "network": "STRING v12 >= 700"}
    versions = {
        "current_ssgsea": Study(project, drug, pre_ss, pre_resp[drug], pat_ss, clinical, treated, feats, "current_ssgsea", notes),
        "current_rank": Study(project, drug, pre_rk, pre_resp[drug], pat_rk, clinical, treated, feats, "current_rank", notes),
    }
    out = RESULTS / tag
    rep = run(versions, prox, out, n_perm=n_perm, n_jobs=jobs, extra_z={"rwr_z": rwr})
    return {"study": tag, "model": model, "targets": ",".join(targets), **versions["current_ssgsea"].summary(),
            "events_treated": rep["verdict"]["events_in_treated"],
            "adj_HR": rep["versions"]["current_ssgsea"]["robust"].get("adj_HR"),
            "adj_p": rep["versions"]["current_ssgsea"]["robust"].get("adj_p"),
            "perm_p": rep["permutation"]["robust_adjHR_empirical_p"],
            "baseline_min_p": rep["permutation"]["baseline_min_p_observed"],
            "baseline_perm_p": rep["permutation"]["baseline_min_p_empirical_p"],
            "interaction_p": rep["versions"]["current_ssgsea"]["robust"].get("interaction_p"),
            "status": rep["verdict"]["signature_status"],
            "robust_pathways": ";".join(rep["robust_pathways"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--perm", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    g = N.load_string()
    rows = []
    for project, drug, source, tissues in STUDIES:
        if args.only and f"{project}_{drug}" not in args.only:
            continue
        rows.append(run_one(project, drug, source, tissues, g, args.perm, args.jobs))
        pd.DataFrame(rows).to_csv(RESULTS / "multicancer_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
