"""Colorectal cancer / 5-fluorouracil: robust analysis across three data versions.

  paper          : the authors' committed ssGSEA matrices (organoid + TCGA HTSeq FPKM-UQ,
                   GENCODE v22) and the 2019 biotab clinical files -> exact reproduction
  current_ssgsea : same organoids, TCGA re-scored from the current GDC STAR FPKM-UQ
                   (GENCODE v36, Xena mirror) with current GDC survival
  current_rank   : current data, cohort-independent single-sample rank scores

Usage: python robust/run_coad_5fu.py [--perm 500]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from obd import RESULTS  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from obd.runner import run  # noqa: E402
from obd.study import Study  # noqa: E402

DRUG = "FLUOROURACIL"


def paper_gmt():
    """The exact Reactome gene sets the authors used for TCGA (restricted to TCGA genes)."""
    gmt = C.DATA.parent / "python/results/COAD/TCGA/reactome.gmt"
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(gmt)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    expr_org, resp = C.coad_organoids()
    prox = N.precomputed_proximity()[DRUG]
    feats = N.proximal([prox])
    drugs = C.tcga_biotab_drugs("COAD")
    treated = drugs.get(DRUG, set())

    # paper version
    org_paper = C.paper_scores("COAD", "organoid")
    paper = Study("COAD", DRUG, org_paper, resp[DRUG], C.paper_scores("COAD", "TCGA"),
                  C.tcga_biotab_clinical("COAD"), treated, feats, label="paper")

    # current GDC data
    tcga = C.tcga_star_fpkm_uq("COAD")
    cli_now = C.tcga_xena_clinical("COAD")
    pat_ss = scoring.cached("COAD_TCGA_star_fpkmuq_reactome_ssgsea", scoring.ssgsea, tcga, paper_gmt())
    current_ss = Study("COAD", DRUG, org_paper, resp[DRUG], pat_ss, cli_now, treated, feats, label="current_ssgsea")

    gs = gene_sets(("REACTOME",))
    org_rank = scoring.cached("COAD_organoid_reactome_rank", scoring.rank_score, expr_org, gs)
    pat_rank = scoring.cached("COAD_TCGA_star_fpkmuq_reactome_rank", scoring.rank_score, tcga, gs)
    current_rank = Study("COAD", DRUG, org_rank, resp[DRUG], pat_rank, cli_now, treated, feats, label="current_rank")

    rep = run({"paper": paper, "current_ssgsea": current_ss, "current_rank": current_rank},
              prox, RESULTS / "COAD_FLUOROURACIL", n_perm=args.perm, n_jobs=args.jobs)
    print((RESULTS / "COAD_FLUOROURACIL" / "REPORT.md").read_text())
    return rep


if __name__ == "__main__":
    main()
