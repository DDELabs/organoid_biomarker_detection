"""Per-patient Reactome rank scores for every TCGA project (TRIAD features).

Downloads each project's current GDC STAR FPKM-UQ from the Xena mirror, keeps one primary
tumour per patient, computes single-sample rank scores (cached in robust/results/cache)
and, with --cleanup, deletes the large expression file afterwards to save disk.

Usage: python robust/prepare_pancancer_scores.py [--projects COAD READ ...] [--cleanup]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from obd import EXTERNAL  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import scoring  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

PROJECTS = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP",
            "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
            "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


def project_scores(project, cleanup=False):
    def compute():
        expr = C.tcga_star_fpkm_uq(project)
        gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
        return scoring.rank_score(expr, {k: v for k, v in gs.items() if len(v) >= 5})

    sc = scoring.cached(f"TCGA_{project}_reactome_rank_all", compute)
    if cleanup:
        f = EXTERNAL / f"TCGA-{project}.star_fpkm-uq.tsv.gz"
        if f.exists():
            f.unlink()
    return sc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--projects", nargs="*", default=PROJECTS)
    ap.add_argument("--cleanup", action="store_true")
    a = ap.parse_args()
    for p in a.projects:
        try:
            s = project_scores(p, a.cleanup)
            print(p, s.shape, flush=True)
        except Exception as exc:
            print(p, "FAILED", exc, flush=True)
