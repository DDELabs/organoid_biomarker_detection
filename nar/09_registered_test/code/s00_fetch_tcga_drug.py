"""Rebuild the TCGA drug-record table and the TCGA-CDR survival table with the existing loader
(robust/obd/tcga_response.py), writing into this folder's _data/ (the shared data/external cache is absent)."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parents[1] / "robust"))
from obd import tcga_response as R  # noqa: E402

R.ROOT = HERE / "_data" / "tcga_response"
R.RAW = R.ROOT / "raw"
R.LONG_PATH = R.ROOT / "drug_response_long.tsv.gz"
R.CDR_PATH = R.ROOT / "survival_cdr.tsv.gz"
R.SUMMARY_PATH = R.ROOT / "SUMMARY.md"
if __name__ == "__main__":
    R.build()
    print(R.load_response_labels().shape, R.load_cdr().shape)
