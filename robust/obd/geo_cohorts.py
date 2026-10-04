"""External colorectal cancer validation cohorts for 5-FU biomarkers.

Files are prepared under data/external/geo/<cohort>/ (see each SOURCE.md):
  expression.tsv.gz  genes x samples, first column 'gene', log2 values
  clinical.tsv       one row per sample: sample, os_months, os_event, rfs_months, rfs_event,
                     stage (I-IV), age, sex (1 male / 0 female), chemo (1/0/NA), regimen,
                     fu_based (1/0/NA), response (CR/PR/SD/PD/NA) + cohort-specific columns
"""
import pandas as pd

from . import EXTERNAL

GEO_DIR = EXTERNAL / "geo"

# cohorts with expression + clinical prepared from verified downloads
AVAILABLE = ["GSE39582", "GSE14333", "GSE28702", "TCGA-READ"]


def load_geo_cohort(gse, canonical=False):
    """Return (expression genes x samples, clinical indexed by sample) for one prepared cohort.

    canonical=True maps symbols to the repo's canonical gene names (obd.cohorts.to_canonical).
    """
    d = GEO_DIR / gse
    if not (d / "expression.tsv.gz").exists():
        raise FileNotFoundError(f"{gse} not prepared under {d}; available: {AVAILABLE}")
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col="gene")
    clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col="sample", na_values=["NA"], keep_default_na=True)
    expr = expr.loc[:, [s for s in expr.columns if s in clin.index]]
    clin = clin.loc[expr.columns]
    if canonical:
        from .cohorts import to_canonical
        expr = to_canonical(expr)
    return expr, clin


def fu_treated(clin, presumed=False):
    """Samples treated with a 5-FU-based regimen (optionally including presumed, e.g. GSE14333)."""
    col = "fu_based_presumed" if presumed and "fu_based_presumed" in clin else "fu_based"
    return clin.index[clin[col] == 1]
