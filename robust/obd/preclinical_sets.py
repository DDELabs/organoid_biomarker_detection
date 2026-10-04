"""Extra pre-clinical pharmacogenomic sets bundled under data/external.

Organoid sets live in data/external/organoids/<name>/ with
  expression.tsv.gz : genes x samples, first column 'gene'
  response.tsv      : long table (sample, drug, response, metric); lower = more sensitive
  SOURCE.md         : provenance, checksums and caveats

Cell-line pre-training data (GDSC1/GDSC2 fitted dose response + CTRPv2 AUC, GDSC
RMA basal expression keyed by COSMIC ID) lives in data/external/gdsc/.

Both loaders return the same shapes as robust.obd.cohorts:
  expression : DataFrame genes (canonical symbols) x samples
  response   : DataFrame samples x drugs (lower = more sensitive)
"""
from functools import lru_cache

import numpy as np
import pandas as pd

from . import EXTERNAL
from .cohorts import to_canonical
from .reference import common_drug_name

ORGANOIDS = EXTERNAL / "organoids"
GDSC = EXTERNAL / "gdsc"


def _available():
    if not ORGANOIDS.exists():
        return []
    return sorted(p.name for p in ORGANOIDS.iterdir()
                  if (p / "expression.tsv.gz").exists() and (p / "response.tsv").exists())


# Organoid sets present on disk (this agent added 'ovarian_vias2023'; other sets
# dropped into data/external/organoids/ in the same layout are picked up too).
AVAILABLE = _available()


def load_organoid_set_extra(name, canonical=True, metric=None):
    """Load one organoid set from data/external/organoids/<name>.

    canonical: map gene symbols to the pipeline's canonical UniProt names.
    metric:    keep only rows with this metric (sets may mix metrics).
    Returns (expression genes x samples, response samples x drugs).
    Replicate measurements of a sample/drug pair are median-aggregated.
    """
    d = ORGANOIDS / name
    if not (d / "expression.tsv.gz").exists():
        raise KeyError(f"unknown organoid set {name!r}; available: {_available()}")
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    expr.columns = expr.columns.astype(str)
    if canonical:
        expr = to_canonical(expr)
    long = pd.read_csv(d / "response.tsv", sep="\t", dtype={"sample": str})
    if metric is not None:
        long = long[long["metric"] == metric]
    long["drug"] = long["drug"].map(common_drug_name)
    resp = long.groupby(["sample", "drug"])["response"].median().unstack()
    keep = [s for s in resp.index if s in expr.columns]
    return expr[keep], resp.loc[keep]


@lru_cache(None)
def _gdsc_expression(canonical=True):
    expr = pd.read_csv(GDSC / "expression.tsv.gz", sep="\t", index_col=0)
    expr.columns = expr.columns.astype(str)
    return to_canonical(expr) if canonical else expr


@lru_cache(None)
def gdsc_response_long():
    """All cell-line response rows (GDSC1, GDSC2, CTRPv2) as a long DataFrame."""
    path = GDSC / "response.tsv.gz"
    if not path.exists():
        path = GDSC / "response.tsv"
    return pd.read_csv(path, sep="\t", low_memory=False,
                       dtype={"sanger_model_id": str, "gdsc_cancer_type": str})


def gdsc_cell_line_map():
    """COSMIC ID <-> DepMap ACH <-> Sanger model ID <-> name <-> tissue table."""
    return pd.read_csv(GDSC / "cell_line_map.tsv", sep="\t", index_col=0)


def load_gdsc(drugs=None, tissues=None, datasets=("GDSC2", "GDSC1"), metric="ln_ic50",
              canonical=True, min_lines=1):
    """Cell lines with GDSC RMA expression and drug response.

    drugs:    iterable of drug names (any synonym; resolved with common_drug_name) or None for all.
    tissues:  iterable of case-insensitive substrings matched against GDSC tissue descriptors,
              TCGA label, GDSC cancer type and DepMap lineage (e.g. ['glioma'], ['GBM', 'LGG'],
              ['stomach'], ['breast']); None = all lines.
    datasets: response sources in order of preference; for each (line, drug) the first dataset
              with a value is used. 'CTRPv2' only has metric='auc'.
    metric:   'ln_ic50' (natural log of IC50 in uM) or 'auc' (fraction of viability AUC for GDSC;
              CTRPv2 AUC is on a different, unnormalised scale). Lower = more sensitive for both.
    Returns (expression genes x COSMIC IDs (str), response COSMIC IDs (str) x drugs).
    """
    if metric not in ("ln_ic50", "auc"):
        raise ValueError("metric must be 'ln_ic50' or 'auc'")
    if isinstance(datasets, str):
        datasets = (datasets,)
    long = gdsc_response_long()
    long = long[long["dataset"].isin(datasets) & long[metric].notna()]
    if drugs is not None:
        want = {common_drug_name(d) for d in drugs}
        long = long[long["drug"].isin(want)]
    if tissues is not None:
        cols = ["tissue", "tissue_sub", "tcga_label", "gdsc_cancer_type", "depmap_lineage"]
        text = long[cols].fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
        hit = np.zeros(len(long), bool)
        for t in tissues:
            hit |= text.str.contains(str(t).lower(), regex=False).to_numpy()
        long = long[hit]
    long = long.assign(cosmic_id=long["cosmic_id"].astype(int).astype(str))
    # duplicate drug IDs inside one dataset (e.g. two cisplatin screens) -> median
    per = long.groupby(["dataset", "cosmic_id", "drug"])[metric].median()
    resp = None
    for ds in datasets:
        if ds not in per.index.get_level_values(0):
            continue
        m = per.loc[ds].unstack()
        resp = m if resp is None else resp.combine_first(m)
    if resp is None:
        resp = pd.DataFrame()
    expr = _gdsc_expression(canonical)
    keep = [c for c in resp.index if c in expr.columns]
    resp = resp.loc[keep]
    resp = resp.loc[:, resp.notna().sum() >= min_lines]
    return expr[keep], resp
