"""Standardised patient-derived organoid pharmacogenomic sets.

Each set lives in data/external/organoids/<name>/ as
  expression.tsv.gz : genes x organoids, first column 'gene' (symbols)
  response.tsv      : long table sample, drug, response, metric, direction
  SOURCE.md         : provenance, sha256, units and caveats

    expr, resp = load_organoid_set("pancreas_tiriac2018")

returns expression (genes x samples) and response (samples x drugs, lower = more
sensitive), restricted by default to organoids present in both.
"""
import pandas as pd

from . import EXTERNAL

ROOT = EXTERNAL / "organoids"

# name -> (cancer, metric, expression units)
SETS = {
    "pancreas_tiriac2018": ("PAAD", "AUC", "log2(FPKM+1), GDC HTSeq, protein-coding"),
    "liver_licob": ("LIHC", "AUC", "log2 scale as distributed by iLICOB"),
}
AVAILABLE = sorted(SETS)


def load_organoid_set(name, overlap_only=True, canonical=False, metric=None):
    """Return (expression genes x samples, response samples x drugs) for a set.

    overlap_only : keep only organoids with both expression and response.
    canonical    : map symbols to the canonical UniProt names used by the
                   rest of the pipeline (cohorts.to_canonical).
    metric       : restrict to one metric if a set ever carries several.
    """
    if name not in SETS:
        raise KeyError(f"unknown organoid set {name!r}; available: {AVAILABLE}")
    d = ROOT / name
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    expr.index.name = "gene"
    long = pd.read_csv(d / "response.tsv", sep="\t")
    if metric is not None:
        long = long[long["metric"] == metric]
    long["sample"] = long["sample"].astype(str).str.strip()
    resp = long.groupby(["sample", "drug"])["response"].median().unstack()
    if canonical:
        from .cohorts import to_canonical
        expr = to_canonical(expr)
    if overlap_only:
        common = [s for s in expr.columns if s in resp.index]
        expr, resp = expr[common], resp.loc[common]
    return expr, resp


def describe():
    """One row per set: organoids with expression+response, drugs, metric."""
    rows = []
    for name in AVAILABLE:
        e, r = load_organoid_set(name)
        rows.append(dict(name=name, cancer=SETS[name][0], organoids=e.shape[1], genes=e.shape[0],
                         drugs=r.shape[1], metric=SETS[name][1], units=SETS[name][2]))
    return pd.DataFrame(rows).set_index("name")


if __name__ == "__main__":
    print(describe().to_string())
