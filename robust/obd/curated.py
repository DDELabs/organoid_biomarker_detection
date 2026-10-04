"""Loaders for the committed curated datasets (data/curated), built by network-enabled
helper sessions: GEO trial cohorts (data/curated/trials) and pre-clinical
pharmacogenomic sets (data/curated/preclinical)."""
import pandas as pd

from . import DATA
from .cohorts import to_canonical

TRIALS = DATA / "curated" / "trials"
PRECLIN = DATA / "curated" / "preclinical"


def trial_catalog():
    return pd.read_csv(TRIALS / "CATALOG_GEO.tsv", sep="\t")


def preclinical_catalog():
    return pd.read_csv(PRECLIN / "CATALOG_PRECLINICAL.tsv", sep="\t")


def load_curated_trial(cid, canonical=True):
    d = TRIALS / cid
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col=0, low_memory=False)
    if canonical:
        expr = to_canonical(expr)
    return expr, clin


def load_curated_preclinical(name, canonical=True):
    d = PRECLIN / name
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    resp = pd.read_csv(d / "response.tsv", sep="\t")
    resp = resp.groupby(["sample", "drug"])["response"].median().unstack()
    if canonical:
        expr = to_canonical(expr)
    return expr, resp
