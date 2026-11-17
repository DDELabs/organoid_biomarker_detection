"""Uniform access to the four drug screens (GDSC1, GDSC2, CTRPv2, PRISM) with expression-based
Reactome rank scores, lineage labels, Sanger model IDs (for CRISPR / mutations) and colon / bladder flags.

Every response matrix is samples x drugs, LOWER = MORE SENSITIVE.
"""
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd

from common import load


@dataclass
class Screen:
    name: str
    resp: pd.DataFrame            # samples x drugs
    scores: pd.DataFrame          # pathways x samples (Reactome rank scores)
    lineage: pd.Series            # sample -> lineage label
    sidm: pd.Series               # sample -> Sanger model id (may be NaN)
    metric: str
    colon: list = field(default_factory=list)
    bladder: list = field(default_factory=list)

    def common(self, drug, subset=None):
        """Samples with both expression scores and a response to `drug`."""
        if drug not in self.resp:
            return []
        s = self.resp[drug].dropna().index.intersection(self.scores.columns)
        if subset is not None:
            s = s.intersection(pd.Index(subset))
        return list(s)


@lru_cache(None)
def models():
    return load("models")


@lru_cache(None)
def gdsc_matrix(dataset, metric="ln_ic50"):
    g = load("gdsc_long")
    g = g[(g["dataset"] == dataset) & g[metric].notna() & g["sanger_model_id"].notna()]
    return g.groupby(["sanger_model_id", "drug"])[metric].median().unstack()


@lru_cache(None)
def screen(name, metric=None):
    m = models()
    if name in ("GDSC1", "GDSC2"):
        metric = metric or "ln_ic50"
        resp = gdsc_matrix(name, metric)
        sc = load("reactome_cmp")
        ids = resp.index.intersection(sc.columns)
        lin = m["cancer_type"].reindex(resp.index)
        sidm = pd.Series(resp.index, index=resp.index)
        colon = [i for i in resp.index if lin.get(i) == "Colorectal Carcinoma"]
        blad = [i for i in resp.index if lin.get(i) == "Bladder Carcinoma"]
        return Screen(name, resp, sc[ids], lin, sidm, metric, colon, blad)
    tag = name.lower()
    long = load(f"{tag}_long")
    samp = load(f"{tag}_samples").set_index("sample")
    resp = long.groupby(["sample", "drug"])["response"].median().unstack()
    sc = load(f"reactome_{tag}")
    lin = samp["lineage"].reindex(resp.index)
    sidm = samp["Sanger_Model_ID"].reindex(resp.index)
    colon = [i for i in resp.index if lin.get(i) == "colorectal"]
    blad = [i for i in resp.index if lin.get(i) == "urinary_tract"]
    return Screen(name, resp, sc[resp.index.intersection(sc.columns)], lin, sidm, "AUC", colon, blad)


SCREENS = ("GDSC1", "GDSC2", "CTRP", "PRISM")


def all_screens():
    return {n: screen(n) for n in SCREENS}


@lru_cache(None)
def crispr():
    return load("crispr_lfc")


@lru_cache(None)
def mutations():
    return load("mutations")


def frozen_signatures():
    from common import ROOT
    f = ROOT / "robust" / "results" / "external_validation" / "frozen_signatures.tsv"
    d = pd.read_csv(f, sep="\t")
    return {k: v.set_index("pathway")["weight"] for k, v in d.groupby("signature")}


def organoid_37():
    """The 37 network-proximal pathways and their organoid weights (robust scorecard, paper data)."""
    from common import ROOT
    t = pd.read_csv(ROOT / "robust" / "results" / "COAD_FLUOROURACIL" / "pathway_scorecard.tsv", sep="\t", index_col=0)
    return t["mean_coef"]
