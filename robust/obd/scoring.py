"""Pathway activity scores per sample."""
import hashlib

import numpy as np
import pandas as pd

from . import RESULTS


def _restrict(expr, gene_sets, min_size):
    genes = set(expr.index)
    gs = {k: [g for g in v if g in genes] for k, v in gene_sets.items()}
    return {k: v for k, v in gs.items() if len(v) >= min_size}


def ssgsea(expr, gene_sets, min_size=2, threads=4):
    """ssGSEA raw enrichment scores (pathways x samples), as in the original pipeline.

    Raw ES is used instead of gseapy's 'normalized' ES: the latter divides by a
    single global range over the cohort, a constant that downstream z-scoring
    removes anyway but that makes files depend on cohort composition.
    """
    import gseapy as gp

    gs = _restrict(expr, gene_sets, min_size)
    res = gp.ssgsea(data=expr, gene_sets=gs, sample_norm_method="rank", permutation_num=0,
                    no_plot=True, outdir=None, min_size=min_size, threads=threads, verbose=False)
    return res.res2d.pivot(index="Term", columns="Name", values="ES").astype(float)


def rank_score(expr, gene_sets, min_size=5):
    """Single-sample rank score (singscore-style, Foroutan 2018), pathways x samples.

    Mean within-sample percentile rank of the set's genes, centred at 0. Depends
    only on the sample itself, so it is invariant to cohort composition and to
    any monotone per-sample transform (FPKM, FPKM-UQ, TPM, log or not).
    """
    gs = _restrict(expr, gene_sets, min_size)
    pct = expr.rank(axis=0).values / (len(expr) + 1)
    idx = {g: i for i, g in enumerate(expr.index)}
    out = np.vstack([pct[[idx[g] for g in genes]].mean(axis=0) - 0.5 for genes in gs.values()])
    return pd.DataFrame(out, index=list(gs), columns=expr.columns)


def cached(name, fn, *args, **kwargs):
    """Cache a scores matrix under robust/results/cache keyed by name."""
    path = RESULTS / "cache" / f"{name}.tsv.gz"
    if path.exists():
        return pd.read_csv(path, sep="\t", index_col=0)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = fn(*args, **kwargs)
    out.to_csv(path, sep="\t")
    return out


def fingerprint(*objs):
    h = hashlib.md5()
    for o in objs:
        h.update(repr(o).encode())
    return h.hexdigest()[:8]
