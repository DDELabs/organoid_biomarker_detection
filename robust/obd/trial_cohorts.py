"""Clinical-trial transcriptomic cohorts with treatment labels, and genetic-interaction networks.

Cohorts live under data/external/trials/<cohort_id>/ (see each SOURCE.md and trials/CATALOG.tsv):
  expression.tsv.gz  genes x samples, first column 'gene' (HUGO symbols), log scale
                     (log2 / log2(x+1) / log2(CPM+1); ENLIGHT microarray sets may hold probe-SUMMED log2
                     values, so standardise genes before pooling)
  clinical.tsv       one row per sample: sample, drug, arm, response (CR/PR/SD/PD, R/NR or pCR/RD),
                     responder (1/0/NA), plus os_/pfs_/rfs_/drfs_ survival columns where available
  SOURCE.md          URLs, sha256, processing, PCA-vs-label QC

CATALOG columns: cohort_id, cancer, drugs, setting, n, responders, n_labelled, has_control_arm, arms,
endpoint, platform, accession, source, qc_label, qc_pc1_auc..qc_pc3_auc, qc_flag (|PC1 AUC-0.5|>0.35),
qc_warn_pc23 (same threshold on PC2/PC3), usable.

Genetic-interaction networks live under data/external/gi_networks/<name>.tsv.gz with columns
gene_a, gene_b, type (SL, ...), score, source (+ network-specific extras); see gi_networks/SOURCE.md.
"""
import pandas as pd

from . import EXTERNAL

TRIALS_DIR = EXTERNAL / "trials"
GI_DIR = EXTERNAL / "gi_networks"


def _read_catalog():
    f = TRIALS_DIR / "CATALOG.tsv"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep="\t").set_index("cohort_id", drop=False)


CATALOG = _read_catalog()


def list_trials(usable_only=False, control_only=False):
    """Cohort ids from CATALOG, optionally restricted to usable / with a control or comparator arm."""
    c = CATALOG
    if usable_only:
        c = c[c["usable"].astype(bool)]
    if control_only:
        c = c[c["has_control_arm"].astype(bool)]
    return list(c.index)


def load_trial(cohort_id, canonical=False, labelled_only=False):
    """Return (expression genes x samples, clinical indexed by sample) for one cohort.

    canonical=True maps symbols to the repo's canonical gene names (obd.cohorts.to_canonical).
    labelled_only=True drops samples with a missing `responder` label.
    """
    d = TRIALS_DIR / cohort_id
    if not (d / "expression.tsv.gz").exists():
        raise FileNotFoundError(f"{cohort_id} not under {TRIALS_DIR}; available: {list(CATALOG.index)}")
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col="gene")
    clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col="sample", na_values=["NA"], keep_default_na=True,
                       dtype={"sample": str}, low_memory=False)
    clin.index = clin.index.astype(str)
    expr.columns = expr.columns.astype(str)
    if labelled_only and "responder" in clin:
        clin = clin[clin["responder"].notna()]
    expr = expr.loc[:, [s for s in clin.index if s in expr.columns]]
    clin = clin.loc[expr.columns]
    if canonical:
        from .cohorts import to_canonical
        expr = to_canonical(expr)
    return expr, clin


def list_gi_networks():
    return sorted(p.name[:-len(".tsv.gz")] for p in GI_DIR.glob("*.tsv.gz"))


def load_gi_network(name, types=None, genes=None):
    """Return a genetic-interaction edge table (gene_a, gene_b, type, score, source, ...).

    types: optional iterable of interaction types to keep (e.g. {"SL"}).
    genes: optional iterable; keep edges with at least one endpoint in it.
    """
    f = GI_DIR / f"{name}.tsv.gz"
    if not f.exists():
        raise FileNotFoundError(f"{name} not under {GI_DIR}; available: {list_gi_networks()}")
    df = pd.read_csv(f, sep="\t", na_values=["NA"], low_memory=False)
    if types is not None:
        df = df[df["type"].isin(set(types))]
    if genes is not None:
        g = set(genes)
        df = df[df["gene_a"].isin(g) | df["gene_b"].isin(g)]
    return df.reset_index(drop=True)


def gi_partners(name, gene, types=None):
    """Partners of `gene` in network `name` (either edge direction)."""
    df = load_gi_network(name, types=types, genes=[gene])
    return sorted(set(df.loc[df.gene_a == gene, "gene_b"]) | set(df.loc[df.gene_b == gene, "gene_a"]))
