"""Shared data access for the BLCA / OV / LUAD / other organoid-only analyses.

Reuses the loaders of robust/obd (Python 3 port of Kong et al. 2020) but points every download and
cache at this folder (nar/02_organoid_only/BLCA_OV_LUAD_OTHER/_data), so nothing outside nar/ is written.
"""
import hashlib
import sys
import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
BASE = HERE.parent                                   # nar/02_organoid_only/BLCA_OV_LUAD_OTHER
PARENT = BASE.parent                                 # nar/02_organoid_only
REPO = BASE.parents[2]
DATA = BASE / "_data"
RAW = DATA / "raw"
CACHE = DATA / "cache"
CACHE.mkdir(parents=True, exist_ok=True)
MOESM3 = DATA / "moesm3"
sys.path.insert(0, str(REPO / "robust"))

import obd  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd.reference import gene_sets, canonical_symbol, common_drug_name, drug_targets, drugbank_targets  # noqa: E402,F401

# redirect the robust loaders' download folder to ours (files already downloaded by download_raw.sh)
C.EXTERNAL = RAW
obd.EXTERNAL = RAW
PRECLIN = REPO / "data" / "curated" / "preclinical"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------ TCGA clinical / treatment
@lru_cache(None)
def treatment_long(project):
    """Long table (patient, drug) from the BCR biotab clinical_drug file (same parser as robust.obd.tcga_response)."""
    f = CACHE / f"treat_{project}.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t")
    rec = R._read_biotab(RAW / f"nationwidechildrens.org_clinical_drug_{project.lower()}.txt", project)
    rows = []
    for p, raw in zip(rec["patient"], rec["drug_raw"]):
        for name, ok in R.normalize_drug(raw):
            rows.append((p, name))
    df = pd.DataFrame(rows, columns=["patient", "drug"]).drop_duplicates()
    # patients listed in the drug file at all (any record, even unresolved) count as "had a drug"
    anyrec = pd.DataFrame({"patient": rec["patient"].unique(), "drug": "__ANY_RECORD__"})
    df = pd.concat([df, anyrec], ignore_index=True)
    df.to_csv(f, sep="\t", index=False)
    return df


def treated(project, drug):
    t = treatment_long(project)
    return set(t.loc[t["drug"] == drug, "patient"])


def any_drug(project):
    t = treatment_long(project)
    return set(t["patient"])


@lru_cache(None)
def survival_current(project):
    """Current GDC survival (Xena GDC hub, OS in months), primary tumour patients."""
    s = C.tcga_xena_clinical(project)
    return s.rename(columns={"months": "T", "event": "E"})


@lru_cache(None)
def survival_biotab(project):
    """Paper-era survival from the BCR biotab clinical_patient file (as python/utilities/parse_patient_clinical_data.py)."""
    df = pd.read_csv(RAW / f"nationwidechildrens.org_clinical_patient_{project.lower()}.txt", sep="\t", dtype=str, skiprows=[1, 2])
    rows = {}
    for _, r in df.iterrows():
        pat, st = r["bcr_patient_barcode"], str(r["vital_status"]).upper()
        lc = r.get("last_contact_days_to", r.get("days_to_last_followup"))
        dd = r.get("death_days_to", r.get("days_to_death"))
        days = lc if st == "ALIVE" else dd if st == "DEAD" else "["
        if "TCGA" not in str(pat) or "[" in str(days) or pd.isna(days):
            continue
        rows[pat] = {"T": float(days) / (365. / 12.), "E": int(st == "DEAD")}
    return pd.DataFrame.from_dict(rows, orient="index")


@lru_cache(None)
def cdr():
    return R._build_cdr(RAW / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv")


# ------------------------------------------------------------------ expression and pathway scores
REACTOME = None


def reactome():
    global REACTOME
    if REACTOME is None:
        REACTOME = dict(gene_sets(("REACTOME",)))
    return REACTOME


def tcga_expr(project):
    return C.tcga_star_fpkm_uq(project)


def ssgsea_cached(name, expr, gs, threads=2):
    """gseapy ssGSEA (raw ES, rank-normalised samples), cached. Same call as robust.obd.scoring.ssgsea."""
    path = CACHE / f"ssgsea_{name}.tsv.gz"
    if path.exists():
        return pd.read_csv(path, sep="\t", index_col=0)
    import gseapy as gp
    res = gp.ssgsea(data=expr, gene_sets=gs, sample_norm_method="rank", permutation_num=0, no_plot=True,
                    outdir=None, min_size=2, threads=threads, verbose=False)
    out = res.res2d.pivot(index="Term", columns="Name", values="ES").astype(float)
    out.to_csv(path, sep="\t")
    return out


def paired_scores(tag, oexpr, texpr):
    """Score organoid and patient samples on the shared gene universe (as python/run_ssGSEA.py)."""
    genes = oexpr.index.intersection(texpr.index)
    gs = {k: [g for g in v if g in genes] for k, v in reactome().items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    return ssgsea_cached(f"{tag}_model", oexpr.loc[genes], gs), ssgsea_cached(f"{tag}_tcga", texpr.loc[genes], gs)


def load_preclinical(name, value="response"):
    """Curated set from data/curated/preclinical (expression genes x samples canonical, long response)."""
    d = PRECLIN / name
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    expr.columns = expr.columns.astype(str)
    expr = C.to_canonical(expr)
    r = pd.read_csv(d / "response.tsv", sep="\t", dtype={"sample": str})
    return expr, r


def vias_ovarian():
    """Vias 2023 HGSOC organoids, rebuilt from the raw repo files exactly as robust/data_builders/build_ov.py."""
    f_e, f_r = CACHE / "vias_expression.tsv.gz", CACHE / "vias_response.tsv"
    if not f_e.exists():
        V = RAW / "vias"
        jb2org = {'JBLAB19902': '119148', 'JBLAB19904': '54327', 'JBLAB19905': '118976', 'JBLAB19906': '119178',
                  'JBLAB19907': '119127', 'JBLAB19917': '23868', 'JBLAB19921': '119058', 'JBLAB19925': '54288',
                  'JBLAB19937': '54059', 'JBLAB19938': '151723', 'JBLAB19939': '151761', 'JBLAB19940': '54276',
                  'JBLAB19941': '32077', 'JBLAB19942': '151773'}
        t = pd.read_csv(V / "20191218_ViasM_BJ_orgaBrs_tpm.csv")
        t = t[t.gene_name.notna()]
        m = t[list(jb2org)].copy()
        m.columns = [jb2org[c] for c in m.columns]
        m.index = t.gene_name
        m = np.log2(m.groupby(level=0).sum() + 1)
        m = m[(m > 0).any(axis=1)]
        m.to_csv(f_e, sep="\t")
        a = pd.read_csv(V / "20200419-AUC-organoids.csv")
        a["sample"] = a["id"].astype(str)
        alias = {'AZD2014': 'VISTUSERTIB', 'AZD5363': 'CAPIVASERTIB', 'AZD2281': 'OLAPARIB', 'AZD6738': 'CERALASERTIB',
                 'AZD1775': 'ADAVOSERTIB', 'Elescamol': 'ELESCLOMOL', 'APR-246': 'EPRENETAPOPT'}
        rows = []
        for c in [c for c in a.columns if c.startswith("auc_ll5.")]:
            d = c.split(".", 1)[1]
            name = common_drug_name(alias.get(d, d))
            for s, v in zip(a["sample"], a[c]):
                if pd.notna(v):
                    rows.append((s, name, float(v), "AUC_LL5_viability"))
        pd.DataFrame(rows, columns=["sample", "drug", "response", "metric"]).to_csv(f_r, sep="\t", index=False)
    expr = C.to_canonical(pd.read_csv(f_e, sep="\t", index_col=0))
    expr.columns = expr.columns.astype(str)
    r = pd.read_csv(f_r, sep="\t", dtype={"sample": str})
    return expr, r


# ------------------------------------------------------------------ network proximity
def paper_proximity():
    """Paper's own z-scores (MOESM3 sheet 'Figure 2a,d' = data/coad_blca_organoid_drugs_zscore_result_reactome.txt)."""
    return N.precomputed_proximity()


@lru_cache(None)
def paper_network():
    """The paper's STRING>700 network (MOESM3 sheet 'Figure 2a-f'), on canonical symbols."""
    import networkx as nx
    e = pd.read_csv(MOESM3 / "Figure_2a_f.tsv", sep="\t")
    a = e["Gene_Symbol_1"].astype(str).map(lambda s: canonical_symbol(s) or s)
    b = e["Gene_Symbol_2"].astype(str).map(lambda s: canonical_symbol(s) or s)
    g = nx.Graph()
    g.add_edges_from(zip(a, b))
    g.remove_edges_from(nx.selfloop_edges(g))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def proximity(drug, targets=None):
    """z-score per Reactome pathway. Paper table when the drug is in it; else closest-distance proximity
    (robust.obd.network.proximity_z, 1000 degree-matched randomisations) on the paper's network."""
    P = paper_proximity()
    if drug in P.columns and targets is None:
        return P[drug].dropna(), "paper table (MOESM3 Fig 2a,d)"
    f = CACHE / f"proximity_papernet_{drug}.tsv"
    if f.exists():
        s = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
        return s, f"computed on paper network, targets {open(str(f) + '.targets').read().strip()}"
    if targets is None:
        targets = drug_targets().get(drug) or drugbank_targets().get(drug)
    if not targets:
        raise KeyError(f"no targets for {drug}")
    pw = {k: v for k, v in reactome().items() if k in P.index}
    z = N.proximity_z(paper_network(), list(targets), pw, n_random=1000)
    z.to_frame("z").to_csv(f, sep="\t")
    open(str(f) + ".targets", "w").write(",".join(targets))
    return z, f"computed on paper network, targets {','.join(targets)}"
