"""Shared loaders for the organoid-only (Kong et al. 2020) strategy on GI lab-model sets.

Reuses the project's existing code in robust/obd (TCGA loader, drug-name normaliser, BCR biotab reader,
ssGSEA wrapper, paper proximity table). All downloads and caches are redirected to
nar/02_organoid_only/_data so nothing outside this folder is written.
"""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                                   # nar/02_organoid_only
REPO = ROOT.parents[1]
DATA = ROOT / "_data"
RAW = DATA / "raw"
SCORES = DATA / "scores"
sys.path.insert(0, str(REPO / "robust"))
warnings.filterwarnings("ignore")

from obd import cohorts as C  # noqa: E402
from obd import scoring  # noqa: E402
from obd import network as N  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd.curated import PRECLIN  # noqa: E402
from obd.reference import common_drug_name, gene_sets, drug_targets, drugbank_targets  # noqa: E402

# redirect every write of the reused modules into this folder
C.EXTERNAL = RAW
scoring.RESULTS = DATA

PAPER_CUTOFF = -1.2816
ALT_CUTOFFS = [-1.182405, -1.2816, -1.6449, -2.1925376]   # user's notebooks: -1.182405 (LUAD), -2.1925376 (BLCA); -1.6449 = 5 %

# set key -> (folder, curated set name, TCGA projects, response column, label)
SETS = {
    "vandewetering": ("COAD_vandewetering", "colorectal_vandewetering2015", ["COAD", "READ"], "ln_ic50",
                      "Colorectal organoids (van de Wetering 2015)"),
    "tiriac": ("PAAD_tiriac", "pancreas_tiriac2018", ["PAAD"], "response", "Pancreatic organoids (Tiriac 2018)"),
    "shi": ("PAAD_shi", "pancreas_shi2022", ["PAAD"], "response", "Pancreatic organoids (Shi 2022)"),
    "ji": ("LIHC_ji", "liver_ji2023", ["LIHC"], "response", "Liver cancer organoids (Ji 2023)"),
    "broutier": ("LIHC_broutier", "liver_broutier2017", ["LIHC"], "response", "Liver cancer organoids (Broutier 2017)"),
    "isella": ("COAD_isella_pdx", "colorectal_pdx_isella2017", ["COAD", "READ"], "response",
               "Colorectal PDX (Isella 2017, cetuximab)"),
}
# lab-model drug name -> (patient drug, proximity column used)
DRUG_MAP = {"SN-38": "IRINOTECAN"}                    # SN-38 is the active metabolite of irinotecan
PROX_PROXY = {"OXALIPLATIN": "CISPLATIN"}             # no protein target of its own: platinum class (as robust/run_survlearn.py)


def set_response(key):
    folder, name, projects, col, label = SETS[key]
    r = pd.read_csv(PRECLIN / name / "response.tsv", sep="\t", dtype={"sample": str})
    r["drug"] = r["drug"].str.upper()
    return r.groupby(["sample", "drug"])[col].median().unstack()


def set_expression(key):
    name = SETS[key][1]
    e = pd.read_csv(PRECLIN / name / "expression.tsv.gz", sep="\t", index_col=0)
    e.columns = e.columns.astype(str)
    return C.to_canonical(e)


def tcga_expr(project):
    return C.tcga_star_fpkm_uq(project)


_TREAT = {}


def treatment_table(project):
    """Long (patient, drug) table from the GDC BCR biotab clinical_drug file (two-pass name resolution
    exactly as robust/obd/tcga_response.build)."""
    if project in _TREAT:
        return _TREAT[project]
    df = R._read_biotab(RAW / f"nationwidechildrens.org_clinical_drug_{project.lower()}.txt", project)
    uniq = df["drug_raw"].unique()
    first = {r: R.normalize_drug(r) for r in uniq}
    freq = pd.Series([n for r in df["drug_raw"] for n, ok in first[r] if ok]).value_counts()
    vocab = tuple(freq[freq >= 3].index)
    rows = []
    for p, raw in zip(df["patient"], df["drug_raw"]):
        for n, ok in R.normalize_drug(raw, fuzzy_vocab=vocab):
            rows.append((p, n))
    out = pd.DataFrame(rows, columns=["patient", "drug"]).drop_duplicates()
    out.attrs["any_record"] = set(df["patient"])
    _TREAT[project] = out
    return out


def cdr():
    return R._build_cdr(RAW / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv")


def proximity(drug):
    """Paper's STRING>700 Reactome proximity z-scores for the drug (or its proxy)."""
    paper = N.precomputed_proximity()
    d = PROX_PROXY.get(drug, drug)
    return paper[d].dropna() if d in paper.columns else None, d


def scores(key, projects):
    """ssGSEA raw ES (pathways x samples) for the lab models and each TCGA project, on the shared gene
    universe (lab models  and  all projects), Reactome sets with >= 5 genes, as robust/run_paper_method_all.py."""
    name = SETS[key][1]
    tag = f"{name}_{'_'.join(projects)}"
    f = SCORES / f"{tag}_models.tsv.gz"
    if f.exists():
        out = {"models": pd.read_csv(f, sep="\t", index_col=0)}
        for p in projects:
            out[p] = pd.read_csv(SCORES / f"{tag}_{p}.tsv.gz", sep="\t", index_col=0)
        return out
    SCORES.mkdir(parents=True, exist_ok=True)
    oexpr = set_expression(key)
    texpr = {p: tcga_expr(p) for p in projects}
    genes = oexpr.index
    for p in projects:
        genes = genes.intersection(texpr[p].index)
    gs = {k: [g for g in v if g in set(genes)] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    out = {"models": scoring.ssgsea(oexpr.loc[genes], gs, threads=2)}
    out["models"].to_csv(f, sep="\t")
    for p in projects:
        out[p] = scoring.ssgsea(texpr[p].loc[genes], gs, threads=2)
        out[p].to_csv(SCORES / f"{tag}_{p}.tsv.gz", sep="\t")
    pd.Series({"genes": len(genes), "pathways": len(gs)}).to_csv(SCORES / f"{tag}_universe.tsv", sep="\t")
    return out


# ------------------------------------------------------------------ authors' committed data (exact reproduction)
def authors_inputs(drug):
    """Authors' committed organoid & TCGA NES, 2019 BCR biotab drug and patient files (as
    robust/run_paper_single_pathway.authors_data, generalised to any COAD organoid drug)."""
    nes = lambda p: pd.read_csv(p, sep="\t", comment="#", index_col=0)  # noqa: E731
    org = nes(REPO / "python/results/COAD/organoid/reactome_ssgsea_result/gseapy.samples.normalized.es.txt")
    pat = nes(REPO / "python/results/COAD/TCGA/reactome_ssgsea_result/gseapy.samples.normalized.es.txt")
    org.columns = [c.upper() for c in org.columns]
    resp = pd.read_csv(REPO / "data/organoid_COAD/drug_response/median_response.txt", sep="\t")
    resp["d"] = resp["Drug Name"].map(common_drug_name)
    resp = resp[resp["d"] == drug]
    resp = pd.Series(dict(zip(resp["Organoid"].str.upper(), resp["median_IC50"])))
    syn = {}
    for line in csv.reader(open(REPO / "data/drugbank_vocabulary.csv", encoding="utf-8", errors="replace")):
        if "DrugBank" in line[0]:
            continue
        for s in line[5].split(" | "):
            s = s.upper()
            if s:
                syn[s] = syn[s.replace("-", "")] = syn[s.replace(" ", "")] = line[2].upper()
    d = pd.read_csv(REPO / "data/TCGA_COAD/clinical_drug/dae92461-d098-40a3-997e-1a2d3e5b8dfe/nationwidechildrens.org_clinical_drug_coad.txt", sep="\t")
    treated, anydrug = set(), set()
    for p, dn in zip(d["bcr_patient_barcode"], d["pharmaceutical_therapy_drug_name"]):
        if "TCGA" not in str(p):
            continue
        anydrug.add(p)
        for t in str(dn).upper().replace("-", "").replace(" ", "").split("AND"):
            if "UNKNOWN" not in t and "NOTAVAILABLE" not in t and syn.get(t, t) == drug:
                treated.add(p)
    c = pd.read_csv(REPO / "data/TCGA_COAD/clinical_patient/4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt", sep="\t")
    sur = {}
    for p, st, lc, dd in zip(c["bcr_patient_barcode"], c["vital_status"], c["last_contact_days_to"], c["death_days_to"]):
        st = str(st).upper()
        days = lc if st == "ALIVE" else dd if st == "DEAD" else None
        if "TCGA" in str(p) and days is not None and "[" not in str(days):
            sur[p] = (float(days) / (365. / 12.), 1 if st == "DEAD" else 0)
    sur = pd.DataFrame(sur, index=["months", "event"]).T
    return org, resp, pat, sur, treated, anydrug


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
