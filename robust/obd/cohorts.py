"""Load organoid and patient cohorts into a common shape.

Every loader returns plain pandas objects:
  expression : DataFrame genes (canonical symbols) x samples
  response   : DataFrame samples x drugs (IC50 / AUC, lower = more sensitive)
  clinical   : DataFrame indexed by 12-char patient barcode with columns
               months, event, stage (I-IV or NaN), age, sex
  treatment  : dict {common drug name: set(patient barcodes)}
"""
import gzip
import urllib.request

import numpy as np
import pandas as pd

from . import DATA, EXTERNAL
from .reference import canonical_symbol, common_drug_name, ensembl_to_symbol

XENA_GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"


# --------------------------------------------------------------------- download
def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print(f"downloading {url}")
        tmp = dest.with_suffix(dest.suffix + ".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)
    return dest


# --------------------------------------------------------------------- organoids
def coad_organoids():
    """van de Wetering 2015 colorectal organoids (GSE64392 median expression, median IC50)."""
    path = DATA / "organoid_COAD/expression/GSE64392/GSE64392_series_matrix.txt/geneID_expression_median.txt"
    expr = pd.read_csv(path, sep="\t", index_col=0).drop(columns="uniprotID")
    ic = {}
    with open(DATA / "organoid_COAD/drug_response/supp_table_S2b.txt") as f:
        for line in f:
            p = line.strip().split("\t")
            if "Organoid" in p[0] or len(p) <= 5:
                continue
            ic.setdefault((p[0].upper(), common_drug_name(p[2])), []).append(float(p[4]))
    resp = pd.Series({k: np.median(v) for k, v in ic.items()}).unstack()
    return expr, resp


def generic_organoids(expr_path, response_path, sample_col="sample", drug_col="drug", value_col="response"):
    """Any organoid set: a genes x samples matrix (symbols in column 1) and a long response table."""
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = to_canonical(expr)
    long = pd.read_csv(response_path, sep="\t")
    long[drug_col] = long[drug_col].map(common_drug_name)
    resp = long.groupby([sample_col, drug_col])[value_col].median().unstack()
    return expr, resp


def to_canonical(expr):
    sym = expr.index.map(canonical_symbol)
    expr = expr[sym.notna()]
    expr.index = sym[sym.notna()]
    return expr.groupby(level=0).median()


# --------------------------------------------------------------------- TCGA
def tcga_star_fpkm_uq(project):
    """Current GDC STAR FPKM-UQ (log2(x+1), GENCODE v36) from the Xena GDC hub.

    Primary tumours only (-01), one aliquot per patient (vial A preferred),
    Ensembl IDs mapped to the canonical symbols used by the original pipeline.
    """
    path = fetch(XENA_GDC + f"TCGA-{project}.star_fpkm-uq.tsv.gz", EXTERNAL / f"TCGA-{project}.star_fpkm-uq.tsv.gz")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.str.split(".").str[0]
    keep = {}
    for c in sorted(c for c in expr.columns if c[13:15] == "01"):
        keep.setdefault(c[:12], c)
    expr = expr[list(keep.values())]
    expr.columns = list(keep.keys())
    e2s = ensembl_to_symbol()
    sym = expr.index.map(lambda e: canonical_symbol(e2s[e]) if e in e2s else None)
    expr = expr[sym.notna()]
    expr.index = sym[sym.notna()]
    return expr.groupby(level=0).median()


def _stage(s):
    s = str(s).upper().replace("STAGE", "").strip()
    for roman in ("IV", "III", "II", "I"):
        if s.startswith(roman):
            return roman
    return np.nan


def tcga_xena_clinical(project):
    """Survival (OS) and covariates from the Xena GDC hub (current GDC release)."""
    sur = pd.read_csv(fetch(XENA_GDC + f"TCGA-{project}.survival.tsv.gz", EXTERNAL / f"TCGA-{project}.survival.tsv.gz"), sep="\t")
    sur = sur[sur["sample"].str[13:15] == "01"].drop_duplicates("_PATIENT").set_index("_PATIENT")
    cli = pd.read_csv(fetch(XENA_GDC + f"TCGA-{project}.clinical.tsv.gz", EXTERNAL / f"TCGA-{project}.clinical.tsv.gz"), sep="\t", low_memory=False)
    cli["patient"] = cli["sample"].str[:12]
    cli = cli.drop_duplicates("patient").set_index("patient")
    out = pd.DataFrame({"months": sur["OS.time"] / (365 / 12), "event": sur["OS"]})
    out["stage"] = cli.reindex(out.index)["ajcc_pathologic_stage.diagnoses"].map(_stage)
    out["age"] = pd.to_numeric(cli.reindex(out.index)["age_at_index.demographic"], errors="coerce")
    out["sex"] = (cli.reindex(out.index)["gender.demographic"].str.lower() == "male").astype(float)
    return out.dropna(subset=["months", "event"])


def tcga_biotab_clinical(project):
    """Survival from the 2019 BCR biotab patient file shipped with the repo (paper-era)."""
    path = next((DATA / f"TCGA_{project}/clinical_patient").glob("*/nationwidechildrens.org_clinical_patient_*.txt"))
    df = pd.read_csv(path, sep="\t")
    rows = {}
    for _, r in df.iterrows():
        pat, status = r["bcr_patient_barcode"], str(r["vital_status"]).upper()
        if "TCGA" not in str(pat):
            continue
        days = r["last_contact_days_to"] if status == "ALIVE" else r["death_days_to"] if status == "DEAD" else "["
        if "[" in str(days):
            continue
        rows[pat] = {"months": float(days) / (365 / 12), "event": int(status == "DEAD"),
                     "stage": _stage(r.get("ajcc_pathologic_tumor_stage")),
                     "age": pd.to_numeric(r.get("age_at_initial_pathologic_diagnosis"), errors="coerce"),
                     "sex": float(str(r.get("gender")).upper() == "MALE")}
    return pd.DataFrame.from_dict(rows, orient="index")


def tcga_biotab_drugs(project):
    """{common drug name: set(patients)} from the BCR biotab clinical_drug file."""
    path = next((DATA / f"TCGA_{project}/clinical_drug").glob("*/nationwidechildrens.org_clinical_drug_*.txt"))
    df = pd.read_csv(path, sep="\t")
    out = {}
    for pat, name in zip(df["bcr_patient_barcode"], df["pharmaceutical_therapy_drug_name"]):
        if "TCGA" not in str(pat):
            continue
        for d in str(name).upper().replace("-", "").replace(" ", "").split("AND"):
            if "UNKNOWN" in d or "NOTAVAILABLE" in d or d == "NAN":
                continue
            out.setdefault(common_drug_name(d), set()).add(pat)
    return out


def paper_scores(cancer, source):
    """Pathway NES matrices committed by the original authors (pathways x samples)."""
    path = DATA.parent / f"python/results/{cancer.upper()}/{source}/reactome_ssgsea_result/gseapy.samples.normalized.es.txt"
    return pd.read_csv(path, sep="\t", comment="#", index_col=0)
