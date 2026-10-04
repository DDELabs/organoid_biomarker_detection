import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
"""Build data/external/geo/<cohort>/{expression.tsv.gz, clinical.tsv} from mirrored sources."""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import rda_eset

RAW = Path("/tmp/obd_raw/geo_raw")
REPO = Path(_REPO)
OUT = REPO / "data/external/geo"

COLS = ["sample", "os_months", "os_event", "rfs_months", "rfs_event", "stage", "age", "sex",
        "chemo", "regimen", "fu_based", "response"]
FU_WORDS = ("5FU", "5-FU", "5 FU", "FLUOROURACIL", "FUFOL", "FOLFOX", "FOLFIRI", "CAPECITABINE",
            "XELODA", "LV5FU", "TEGAFUR", "UFT", "S-1", "XELOX", "CAPOX", "FLOX")
ROMAN = {0: np.nan, 1: "I", 2: "II", 3: "III", 4: "IV"}
DUKES = {"A": "I", "B": "II", "C": "III", "D": "IV"}


def is_fu(regimen):
    if regimen is None or (isinstance(regimen, float) and np.isnan(regimen)):
        return np.nan
    r = str(regimen).upper()
    return int(any(w in r for w in FU_WORDS))


def write(name, expr, clin):
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    expr = expr.loc[:, [s for s in expr.columns if s in set(clin["sample"])]]
    clin = clin[clin["sample"].isin(expr.columns)].copy()
    expr.index.name = "gene"
    expr.round(4).to_csv(d / "expression.tsv.gz", sep="\t", compression="gzip")
    clin = clin[COLS + [c for c in clin.columns if c not in COLS]]
    clin.to_csv(d / "clinical.tsv", sep="\t", index=False, na_rep="NA")
    print(f"{name}: expression {expr.shape}, clinical {clin.shape}")


def gpl570_map():
    ann = pd.read_csv(RAW / "GPL570/GPL570-55999.txt", sep="\t", comment="#", usecols=["ID", "Gene Symbol"], dtype=str)
    ann = ann.dropna()
    ann = ann[~ann["Gene Symbol"].str.contains("///")]  # unique-gene probes only
    return dict(zip(ann["ID"], ann["Gene Symbol"].str.strip()))


def probes_to_genes(expr, mapping):
    sym = expr.index.map(mapping)
    expr = expr[sym.notna()]
    expr.index = sym[sym.notna()]
    return expr.groupby(level=0).median()


# ------------------------------------------------------------------ GSE39582
def gse39582():
    raw = RAW / "GSE39582"
    expr = pd.read_csv(raw / "poshine_GSE39582.frma.tsv", sep="\t", index_col=0)
    expr = probes_to_genes(expr, gpl570_map())

    j = pd.read_csv(raw / "johnny.csv")  # GEO characteristics excerpt (566 tumours)
    j.columns = ["sample", "gender", "age", "tnm_stage", "rfs_event", "rfs_months", "os_event", "os_months",
                 "tnm_t", "tnm_n", "tnm_m", "tumor_location", "chemotherapy_adjuvant", "chemotherapy_adjuvant_type"]
    for c in ("tumor_location", "chemotherapy_adjuvant", "chemotherapy_adjuvant_type"):
        j[c] = j[c].str.split(": ", n=1).str[1].replace("NA", np.nan)
    x = pd.read_csv(raw / "xmuyu_clin.txt", sep="\t").rename(columns={"Patient_ID": "sample"})
    x = x[["sample", "MMR_status", "CIMP_status", "CIN_status", "CIT_molecularsubtype", "KRAS_mutation",
           "BRAF_mutation", "TP53_mutation", "Overall_time", "Overall_event", "RFS_time", "RFS_event", "Adjuvant_chemotherapy"]]
    c = j.merge(x, on="sample", how="left")
    # consistency checks between the two independent mirrors
    for a, b in (("os_months", "Overall_time"), ("os_event", "Overall_event"), ("rfs_months", "RFS_time"), ("rfs_event", "RFS_event")):
        m = c[a].notna() & c[b].notna()
        assert (np.abs(c.loc[m, a] - c.loc[m, b]) < 1e-6).all(), a
    assert (c["chemotherapy_adjuvant"].fillna("NA") == c["Adjuvant_chemotherapy"].fillna("NA")).all()
    c = c.drop(columns=["Overall_time", "Overall_event", "RFS_time", "RFS_event", "Adjuvant_chemotherapy"])
    c["stage"] = pd.to_numeric(c["tnm_stage"], errors="coerce").map(ROMAN)
    c["sex"] = c["gender"].str.upper().map({"MALE": 1, "FEMALE": 0})
    c["chemo"] = c["chemotherapy_adjuvant"].map({"Y": 1, "N": 0})
    reg = c["chemotherapy_adjuvant_type"].copy()
    reg[(c["chemo"] == 1) & reg.isna()] = "unspecified"
    reg[c["chemo"] == 0] = "none"
    c["regimen"] = reg
    c["fu_based"] = np.where(c["chemo"] == 0, 0, c["regimen"].map(lambda r: np.nan if r in ("unspecified", "other", None) or pd.isna(r) else is_fu(r)))
    c["response"] = np.nan
    c["setting"] = "adjuvant"
    write("GSE39582", expr, c)


# ------------------------------------------------------------------ curatedCRCData helpers
def curated(gse):
    _, expr, p, feat = rda_eset.read_eset(RAW / f"curatedCRC/{gse}_eset.rda")
    return expr, p


def gse14333():
    expr, p = curated("GSE14333")
    u = p["uncurated_author_metadata"].str.split("///", expand=True)
    meta = {}
    for col in u.columns:
        kv = u[col].str.replace(r"^X\d+: ", "", regex=True).str.split(": ", n=1)
        key = kv.str[0].iloc[0]
        meta[key] = kv.str[1].str.strip().replace({"": np.nan, "NA": np.nan})
    m = pd.DataFrame(meta, index=p.index)
    c = pd.DataFrame({"sample": p.index})
    c["os_months"] = np.nan
    c["os_event"] = np.nan
    c["rfs_months"] = pd.to_numeric(m["DFS_Time"], errors="coerce").values
    # authors coded DFS_Cens 1 = censored, 0 = event (see curatedCRCData GSE14333_curation.r)
    c["rfs_event"] = pd.to_numeric(m["DFS_Cens"], errors="coerce").map({1: 0, 0: 1}).values
    c["stage"] = m["DukesStage"].map(DUKES).values
    c["age"] = pd.to_numeric(m["Age_Diag"], errors="coerce").values
    c["sex"] = m["Gender"].map({"M": 1, "F": 0}).values
    c["chemo"] = m["AdjCTX"].map({"Y": 1, "N": 0}).values
    c["regimen"] = np.where(c["chemo"] == 1, "adjuvant_unspecified", np.where(c["chemo"] == 0, "none", None))
    c["fu_based"] = np.where(c["chemo"] == 0, 0, np.nan)
    c["fu_based_presumed"] = np.where(c["chemo"] == 1, 1, c["fu_based"])
    c["response"] = np.nan
    c["setting"] = "adjuvant"
    c["dukes_stage"] = m["DukesStage"].values
    c["location"] = m["Location"].values
    c["adjuvant_radiotherapy"] = m["AdjXRT"].map({"Y": 1, "N": 0}).values
    write("GSE14333", expr, c)


def gse28702():
    expr, p = curated("GSE28702")
    meta = p["uncurated_author_metadata"]
    get = lambda key: meta.str.extract(rf"{re.escape(key)}: ([^/]+)")[0].str.strip()
    c = pd.DataFrame({"sample": p.index})
    for k in ("os_months", "os_event", "rfs_months", "rfs_event", "stage"):
        c[k] = np.nan
    c["age"] = pd.to_numeric(p["age_at_initial_pathologic_diagnosis"], errors="coerce").values
    c["sex"] = p["gender"].map({"m": 1, "f": 0}).values
    c["chemo"] = 1
    c["regimen"] = "mFOLFOX6"
    c["fu_based"] = 1
    c["response"] = np.nan  # only responder / non-responder published (responder = CR+PR, RECIST)
    c["responder"] = p["drug_response"].map({"y": 1, "n": 0}).values
    c["setting"] = "first-line metastatic (unresectable)"
    c["sample_type"] = p["sample_type"].values
    c["lesion"] = get("characteristics_ch1: lesion").values
    c["location"] = get("characteristics_ch1.1: location").values
    c["set"] = get("characteristics_ch1.4: type").values
    c["title"] = meta.str.extract(r"title: ([^/]+)")[0].str.strip().values
    c["grade"] = get("Grade").values
    write("GSE28702", expr, c)


# ------------------------------------------------------------------ TCGA-READ
def tcga_read():
    raw = RAW / "TCGA-READ"
    sys.path.insert(0, str(REPO / "robust"))
    from obd.reference import ensembl_to_symbol
    expr = pd.read_csv(raw / "TCGA-READ.star_fpkm-uq.tsv.gz", sep="\t", index_col=0)
    expr.index = expr.index.str.split(".").str[0]
    keep = {}
    for col in sorted(col for col in expr.columns if col[13:15] == "01"):
        keep.setdefault(col[:12], col)
    expr = expr[list(keep.values())]
    expr.columns = list(keep.keys())
    e2s = ensembl_to_symbol()
    sym = expr.index.map(lambda e: e2s.get(e))
    expr = expr[sym.notna()]
    expr.index = sym[sym.notna()]
    expr = expr.groupby(level=0).median()

    sur = pd.read_csv(raw / "TCGA-READ.survival.tsv.gz", sep="\t")
    sur = sur[sur["sample"].str[13:15] == "01"].drop_duplicates("_PATIENT").set_index("_PATIENT")
    cli = pd.read_csv(raw / "TCGA-READ.clinical.tsv.gz", sep="\t", low_memory=False)
    cli["patient"] = cli["sample"].str[:12]
    cli = cli.drop_duplicates("patient").set_index("patient")
    drug = pd.read_csv(raw / "nationwidechildrens.org_clinical_drug_read.txt", sep="\t", skiprows=[1, 2])
    drug = drug[drug["bcr_patient_barcode"].astype(str).str.startswith("TCGA")]
    drug = drug[drug["pharmaceutical_therapy_type"].astype(str).str.contains("Chemotherapy", case=False)]
    resp_map = {"Complete Response": "CR", "Partial Response": "PR", "Stable Disease": "SD",
                "Clinical Progressive Disease": "PD", "Progressive Disease": "PD"}
    agg = {}
    for pat, g in drug.groupby("bcr_patient_barcode"):
        names = sorted({str(n).strip() for n in g["pharmaceutical_therapy_drug_name"] if "[" not in str(n)})
        settings = sorted({str(s) for s in g["therapy_regimen"] if "[" not in str(s)})
        resp = [resp_map.get(str(r)) for r in g["treatment_best_response"]]
        resp = [r for r in resp if r]
        order = ["CR", "PR", "SD", "PD"]
        agg[pat] = {"regimen": "+".join(names) if names else "unspecified",
                    "therapy_setting": ";".join(settings) if settings else np.nan,
                    "response": min(resp, key=order.index) if resp else np.nan}
    agg = pd.DataFrame.from_dict(agg, orient="index")

    pats = list(expr.columns)
    c = pd.DataFrame({"sample": pats})
    s = sur.reindex(pats)
    k = cli.reindex(pats)
    c["os_months"] = (s["OS.time"] / (365.25 / 12)).values
    c["os_event"] = s["OS"].values
    c["rfs_months"] = np.nan
    c["rfs_event"] = np.nan
    stage = k["ajcc_pathologic_stage.diagnoses"].fillna("").astype(str).str.upper().str.replace("STAGE", "").str.strip()
    c["stage"] = stage.map(lambda v: next((r for r in ("IV", "III", "II", "I") if str(v).startswith(r)), np.nan)).values
    c["age"] = pd.to_numeric(k["age_at_index.demographic"], errors="coerce").values
    c["sex"] = k["gender.demographic"].str.lower().map({"male": 1, "female": 0}).values
    a = agg.reindex(pats)
    c["chemo"] = np.where(a["regimen"].notna(), 1, np.nan)  # absence of a drug record is not evidence of no chemo
    c["regimen"] = a["regimen"].values
    c["fu_based"] = [is_fu(r) if r != "unspecified" and isinstance(r, str) else np.nan for r in a["regimen"]]
    c["response"] = a["response"].values
    c["setting"] = a["therapy_setting"].values
    write("TCGA-READ", expr, c)


if __name__ == "__main__":
    for f in sys.argv[1:] or ["gse39582", "gse14333", "gse28702", "tcga_read"]:
        globals()[f]()
