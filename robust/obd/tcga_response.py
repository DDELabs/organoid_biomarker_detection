"""Pan-cancer TCGA drug-response labels (BCR biotab clinical_drug) and TCGA-CDR survival.

Resource lives in data/external/tcga_response/:
  drug_response_long.tsv.gz  one row per (patient, drug, regimen line / drug record)
  survival_cdr.tsv.gz        TCGA-CDR (Liu et al., Cell 2018) per patient
  SUMMARY.md                 counts per drug / project and source URLs + sha256
  raw/                       verbatim downloads used by build()

Public API
  load_response_labels(drugs=None, projects=None, labelled_only=False) -> long DataFrame
  load_cdr() -> DataFrame indexed by 12-char patient barcode
  patient_labels(drug, projects=None, how="first") -> one label per patient for a drug
  drug_label_matrix(min_patients=30) -> per-drug summary (patients, responders, projects)
  build() -> (re)download raw tables and rewrite the resource

Response semantics: `treatment_best_response` in the biotab drug form is the best
response to that regimen as reported by the tissue source site. It is NOT centrally
reviewed RECIST; for adjuvant regimens (no measurable disease) "Complete Response"
usually means "no recurrence during/after treatment". `regimen_indication`
(ADJUVANT / PROGRESSION / ...) is never filled on rows that carry a response, so the
column `setting_proxy` (post_progression / early / late_no_event / unknown; derived from
drug start day vs. the TCGA-CDR PFI event) is provided to restrict to the
recurrent/progressive setting when a tumour-shrinkage label is wanted.
"""
import difflib
import hashlib
import re
import urllib.request
from functools import lru_cache

import numpy as np
import pandas as pd

from . import EXTERNAL
from .reference import BRANDS, common_drug_name, drug_synonyms

ROOT = EXTERNAL / "tcga_response"
RAW = ROOT / "raw"
LONG_PATH = ROOT / "drug_response_long.tsv.gz"
CDR_PATH = ROOT / "survival_cdr.tsv.gz"
SUMMARY_PATH = ROOT / "SUMMARY.md"

PROJECTS = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
            "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
            "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

KEMPLAB = ("https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/"
           "_data_/input/TCGA/nationwidechildrens.org_clinical_drug_{p}.txt")
# kemplab mirror lacks ACC (and LAML: TCGA never released a LAML clinical_drug biotab table).
MIRRORS = {"ACC": "https://raw.githubusercontent.com/nikcheerla/mirnanalyze/master/Scripts/cancer_clinical_data/"
                  "acc/nationwidechildrens.org_clinical_drug_acc.txt"}
# Ding, Zu & Gu, Bioinformatics 2016 curated pan-cancer table (as redistributed by VAEN, Jia et al. 2021)
DING2016 = "https://raw.githubusercontent.com/bsml320/VAEN/master/DATA/response/drug_response.txt"
CDR_URL = ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
           "Survival_SupplementalTable_S1_20171025_xena_sp")

RECIST = {"COMPLETE RESPONSE": "CR", "PARTIAL RESPONSE": "PR", "STABLE DISEASE": "SD",
          "CLINICAL PROGRESSIVE DISEASE": "PD", "RADIOGRAPHIC PROGRESSIVE DISEASE": "PD",
          "PROGRESSIVE DISEASE": "PD"}

# Brand names / abbreviations / misspellings seen in TCGA clinical_drug files that neither
# reference.BRANDS nor the DrugBank vocabulary resolve. Keys: upper case.
EXTRA_BRANDS = {
    "LUPRON": "LEUPROLIDE", "ELIGARD": "LEUPROLIDE", "LEUPORELINE": "LEUPROLIDE", "CASODEX": "BICALUTAMIDE",
    "ZOLADEX": "GOSERELIN", "GOSERELINE": "GOSERELIN", "TRELSTAR": "TRIPTORELIN", "LHRH AGONIST": "LHRH AGONIST",
    "CELEBREX": "CELECOXIB", "GLIADEL": "CARMUSTINE", "GLIADEL WAFER": "CARMUSTINE",
    "GLIADEL WAFERS": "CARMUSTINE", "CARMUSTIN": "CARMUSTINE", "LOMUSTIN": "LOMUSTINE", "CCNG": "LOMUSTINE",
    "SYNTHROID": "LEVOTHYROXINE", "LEVOXYL": "LEVOTHYROXINE", "CYTOMEL": "LIOTHYRONINE",
    "MEGACE": "MEGESTROL ACETATE", "SUTENT": "SUNITINIB", "ZOMETA": "ZOLEDRONIC ACID",
    "GLEEVEC": "IMATINIB", "GLEEVAC": "IMATINIB", "NEULASTA": "PEGFILGRASTIM", "PEG G-CSF": "PEGFILGRASTIM",
    "G-CSF": "FILGRASTIM", "NEUPOGEN": "FILGRASTIM", "LEUKINE": "SARGRAMOSTIM",
    "HEXALEN": "ALTRETAMINE", "HEXALIN": "ALTRETAMINE", "HEXAMETHYLMELAMINE": "ALTRETAMINE",
    "IL-2": "ALDESLEUKIN", "IL2": "ALDESLEUKIN", "INTERLEUKIN-2": "ALDESLEUKIN", "INTERLEUKIN 2": "ALDESLEUKIN",
    "PROLEUKIN": "ALDESLEUKIN", "TEMADOR": "TEMOZOLOMIDE", "TEMODAL": "TEMOZOLOMIDE",
    "VECTIBIX": "PANITUMUMAB", "ONCOVIN": "VINCRISTINE", "HYDROXYDAUNOMYCIN": "DOXORUBICIN",
    "DOXORUBICINA": "DOXORUBICIN", "DOXORUBICINE": "DOXORUBICIN", "ADRIAMICIN": "DOXORUBICIN",
    "ADRIMYCIN": "DOXORUBICIN", "ADRIMICIN": "DOXORUBICIN", "ADRIAMYICIN": "DOXORUBICIN",
    "ACCUTANE": "ISOTRETINOIN", "CIS RETINOIC ACID": "ISOTRETINOIN", "CIS-RETINOIC ACID": "ISOTRETINOIN",
    "13-CIS RETINOIC ACID": "ISOTRETINOIN", "AFINITOR": "EVEROLIMUS", "RAD001": "EVEROLIMUS",
    "YERVOY": "IPILIMUMAB", "TS-1": "TEGAFUR", "S-1": "TEGAFUR", "FARESTON": "TOREMIFENE",
    "INTRON A": "INTERFERON ALFA-2B", "LAFERON": "INTERFERON ALFA-2B", "SYLATRON": "PEGINTERFERON ALFA-2B",
    "PD 0332991": "PALBOCICLIB", "PD0332991": "PALBOCICLIB", "ALOXI": "PALONOSETRON", "XGEVA": "DENOSUMAB",
    "PROLIA": "DENOSUMAB", "RITUXAN": "RITUXIMAB", "TORISEL": "TEMSIROLIMUS", "VOTRIENT": "PAZOPANIB",
    "AZD6244": "SELUMETINIB", "ABT-888": "VELIPARIB", "XL 184": "CABOZANTINIB", "XL184": "CABOZANTINIB",
    "ARQ-197": "TIVANTINIB", "XYOTAX": "PACLITAXEL POLIGLUMEX", "ET-743": "TRABECTEDIN",
    "YONDELIS": "TRABECTEDIN", "DECADRON": "DEXAMETHASONE", "REVLIMID": "LENALIDOMIDE", "ARA-C": "CYTARABINE",
    "PREDISONE": "PREDNISONE", "NOVADEX": "TAMOXIFEN", "TAMOXIPHENE": "TAMOXIFEN", "TAMOXIPHEN": "TAMOXIFEN",
    "TAMOXITEN": "TAMOXIFEN", "PHOTOFRIN": "PORFIMER", "PHOTOFIN": "PORFIMER", "PORFIMER SODIUM": "PORFIMER",
    "AMITOSTINE": "AMIFOSTINE", "AMITOSTIN": "AMIFOSTINE", "ETHYOL": "AMIFOSTINE",
    "6 THIGUANINE": "TIOGUANINE", "6-THIOGUANINE": "TIOGUANINE", "06-BG": "6-O-BENZYLGUANINE",
    "06BG": "6-O-BENZYLGUANINE", "06GB": "6-O-BENZYLGUANINE", "O6B6": "6-O-BENZYLGUANINE",
    "O6BG": "6-O-BENZYLGUANINE", "O6-BG": "6-O-BENZYLGUANINE", "CALCIUM FOLIATUM": "LEUCOVORIN",
    "FOLINIC ACID": "LEUCOVORIN", "LEVCOVORIN": "LEUCOVORIN", "BCG": "BCG VACCINE",
    "BACILLUS CALMETTE-GUERIN": "BCG VACCINE", "TICE BCG": "BCG VACCINE", "THERACYS": "BCG VACCINE",
    "CAPECYTABINUM": "CAPECITABINE", "CAPECETABINE": "CAPECITABINE", "METOTREKSAT": "METHOTREXATE",
    "METHOTREXATUM": "METHOTREXATE", "FLUOROURACILLUM": "FLUOROURACIL", "5 FLUOROURACIL": "FLUOROURACIL",
    "5-FLOUROURACIL": "FLUOROURACIL", "5-FLUROURACIL": "FLUOROURACIL", "5-FLUORUORACIL": "FLUOROURACIL",
    "FLOUROURACIL": "FLUOROURACIL", "CARBO": "CARBOPLATIN", "CARBPLATIN": "CARBOPLATIN",
    "PACILTAXEL": "PACLITAXEL", "PACILTAXLE": "PACLITAXEL", "PACITAXEL": "PACLITAXEL", "PACITAXOL": "PACLITAXEL",
    "PACLITAXOL": "PACLITAXEL", "NAB-PACLITAXEL": "PACLITAXEL", "ALBUMIN-BOUND PACLITAXEL": "PACLITAXEL",
    "DOXETAXEL": "DOCETAXEL", "TRUSTUZUMAB": "TRASTUZUMAB", "TOPTECAN": "TOPOTECAN", "TOPETECAN": "TOPOTECAN",
    "TOPECAN": "TOPOTECAN", "GEMCITIBINE": "GEMCITABINE", "GEMICITABINE": "GEMCITABINE",
    "CYCLOPHOSPHANE": "CYCLOPHOSPHAMIDE", "CYCLOPHOSPHAMID": "CYCLOPHOSPHAMIDE", "CYTOXEN": "CYCLOPHOSPHAMIDE",
    "CYOTXAN": "CYCLOPHOSPHAMIDE", "CUCLOPHOSPHAMIDE": "CYCLOPHOSPHAMIDE",
    "CYCLOPHASPHAMIDE": "CYCLOPHOSPHAMIDE", "CYCLOPHOSPAMIDE": "CYCLOPHOSPHAMIDE",
    "HYDROXUREA": "HYDROXYUREA", "ETOPSIDE": "ETOPOSIDE", "IRINTOTECAN": "IRINOTECAN",
    "TEMOZOLOMODE": "TEMOZOLOMIDE", "TEMOZLOMIDE": "TEMOZOLOMIDE", "ANASTRAZOLE": "ANASTROZOLE",
    "ANASTRAZOLUM": "ANASTROZOLE", "ARMIDEX": "ANASTROZOLE", "CISPLATNIN": "CISPLATIN",
    "DOXOBUBICIN": "DOXORUBICIN", "EPIRUBICOIN": "EPIRUBICIN", "SAHA": "VORINOSTAT", "ZACTIMA": "VANDETANIB",
    "AT-101": "AT-101", "AT 101": "AT-101", "CI 980": "CI-980", "OVAREX": "OREGOVOMAB",
    "TIVANTINIB (ARQ-197)": "TIVANTINIB", "KEYTRUDA": "PEMBROLIZUMAB", "OPDIVO": "NIVOLUMAB",
    "ZELBORAF": "VEMURAFENIB", "TAFINLAR": "DABRAFENIB", "MEKINIST": "TRAMETINIB", "TYKERB": "LAPATINIB",
    "PERJETA": "PERTUZUMAB", "KADCYLA": "TRASTUZUMAB EMTANSINE", "IBRANCE": "PALBOCICLIB",
    "ZYTIGA": "ABIRATERONE", "XTANDI": "ENZALUTAMIDE", "FIRMAGON": "DEGARELIX", "INLYTA": "AXITINIB",
    "STIVARGA": "REGORAFENIB", "ZALTRAP": "AFLIBERCEPT", "LYNPARZA": "OLAPARIB", "TEMSIROLIMUS": "TEMSIROLIMUS",
    "PROVERA": "MEDROXYPROGESTERONE ACETATE", "DEPO-PROVERA": "MEDROXYPROGESTERONE ACETATE",
    "MITOMYCIN C": "MITOMYCIN", "MUTAMYCIN": "MITOMYCIN", "VALRUBICIN": "VALRUBICIN", "VALSTAR": "VALRUBICIN",
    "IFEX": "IFOSFAMIDE", "MESNEX": "MESNA", "ZOFRAN": "ONDANSETRON", "EMEND": "APREPITANT",
    "PROCRIT": "EPOETIN ALFA", "ARANESP": "DARBEPOETIN ALFA", "AREDIA": "PAMIDRONATE",
    "THALOMID": "THALIDOMIDE", "VELCADE": "BORTEZOMIB", "TARGRETIN": "BEXAROTENE",
    "DTIC": "DACARBAZINE", "DTIC-DOME": "DACARBAZINE", "MATULANE": "PROCARBAZINE", "BLENOXANE": "BLEOMYCIN",
    "VELBAN": "VINBLASTINE", "MUSTARGEN": "MECHLORETHAMINE", "ALKERAN": "MELPHALAN", "LEUKERAN": "CHLORAMBUCIL",
    "CEENU": "LOMUSTINE", "GLEOSTINE": "LOMUSTINE", "BICNU": "CARMUSTINE", "MITOXANTRONE HCL": "MITOXANTRONE",
    "MESNA": "MESNA", "INTERFERON": "INTERFERON ALFA-2B", "INTERLEUKIN": "ALDESLEUKIN", "TARVECA": "ERLOTINIB",
    "TANCEVA": "ERLOTINIB", "OSI-774": "ERLOTINIB", "OS1-774": "ERLOTINIB", "LOMUSTINE CCNU": "LOMUSTINE",
    "IRINTOCEAN": "IRINOTECAN", "BEVACOZIMAB": "BEVACIZUMAB", "PS341": "BORTEZOMIB", "PS-341": "BORTEZOMIB",
    "CELBREX": "CELECOXIB", "SARASAR": "LONAFARNIB", "TIPFARNIB": "TIPIFARNIB", "ZARNESTRA": "TIPIFARNIB",
    "HYDROCHOROQUINE": "HYDROXYCHLOROQUINE", "DOXORUBICIN HCL LIPOSOMAL": "DOXORUBICIN",
    "METILPREDNISONA": "METHYLPREDNISOLONE", "LY317615": "ENZASTAURIN", "ABT-888 PARP INHIBITOR": "VELIPARIB",
    "PLEXXIKON PLX3397": "PEXIDARTINIB", "PLX3397": "PEXIDARTINIB", "ABTC 0603 HYDROXYCHLOROQUINE": "HYDROXYCHLOROQUINE",
    "AMG 655": "CONATUMUMAB", "INTERFERON ALPHA": "INTERFERON ALFA-2B", "INTERFERON ALFA": "INTERFERON ALFA-2B",
    "INTERFERON-ALFA": "INTERFERON ALFA-2B", "INTERFERON-ALPHA": "INTERFERON ALFA-2B",
    "INTERFERON ALFA-2B": "INTERFERON ALFA-2B", "STUDY DRUG AMG 655": "CONATUMUMAB",
}

# Regimen acronyms -> component drugs (only unambiguous ones; "TC" is left unresolved)
REGIMENS = {
    "FOLFOX": ["FLUOROURACIL", "LEUCOVORIN", "OXALIPLATIN"], "FOLFOX4": ["FLUOROURACIL", "LEUCOVORIN", "OXALIPLATIN"],
    "FOLFOX6": ["FLUOROURACIL", "LEUCOVORIN", "OXALIPLATIN"], "MFOLFOX6": ["FLUOROURACIL", "LEUCOVORIN", "OXALIPLATIN"],
    "FOLFIRI": ["FLUOROURACIL", "LEUCOVORIN", "IRINOTECAN"],
    "FOLFIRINOX": ["FLUOROURACIL", "LEUCOVORIN", "IRINOTECAN", "OXALIPLATIN"],
    "XELOX": ["CAPECITABINE", "OXALIPLATIN"], "CAPOX": ["CAPECITABINE", "OXALIPLATIN"],
    "CHOP": ["CYCLOPHOSPHAMIDE", "DOXORUBICIN", "VINCRISTINE", "PREDNISONE"],
    "R-CHOP": ["RITUXIMAB", "CYCLOPHOSPHAMIDE", "DOXORUBICIN", "VINCRISTINE", "PREDNISONE"],
    "EPOCH": ["ETOPOSIDE", "PREDNISONE", "VINCRISTINE", "CYCLOPHOSPHAMIDE", "DOXORUBICIN"],
    "ICE": ["IFOSFAMIDE", "CARBOPLATIN", "ETOPOSIDE"],
    "CVAD": ["CYCLOPHOSPHAMIDE", "VINCRISTINE", "DOXORUBICIN", "DEXAMETHASONE"],
    "TCH": ["DOCETAXEL", "CARBOPLATIN", "TRASTUZUMAB"], "TAC": ["DOCETAXEL", "DOXORUBICIN", "CYCLOPHOSPHAMIDE"],
    "AC": ["DOXORUBICIN", "CYCLOPHOSPHAMIDE"], "CMF": ["CYCLOPHOSPHAMIDE", "METHOTREXATE", "FLUOROURACIL"],
    "FEC": ["FLUOROURACIL", "EPIRUBICIN", "CYCLOPHOSPHAMIDE"], "CAF": ["CYCLOPHOSPHAMIDE", "DOXORUBICIN", "FLUOROURACIL"],
    "FAC": ["FLUOROURACIL", "DOXORUBICIN", "CYCLOPHOSPHAMIDE"], "MVAC": ["METHOTREXATE", "VINBLASTINE", "DOXORUBICIN", "CISPLATIN"],
    "M-VAC": ["METHOTREXATE", "VINBLASTINE", "DOXORUBICIN", "CISPLATIN"], "GC": ["GEMCITABINE", "CISPLATIN"],
    "BEP": ["BLEOMYCIN", "ETOPOSIDE", "CISPLATIN"], "EP": ["ETOPOSIDE", "CISPLATIN"],
    "PCV": ["PROCARBAZINE", "LOMUSTINE", "VINCRISTINE"], "ECF": ["EPIRUBICIN", "CISPLATIN", "FLUOROURACIL"],
    "ECX": ["EPIRUBICIN", "CISPLATIN", "CAPECITABINE"], "EOX": ["EPIRUBICIN", "OXALIPLATIN", "CAPECITABINE"],
    "VIP": ["ETOPOSIDE", "IFOSFAMIDE", "CISPLATIN"], "TIP": ["PACLITAXEL", "IFOSFAMIDE", "CISPLATIN"],
}

# Non-specific entries that name no drug; such records are dropped from the long table.
NONSPECIFIC = {"", "NAN", "[NOT AVAILABLE]", "[UNKNOWN]", "[NOT APPLICABLE]", "[DISCREPANCY]", "[NOT EVALUATED]",
               "UNKNOWN", "CHEMO, NOS", "CHEMO NOS", "CHEMOTHERAPY", "HORMONE, NOS", "HORMONE NOS",
               "NOT OTHERWISE SPECIFIED", "CLINICAL TRIAL", "STUDY DRUG", "PLACEBO", "NOS", "OTHER", "TAXANE",
               "NONE", "UNKNOWN CHEMOTHERAPY", "UNKNOWN CHEMO", "CHEMO", "TC", "O-ICE", "IND", "VS", "VERSUS"}

_NOISE = [r"\b(?:VS\.?|VERSUS)\s.*$", r"^[A-Z]{2,5}\.?\d+\s+-\s+", r"\bSTUDY DRUG\b", r"\bOR PLACEBO\b", r"/?\s*PLACEBO\b.*$", r"\bPROVIDED BY STUDY\b", r"\+/-",
          r",?\s*INTRAVESIC(?:U|A)LAR\b", r"\bINTRATHECAL\b", r"\bHIGH DOSE\b", r"-XRT\b", r"^C\d+\s+",
          r"\bAROMATASE\b(?=\s+\w)", r"(?<=[A-Z]{4})-\d$"]
_SALTS = [" HYDROCHLORIDE", " HCL", " ACETATE", " TARTRATE", " CITRATE", " DISODIUM", " SODIUM", " MESYLATE",
          " MALATE", " TOSYLATE", " SULFATE", " CALCIUM", " LIPOSOMAL", " LIPOSOME", " PEGYLATED"]


@lru_cache(None)
def _vocab():
    return set(drug_synonyms().values())


def _resolve_token(tok):
    """Resolve one drug token -> (name, resolved flag)."""
    tok = re.sub(r"\s+", " ", tok.strip(" .;:-")).upper()
    if not tok or tok in NONSPECIFIC:
        return None, False
    vocab = _vocab()
    cands = [tok, re.sub(r"^LIPOSOMAL\s+", "", tok)]
    for s in _SALTS:
        if tok.endswith(s):
            cands.append(tok[: -len(s)].strip())
    if tok.endswith("UM") and len(tok) > 7:  # latinised names (CISPLATINUM, OXALIPLATINUM, LETROZOLUM ...)
        cands += [tok[:-2], tok[:-2] + "E", tok[:-3] + "E"]
    for c in cands:
        for key in (c, c.replace(" ", ""), c.replace("-", "")):
            if key in EXTRA_BRANDS:
                n = EXTRA_BRANDS[key]
                return common_drug_name(n) if common_drug_name(n) in vocab else n, True
        n = common_drug_name(c)
        if n in vocab:
            return n, True
        if c.replace(" ", "").replace("-", "") in BRANDS:
            return common_drug_name(BRANDS[c.replace(" ", "").replace("-", "")]), True
    return tok, False


def normalize_drug(raw, fuzzy_vocab=()):
    """Raw clinical_drug name -> list of (common name, resolved flag). Splits combinations
    ('A+B', 'A and B', 'A/B', 'A, B'), expands regimen acronyms, handles 'Brand (generic)'."""
    s = str(raw).upper().strip()
    if s in NONSPECIFIC or s.startswith("[") or s == "NAN":
        return []
    for pat in _NOISE:
        s = re.sub(pat, " ", s).strip()
    out = []
    # 'X (Y)': resolve each half separately and keep the resolvable one
    parts = re.split(r"\s*(?:\+|,|/|&|;|\bAND\b|\bTHEN\b|\bPLUS\b)\s*", s)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if part in REGIMENS or part.replace(" ", "") in REGIMENS:
            out += [(common_drug_name(d), True) for d in REGIMENS.get(part, REGIMENS.get(part.replace(" ", "")))]
            continue
        m = re.match(r"^(.*?)\s*\((.*?)\)?\s*$", part)
        halves = [m.group(1), m.group(2)] if m and "(" in part else [part]
        best = None
        for h in halves:
            if h and (h in REGIMENS):
                best = [(common_drug_name(d), True) for d in REGIMENS[h]]
                break
            n, ok = _resolve_token(h or "")
            if n is None:
                continue
            if ok:
                best = [(n, True)]
                break
            if best is None:
                best = [(n, False)]
        if best and not best[0][1] and fuzzy_vocab and len(best[0][0]) >= 6:
            hit = difflib.get_close_matches(best[0][0], list(fuzzy_vocab), n=1, cutoff=0.86)
            if hit:
                best = [(hit[0], True)]
        if best:
            out += best
    seen, dedup = set(), []
    for n, ok in out:
        if n not in seen:
            seen.add(n)
            dedup.append((n, ok))
    return dedup


def harmonise_recist(value):
    return RECIST.get(str(value).strip().upper(), np.nan)


# --------------------------------------------------------------------- build
def _download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)
    return dest


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_raw():
    """Download all raw tables into RAW; returns {name: (url, path)} for files that exist."""
    got = {}
    for p in PROJECTS:
        url = MIRRORS.get(p, KEMPLAB.format(p=p.lower()))
        dest = RAW / f"nationwidechildrens.org_clinical_drug_{p.lower()}.txt"
        try:
            _download(url, dest)
            got[p] = (url, dest)
        except Exception as e:  # LAML has no clinical_drug table anywhere
            print(f"{p}: no clinical_drug table ({e})")
    for name, url, dest in [("DING2016", DING2016, RAW / "ding2016_vaen_drug_response.txt"),
                            ("CDR", CDR_URL, RAW / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv")]:
        _download(url, dest)
        got[name] = (url, dest)
    return got


def _clean(v):
    s = str(v).strip()
    return np.nan if s in ("", "nan") or s.startswith("[") else s


def _read_biotab(path, project):
    df = pd.read_csv(path, sep="\t", dtype=str, skiprows=[1, 2])
    df = df[df["bcr_patient_barcode"].astype(str).str.startswith("TCGA")]
    col = lambda *names: next((df[n] for n in names if n in df.columns), pd.Series(np.nan, index=df.index))  # noqa: E731
    ind = col("regimen_indication").map(_clean)
    ind = ind.fillna(col("therapy_regimen").map(_clean))
    return pd.DataFrame({
        "project": project, "patient": df["bcr_patient_barcode"].str[:12], "bcr_drug_barcode": df["bcr_drug_barcode"],
        "drug_raw": col("pharmaceutical_therapy_drug_name").astype(str).str.strip(),
        "therapy_type": col("pharmaceutical_therapy_type").map(_clean),
        "start_days": pd.to_numeric(col("pharmaceutical_tx_started_days_to"), errors="coerce"),
        "end_days": pd.to_numeric(col("pharmaceutical_tx_ended_days_to"), errors="coerce"),
        "therapy_ongoing": col("pharmaceutical_tx_ongoing_indicator").map(_clean),
        "best_response": col("treatment_best_response").map(_clean),
        "regimen_indication": ind.astype(object).map(lambda v: v.upper() if isinstance(v, str) else np.nan),
        "regimen_number": pd.to_numeric(col("regimen_number"), errors="coerce"),
        "on_clinical_trial": col("tx_on_clinical_trial").map(_clean),
        "cycles": pd.to_numeric(col("pharma_adjuvant_cycles_count", "number_cycles"), errors="coerce"),
    })


def build():
    """Download sources and write drug_response_long.tsv.gz, survival_cdr.tsv.gz and SUMMARY.md."""
    got = fetch_raw()
    recs = pd.concat([_read_biotab(path, p) for p, (_, path) in got.items() if p in PROJECTS], ignore_index=True)

    # two passes: exact resolution, then fuzzy matching of misspellings against frequent resolved names
    uniq = recs["drug_raw"].unique()
    first = {r: normalize_drug(r) for r in uniq}
    freq = pd.Series([n for r in recs["drug_raw"] for n, ok in first[r] if ok]).value_counts()
    vocab = tuple(freq[freq >= 3].index)
    mapping = {r: normalize_drug(r, fuzzy_vocab=vocab) for r in uniq}

    rows = []
    for rec in recs.itertuples(index=False):
        for name, ok in mapping[rec.drug_raw]:
            d = rec._asdict()
            d.update(drug=name, drug_resolved=ok)
            rows.append(d)
    long = pd.DataFrame(rows)
    long["recist"] = long["best_response"].map(harmonise_recist)
    long["responder"] = long["recist"].map({"CR": 1, "PR": 1, "SD": 0, "PD": 0}).astype("Int64")
    long = long.drop_duplicates(["patient", "drug", "start_days", "end_days", "best_response", "regimen_number"])
    key = long["start_days"].astype(str)
    grp = long.groupby(["patient", key])["drug"].transform("nunique")
    # drugs without a start day only share a regimen with drugs from the same biotab record
    grp_rec = long.groupby("bcr_drug_barcode")["drug"].transform("nunique")
    long["n_drugs_in_regimen"] = np.where(long["start_days"].notna(), grp, grp_rec).astype(int)

    # cross-reference: Ding, Zu & Gu 2016 curated label for the same (patient, drug)
    ding = pd.read_csv(got["DING2016"][1], sep="\t", dtype=str, encoding="latin-1")
    ding["drug"] = [(normalize_drug(x, vocab) or [(str(x).upper(), False)])[0][0] for x in ding["drug.name"]]
    ding["patient"] = ding["patient.arr"].str[:12]
    ding = ding.drop_duplicates(["patient", "drug"]).set_index(["patient", "drug"])["response"]
    idx = pd.MultiIndex.from_frame(long[["patient", "drug"]])
    long["in_ding2016"] = idx.isin(ding.index).astype(int)

    # treatment-setting proxy (regimen_indication is never filled where a response is recorded):
    # post_progression = drug started on/after the CDR progression-free-interval event (measurable
    # disease likely, response ~ RECIST); early = started <= 180 d from diagnosis without prior PFI
    # event (adjuvant / neoadjuvant / first-line; CR often means "no recurrence"); late_no_event otherwise.
    cdr = _build_cdr(got["CDR"][1])
    pfi = long["patient"].map(cdr["PFI"])
    pfit = long["patient"].map(cdr["PFI.time"])
    post = (pfi == 1) & pfit.notna() & long["start_days"].notna() & (pfit <= long["start_days"] + 14)
    early = long["start_days"].notna() & (long["start_days"] <= 180) & ~post
    setting = np.select([post, early, long["start_days"].notna()], ["post_progression", "early", "late_no_event"],
                        default="unknown")
    long["setting_proxy"] = np.where(long["regimen_indication"].isin(["PROGRESSION", "RECURRENCE"]),
                                     "post_progression", setting)
    long.loc[long["regimen_indication"].eq("ADJUVANT"), "setting_proxy"] = "early"

    cols = ["project", "patient", "drug", "drug_raw", "drug_resolved", "therapy_type", "regimen_indication",
            "regimen_number", "start_days", "end_days", "therapy_ongoing", "cycles", "on_clinical_trial",
            "best_response", "recist", "responder", "n_drugs_in_regimen", "setting_proxy", "in_ding2016",
            "bcr_drug_barcode"]
    long = long[cols].sort_values(["project", "patient", "start_days", "drug"]).reset_index(drop=True)
    ROOT.mkdir(parents=True, exist_ok=True)
    long.to_csv(LONG_PATH, sep="\t", index=False, compression="gzip")

    cdr.to_csv(CDR_PATH, sep="\t", compression="gzip")
    _write_summary(long, cdr, ding, got)
    load_response_labels.cache_clear()
    load_cdr.cache_clear()
    return long


def _build_cdr(path):
    raw = pd.read_csv(path, sep="\t", dtype=str)
    raw = raw.drop_duplicates("_PATIENT").set_index("_PATIENT")
    stage = raw["ajcc_pathologic_tumor_stage"].map(_clean).fillna(raw["clinical_stage"].map(_clean))
    num = lambda c: pd.to_numeric(raw[c], errors="coerce")  # noqa: E731
    out = pd.DataFrame({
        "project": raw["cancer type abbreviation"], "age": num("age_at_initial_pathologic_diagnosis"),
        "gender": raw["gender"].map(_clean), "race": raw["race"].map(_clean), "stage": stage,
        "histology": raw["histological_type"].map(_clean), "grade": raw["histological_grade"].map(_clean),
        "vital_status": raw["vital_status"].map(_clean), "tumor_status": raw["tumor_status"].map(_clean),
        "treatment_outcome_first_course": raw["treatment_outcome_first_course"].map(_clean),
        "residual_tumor": raw["residual_tumor"].map(_clean),
        "OS": num("OS"), "OS.time": num("OS.time"), "DSS": num("DSS"), "DSS.time": num("DSS.time"),
        "DFI": num("DFI"), "DFI.time": num("DFI.time"), "PFI": num("PFI"), "PFI.time": num("PFI.time"),
        "redaction": raw["Redaction"].map(_clean),
    })
    out.index.name = "patient"
    return out


def _write_summary(long, cdr, ding, got):
    lab = long.dropna(subset=["recist"])
    summ = drug_label_matrix(min_patients=1, df=long).head(40)
    adj = lab["setting_proxy"].eq("early")
    post = lab["setting_proxy"].eq("post_progression")
    per_proj = pd.DataFrame({
        "drug_records": long.groupby("project").size(),
        "patients_treated": long.groupby("project")["patient"].nunique(),
        "labelled_records": lab.groupby("project").size(),
        "labelled_patients": lab.groupby("project")["patient"].nunique(),
        "responder_records": lab.groupby("project")["responder"].sum(),
        "cdr_patients": cdr.groupby("project").size(),
    }).fillna(0).astype(int)
    per_proj.loc["LAML"] = [0, 0, 0, 0, 0, int((cdr["project"] == "LAML").sum())]
    per_proj = per_proj.sort_index()
    # agreement with Ding et al. on shared (patient, drug) pairs
    m = lab.drop_duplicates(["patient", "drug"]).set_index(["patient", "drug"])["recist"]
    shared = m.index.intersection(ding.index)
    agree = (m.loc[shared].values == ding.loc[shared].map(harmonise_recist).values).mean() if len(shared) else np.nan

    def table(df):
        cols = [str(c) for c in df.columns]
        lines = ["| " + " | ".join([df.index.name or ""] + cols) + " |", "|" + "---|" * (len(cols) + 1)]
        for i, r in df.iterrows():
            lines.append("| " + " | ".join([str(i)] + [str(v) for v in r.values]) + " |")
        return "\n".join(lines)

    top = summ[["n_patients", "n_responders", "responder_rate", "n_patients_post_progression",
                "responder_rate_post_progression", "n_projects", "projects"]].copy()
    top["responder_rate"] = top["responder_rate"].round(2)
    top["responder_rate_post_progression"] = top["responder_rate_post_progression"].round(2)
    per_proj.index.name = "project"
    src = ["| file | url | sha256 |", "|---|---|---|"]
    for name, (url, path) in sorted(got.items()):
        src.append(f"| raw/{path.name} | {url} | `{_sha256(path)}` |")
    for path in (LONG_PATH, CDR_PATH):
        src.append(f"| {path.name} | derived (obd.tcga_response.build) | `{_sha256(path)}` |")
    text = f"""# TCGA pan-cancer drug-response labels

Built by `robust/obd/tcga_response.py::build()`; load with `load_response_labels()`, `load_cdr()`,
`drug_label_matrix()`.

* Drug records (rows, one per patient x drug x regimen record): **{len(long):,}**
  ({long['patient'].nunique():,} patients, {long['project'].nunique()} projects, {long['drug'].nunique()} distinct drug names,
  {long['drug_resolved'].mean():.1%} of rows resolved to a DrugBank/curated name)
* Rows with a RECIST-like label (CR/PR/SD/PD): **{len(lab):,}** ({lab['patient'].nunique():,} patients);
  responders (CR/PR) {int(lab['responder'].sum()):,} ({lab['responder'].mean():.1%});
  setting_proxy of labelled rows: {adj.mean():.1%} early (started <=180 d after diagnosis, no prior
  progression; adjuvant/neoadjuvant/first-line), {post.mean():.1%} post_progression (started after the CDR PFI event).
* `recist` harmonised from `treatment_best_response`: Complete Response->CR, Partial Response->PR,
  Stable Disease->SD, Clinical/Radiographic Progressive Disease->PD; `responder` = 1 for CR/PR, 0 for SD/PD.
* `n_drugs_in_regimen`: distinct drugs with the same patient and start day (same biotab record if start missing).
* `regimen_indication`: `regimen_indication` or `therapy_regimen` column (ADJUVANT, PROGRESSION, RECURRENCE,
  PALLIATIVE, PRIMARY, OTHER...); absent for {long['regimen_indication'].isna().mean():.0%} of rows and, in these
  biotab versions, NEVER present on a row that carries a response (the form versions with `therapy_regimen`
  dropped `treatment_best_response`). Use `setting_proxy` instead.
* `setting_proxy`: post_progression (start >= CDR PFI event time - 14 d, or indication PROGRESSION/RECURRENCE),
  early (start <= 180 d from diagnosis, or indication ADJUVANT), late_no_event, unknown (no start day).
* `in_ding2016`: (patient, drug) also present in the Ding, Zu & Gu 2016 curated table ({len(ding):,} pairs);
  RECIST agreement on {len(shared):,} shared labelled pairs: {agree:.1%}.
* LAML: TCGA released no clinical_drug biotab table (no labels). DLBC's table has no response field.
* survival_cdr.tsv.gz: TCGA-CDR (Liu et al. Cell 2018) for {len(cdr):,} patients.

## Top 40 drugs by patients with a RECIST label

Patient-level label = earliest labelled record of the drug for that patient. `post_progression` restricts to
records with setting_proxy == post_progression (treatment of measurable recurrent/progressive disease).

{table(top)}

## Per project

{table(per_proj)}

## SOURCE

Mirrors of the TCGA BCR biotab `nationwidechildrens.org_clinical_drug_<project>.txt` (GDC legacy archive,
data frozen ~2016; GDC API blocked from the build environment). kemplab/FBA-pipeline covers 31 projects,
nikcheerla/mirnanalyze provides ACC. Shicheng-Guo/HowtoBook `TCGA/drug_response/mutation/tcga.TCGA-*.drug.txt`
(TCGAbiolinks-prepared) was checked as a cross-mirror (same records, ~1% fewer) but not used.

{chr(10).join(src)}
"""
    SUMMARY_PATH.write_text(text)


# --------------------------------------------------------------------- loaders
@lru_cache(None)
def _long():
    df = pd.read_csv(LONG_PATH, sep="\t", low_memory=False)
    df["responder"] = df["responder"].astype("Int64")
    return df


def load_response_labels(drugs=None, projects=None, labelled_only=False):
    """Long table of TCGA drug records. `drugs`: names (any synonym/brand), `projects`: e.g. ['BLCA']."""
    df = _long()
    if drugs is not None:
        drugs = [drugs] if isinstance(drugs, str) else drugs
        want = {common_drug_name(EXTRA_BRANDS.get(str(d).upper(), d)) for d in drugs} | {str(d).upper() for d in drugs}
        df = df[df["drug"].isin(want)]
    if projects is not None:
        projects = [projects] if isinstance(projects, str) else projects
        df = df[df["project"].isin([p.upper().replace("TCGA-", "") for p in projects])]
    if labelled_only:
        df = df[df["recist"].notna()]
    return df.copy()


load_response_labels.cache_clear = _long.cache_clear


@lru_cache(None)
def _cdr():
    return pd.read_csv(CDR_PATH, sep="\t", index_col="patient")


def load_cdr():
    """TCGA-CDR per patient (OS, OS.time, DSS, DFI, PFI (+.time, days), stage, age, gender, histology...)."""
    return _cdr().copy()


load_cdr.cache_clear = _cdr.cache_clear


def patient_labels(drug, projects=None, how="first", settings=None):
    """One label per patient for `drug`: how='first' (earliest labelled record), 'best' (any CR/PR) or
    'worst'. `settings`: restrict to setting_proxy values, e.g. 'post_progression'. Returns DataFrame indexed by patient with project, recist, responder, start_days, n_records."""
    df = load_response_labels(drug, projects, labelled_only=True)
    if settings is not None:
        df = df[df["setting_proxy"].isin([settings] if isinstance(settings, str) else settings)]
    order = {"CR": 0, "PR": 1, "SD": 2, "PD": 3}
    df = df.assign(_rank=df["recist"].map(order))
    if how == "first":
        df = df.sort_values(["start_days", "_rank"], na_position="last")
    elif how == "best":
        df = df.sort_values("_rank")
    elif how == "worst":
        df = df.sort_values("_rank", ascending=False)
    else:
        raise ValueError(how)
    n = df.groupby("patient").size().rename("n_records")
    out = df.drop_duplicates("patient").set_index("patient")
    return out[["project", "drug", "recist", "responder", "start_days", "setting_proxy",
                "n_drugs_in_regimen"]].join(n)


def drug_label_matrix(min_patients=30, df=None):
    """Per-drug summary of patient-level labels (earliest labelled record per patient & drug):
    n_patients, n_responders, responder_rate, the same restricted to setting_proxy == post_progression,
    n_patients_early / responder_rate_early, n_projects,
    projects ('BLCA:40,...'). Only drugs with >= min_patients labelled patients."""
    df = (_long() if df is None else df)
    lab = df[df["recist"].notna()].sort_values(["start_days"], na_position="last")
    first = lab.drop_duplicates(["patient", "drug"])
    sub = {k: lab[lab["setting_proxy"].eq(k)].drop_duplicates(["patient", "drug"]) for k in ("post_progression", "early")}
    g = first.groupby("drug")
    out = pd.DataFrame({"n_patients": g.size(), "n_responders": g["responder"].sum().astype(int)})
    out["responder_rate"] = out["n_responders"] / out["n_patients"]
    for k, d in sub.items():
        gn = d.groupby("drug")
        out[f"n_patients_{k}"] = gn.size().reindex(out.index).fillna(0).astype(int)
        out[f"responder_rate_{k}"] = (gn["responder"].sum().reindex(out.index) / out[f"n_patients_{k}"]).astype(float)
    out["n_projects"] = g["project"].nunique()
    out["projects"] = g["project"].apply(
        lambda s: ",".join(f"{k}:{v}" for k, v in s.value_counts().items()))
    out.index.name = "drug"
    return out[out["n_patients"] >= min_patients].sort_values("n_patients", ascending=False)


if __name__ == "__main__":
    build()
