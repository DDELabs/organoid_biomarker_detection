"""Gene ID mapping, drug synonyms, drug targets and gene sets."""
import csv
from functools import lru_cache

import pandas as pd

from . import DATA


@lru_cache(None)
def uniprot_maps():
    """Return (gene symbol -> UniProt, UniProt -> canonical symbol)."""
    df = pd.read_csv(DATA / "uniprot_homoSapiens_multipleGeneName_20180802.tab", sep="\t")
    g2u, u2g = {}, {}
    for uniprot, names in zip(df["Entry"], df["Gene names"]):
        if pd.isnull(names):
            continue
        names = names.split()
        u2g[uniprot] = names[0]
        for g in names:
            g2u[g] = uniprot
    return g2u, u2g


def canonical_symbol(symbol):
    """Map any alias to the canonical UniProt gene name used by the original pipeline."""
    g2u, u2g = uniprot_maps()
    u = g2u.get(symbol)
    return u2g.get(u) if u else None


@lru_cache(None)
def ensembl_to_symbol():
    df = pd.read_csv(DATA / "2017_07_31_biomart_protein_coding_genes.txt", sep="\t")
    return dict(zip(df["Gene stable ID"], df["Gene name"]))


@lru_cache(None)
def drug_synonyms():
    """{synonym (upper, also without '-' / ' ') : DrugBank common name}."""
    out = {}
    with open(DATA / "drugbank_vocabulary.csv") as f:
        for line in csv.reader(f):
            if line[0].startswith("DrugBank"):
                continue
            name = line[2].upper()
            out[name] = name
            out[name.replace("-", "")] = name
            out[name.replace(" ", "")] = name
            for s in line[5].split(" | "):
                s = s.upper()
                if s:
                    out[s] = out[s.replace("-", "")] = out[s.replace(" ", "")] = name
    return out


# Brand names, abbreviations and common misspellings in TCGA clinical_drug files
# that the DrugBank vocabulary (synonyms only, no brands) does not resolve.
BRANDS = {
    "GEMZAR": "GEMCITABINE", "TEMODAR": "TEMOZOLOMIDE", "TEMODOR": "TEMOZOLOMIDE", "TEMOZOLAMIDE": "TEMOZOLOMIDE",
    "TAXOL": "PACLITAXEL", "TAXOTERE": "DOCETAXEL", "ABRAXANE": "PACLITAXEL", "CYTOXAN": "CYCLOPHOSPHAMIDE",
    "ADRIAMYCIN": "DOXORUBICIN", "DOXIL": "DOXORUBICIN", "ALIMTA": "PEMETREXED", "NAVELBINE": "VINORELBINE",
    "NEXAVAR": "SORAFENIB", "NAXAVAR": "SORAFENIB", "XELODA": "CAPECITABINE", "ARIMIDEX": "ANASTROZOLE",
    "FEMARA": "LETROZOLE", "AROMASIN": "EXEMESTANE", "HERCEPTIN": "TRASTUZUMAB", "AVASTIN": "BEVACIZUMAB",
    "ELOXATIN": "OXALIPLATIN", "CAMPTOSAR": "IRINOTECAN", "CPT11": "IRINOTECAN", "5FU": "FLUOROURACIL",
    "5FLUOROURACIL": "FLUOROURACIL", "PARAPLATIN": "CARBOPLATIN", "PLATINOL": "CISPLATIN", "CCNU": "LOMUSTINE",
    "BCNU": "CARMUSTINE", "GLIADEL": "CARMUSTINE", "TARCEVA": "ERLOTINIB", "IRESSA": "GEFITINIB",
    "HYCAMTIN": "TOPOTECAN", "ELLENCE": "EPIRUBICIN", "VEPESID": "ETOPOSIDE", "VP16": "ETOPOSIDE",
    "ERBITUX": "CETUXIMAB", "GEMCITABINEHCL": "GEMCITABINE", "GEMCITABINE HCL": "GEMCITABINE", "NOLVADEX": "TAMOXIFEN", "FASLODEX": "FULVESTRANT", "LEUCOVORINCALCIUM": "LEUCOVORIN",
}


def common_drug_name(name):
    name = str(name).upper().strip()
    for key in (name, name.replace("-", "").replace(" ", "")):
        if key in BRANDS:
            name = BRANDS[key]
            break
    syn = drug_synonyms()
    for key in (name, name.replace("-", ""), name.replace(" ", ""), name.replace("-", "").replace(" ", "")):
        if key in syn:
            return syn[key]
    return name


@lru_cache(None)
def drug_targets():
    """Curated targets of the organoid drugs ({drug: [symbols]}), plus DrugBank targets."""
    out = {}
    df = pd.read_csv(DATA / "drug_drugTarget.txt", sep="\t")
    for drugs, targets in zip(df["drugs"], df["targets"]):
        for d in drugs.split(","):
            out[common_drug_name(d)] = sorted(set(t.strip() for t in targets.split(",")))
    return out


@lru_cache(None)
def drugbank_targets():
    """All human DrugBank targets ({common drug name: [symbols]})."""
    _, u2g = uniprot_maps()
    ann = {}
    with open(DATA / "drugbank_vocabulary.csv") as f:
        for line in csv.reader(f):
            if not line[0].startswith("DrugBank"):
                ann[line[0]] = line[2].upper()
    out = {}
    with open(DATA / "all.csv") as f:
        for line in csv.reader(f):
            if line[0] == "ID" or line[11] != "Humans" and line[11] != "Human":
                continue
            gene = u2g.get(line[5])
            if not gene:
                continue
            for dbid in line[12].split("; "):
                if dbid in ann:
                    out.setdefault(ann[dbid], set()).add(gene)
    return {k: sorted(v) for k, v in out.items()}


@lru_cache(None)
def gene_sets(collections=("REACTOME",)):
    """Gene sets from the bundled MSigDB v6.1 file, symbols mapped to canonical names.

    collections: prefixes such as REACTOME, KEGG, BIOCARTA, PID, HALLMARK.
    """
    out = {}
    with open(DATA / "msigdb.v6.1.symbols.gmt.txt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].split("_")[0] not in collections:
                continue
            genes = {canonical_symbol(g) for g in parts[2:]}
            genes.discard(None)
            out[parts[0]] = sorted(genes)
    extra = DATA / "external" / "h.all.symbols.gmt"
    if "HALLMARK" in collections and extra.exists():
        for line in open(extra):
            parts = line.rstrip("\n").split("\t")
            genes = {canonical_symbol(g) for g in parts[2:]}
            genes.discard(None)
            out[parts[0]] = sorted(genes)
    return out
