#!/usr/bin/env python3
"""Build the curated pre-clinical pharmacogenomic sets for ATLAS.

Every set is written to data/curated/preclinical/<SET>/ as
  expression.tsv.gz  genes (HGNC symbol, protein-coding) x models, log2 scale
  response.tsv       long table: sample, drug, response, metric, direction (+ extra columns)
  SOURCE.md          provenance, sha256 of every raw input, derivation, ID matching
and summarised in CATALOG_PRECLINICAL.tsv.

Raw downloads are cached in --raw (default $ATLAS_RAW or ~/.cache/atlas_preclinical_raw);
nothing outside data/curated/preclinical is written.

    python3 data/curated/preclinical/build_preclinical.py            # all sets
    python3 data/curated/preclinical/build_preclinical.py bladder_lee2018 prism_repurposing
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = Path(os.environ.get("ATLAS_RAW", Path.home() / ".cache" / "atlas_preclinical_raw"))

FIG = "https://ndownloader.figshare.com/files/{}"
GEO_COUNTS = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file&acc={acc}&file={file}"
EPMC_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/supplementaryFiles"
SYNAPSE = "https://repo-prod.prod.sagebase.org/repo/v1"

# raw inputs: key -> (url, local file name)
SOURCES = {
    "gene_info": ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                  "Homo_sapiens.gene_info.gz"),
    # PRISM Repurposing 19Q4 secondary screen (Corsello et al. 2020, figshare 9393293)
    "prism_curves": (FIG.format(20237739), "secondary-screen-dose-response-curve-parameters.csv"),
    # DepMap Public 21Q4 (figshare 16924132): CCLE RNA-seq log2(TPM+1), protein-coding; sample_info
    "ccle_expr": (FIG.format(31315882), "CCLE_expression_21Q4.csv"),
    "ccle_info": (FIG.format(31316011), "sample_info_21Q4.csv"),
    # Lee et al. 2018 bladder PDOs: NCBI-recomputed RNA-seq TPM for GSE103990
    "lee_tpm": (GEO_COUNTS.format(acc="GSE103990", file="GSE103990_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
                "GSE103990_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
    "lee_gsm": ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103990&targ=gsm&form=text&view=brief",
                "GSE103990_gsm_brief.txt"),
    # CoderData 2.1.0 (figshare): bladder PDO, pancreatic PDO (Tiriac 2018), sarcoma PDO, CTRPv2
    "cd_bladder_exp": (FIG.format(53956220), "coderdata_bladderpdo_experiments.tsv.gz"),
    "cd_bladder_samples": (FIG.format(53956226), "coderdata_bladderpdo_samples.csv"),
    "cd_bladder_drugs": (FIG.format(53956214), "coderdata_bladderpdo_drugs.tsv.gz"),
    "cd_panc_exp": (FIG.format(53779676), "coderdata_pancpdo_experiments.tsv.gz"),
    "cd_panc_samples": (FIG.format(53779637), "coderdata_pancpdo_samples.csv"),
    "cd_panc_drugs": (FIG.format(53779610), "coderdata_pancpdo_drugs.tsv.gz"),
    "cd_panc_tx": (FIG.format(53779682), "coderdata_pancpdo_transcriptomics.csv.gz"),
    "cd_sarc_exp": (FIG.format(53956238), "coderdata_sarcpdo_experiments.tsv.gz"),
    "cd_sarc_samples": (FIG.format(53956244), "coderdata_sarcpdo_samples.csv"),
    "cd_sarc_drugs": (FIG.format(53956235), "coderdata_sarcpdo_drugs.tsv.gz"),
    "cd_sarc_tx": (FIG.format(53956247), "coderdata_sarcpdo_transcriptomics.csv.gz"),
    "cd_ctrp_exp": (FIG.format(53779400), "coderdata_ctrpv2_experiments.tsv.gz"),
    "cd_ctrp_samples": (FIG.format(53779412), "coderdata_ctrpv2_samples.csv"),
    "cd_ctrp_drugs": (FIG.format(53779397), "coderdata_ctrpv2_drugs.tsv.gz"),
    # CoderData 2.2.x full release (figshare article 29923646 v4, 2.8 GB zip): liver PDOs (Ji 2023)
    "cd22_zip": ("https://ndownloader.figshare.com/articles/29923646/versions/4",
                 "coderdata_v2.2_article29923646_v4.zip"),
    # Shi et al. 2022 Nat Commun pancreatic PDOs
    "shi_fpkm": ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194249/suppl/GSE194249_PDPCOs_FPKM.txt.gz",
                 "GSE194249_PDPCOs_FPKM.txt.gz"),
    "shi_supp": (EPMC_SUPP.format("PMC9023604"), "PMC9023604_supplementary.zip"),
    # Broutier et al. 2017 Nat Med liver tumouroids
    "broutier_supp": (EPMC_SUPP.format("PMC5722201"), "PMC5722201_supplementary.zip"),
}

DIRECTION = "lower = more sensitive"


# ----------------------------------------------------------------------------- utilities
def fetch(key):
    url, name = SOURCES[key]
    path = RAW / name
    if not path.exists() or path.stat().st_size == 0:
        RAW.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".part")
        subprocess.run(["curl", "-sSfL", "--retry", "4", "-m", "1800", "-o", str(tmp), url], check=True)
        tmp.rename(path)
    return path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def provenance(keys):
    rows = []
    for k in keys:
        p = fetch(k)
        rows.append(f"| `{SOURCES[k][1]}` | {SOURCES[k][0]} | `{sha256(p)}` |")
    return "| file | URL | sha256 |\n|---|---|---|\n" + "\n".join(rows)


_GENES = {}


def gene_tables():
    """NCBI gene_info: Entrez -> HGNC symbol, symbol/synonym -> symbol, Ensembl -> symbol (protein-coding)."""
    if _GENES:
        return _GENES
    gi = pd.read_csv(fetch("gene_info"), sep="\t", dtype=str)
    gi = gi[gi["type_of_gene"] == "protein-coding"].copy()
    gi["sym"] = np.where(gi["Symbol_from_nomenclature_authority"] != "-",
                         gi["Symbol_from_nomenclature_authority"], gi["Symbol"])
    entrez = dict(zip(gi["GeneID"], gi["sym"]))
    current = set(gi["sym"])
    syn = {}
    for s, alts in zip(gi["sym"], gi["Synonyms"]):
        for a in str(alts).split("|"):
            if a not in ("-", "") and a not in current:
                syn.setdefault(a.upper(), set()).add(s)
    symbol = {s.upper(): s for s in current}
    symbol.update({a: next(iter(v)) for a, v in syn.items() if len(v) == 1 and a not in symbol})
    ens = {}
    for s, x in zip(gi["sym"], gi["dbXrefs"]):
        for m in re.findall(r"Ensembl:(ENSG\d+)", str(x)):
            ens[m] = s
    _GENES.update(entrez=entrez, symbol=symbol, ensembl=ens)
    return _GENES


def collapse(expr, mapping):
    """Map index through `mapping` (drop unmapped); duplicated symbols -> row with the highest mean."""
    expr = expr.copy()
    expr.index = [mapping.get(str(i)) for i in expr.index]
    expr = expr[expr.index.notna()]
    mean = expr.mean(axis=1)
    expr = expr.assign(_m=mean.values).sort_values("_m", ascending=False)
    expr = expr[~expr.index.duplicated()].drop(columns="_m")
    expr.index.name = "gene"
    return expr.sort_index()


def map_symbols(expr):
    g = gene_tables()["symbol"]
    return collapse(expr, {i: g.get(str(i).upper()) for i in expr.index})


def write_set(name, expr, resp, source_md, catalog, extra_files=None, decimals=3):
    out = HERE / name
    out.mkdir(parents=True, exist_ok=True)
    expr = expr.loc[:, sorted(expr.columns)].astype(float).round(decimals)
    with open(out / "expression.tsv.gz", "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0, filename="") as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="") as f:      # mtime=0: byte-reproducible output
        expr.to_csv(f, sep="\t")
    resp = resp.copy()
    resp["direction"] = DIRECTION
    cols = ["sample", "drug", "response", "metric", "direction"]
    resp = resp[cols + [c for c in resp.columns if c not in cols]]
    resp = resp.sort_values(["metric", "drug", "sample"])
    resp.to_csv(out / "response.tsv", sep="\t", index=False, float_format="%.6g")
    for fname, df in (extra_files or {}).items():
        df.to_csv(out / fname, sep="\t", index=False)
    both = sorted(set(resp["sample"]) & set(expr.columns))
    stats = {
        "set": name,
        "n_expr_models": expr.shape[1],
        "n_resp_models": resp["sample"].nunique(),
        "n_both": len(both),
        "n_genes": expr.shape[0],
        "n_drugs": resp["drug"].nunique(),
        "metrics": ",".join(sorted(resp["metric"].unique())),
        "n_response_rows": len(resp),
    }
    stats.update(catalog)
    (out / "SOURCE.md").write_text(source_md.format(**stats).rstrip() + "\n")
    for f in out.iterdir():
        mb = f.stat().st_size / 1e6
        if mb > 90:
            raise RuntimeError(f"{f} is {mb:.1f} MB (> 90 MB limit)")
    print(f"[{name}] expr {expr.shape}, resp rows {len(resp)}, both {len(both)}, drugs {stats['n_drugs']}")
    return stats


# ----------------------------------------------------------------------------- drug names
_DB = {}


def drugbank():
    """DrugBank vocabulary (repo data/): InChIKey -> name, InChIKey skeleton -> name, synonym -> name."""
    if _DB:
        return _DB
    v = pd.read_csv(REPO / "data" / "drugbank_vocabulary.csv", dtype=str)
    ik, sk, syn = {}, {}, {}
    for name, key, syns in zip(v["Common name"], v["Standard InChI Key"], v["Synonyms"]):
        n = name.upper()
        if isinstance(key, str):
            ik.setdefault(key, n)
            sk.setdefault(key[:14], n)
        syn.setdefault(norm(name), n)
        for s in str(syns).split(" | "):
            if s and s != "nan":
                syn.setdefault(norm(s), n)
    _DB.update(ik=ik, sk=sk, syn=syn, names={norm(n) for n in v["Common name"]})
    return _DB


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# Code names / misspellings seen in the sources -> generic name.
ALIASES = {
    "5FLUOROURACIL": "FLUOROURACIL", "5FU": "FLUOROURACIL", "GEMCITIBINE": "GEMCITABINE",
    "GEM": "GEMCITABINE", "PTX": "PACLITAXEL", "OXA": "OXALIPLATIN", "IRI": "IRINOTECAN",
    "1OHP": "OXALIPLATIN", "DASATANIB": "DASATINIB", "CH5424802": "ALECTINIB",
    "PD0332991": "PALBOCICLIB", "BIRB0796": "DORAMAPIMOD", "NUTLIN3A": "NUTLIN-3A",
    "EMD1214063": "TEPOTINIB", "LGK974": "WNT-974", "AZD8931": "SAPITINIB", "PF05212384": "GEDATOLISIB",
    "JNJ42756493": "ERDAFITINIB", "PF477736": "PF-477736", "SIROLIMUS": "SIROLIMUS",
    "RAPAMYCIN": "SIROLIMUS", "MITOMYCINC": "MITOMYCIN", "MITOMYCIN": "MITOMYCIN", "FK866": "FK-866", "GDC0032": "TASELISIB", "HKI272": "NERATINIB", "ABT263": "NAVITOCLAX",
    "LBH589": "PANOBINOSTAT", "SAHA": "VORINOSTAT", "BI2536": "BI-2536", "SCH772984": "SCH-772984",
    "SN38": "SN-38", "7ETHYL10HYDROXYCAMPTOTHECIN": "SN-38", "AZD6738": "CERALASERTIB", "AT406XEVINAPANT": "XEVINAPANT", "AT406": "XEVINAPANT", "AG881": "VORASIDENIB", "LOXO101": "LAROTRECTINIB",
    "NOFINDER": "GSK-1838705A",  # CoderData synonym placeholder for GSK1838705A (PubChem 25113169)
    "AURORAAINHIBITORI": "AURORA A INHIBITOR I", "KU55933": "KU-55933", "OSI027": "OSI-027", "MK2206": "MK-2206",
    "PD173074": "PD-173074", "BIBR1532": "BIBR-1532", "LY2109761": "LY-2109761", "GSK126": "GSK-126",
    "PAC1": "PAC-1", "XAV939": "XAV-939", "PLX4720": "PLX-4720", "AZD8055": "AZD-8055", "JQ1": "JQ1",
    "TALAZOPARIB": "TALAZOPARIB", "RAD001": "EVEROLIMUS", "CHIR258": "DOVITINIB", "PHA739358": "DANUSERTIB",
}


SALT = r"\s+(\d?HCL|HYDROCHLORIDE|FREE BASE|L-\s*-?\s*TARTARIC ACID|TARTRATE|MESYLATE|SODIUM|MALEATE|2HCL)$"


def generic(name):
    """Best generic UPPERCASE name for a free-text drug name."""
    raw = re.sub(SALT, "", str(name).strip(), flags=re.I).strip()
    base = re.sub(r"\s*\(.*?\)\s*", " ", raw).strip()          # drop "(LBH589)" etc.
    base = re.sub(SALT, "", base, flags=re.I).strip()
    paren = re.findall(r"\((.*?)\)", raw)
    db = drugbank()["syn"]
    for cand in [raw, base] + paren:
        k = norm(cand)
        if k in ALIASES:
            return ALIASES[k]
    for cand in [raw, base]:
        k = norm(cand)
        if k in db:
            return db[k]
    for cand in paren:
        k = norm(cand)
        if k in db:
            return db[k]
    return base.upper()


JUNK = re.compile(r"^(schembl|chembl|akos|hms\d|ncgc|bdbm|q\d|sr-|ex-a|cs-|hy-|bcp|mfcd|en300|z\d|dtx|sdccg|ccg-|"
                  r"gtpl|nsc|d\d|unii|smr|ab\d|f\d|ac1|ks-|da-|db-|as-|sw\d|glxc|mls|cid|chebi|brd-|kbio|spectrum|"
                  r"tox21|cas-|bspbio|pubchem|zinc|nci|s\d|mcule|stk|vu\d|db\d|fg\d|bp-|ms-|orb\d|[0-9]+-[0-9]+-[0-9]$)")


def coderdata_names(drugs, override=None, brd=None):
    """improve_drug_id -> generic name.

    Order: explicit override; DrugBank common name via InChIKey (full key, then the 14-char skeleton) unless that
    name is an IUPAC-like string; CoderData's source-dataset name (the synonym it suffixes with '?'); any synonym
    known to DrugBank; else the shortest readable synonym. ALIASES is applied to the result.
    brd: optional {BRD-K######## (13 chars): name} used right after the InChIKey step (Broad compound IDs).
    """
    db = drugbank()
    override = override or {}
    out = {}
    for did, g in drugs.groupby("improve_drug_id"):
        if did in override:
            out[did] = override[did]
            continue
        syns = g["chem_name"].astype(str).tolist()
        key = g["InChIKey"].dropna().astype(str)
        name = None
        if len(key):
            k = key.iloc[0]
            name = db["ik"].get(k) or db["sk"].get(k[:14])
            if name and (len(name) > 30 or re.search(r"[\[\]{}(),]", name)):
                name = None
        if name is None and brd:
            for x in syns:
                if x.upper().startswith("BRD-") and x.upper()[:13] in brd:
                    name = brd[x.upper()[:13]]
                    break
        if name is None:
            src = [x[:-1].strip() for x in syns if x.endswith("?")]
            if src:
                name = generic(src[0])
        if name is None:
            for x in syns:
                n = norm(x)
                base = norm(re.sub(SALT, "", re.sub(r"\s*\(.*?\)\s*", " ", x).strip(), flags=re.I))
                if n in ALIASES or n in db["syn"]:
                    name = ALIASES.get(n) or db["syn"][n]
                    break
                if base in db["names"]:
                    name = db["syn"][base]
                    break
        if name is None:
            cands = [x for x in syns
                     if not JUNK.match(x) and len(x) <= 25 and not re.search(r"[\[\]\(\),;:=]|\s.*\s", x)]
            name = (sorted(cands, key=lambda x: (len(x) < 5, len(x)))[0] if cands else syns[0]).upper()
        out[did] = ALIASES.get(norm(name), name.upper())
    return out


# ----------------------------------------------------------------------------- sets
def build_prism():
    curves = pd.read_csv(fetch("prism_curves"), low_memory=False)
    info = pd.read_csv(fetch("ccle_info"))
    n0 = len(curves)
    curves = curves[curves["passed_str_profiling"].astype(str).str.upper() == "TRUE"]
    curves = curves[curves["auc"].notna() & curves["depmap_id"].notna() & curves["name"].notna()]
    # one value per line x compound: prefer MTS010 > MTS006 > MTS005 > HTS002 (DepMap secondary AUC convention)
    prio = {"MTS010": 0, "MTS006": 1, "MTS005": 2, "HTS002": 3}
    curves["prio"] = curves["screen_id"].map(prio).fillna(9)
    curves = curves.sort_values("prio").drop_duplicates(["depmap_id", "broad_id"])
    curves["drug"] = curves["name"].str.upper().str.strip()
    resp = (curves.groupby(["depmap_id", "drug"])
            .agg(response=("auc", "median"), screen_id=("screen_id", lambda s: ",".join(sorted(set(s)))),
                 broad_id=("broad_id", lambda s: ",".join(sorted(set(s)))), n_compounds=("broad_id", "nunique"))
            .reset_index().rename(columns={"depmap_id": "sample"}))
    resp["metric"] = "AUC"

    expr = pd.read_csv(fetch("ccle_expr"), index_col=0).T
    expr.index = [re.search(r"\((\d+)\)", c).group(1) for c in expr.index]
    keep = sorted(set(expr.columns) & set(resp["sample"]))
    expr = collapse(expr[keep], gene_tables()["entrez"])
    samples = (info[info["DepMap_ID"].isin(set(resp["sample"]))]
               .rename(columns={"DepMap_ID": "sample"}).sort_values("sample"))
    samples = samples[["sample", "CCLE_Name", "stripped_cell_line_name", "lineage", "lineage_subtype",
                       "primary_disease", "Subtype", "COSMICID", "Sanger_Model_ID"]]
    samples["has_expression"] = samples["sample"].isin(expr.columns)
    resp = resp.merge(samples[["sample", "lineage"]], on="sample", how="left")
    md = f"""# PRISM Repurposing secondary screen + CCLE RNA-seq

**Paper**: Corsello SM et al. *Discovering the anticancer potential of non-oncology drugs by systematic viability
profiling.* Nat Cancer 2020;1:235-248. doi:10.1038/s43018-019-0018-6. Expression: Ghandi M et al. Nature 2019
(CCLE) / DepMap Public 21Q4.

**Accessions / URLs**: PRISM Repurposing 19Q4 secondary screen (figshare article 9393293, file 20237739);
DepMap Public 21Q4 (figshare article 16924132: `CCLE_expression.csv` file 31315882, `sample_info.csv` file 31316011).
The DepMap portal itself sits behind a Cloudflare challenge here, so the 21Q4 figshare release
was used instead of the newer `OmicsExpressionProteinCodingGenesTPMLogp1.csv`.

{provenance(["prism_curves", "ccle_expr", "ccle_info", "gene_info"])}

## Derivation
* Response: `secondary-screen-dose-response-curve-parameters.csv` ({n0} curves) -> kept curves with
  `passed_str_profiling == TRUE` and a fitted `auc`. For a cell line x compound (broad_id) measured in several
  screens one curve is kept, preferring MTS010 > MTS006 > MTS005 > HTS002 (as in DepMap's secondary AUC matrix).
  Compound `name` upper-cased = `drug`; if several broad_ids share a name, the median AUC is reported
  (`n_compounds`, `broad_id`, `screen_id` columns record this).
* `metric` = AUC of the fitted 4-parameter log-logistic viability curve over the screened dose range
  (8 doses, 3-fold, max 10 uM; the curve is normalised so AUC = 1 is no effect; values can exceed 1).
  **Lower = more sensitive** (no conversion needed). `lineage` column copied from DepMap sample_info.
* Expression: `CCLE_expression.csv` (RSEM log2(TPM+1), protein-coding, columns `SYMBOL (ENTREZ)`) -> Entrez IDs
  mapped to current HGNC symbols (NCBI gene_info, protein-coding only), restricted to cell lines with PRISM
  secondary data. Units: **log2(TPM + 1)**.
* `samples.tsv`: DepMap annotations (CCLE name, lineage, subtype, COSMIC/Sanger IDs) for every PRISM line.

## Sample-ID matching
Both tables use DepMap IDs (ACH-xxxxxx). PRISM lines (STR-passed): {{n_resp_models}};
with CCLE RNA-seq: **{{n_both}}** (expression columns = {{n_expr_models}}). Genes: {{n_genes}}. Drugs: {{n_drugs}}.
"""
    return write_set("prism_repurposing", expr, resp, md,
                     {"type": "cell line", "tissue": "pan-cancer", "source": "figshare (Broad/DepMap)"},
                     extra_files={"samples.tsv": samples})


def synapse_drug_files(folder="syn64765430"):
    """Names of the per-drug files in the Lee bladder PDO Synapse folder (anonymous metadata read)."""
    import urllib.request
    names, token = [], None
    while True:
        body = {"parentId": folder, "includeTypes": ["file"]}
        if token:
            body["nextPageToken"] = token
        req = urllib.request.Request(SYNAPSE + "/entity/children", data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
        r = json.load(urllib.request.urlopen(req, timeout=120))
        names += [(p["id"], p["name"]) for p in r["page"]]
        token = r.get("nextPageToken")
        if not token:
            return names


def build_bladder():
    # ---- response (CoderData refit of the Lee 2018 raw dose-response data on Synapse)
    exp = pd.read_csv(fetch("cd_bladder_exp"), sep="\t")
    smp = pd.read_csv(fetch("cd_bladder_samples"))
    drugs = pd.read_csv(fetch("cd_bladder_drugs"), sep="\t")
    files = synapse_drug_files()
    (RAW / "syn64765430_children.json").write_text(json.dumps(files, indent=0))
    src_names = [re.match(r"\d+\) (.*)\.txt", n).group(1) for _, n in files if re.match(r"\d+\) .*\.txt", n)]
    used = set(exp["improve_drug_id"])
    drugs["n"] = drugs["chem_name"].map(norm)
    override, unmatched = {}, []
    extra = {"NUTLIN3A": {"NUTLIN3", "NUTLIN3A"}, "SIROLIMUS": {"RAPAMYCIN"}, "MITOMYCINC": {"MITOMYCIN"}}
    for nm in src_names:
        keys = {norm(nm), norm(re.sub(r"\(.*\)", "", nm))} | {norm(x) for x in re.findall(r"\((.*?)\)", nm)}
        keys |= extra.get(norm(nm), set())
        hit = drugs[drugs["n"].isin(keys) & drugs["improve_drug_id"].isin(used)]["improve_drug_id"].unique()
        if len(hit) == 1:
            override[hit[0]] = generic(nm)
        else:
            unmatched.append(nm)
    if unmatched or len(override) != len(used):
        raise RuntimeError(f"bladder drug mapping incomplete: {unmatched}")
    sid = smp.drop_duplicates("improve_sample_id").set_index("improve_sample_id")
    exp = exp[exp["dose_response_metric"] == "fit_auc"].copy()
    exp["screen"] = exp["improve_sample_id"].map(sid["other_id"])
    exp["line"] = exp["improve_sample_id"].map(sid["common_name"]).str.replace(".", "_", regex=False)
    # organoid screens only: '<line>_Organoid_P<n>' and '<line>_Parental[_k]' (parental organoid line);
    # xenograft-derived organoids ('XenoOrganoid', 'Xenograft') are excluded.
    is_org = exp["screen"].str.contains(r"_Organoid_P\d+$|_Parental(?:_\d+)?$", regex=True)
    n_all_screens = exp["screen"].nunique()
    exp = exp[is_org]
    exp["drug"] = exp["improve_drug_id"].map(override)
    resp = (exp.groupby(["line", "drug"])
            .agg(response=("dose_response_value", "median"), n_screens=("screen", "nunique"),
                 screens=("screen", lambda s: ",".join(sorted(s))))
            .reset_index().rename(columns={"line": "sample"}))
    resp["metric"] = "AUC"

    # ---- expression (GEO GSE103990, NCBI-recomputed TPM, organoid samples)
    txt = fetch("lee_gsm").read_text()
    acc = re.findall(r"!Sample_geo_accession = (\S+)", txt)
    tit = re.findall(r"!Sample_title = (.+)", txt)
    title = dict(zip(acc, [t.strip() for t in tit]))
    tpm = pd.read_csv(fetch("lee_tpm"), sep="\t", index_col=0)
    tpm.index = tpm.index.astype(str)
    org = {g: t for g, t in title.items() if "_org" in t and g in tpm.columns}
    lines = {g: re.sub(r"_orgP\d+$", "", t).replace(".", "_") for g, t in org.items()}
    lg = np.log2(tpm[list(org)] + 1)
    expr = lg.T.groupby(pd.Series(lines)).mean().T              # mean over passages (log scale)
    expr = collapse(expr, gene_tables()["entrez"])
    gsm_tab = pd.DataFrame({"gsm": list(org), "title": [org[g] for g in org], "sample": [lines[g] for g in org]})
    both = sorted(set(expr.columns) & set(resp["sample"]))
    md = f"""# Bladder cancer PDOs (Lee et al. 2018)

**Paper**: Lee SH, Hu W, Matulay JT, et al. *Tumor Evolution and Drug Response in Patient-Derived Organoid Models
of Bladder Cancer.* Cell 2018;173(2):515-528.e17. doi:10.1016/j.cell.2018.03.017 (PMID 29625057).

**Accessions**: RNA-seq GEO **GSE103990** (SRA SRP118077). Drug screen raw dose-response files: Synapse
folder **syn64765430** ("Lee Bladder PDO Datasets", deposited by CoderData), curve-refit and released in
**CoderData 2.1.0** (figshare files 53956220 experiments, 53956226 samples, 53956214 drugs).

{provenance(["lee_tpm", "lee_gsm", "cd_bladder_exp", "cd_bladder_samples", "cd_bladder_drugs", "gene_info"])}

## Access notes
* The author-supplied `GSE103990_Normalized_counts.txt.gz` (DESeq2 VST) answers HTTP 403 at GEO. NCBI's
  uniform re-quantification of the SRA runs (GRCh38.p13, `..._norm_counts_TPM_...`) was used instead.
* The Cell supplementary tables (PMC5890941) are not open access, and PMC is not reachable from the build host.
  Synapse file *contents* need a login, but folder metadata is public; the per-drug file names
  (`1) Gemcitabine.txt` ... `50) GSK126.txt`) were read anonymously to assign generic drug names to CoderData's
  `improve_drug_id`s (50/50 matched 1:1 by synonym).

## Derivation
* Response: CoderData `dose_response_metric == fit_auc`: the AUC of a fitted Hill curve on fraction viability over
  the tested dose range (6-day assay), 0-1 scale. **Lower = more sensitive**. CoderData labels screens as
  `<line>_Organoid_P<passage>`, `<line>_Parental`, `<line>_XenoOrganoid_P<n>` and `<line>_Xenograft`.
  Only organoid screens (`Organoid_P*`, and `Parental*`, which we take to be the parental organoid line before
  xenografting; CoderData does not define it) are kept ({n_all_screens} screens in total before
  filtering), and per line the **median** AUC over passages is reported (`n_screens`, `screens` columns).
  SCBO-3.2 / SCBO-11.2 etc. (organoids from recurrences) are separate lines written as `SCBO-3_2`, `SCBO-11_2`.
* Expression: NCBI TPM -> log2(TPM + 1); organoid samples only (GEO titles `SCBO-x_orgP<n>`, tumour tissue dropped);
  multiple passages of one line averaged on the log scale (`gsm_map.tsv` lists the GSMs per line). Entrez IDs ->
  HGNC symbols, protein-coding only. Units: **log2(TPM + 1)**.

## Sample-ID matching
Line IDs (`SCBO-5`, `SCBO-3_2`, ...) shared by both tables. Expression lines: {{n_expr_models}};
drug-screened lines: {{n_resp_models}}; **overlap n = {{n_both}}** ({", ".join(both)}).
Genes: {{n_genes}}; drugs: {{n_drugs}}. Note that the screened passage is not always the sequenced passage.
"""
    return write_set("bladder_lee2018", expr, resp, md,
                     {"type": "organoid", "tissue": "bladder", "source": "GEO + Synapse/CoderData"},
                     extra_files={"gsm_map.tsv": gsm_tab})


def build_pancreas_tiriac():
    exp = pd.read_csv(fetch("cd_panc_exp"), sep="\t")
    smp = pd.read_csv(fetch("cd_panc_samples"))
    drugs = pd.read_csv(fetch("cd_panc_drugs"), sep="\t")
    tx = pd.read_csv(fetch("cd_panc_tx"))
    name = smp[smp["other_id_source"] == "experimentId"].drop_duplicates("improve_sample_id") \
        .set_index("improve_sample_id")["other_id"]
    dn = coderdata_names(drugs)
    exp = exp[exp["dose_response_metric"] == "fit_auc"].copy()
    exp["sample"] = exp["improve_sample_id"].map(name)
    exp["drug"] = exp["improve_drug_id"].map(dn)
    resp = (exp.groupby(["sample", "drug"]).agg(response=("dose_response_value", "median"),
                                                n=("dose_response_value", "size")).reset_index())
    resp["metric"] = "AUC"
    mat = tx.pivot_table(index="entrez_id", columns="improve_sample_id", values="transcriptomics", aggfunc="mean")
    mat.index = mat.index.astype("int64").astype(str)
    mat = mat[[c for c in mat.columns if c in name.index]]
    mat.columns = [name[c] for c in mat.columns]
    mat = mat.T.groupby(level=0).mean().T
    expr = collapse(np.log2(mat + 1), gene_tables()["entrez"])
    md = f"""# Pancreatic cancer PDOs (Tiriac et al. 2018; HCMI)

**Paper**: Tiriac H, Belleau P, Engle DD, et al. *Organoid Profiling Identifies Common Responders to Chemotherapy
in Pancreatic Cancer.* Cancer Discov 2018;8(9):1112-1129. doi:10.1158/2159-8290.CD-18-0349.

**Accessions**: dose-response from Tiriac 2018 (AACR figshare 39996295, refit by CoderData); RNA-seq of the same
organoids from the NCI Human Cancer Models Initiative (GDC, HCMI). Both harmonised in **CoderData 2.1.0**
(`pancpdo_*`, figshare files 53779676, 53779637, 53779610, 53779682).

{provenance(["cd_panc_exp", "cd_panc_samples", "cd_panc_drugs", "cd_panc_tx", "gene_info"])}

## Derivation
* Response: CoderData `fit_auc` (AUC of a fitted Hill curve on fraction viability, 0-1, 120 h assay).
  **Lower = more sensitive**. Drug names from DrugBank via InChIKey (SN-38 kept as SN-38, the active irinotecan
  metabolite; `1-ohp` = oxaliplatin).
* Expression: CoderData transcriptomics (GDC STAR-count TPM of the HCMI organoid RNA-seq) -> log2(TPM + 1);
  Entrez -> HGNC symbols, protein-coding only. Units: **log2(TPM + 1)**.

## Sample-ID matching
Models are named by the Tiriac organoid IDs (`hF32`, `hM1A`, ...), taken from CoderData's `experimentId` other_id.
Drug-screened organoids: {{n_resp_models}}; with RNA-seq: {{n_expr_models}}; **overlap n = {{n_both}}**.
Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    return write_set("pancreas_tiriac2018", expr, resp, md,
                     {"type": "organoid", "tissue": "pancreas", "source": "CoderData (AACR figshare + GDC/HCMI)"})


def build_sarcoma():
    exp = pd.read_csv(fetch("cd_sarc_exp"), sep="\t")
    smp = pd.read_csv(fetch("cd_sarc_samples"))
    drugs = pd.read_csv(fetch("cd_sarc_drugs"), sep="\t")
    tx = pd.read_csv(fetch("cd_sarc_tx"))
    s = smp.drop_duplicates("improve_sample_id").set_index("improve_sample_id")
    dn = coderdata_names(drugs)
    exp = exp.copy()
    exp["sample"] = exp["improve_sample_id"].map(s["common_name"])
    exp["drug"] = exp["improve_drug_id"].map(dn)
    resp = exp.groupby(["sample", "drug"]).agg(response=("dose_response_value", "median")).reset_index()
    resp["metric"] = "viability_pct"
    org = s[s["model_type"] == "organoid"]
    mat = tx[tx["improve_sample_id"].isin(org.index)].pivot_table(
        index="entrez_id", columns="improve_sample_id", values="transcriptomics", aggfunc="mean")
    mat.index = mat.index.astype("int64").astype(str)
    mat.columns = [org.loc[c, "common_name"] for c in mat.columns]
    expr = collapse(np.log2(mat + 1), gene_tables()["entrez"])
    ann = s.reset_index()[["common_name", "cancer_type"]].drop_duplicates("common_name") \
        .rename(columns={"common_name": "sample"})
    md = f"""# Sarcoma PDOs (Al Shihabi et al. 2024)

**Paper**: Al Shihabi A, Tebon PJ, Nguyen HTL, et al. *The landscape of drug sensitivity and resistance in
sarcoma.* Cell Stem Cell 2024;31(10):1524-1542.e4. doi:10.1016/j.stem.2024.08.010.

**Accessions**: Synapse syn61892224 (drug table), syn61894695/syn61894699 (omics); harmonised in
**CoderData 2.1.0** (`sarcpdo_*`, figshare files 53956238, 53956244, 53956235, 53956247).

{provenance(["cd_sarc_exp", "cd_sarc_samples", "cd_sarc_drugs", "cd_sarc_tx", "gene_info"])}

## Derivation and caveats
* Response: the published per-drug **viability score** (`Viability_Score` in syn61892224; CoderData calls it
  `published_auc` but notes that full dose-response data were not released). It is % viability relative to
  vehicle, so `metric = viability_pct`, *not* an AUC. Higher = more viable, so **lower = more sensitive** without
  conversion. Values > 100 occur (growth above control).
  CoderData attaches these scores to the `_Tumor` sample IDs, but the screens were run on the PDOs of the
  same patient. They are keyed here by patient/model ID (`SARC0065`, `SARC0139_1`, ...).
* Expression: CoderData transcriptomics (TPM) of the **organoid** samples only -> log2(TPM + 1); Entrez -> HGNC
  symbols, protein-coding. Units: **log2(TPM + 1)**. Histology per model in `samples.tsv`.
* Drug names: DrugBank via InChIKey, then synonyms.

## Sample-ID matching
Organoid RNA-seq models: {{n_expr_models}}; screened models: {{n_resp_models}}; **overlap n = {{n_both}}**.
Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    return write_set("sarcoma_alshihabi2024", expr, resp, md,
                     {"type": "organoid", "tissue": "sarcoma (mixed histologies)", "source": "Synapse/CoderData"},
                     extra_files={"samples.tsv": ann})


def build_shi():
    zf = zipfile.ZipFile(fetch("shi_supp"))
    auc = pd.read_excel(io.BytesIO(zf.read("41467_2022_29857_MOESM11_ESM.xlsx")), header=1)
    cat = pd.read_excel(io.BytesIO(zf.read("41467_2022_29857_MOESM10_ESM.xlsx")), sheet_name=None, header=1)
    names = {}
    for sheet in cat.values():
        for c, p in zip(sheet["Catalog Number"], sheet["Product Name"]):
            if isinstance(c, str):
                names[c.strip()] = generic(p)
    chemo = {"GEM": "GEMCITABINE", "5-FU": "FLUOROURACIL", "PTX": "PACLITAXEL", "OXA": "OXALIPLATIN",
             "IRI": "IRINOTECAN"}
    auc = auc.rename(columns={"Drug list": "code"}).dropna(subset=["code"])
    auc["code"] = auc["code"].astype(str).str.strip()
    missing = [c for c in auc["code"] if c not in names and c not in chemo]
    if missing:
        raise RuntimeError(f"Shi drug codes without a name: {missing}")
    long = auc.melt(id_vars="code", var_name="sample", value_name="response").dropna()
    long["drug"] = long["code"].map(lambda c: chemo.get(c, names.get(c)))
    long["catalog"] = long["code"]
    long["metric"] = "AUC"
    resp = (long.groupby(["sample", "drug", "metric"])
            .agg(response=("response", "median"), catalog=("catalog", lambda c: ",".join(sorted(c))))
            .reset_index())
    fp = pd.read_csv(fetch("shi_fpkm"), sep="\t", index_col=0)
    fp = fp.T.groupby(level=0).mean().T if fp.columns.duplicated().any() else fp
    expr = map_symbols(np.log2(fp.astype(float) + 1))
    md = f"""# Pancreatic cancer PDOs (Shi et al. 2022)

**Paper**: Shi X, Li Y, Yuan Q, et al. *Integrated profiling of human pancreatic cancer organoids reveals chromatin
accessibility features associated with drug sensitivity.* Nat Commun 2022;13:2169.
doi:10.1038/s41467-022-29857-6 (PMC9023604, open access).

**Accessions**: RNA-seq GEO **GSE194249** (`GSE194249_PDPCOs_FPKM.txt.gz`); drug screen = Supplementary Data 8
(normalised AUC, 39 exocrine PDPCOs x 59 compounds + 5 chemotherapies) with compound list in Supplementary Data 7,
fetched as the Europe PMC supplementary bundle.

{provenance(["shi_fpkm", "shi_supp"])}

## Derivation
* Response: Supplementary Data 8 "Normalized AUC" (area under the viability curve, normalised to the dose range;
  0-1). **Lower = more sensitive** (no conversion). Rows are Selleck catalogue numbers mapped to names with
  Supplementary Data 7 (`catalog` column kept). Chemotherapy abbreviations: GEM gemcitabine, 5-FU fluorouracil,
  PTX paclitaxel, OXA oxaliplatin, IRI irinotecan.
* Expression: FPKM -> log2(FPKM + 1); GEO gene symbols -> current HGNC symbols (NCBI gene_info symbol or unique
  synonym), protein-coding only; duplicates -> highest-mean row. Units: **log2(FPKM + 1)**.

## Sample-ID matching
Organoid IDs (`CAS-DAC-1`, `CAS-IPMN-1`, `CAS-NEN-*` ...) identical in both tables. Expression models:
{{n_expr_models}} (includes neuroendocrine/normal lines without drug data); screened: {{n_resp_models}};
**overlap n = {{n_both}}**. Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    return write_set("pancreas_shi2022", expr, resp, md,
                     {"type": "organoid", "tissue": "pancreas", "source": "GEO + Nat Commun supplement"})


def build_broutier():
    zf = zipfile.ZipFile(fetch("broutier_supp"))
    raw = pd.read_excel(io.BytesIO(zf.read("NIHMS74480-supplement-Supplementary_Dataset_5.xlsx")),
                        sheet_name="3_Raw data")
    raw = raw.rename(columns=lambda c: str(c).strip())
    raw["drug"] = raw["Drug"].map(generic)
    raw["sample"] = raw["Organoid"].str.replace("-", "", regex=False)
    rows = []
    for metric, col in (("lnIC50", "IC50"), ("AUC", "AUC")):
        g = raw.groupby(["sample", "drug"]).agg(response=(col, "median"), n=(col, "size"),
                                                 max_conc_uM=("maxc", "first"), gdsc_drug_id=("DRUG_ID", "first"))
        g = g.reset_index()
        g["metric"] = metric
        rows.append(g)
    resp = pd.concat(rows, ignore_index=True)
    ex = pd.read_excel(io.BytesIO(zf.read("NIHMS74480-supplement-Supplementary_Dataset_1.xlsx")),
                       sheet_name="S1_RPKM_values", header=2)
    ex = ex.set_index("Ensembl ID")
    tum = [c for c in ex.columns if re.match(r"^(HCC|CHC|CC)\d+_O[a-c]?$", str(c))]
    lines = {c: re.sub(r"_O[a-c]?$", "", c) for c in tum}
    lg = np.log2(ex[tum].astype(float) + 1)
    lg = lg.T.groupby(pd.Series(lines)).mean().T
    lg.index = lg.index.astype(str).str.replace(r"\..*$", "", regex=True)
    expr = collapse(lg, gene_tables()["ensembl"])
    md = f"""# Primary liver cancer tumouroids (Broutier et al. 2017)

**Paper**: Broutier L, Mastrogiovanni G, Verstegen MM, et al. *Human primary liver cancer-derived organoid cultures
for disease modeling and drug screening.* Nat Med 2017;23(12):1424-1435. doi:10.1038/nm.4438 (PMC5722201).

**Accessions**: RNA-seq GEO GSE84073 (here taken from the authors' processed RPKM table, Supplementary Dataset 1,
sheet `S1_RPKM_values`); drug screen = Supplementary Dataset 5, sheet `3_Raw data` (29 compounds, 5 tumouroid lines,
4 replicate fits each; GDSC pipeline). Both from the Europe PMC supplementary bundle.

{provenance(["broutier_supp", "gene_info"])}

## Derivation
* Response: two metrics per line x drug, median over the 4 replicate fits (`n`):
  `lnIC50` = column `IC50` (natural log of IC50 in uM; the sheet defines IC50 (uM) = EXP(IC50)) and
  `AUC` = fitted-curve AUC (fraction viability, 0-1). **Lower = more sensitive** for both. IC50s extrapolated beyond
  `max_conc_uM` are kept as published. Drug names normalised (Gemcitibine -> GEMCITABINE, Dasatanib -> DASATINIB,
  CH5424802 -> ALECTINIB, PD-0332991 -> PALBOCICLIB, BIRB 0796 -> DORAMAPIMOD, AZD8931 -> SAPITINIB,
  EMD-1214063 -> TEPOTINIB, LGK974 -> WNT-974).
* Expression: RPKM -> log2(RPKM + 1); tumouroid columns only (`HCC1_O`, `HCC3_O`, `CHC1_O[a,b]`, `CHC2_O`,
  `CC1_O[a-c]`, `CC2_O`, `CC3_O`); technical/clone replicates averaged on the log scale. Ensembl IDs -> HGNC symbols,
  protein-coding. Units: **log2(RPKM + 1)**.

## Sample-ID matching
Drug-sheet IDs `HCC-1` <-> expression `HCC1_O` (hyphen removed -> `HCC1`). Expression lines: {{n_expr_models}};
screened lines: {{n_resp_models}}; **overlap n = {{n_both}}**. Genes: {{n_genes}}; drugs: {{n_drugs}}.
Very small n: useful for sanity checks / pooled analyses only.
"""
    return write_set("liver_broutier2017", expr, resp, md,
                     {"type": "organoid", "tissue": "liver (HCC, CHC, CC)", "source": "Nat Med supplement"})


def build_ctrp():
    exp = pd.read_csv(fetch("cd_ctrp_exp"), sep="\t")
    smp = pd.read_csv(fetch("cd_ctrp_samples"))
    drugs = pd.read_csv(fetch("cd_ctrp_drugs"), sep="\t")
    info = pd.read_csv(fetch("ccle_info"))
    pr = pd.read_csv(fetch("prism_curves"), usecols=["broad_id", "name"]).dropna().drop_duplicates()
    brd = {b[:13]: n.upper() for b, n in zip(pr["broad_id"], pr["name"])}
    ach = smp[smp["other_id_source"] == "DepMap"].drop_duplicates("improve_sample_id") \
        .set_index("improve_sample_id")["other_id"]
    dn = coderdata_names(drugs, brd=brd)
    exp = exp[exp["dose_response_metric"] == "fit_auc"].copy()
    exp["sample"] = exp["improve_sample_id"].map(ach)
    exp = exp[exp["sample"].notna()]
    exp["drug"] = exp["improve_drug_id"].map(dn)
    resp = (exp.groupby(["sample", "drug"]).agg(response=("dose_response_value", "median"),
                                                n_compounds=("improve_drug_id", "nunique")).reset_index())
    resp["metric"] = "AUC"
    expr = pd.read_csv(fetch("ccle_expr"), index_col=0).T
    expr.index = [re.search(r"\((\d+)\)", c).group(1) for c in expr.index]
    keep = sorted(set(expr.columns) & set(resp["sample"]))
    expr = collapse(expr[keep], gene_tables()["entrez"])
    samples = info[info["DepMap_ID"].isin(set(resp["sample"]))][
        ["DepMap_ID", "CCLE_Name", "stripped_cell_line_name", "lineage", "lineage_subtype", "primary_disease",
         "COSMICID", "Sanger_Model_ID"]].rename(columns={"DepMap_ID": "sample"})
    resp = resp.merge(samples[["sample", "lineage"]], on="sample", how="left")
    dmap = pd.DataFrame(sorted(dn.items()), columns=["improve_drug_id", "drug"])
    md = f"""# CTRPv2 (CoderData refit) + CCLE RNA-seq

**Papers**: Seashore-Ludlow B et al. Cancer Discov 2015;5:1210 and Rees MG et al. Nat Chem Biol 2016;12:109 (CTRPv2);
curves refit with PharmacoGx and harmonised by CoderData 2.1.0. Expression: DepMap Public 21Q4 CCLE RNA-seq.

**Overlap with the repo**: `robust/obd/preclinical_sets.py` already pairs CTRPv2 AUC with **GDSC RMA microarray**
expression keyed by COSMIC ID (data/external/gdsc). This set is complementary: CoderData's uniform refit
(`fit_auc` on a 0-1 scale rather than CTRP's unnormalised area) paired with **CCLE RNA-seq** keyed by DepMap ID.
Use one or the other, not both, in a pooled analysis. The original CTD2 portal
(ctd2-data.nci.nih.gov) is blocked by the build host's egress policy.

{provenance(["cd_ctrp_exp", "cd_ctrp_samples", "cd_ctrp_drugs", "ccle_expr", "ccle_info", "gene_info"])}

## Derivation
* Response: CoderData `fit_auc` (Hill-curve AUC on fraction viability over the tested range, 0-1, 72 h).
  **Lower = more sensitive**. Drug names: DrugBank via InChIKey; else the PRISM Repurposing name of the same
  Broad compound (BRD-K ID among CoderData's synonyms); else CoderData's source name / a DrugBank synonym / the
  shortest readable CoderData synonym; a few tool compounds keep IUPAC-like or BRD-ID names (`drug_map.tsv` lists improve_drug_id -> name). If two improve_drug_ids share a name, the
  median is used (`n_compounds`).
* Expression: as for `prism_repurposing` (CCLE 21Q4 log2(TPM + 1), protein-coding, Entrez -> HGNC), restricted to
  CTRPv2 lines.

## Sample-ID matching
CoderData sample -> DepMap ID (other_id_source == DepMap). Screened lines with a DepMap ID: {{n_resp_models}};
with RNA-seq: **{{n_both}}**. Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    return write_set("ctrpv2_ccle", expr, resp, md,
                     {"type": "cell line", "tissue": "pan-cancer", "source": "CoderData (CTRPv2) + figshare (DepMap)"},
                     extra_files={"samples.tsv": samples, "drug_map.tsv": dmap})


def build_liver_ji():
    zf = zipfile.ZipFile(fetch("cd22_zip"))
    rd = lambda n, **kw: pd.read_csv(io.BytesIO(zf.read(n)), compression="gzip" if n.endswith(".gz") else None, **kw)
    exp = rd("liver_experiments.tsv.gz", sep="\t")
    smp = rd("liver_samples.csv")
    drugs = rd("liver_drugs.tsv.gz", sep="\t")
    tx = rd("liver_transcriptomics.csv.gz")
    member_sha = {n: hashlib.sha256(zf.read(n)).hexdigest() for n in
                  ("liver_experiments.tsv.gz", "liver_samples.csv", "liver_drugs.tsv.gz", "liver_transcriptomics.csv.gz")}
    s = smp.drop_duplicates("improve_sample_id").set_index("improve_sample_id")
    dn = coderdata_names(drugs)
    e = exp.pivot_table(index=["improve_sample_id", "improve_drug_id"], columns="dose_response_metric",
                        values="dose_response_value", aggfunc="first").reset_index()
    e["sample"] = e["improve_sample_id"].map(s["common_name"])
    e["drug"] = e["improve_drug_id"].map(dn)
    resp = (e.groupby(["sample", "drug"])
            .agg(response=("fit_auc", "median"), fit_r2=("fit_r2", "median"), n_compounds=("improve_drug_id", "nunique"))
            .reset_index())
    resp["metric"] = "AUC"
    mat = tx.pivot_table(index="entrez_id", columns="improve_sample_id", values="transcriptomics", aggfunc="mean")
    mat.index = mat.index.astype("int64").astype(str)
    mat.columns = [s.loc[c, "common_name"] for c in mat.columns]
    expr = collapse(np.log2(mat + 1), gene_tables()["entrez"])
    ann = s.reset_index()[["common_name", "cancer_type"]].rename(columns={"common_name": "sample"})
    members = "\n".join(f"| `{k}` (zip member) | | `{v}` |" for k, v in member_sha.items())
    md = f"""# Primary liver cancer PDOs (Ji et al. 2023)

**Paper**: Ji S, Feng L, Fu Z, et al. *Pharmaco-proteogenomic characterization of liver cancer organoids for
precision oncology.* Sci Transl Med 2023;15(706):eadg3358. doi:10.1126/scitranslmed.adg3358 (PMID 37494474).

**Accessions**: data deposited on Synapse by CoderData (syn66401300-syn66401303, syn66593307; login needed for
file contents), harmonised in **CoderData 2.2.x** and taken from the full release zip of figshare article
29923646 (version 4). Individual file IDs of that release are not listable here (api.figshare.com is blocked by the
egress policy), but the whole-article download works.

{{provenance}}
{members}

## Derivation
* Response: CoderData `fit_auc` (AUC of a fitted Hill curve on fraction viability over the tested range, 0-1,
  72 h). **Lower = more sensitive**. The fit R^2 is kept as `fit_r2` (median 0.68; filter on it if needed).
  Drug names: DrugBank via InChIKey, else CoderData's source name / synonyms.
* Expression: CoderData transcriptomics (Synapse RNA-seq table, TPM-scale: ~1.07e6 per sample over all genes) ->
  log2(x + 1); Entrez -> HGNC symbols, protein-coding only. Units: **log2(TPM + 1)** (TPM-scale input).
  Histology per model (HCC, ICC, combined HCC-CC, hepatoblastoma) in `samples.tsv`.

## Sample-ID matching
Organoid IDs (`HCCO4`, `ICCO1`, `CHCO1`, `HBO1` ...) identical in both tables. Expression models: {{n_expr_models}};
screened: {{n_resp_models}}; **overlap n = {{n_both}}**. Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    md = md.replace("{provenance}", provenance(["cd22_zip", "gene_info"]))
    return write_set("liver_ji2023", expr, resp, md,
                     {"type": "organoid", "tissue": "liver (HCC, ICC, CHC, HB)", "source": "Synapse/CoderData 2.2"},
                     extra_files={"samples.tsv": ann})


BUILDERS = {
    "bladder_lee2018": build_bladder,
    "prism_repurposing": build_prism,
    "pancreas_tiriac2018": build_pancreas_tiriac,
    "pancreas_shi2022": build_shi,
    "liver_broutier2017": build_broutier,
    "liver_ji2023": build_liver_ji,
    "sarcoma_alshihabi2024": build_sarcoma,
    "ctrpv2_ccle": build_ctrp,
}


def write_catalog():
    rows = []
    for d in sorted(HERE.iterdir()):
        if not (d / "response.tsv").exists():
            continue
        r = pd.read_csv(d / "response.tsv", sep="\t", usecols=["sample", "drug", "metric"])
        e = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0, nrows=1)
        both = sorted(set(r["sample"]) & set(e.columns))
        rows.append(dict(set=d.name, n_models_both=len(both), n_expr_models=e.shape[1],
                         n_resp_models=r["sample"].nunique(), n_drugs=r["drug"].nunique(),
                         metrics=",".join(sorted(r["metric"].unique())), **META.get(d.name, {})))
    cat = pd.DataFrame(rows)
    cat.to_csv(HERE / "CATALOG_PRECLINICAL.tsv", sep="\t", index=False)
    print(cat.to_string(index=False))


META = {
    "bladder_lee2018": dict(model="organoid", tissue="bladder", expression_units="log2(TPM+1)",
                            reference="Lee 2018 Cell 173:515", accession="GSE103990; syn64765430; CoderData 2.1.0"),
    "prism_repurposing": dict(model="cell line", tissue="pan-cancer", expression_units="log2(TPM+1)",
                              reference="Corsello 2020 Nat Cancer 1:235", accession="figshare 20237739; DepMap 21Q4"),
    "pancreas_tiriac2018": dict(model="organoid", tissue="pancreas", expression_units="log2(TPM+1)",
                                reference="Tiriac 2018 Cancer Discov 8:1112", accession="CoderData 2.1.0 pancpdo; HCMI"),
    "pancreas_shi2022": dict(model="organoid", tissue="pancreas", expression_units="log2(FPKM+1)",
                             reference="Shi 2022 Nat Commun 13:2169", accession="GSE194249; PMC9023604"),
    "liver_broutier2017": dict(model="organoid", tissue="liver", expression_units="log2(RPKM+1)",
                               reference="Broutier 2017 Nat Med 23:1424", accession="GSE84073; PMC5722201"),
    "liver_ji2023": dict(model="organoid", tissue="liver", expression_units="log2(TPM+1)",
                         reference="Ji 2023 Sci Transl Med 15:eadg3358", accession="CoderData 2.2 (figshare 29923646)"),
    "sarcoma_alshihabi2024": dict(model="organoid", tissue="sarcoma", expression_units="log2(TPM+1)",
                                  reference="Al Shihabi 2024 Cell Stem Cell 31:1524",
                                  accession="syn61892224; CoderData 2.1.0 sarcpdo"),
    "ctrpv2_ccle": dict(model="cell line", tissue="pan-cancer", expression_units="log2(TPM+1)",
                        reference="Seashore-Ludlow 2015 Cancer Discov 5:1210",
                        accession="CoderData 2.1.0 ctrpv2; DepMap 21Q4"),
}


def main():
    global RAW
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="*", help=f"subset of {list(BUILDERS)}")
    ap.add_argument("--raw", default=str(RAW), help="download cache directory")
    a = ap.parse_args()
    RAW = Path(a.raw)
    for s in a.sets or BUILDERS:
        BUILDERS[s]()
    write_catalog()


if __name__ == "__main__":
    sys.exit(main())
