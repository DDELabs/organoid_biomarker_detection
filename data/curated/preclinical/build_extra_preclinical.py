#!/usr/bin/env python3
"""Build additional curated pre-clinical sets (patient-derived xenografts / organoids / primary cultures).

Re-uses the helpers of build_preclinical.py (download cache, HGNC mapping, drug-name harmonisation,
byte-reproducible writer) and writes, per set, data/curated/preclinical/<SET>/
  expression.tsv.gz      genes (HGNC symbol, protein-coding) x models, log2 scale
  response.tsv           sample, drug, response, metric, direction, n_screens, screens (+ extra columns)
  patient_response.tsv   (only when the matched patient's clinical response is public)
  mrecist.tsv            (PDXE only)
  SOURCE.md
then regenerates CATALOG_PRECLINICAL.tsv for all sets.

    python3 data/curated/preclinical/build_extra_preclinical.py                 # all extra sets
    python3 data/curated/preclinical/build_extra_preclinical.py pdx_novartis_gao2015 --raw /tmp/cache
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_preclinical as bp  # noqa: E402

SPRINGER = "https://static-content.springer.com/esm/{}"
bp.SOURCES.update({
    # Gao et al. 2015 Nat Med 21:1318 (Novartis PDX Encyclopedia) - Supplementary Table S1 (nm.3954-S2.xlsx)
    "pdxe_xlsx": (SPRINGER.format("art%3A10.1038%2Fnm.3954/MediaObjects/41591_2015_BFnm3954_MOESM10_ESM.xlsx"),
                  "Gao2015_nm3954_TableS1.xlsx"),
    # Isella et al. 2017 Nat Commun 8:15107 (Candiolo CRC liver-metastasis PDX): Supplementary Data 4 + GSE76402
    "isella_sd4": ("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5499209/supplementaryFiles",
                   "PMC5499209_supplementary.zip"),
    "gse76402": ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76402/matrix/GSE76402_series_matrix.txt.gz",
                 "GSE76402_series_matrix.txt.gz"),
    "gpl10558": ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
                 "GPL10558.annot.gz"),
    # van de Wetering et al. 2015 Cell 161:933 (CRC living organoid biobank): Table S2 (drug screen) + GSE64392
    "vdw_s2": ("https://ars.els-cdn.com/content/image/1-s2.0-S0092867415003736-mmc3.xlsx", "vdWetering2015_mmc3.xlsx"),
    "gse64392": ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64392/matrix/GSE64392_series_matrix.txt.gz",
                 "GSE64392_series_matrix.txt.gz"),
    "gpl16686": ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16686&targ=self&view=data&form=text",
                 "GPL16686_full.txt"),
    "hgnc": ("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
             "hgnc_complete_set.txt"),
})

LOWER = "lower = more sensitive"
HIGHER = "higher = more sensitive"


def write(name, expr, resp, md, catalog, direction=LOWER, extra_files=None):
    """bp.write_set with a per-set response direction; ensures n_screens / screens columns."""
    resp = resp.copy()
    if "n_screens" not in resp:
        resp["n_screens"] = 1
    if "screens" not in resp:
        resp["screens"] = resp["sample"]
    bp.DIRECTION = direction
    try:
        return bp.write_set(name, expr, resp, md, catalog, extra_files=extra_files)
    finally:
        bp.DIRECTION = LOWER


# ============================================================================ Novartis PDXE (Gao 2015)
# Novartis compound codes -> generic names (INN where one exists; otherwise the code is kept).
PDXE_DRUGS = {
    "5FU": "FLUOROURACIL", "BGJ398": "INFIGRATINIB", "BKM120": "BUPARLISIB", "BYL719": "ALPELISIB",
    "CGM097": "CGM097", "CKX620": "CKX620", "CLR457": "CLR457", "HDM201": "SIREMADLIN", "HSP990": "HSP990",
    "INC280": "CAPMATINIB", "INC424": "RUXOLITINIB", "LCL161": "LCL161", "LDE225": "SONIDEGIB",
    "LDK378": "CERITINIB", "LEE011": "RIBOCICLIB", "LFA102": "LFA102", "LFW527": "LFW527", "LGH447": "PIM447",
    "LGW813": "LGW813", "LJC049": "LJC049", "LJM716": "ELGEMTUMAB", "LKA136": "LKA136", "LLM871": "LLM871",
    "TAS266": "TAS266", "WNT974": "WNT-974", "abraxane": "NAB-PACLITAXEL", "binimetinib": "BINIMETINIB",
    "binimetinib-3.5mpk": "BINIMETINIB", "cetuximab": "CETUXIMAB", "dacarbazine": "DACARBAZINE",
    "encorafenib": "ENCORAFENIB", "erlotinib": "ERLOTINIB", "everolimus": "EVEROLIMUS",
    'figitumumab"': "FIGITUMUMAB", "gemcitabine-50mpk": "GEMCITABINE", "gemcitabine": "GEMCITABINE", "paclitaxel": "PACLITAXEL",
    "tamoxifen": "TAMOXIFEN", "trametinib": "TRAMETINIB", "trastuzumab": "TRASTUZUMAB",
}
PDXE_TUMOUR = {"BRCA": "breast", "CRC": "colorectal", "NSCLC": "lung (NSCLC)", "GC": "gastric",
               "PDAC": "pancreas", "CM": "melanoma"}


def mrecist(best, best_avg):
    """mRECIST call of Gao et al. 2015 (Online Methods), on % tumour-volume change from day 0."""
    if best < -95 and best_avg < -40:
        return "CR"
    if best < -50 and best_avg < -20:
        return "PR"
    if best < 35 and best_avg < 30:
        return "SD"
    return "PD"


def build_pdxe():
    path = bp.fetch("pdxe_xlsx")
    fpkm = pd.read_excel(path, sheet_name="RNAseq_fpkm", index_col=0)
    met = pd.read_excel(path, sheet_name="PCT curve metrics")
    raw = pd.read_excel(path, sheet_name="PCT raw data", usecols=["Model", "Tumor Type"]).drop_duplicates()
    ttype = raw.drop_duplicates("Model").set_index("Model")["Tumor Type"]

    met = met[met["Treatment"] != "untreated"].copy()
    met["drug"] = [" + ".join(PDXE_DRUGS[p.strip()] for p in t.split(" + ")) for t in met["Treatment"]]
    met["mrecist"] = [mrecist(b, a) for b, a in zip(met["BestResponse"], met["BestAvgResponse"])]
    met["responder"] = met["mrecist"].isin(["CR", "PR"]).astype(int)
    met["tumour_type"] = met["Model"].map(ttype)
    met["paper_category"] = met["ResponseCategory"]
    resp_rate = met["responder"].mean()
    agree = (met["mrecist"] == met["paper_category"].str.split("-").str[0]).mean()

    mr = met.rename(columns={"Model": "sample", "BestAvgResponse": "best_avg_response",
                             "BestResponse": "best_response", "Treatment": "treatment"})
    mr = mr[["sample", "drug", "best_avg_response", "mrecist", "responder", "best_response", "treatment",
             "tumour_type", "paper_category", "TimeToDouble", "Day_Last"]].rename(
        columns={"TimeToDouble": "time_to_double", "Day_Last": "day_last"})
    mr = mr.sort_values(["drug", "sample"]).round({"best_avg_response": 3, "best_response": 3})

    # response.tsv: BestAvgResponse. The low-dose arms (binimetinib-3.5mpk) are kept only in mrecist.tsv;
    # gemcitabine-50mpk is the only gemcitabine arm and is kept as GEMCITABINE.
    r = met[met["Treatment"] != "binimetinib-3.5mpk"]
    resp = pd.DataFrame({"sample": r["Model"], "drug": r["drug"], "response": r["BestAvgResponse"],
                         "metric": "BestAvgResponse", "n_screens": 1, "screens": r["Model"] + "|" + r["Treatment"],
                         "mrecist": r["mrecist"], "treatment_type": r["Treatment type"],
                         "tumour_type": r["tumour_type"]})
    assert not resp.duplicated(["sample", "drug"]).any()

    expr = bp.map_symbols(np.log2(fpkm.astype(float) + 1))
    samples = pd.DataFrame({"sample": sorted(set(expr.columns) | set(ttype.index))})
    samples["tumour_type"] = samples["sample"].map(ttype)
    samples["cancer"] = samples["tumour_type"].map(PDXE_TUMOUR)
    samples["has_expression"] = samples["sample"].isin(expr.columns).astype(int)
    samples["has_pct"] = samples["sample"].isin(resp["sample"]).astype(int)

    both = set(expr.columns) & set(resp["sample"])
    per_type = resp[resp["sample"].isin(both)].drop_duplicates("sample")["tumour_type"].value_counts()
    cnt = mr["mrecist"].value_counts()
    md = f"""# Novartis PDX Encyclopedia (Gao et al. 2015)

**Paper**: Gao H, Korn JM, Ferretti S, et al. *High-throughput screening using patient-derived tumor xenografts to
predict clinical trial drug response.* Nat Med 2015;21(11):1318-1325. doi:10.1038/nm.3954 (PMID 26479923).

**Accession**: Supplementary Table S1 (`nm.3954-S2.xlsx`, Springer static content), sheets `RNAseq_fpkm`,
`PCT curve metrics`, `PCT raw data`. Raw sequencing is not used.

{bp.provenance(["pdxe_xlsx", "gene_info"])}

## Design
PDX clinical trial (PCT), 1 x 1 x 1 design: one mouse per PDX model x treatment, tumour volume followed for
~3 weeks or more. {{n_resp_models}} models are treated; tumour types of the models with both data types:
{", ".join(f"{k} {v}" for k, v in per_type.items())} (BRCA breast, CRC colorectal, GC gastric, NSCLC lung,
PDAC pancreas, CM cutaneous melanoma).

## Derivation
* Response (`response.tsv`): **BestAvgResponse** from the paper (minimum over t >= 10 d of the running mean of the
  % tumour-volume change from day 0). **Lower = more sensitive** (negative = shrinkage). One row per model x
  treatment arm; the `untreated` arm and the low-dose `binimetinib-3.5mpk` arm are excluded (the latter is kept in
  `mrecist.tsv`). Combination arms are kept with the components joined by ` + ` (e.g. `ALPELISIB + ELGEMTUMAB`),
  `treatment_type` = single / combo. Novartis codes were mapped to INNs where one exists
  (BYL719 alpelisib, BKM120 buparlisib, LEE011 ribociclib, INC280 capmatinib, INC424 ruxolitinib, LDK378 ceritinib,
  LDE225 sonidegib, BGJ398 infigratinib, HDM201 siremadlin, LGH447 PIM447, LJM716 elgemtumab, WNT974,
  abraxane nab-paclitaxel); research codes without an INN are kept (CGM097, CLR457, CKX620, HSP990, LCL161,
  LFA102, LFW527, LGW813, LJC049, LKA136, LLM871, TAS266).
* `mrecist.tsv`: mRECIST class recomputed from BestResponse / BestAvgResponse with the thresholds of the paper's
  Online Methods: CR if BestResponse < -95% and BestAvgResponse < -40%; PR if BestResponse < -50% and
  BestAvgResponse < -20%; SD if BestResponse < 35% and BestAvgResponse < 30%; otherwise PD. `responder` = 1 for
  CR/PR. Counts: {", ".join(f"{k} {cnt.get(k, 0)}" for k in ["CR", "PR", "SD", "PD"])}. The paper's own
  `ResponseCategory` (which also tracks progression after the best response, e.g. `SD-->PD`) is kept as
  `paper_category`; its first class agrees with the recomputed call in {agree:.1%} of arms.
* Expression: `RNAseq_fpkm` (FPKM, Novartis pipeline) -> log2(FPKM + 1); symbols re-mapped to current HGNC
  protein-coding symbols. Units: **log2(FPKM + 1)**. Mouse reads were removed by the authors.

## Caveats
* No matched patient clinical response is public for PDXE (there is no `patient_response.tsv`).
* Single animal per arm: noisy. Treat the binary responder call as the most robust endpoint.
* Response rates are low (CR+PR {resp_rate:.0%} of arms), and most arms are targeted agents.

## Sample-ID matching
Model IDs (`X-1004` ...) shared by both sheets. Expression models: {{n_expr_models}}; treated models:
{{n_resp_models}}; **overlap n = {{n_both}}**. Genes: {{n_genes}}; drugs/regimens: {{n_drugs}}.
"""
    return write("pdx_novartis_gao2015", expr, resp, md,
                 {"type": "PDX", "tissue": "pan-cancer", "source": "Nat Med supplement"},
                 extra_files={"mrecist.tsv": mr, "samples.tsv": samples})


# ============================================================================ helpers for GEO arrays
def series_matrix(path):
    """GEO series matrix -> (data frame probes x GSM, per-sample annotation frame)."""
    import gzip
    meta, rows = {}, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t")
                v = [x.strip('"') for x in v]
                if k == "!Sample_characteristics_ch1":
                    key = v[0].split(":")[0].strip()
                    meta[key] = [x.split(":", 1)[1].strip() if ":" in x else x for x in v]
                else:
                    meta.setdefault(k[8:], v)
        data = pd.read_csv(f, sep="\t", index_col=0, comment="!", low_memory=False)
    data.index = data.index.astype(str)
    ann = pd.DataFrame(meta)
    ann.index = ann["geo_accession"]
    return data, ann


def probes_to_symbols(data, annot_path, col="Gene symbol"):
    """Probe x sample (log scale) -> HGNC symbol x sample (probe with the highest mean per symbol)."""
    import gzip
    with gzip.open(annot_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        a = pd.read_csv(f, sep="\t", low_memory=False, dtype=str, comment=None)
    a = a[a["ID"].notna() & ~a["ID"].str.startswith(("!", "^"))]
    sym = a.set_index("ID")[col].dropna()
    sym = sym[~sym.str.contains("///")]
    g = bp.gene_tables()["symbol"]
    return bp.collapse(data, {p: g.get(s.upper()) for p, s in sym.items()})


def refseq_to_symbol():
    """RefSeq transcript accession (no version) -> HGNC symbol, from the HGNC complete set (refseq + MANE)."""
    h = pd.read_csv(bp.fetch("hgnc"), sep="\t", dtype=str, low_memory=False)
    m = {}
    for sym, ref, mane in zip(h["symbol"], h["refseq_accession"], h["mane_select"]):
        for a in str(ref).split("|") + str(mane).split("|"):
            if a[:2] in ("NM", "NR", "XM", "XR"):
                m[a.split(".")[0]] = sym
    return m


# ============================================================================ CRC PDO (van de Wetering 2015)
VDW_DRUGS = {
    "(5Z)-7-Oxozeaenol": "5Z-7-OXOZEAENOL", "17-AAG": "TANESPIMYCIN", "5-Fluorouracil": "FLUOROURACIL",
    "681640": "681640", "ABT-263": "NAVITOCLAX", "AICAR": "ACADESINE", "AMG-706": "MOTESANIB",
    "AZD-2281": "OLAPARIB", "AZD6482": "AZD-6482", "AZD8055": "AZD-8055", "AZD8931": "SAPITINIB",
    "BIBW2992/Afatanib": "AFATINIB", "BIRB 0796": "DORAMAPIMOD", "BMN-673": "TALAZOPARIB",
    "BMS-708163": "AVAGACESTAT", "BYL719": "ALPELISIB", "CEP-701": "LESTAURTINIB", "Dasatanib": "DASATINIB",
    "GDC-0449": "VISMODEGIB", "GDC0941": "PICTILISIB", "GSK1120212/Tramatenib": "TRAMETINIB",
    "GSK2118436/Dabrafenib": "DABRAFENIB", "Gemcitibine": "GEMCITABINE", "INCB-18424/Ruxolitinib": "RUXOLITINIB",
    "JNJ-26854165": "SERDEMETAN", "LY317615": "ENZASTAURIN", "MLN8237": "ALISERTIB", "NVP-BEZ235": "DACTOLISIB",
    "Nutlin-3a": "NUTLIN-3A", "OSI-906": "LINSITINIB", "Obatoclax Mesylate": "OBATOCLAX",
    "PD-0332991": "PALBOCICLIB", "PF-02341066": "CRIZOTINIB", "PF477736": "PF-477736", "PLX4720": "PLX-4720",
    "SCH772984": "SCH-772984", "XAV 939": "XAV-939", "FK866": "FK-866",
}

def build_vdw():
    import io
    d = pd.read_excel(bp.fetch("vdw_s2"), sheet_name="Supplemental Table S2b").iloc[:, :7]
    d = d.dropna(subset=["Organoid", "AUC"])
    d["sample"] = d["Organoid"].astype(str).str.upper().str.replace(r"^P(\d+)A$", r"P\1TA", regex=True) \
        .str.replace(r"^P(\d+)B$", r"P\1TB", regex=True)
    d["drug"] = [VDW_DRUGS.get(str(n).strip()) or bp.generic(n) for n in d["Drug Name"]]
    d["rep"] = d.groupby(["sample", "Drug_ID"]).cumcount() + 1
    resp = (d.groupby(["sample", "drug"])
            .agg(response=("AUC", "median"), ln_ic50=("IC50", "median"), max_conc_uM=("Max_conc", "first"),
                 n_screens=("AUC", "size"), gdsc_drug_id=("Drug_ID", lambda x: ",".join(sorted({str(i) for i in x}))))
            .reset_index())
    resp["screens"] = resp["sample"] + "x" + resp["n_screens"].astype(str)
    resp["metric"] = "AUC"
    data, ann = series_matrix(bp.fetch("gse64392"))
    tum = ann[ann["title"].str.contains(r"t[ab]?$")]
    lines = tum["title"].str.upper().str.replace(r"T$", "", regex=True).str.replace(r"^(P\d+)T([AB])$", r"\1T\2", regex=True)
    lines = lines.where(~tum["title"].str.contains(r"t[ab]$"), tum["title"].str.upper())
    txt = bp.fetch("gpl16686").read_text().splitlines()
    i = next(k for k, l in enumerate(txt) if l.startswith("ID\t"))
    gpl = pd.read_csv(io.StringIO("\n".join(l for l in txt[i:] if not l.startswith("!"))), sep="\t", dtype=str)
    gpl = gpl[gpl["GB_ACC"].notna()]
    r2s = refseq_to_symbol()
    probe_sym = {p: r2s.get(a.split(".")[0]) for p, a in zip(gpl["ID"], gpl["GB_ACC"])}
    x = data[tum.index].astype(float)
    x.columns = lines.values
    g = bp.gene_tables()["symbol"]
    expr = bp.collapse(x, {p: g.get(str(s).upper()) for p, s in probe_sym.items() if s})
    gsm = tum[["geo_accession", "title"]].assign(sample=lines.values)
    md = f"""# Colorectal cancer PDOs - living organoid biobank (van de Wetering et al. 2015)

**Paper**: van de Wetering M, Francies HE, Francis JM, et al. *Prospective derivation of a living organoid biobank
of colorectal cancer patients.* Cell 2015;161(4):933-945. doi:10.1016/j.cell.2015.03.053 (PMID 25957691).

**Accessions**: drug screen = Table S2 (`mmc3.xlsx`, Elsevier supplementary-content CDN), sheet S2b (GDSC-style
screen at the Sanger Institute, 83 compounds, technical/biological replicates per organoid); expression = GEO
**GSE64392** (Affymetrix HuGene 2.0 ST, RMA-sketch, tumour and matched normal organoids).

{bp.provenance(["vdw_s2", "gse64392", "gpl16686", "hgnc", "gene_info"])}

## Derivation
* Response: **AUC** of the fitted dose-response curve as published (fraction of the tested range, 0-1; screen
  used a 5-dose, 256-fold range up to `max_conc_uM`). Median over replicate screens (`n_screens`). **Lower = more
  sensitive**. The published ln(IC50 / uM) median is kept as `ln_ic50` (not used by the loader). Drug names ->
  generic via DrugBank vocabulary (e.g. `Nutlin-3a (-)` -> NUTLIN-3A); research codes are kept.
  This replaces the degenerate CoderData 2.2 colorectal refit noted in REPORT_PRECLINICAL.md.
* Expression: GEO series-matrix RMA values (already log2), **tumour organoids only** (`p<n>t`, `p19ta/b`,
  `p24ta/b`; normal organoids `p<n>n` dropped). Transcript clusters -> RefSeq (`GB_ACC` of GPL16686) -> HGNC
  symbol (HGNC complete set, RefSeq/MANE columns), highest-mean cluster per gene, protein-coding only. Units:
  **log2 RMA**. Only ~13.5k genes map (GPL16686 lacks a gene-symbol column).

## Sample-ID matching
Drug-screen IDs `P10`, `P19a`, `P24a` = GEO titles `p10t`, `p19ta`, `p24ta` (written `P10`, `P19TA`, `P24TA`).
Expression organoids: {{n_expr_models}}; screened: {{n_resp_models}}; **overlap n = {{n_both}}**.
Genes: {{n_genes}}; drugs: {{n_drugs}}. No matched patient clinical response is public.
"""
    return write("colorectal_vandewetering2015", expr, resp, md,
                 {"type": "organoid", "tissue": "colorectal", "source": "Cell supplement + GEO"},
                 extra_files={"gsm_map.tsv": gsm})


# ============================================================================ CRC PDX cetuximab (Isella 2017)
def build_isella():
    import zipfile
    zf = zipfile.ZipFile(bp.fetch("isella_sd4"))
    sd4 = pd.read_excel(zf.open("ncomms15107-s5.xlsx"), header=None).iloc[4:, 1:16]
    sd4.columns = ["barcode", "sample", "msi", "kras", "nras", "braf", "fgfr1", "pdgfra", "map2k1", "erbb2_mut",
                   "erbb2_amp", "met_amp", "vol_change_3w", "vol_change_6w", "response_class"]
    sd4 = sd4.dropna(subset=["barcode"])
    data, ann = series_matrix(bp.fetch("gse76402"))
    ann["sample"] = ann["description"].map(sd4.set_index("barcode")["sample"])
    ann["sample"] = ann["sample"].fillna(ann["case_unique_id"] + "LM")
    lg = np.log2(data.astype(float).clip(lower=1))
    lg = lg.T.groupby(ann.loc[lg.columns, "sample"]).mean().T          # mean over arrays (regions / replicates)
    expr = probes_to_symbols(lg, bp.fetch("gpl10558"))
    r = sd4.dropna(subset=["vol_change_3w"]).drop_duplicates("sample")
    nbar = sd4.groupby("sample")["barcode"].apply(lambda s: ",".join(sorted(s)))
    resp = pd.DataFrame({"sample": r["sample"], "drug": "CETUXIMAB",
                         "response": r["vol_change_3w"].astype(float), "metric": "tumour_volume_change_3w",
                         "n_screens": 1, "screens": r["sample"],
                         "vol_change_6w": pd.to_numeric(r["vol_change_6w"], errors="coerce"),
                         "response_class": r["response_class"],
                         "kras_mut": r["kras"], "nras_mut": r["nras"], "braf_mut": r["braf"]})
    resp.loc[resp["vol_change_6w"] == 0, "vol_change_6w"] = np.nan      # 0 = not done in the source table
    gsm = ann[["geo_accession", "title", "description", "sample"]].rename(
        columns={"geo_accession": "gsm", "description": "array_barcode"})
    cls = resp["response_class"].value_counts()
    md = f"""# Colorectal cancer liver-metastasis PDX, cetuximab (Isella et al. 2017 / Bertotti et al.)

**Papers**: Isella C, Brundu F, Bellomo SE, et al. *Selective analysis of cancer-cell intrinsic transcriptional traits
defines novel clinically relevant subtypes of colorectal cancer.* Nat Commun 2017;8:15107. doi:10.1038/ncomms15107
(PMID 28561063). PDX cetuximab trials: Bertotti A, et al. Cancer Discov 2011;1:508 and Nature 2015;526:263.

**Accessions**: expression GEO **GSE76402** (Illumina HumanHT-12 v4, GPL10558, 529 arrays of 244 PDX models,
lumi/loess-normalised, human-specific probes); response: Supplementary Data 4 (`ncomms15107-s5.xlsx`, Europe PMC
open-access supplement of PMC5499209), keyed by array barcode (= GEO `Sample_description`).

{bp.provenance(["gse76402", "gpl10558", "isella_sd4", "gene_info"])}

## Derivation
* Response: **tumour volume change after 3 weeks of cetuximab** (20 mg/kg twice weekly), as a fraction of the
  volume at treatment start (-0.5 = 50% shrinkage; the paper's PR is < -0.5, PD > +0.35). **Lower = more
  sensitive**. The 6-week value (`vol_change_6w`, when available) and the paper's response class
  (`response_class`: {", ".join(f"{k} {v}" for k, v in cls.items())}) are kept. PR and SD were the paper's
  "cetuximab-sensitive". KRAS/NRAS/BRAF status from the same table is kept.
* Expression: series-matrix values (linear) -> log2; arrays of the same PDX model (regions A/B, replicate
  hybridisations) averaged on the log scale (`gsm_map.tsv`); probes -> HGNC symbols via the GPL10558 GEO
  annotation, highest-mean probe per symbol, protein-coding only. Units: **log2 (lumi-normalised intensity)**.

## Caveats
* Single agent (cetuximab, anti-EGFR antibody) only. Mostly KRAS-wild-type selected for EGFR biology; response is
  strongly driven by RAS/BRAF status.
* No matched patient clinical response is public for these models (no `patient_response.tsv`).

## Sample-ID matching
PDX model IDs (`CRC0014LM` ...) from Supplementary Data 4; GEO arrays matched by barcode. Expression models:
{{n_expr_models}}; models with a cetuximab response: {{n_resp_models}}; **overlap n = {{n_both}}**.
Genes: {{n_genes}}; drugs: {{n_drugs}}.
"""
    return write("colorectal_pdx_isella2017", expr, resp, md,
                 {"type": "PDX", "tissue": "colorectal (liver metastasis)", "source": "GEO + Nat Commun supplement"},
                 extra_files={"gsm_map.tsv": gsm})


BUILDERS = {
    "pdx_novartis_gao2015": build_pdxe,
    "colorectal_pdx_isella2017": build_isella,
    "colorectal_vandewetering2015": build_vdw,
}

bp.META.update({
    "pdx_novartis_gao2015": dict(model="PDX", tissue="pan-cancer (breast, CRC, gastric, NSCLC, pancreas, melanoma)",
                                 expression_units="log2(FPKM+1)", reference="Gao 2015 Nat Med 21:1318",
                                 accession="nm.3954 Supplementary Table S1"),
    "colorectal_pdx_isella2017": dict(model="PDX", tissue="colorectal (liver metastasis)",
                                      expression_units="log2(lumi intensity)", reference="Isella 2017 Nat Commun 8:15107",
                                      accession="GSE76402; PMC5499209 Suppl. Data 4"),
    "colorectal_vandewetering2015": dict(model="organoid", tissue="colorectal", expression_units="log2(RMA)",
                                         reference="van de Wetering 2015 Cell 161:933",
                                         accession="GSE64392; Cell Table S2 (mmc3)"),
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="*", help=f"subset of {list(BUILDERS)}")
    ap.add_argument("--raw", default=str(bp.RAW), help="download cache directory")
    a = ap.parse_args()
    bp.RAW = Path(a.raw)
    for s in a.sets or BUILDERS:
        BUILDERS[s]()
    bp.write_catalog()


if __name__ == "__main__":
    sys.exit(main())
