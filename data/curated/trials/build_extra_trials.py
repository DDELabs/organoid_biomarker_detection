#!/usr/bin/env python3
"""Extra public treated cohorts for cancers with few TCGA drug-response labels
(liver, kidney, glioma; plus melanoma anti-PD-1).

Writes data/curated/trials/<COHORT_ID>/{expression.tsv.gz, clinical.tsv, SOURCE.md} in the same
format as build_geo_trials.py (whose helpers it reuses) and appends / replaces the cohort rows in
CATALOG_GEO.tsv.  Differences from build_geo_trials.py:
  * microarray probes are collapsed to genes by keeping the probe with the highest mean
    expression ("max mean"), not the median over probes;
  * RNA-seq cohorts are log2(x + 1) and restricted to symbols resolvable to UniProt gene names
    (data/uniprot_homoSapiens_multipleGeneName_20180802.tab, i.e. what obd.to_canonical keeps);
  * clinical.tsv has a per-sample `cancer` column.

Cohorts
  GSE104580            HCC, TACE (pre-treatment biopsies), responder per GEO.            GPL570
  GSE140901            advanced HCC, anti-PD-1/PD-L1 ICI (Hsu 2021), n=24.               NanoString IO360
  EMTAB3267            metastatic ccRCC, first-line sunitinib (Beuselinck 2015).         HuGene 1.0 ST CEL -> RMA
  BRAUN2020_CHECKMATE  advanced ccRCC, nivolumab vs everolimus (CheckMate 009/010/025),  RNA-seq (paper suppl.)
  CGGA693              glioma WHO II-IV, TMZ-treated vs not (treatment-interaction, OS only), RNA-seq FPKM
  GSE7696              primary GBM, RT+TMZ vs RT (EORTC/NCIC trial), OS only.             GPL570
  GSE78220             metastatic melanoma, pembrolizumab/anti-PD-1 (Hugo 2016).          RNA-seq FPKM
  GSE91061             metastatic melanoma, nivolumab, pre-treatment (Riaz 2017).         RNA-seq FPKM
  GSE67501             metastatic RCC, nivolumab (Ascierto 2016), n=11.                    GPL14951

Usage:  GEO_CACHE=/path/to/cache python build_extra_trials.py [COHORT_ID ...]   (default: all)
Requirements: pandas, numpy, scikit-learn, openpyxl, curl.  EMTAB3267 additionally needs an R with
Bioconductor `oligo` + `pd.hugene.1.0.st.v1` ($RSCRIPT, default `Rscript`) to read the CEL files;
RMA (background, quantile normalisation, median polish) is then done in numpy because the
preprocessCore threaded C routines fail (pthread_create EINVAL) in some containers.
"""
import gzip
import os
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_geo_trials as B  # noqa: E402

CACHE = B.CACHE
REPO = B.REPO
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
MAX_MB = 40
# QC rule requested for these cohorts: flag if any PC1-3 AUC vs responder is < 0.3 or > 0.7
# (stricter than the 0.15/0.85 rule used by build_geo_trials.py)
B.QC_THRESH = 0.2
num = B.num


# ------------------------------------------------------------------------------------------ helpers
def collapse_maxmean(expr, mapping):
    """One row per gene symbol: the probe with the highest mean expression."""
    e = expr.loc[expr.index.intersection(mapping.index)]
    sym = mapping.loc[e.index]
    mean = e.mean(axis=1)
    order = pd.DataFrame({"sym": sym.values, "mean": mean.values}, index=e.index)
    best = order.sort_values("mean", ascending=False).groupby("sym", sort=False).head(1)
    out = e.loc[best.index]
    out.index = best["sym"].values
    out.index.name = "gene"
    return out.sort_index()


def uniprot_symbols():
    df = pd.read_csv(REPO / "data" / "uniprot_homoSapiens_multipleGeneName_20180802.tab", sep="\t")
    return {g for names in df["Gene names"].dropna() for g in names.split()}


def rnaseq_genes(expr):
    """Keep symbols resolvable to UniProt gene names; duplicates -> max mean."""
    keep = uniprot_symbols()
    e = expr[expr.index.isin(keep)]
    if e.index.duplicated().any():
        e = e.assign(_m=e.mean(axis=1)).sort_values("_m", ascending=False)
        e = e[~e.index.duplicated()].drop(columns="_m")
    e.index.name = "gene"
    return e.sort_index()


def microarray_maxmean(gse, gpl, name=None):
    p = B.series_matrix(gse, name)
    meta, raw = B.read_series_matrix(p)
    mapping, gp = B.probe_to_symbol(gpl)
    raw, lognote = B.to_log2(raw)
    expr = collapse_maxmean(raw, mapping)
    notes = [f"Series matrix {p.name} ({raw.shape[0]} probes); {lognote}.",
             f"Probes mapped with {gp.name}; probes with no / multiple symbols dropped; "
             f"probe with the highest mean kept per symbol -> {expr.shape[0]} genes."]
    return meta, expr, [p, gp], notes


def fetch(url, rel):
    return B.fetch(url, CACHE / rel)


def ensure_clin(c, cancer):
    c = c.copy()
    c["cancer"] = cancer if not isinstance(cancer, pd.Series) else cancer.values
    for col in ["responder", "os_months", "os_event", "pfs_months", "pfs_event"]:
        if col not in c:
            c[col] = np.nan
    return c


# -------------------------------------------------------------------------------------- RMA (numpy)
def _density_mode(x, n=16384):
    """Mode of a Gaussian KDE (bw.nrd0) evaluated on an n-point grid (as affy max.density)."""
    x = x[np.isfinite(x)]
    sd = x.std(ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo_ = min(sd, iqr / 1.34) if iqr > 0 else sd
    bw = 0.9 * lo_ * len(x) ** (-0.2)
    lo, hi = x.min() - 3 * bw, x.max() + 3 * bw
    hist, edges = np.histogram(x, bins=n, range=(lo, hi))
    step = edges[1] - edges[0]
    k = np.arange(-int(4 * bw / step) - 1, int(4 * bw / step) + 2) * step
    dens = np.convolve(hist, np.exp(-0.5 * (k / bw) ** 2), mode="same")
    return edges[np.argmax(dens)] + step / 2


def rma_bg(pm):
    """RMA convolution background correction for one array (affy::bg.adjust / bg.parameters)."""
    from math import sqrt
    from scipy.stats import norm
    mu = _density_mode(pm)
    mu = _density_mode(pm[pm < mu])
    bgd = pm[pm < mu] - mu
    sigma = sqrt((bgd ** 2).sum() / (len(bgd) - 1)) * sqrt(2)
    alpha = 1.0 / _density_mode(pm[pm > mu] - mu)
    a = pm - mu - alpha * sigma ** 2
    return a + sigma * norm.pdf(a / sigma) / norm.cdf(a / sigma)


def quantile_normalise(m):
    srt = np.sort(m, axis=0).mean(axis=1)
    ranks = m.argsort(axis=0).argsort(axis=0)
    return srt[ranks]


def median_polish_sets(y, sets, max_iter=10, eps=0.01):
    """Tukey median polish per probeset on log2 data (rows probes, cols arrays); returns sets x arrays."""
    order = np.argsort(sets, kind="stable")
    y, sets = y[order], sets[order]
    uniq, start = np.unique(sets, return_index=True)
    bounds = list(start) + [len(sets)]
    out = np.empty((len(uniq), y.shape[1]))
    for i in range(len(uniq)):
        z = y[bounds[i]:bounds[i + 1]].copy()
        t = 0.0
        r = np.zeros(z.shape[0]); c = np.zeros(z.shape[1])
        old = np.inf
        for _ in range(max_iter):
            rd = np.median(z, axis=1); z -= rd[:, None]; r += rd
            cd = np.median(z, axis=0); z -= cd[None, :]; c += cd
            s = np.abs(z).sum()
            if old == 0 or abs(s - old) < eps * old:
                break
            old = s
        out[i] = np.median(r) + c  # overall + column effects (affy convention)
    return pd.DataFrame(out, index=uniq)


R_PM = r'''
args <- commandArgs(TRUE)
suppressMessages(library(oligo))
f <- list.files(args[1], pattern="\\.CEL$", full.names=TRUE)
raw <- read.celfiles(f)
pi <- getProbeInfo(raw, target="core", field="fid", sortBy="none")
pm <- exprs(raw)[pi$fid, , drop=FALSE]
colnames(pm) <- sub("\\.CEL$", "", basename(colnames(pm)))
write.table(data.frame(fsetid=pi$man_fsetid, pm, check.names=FALSE), gzfile(args[2]),
            sep="\t", quote=FALSE, row.names=FALSE)
'''


# ------------------------------------------------------------------------------------------ cohorts
def GSE104580():
    meta, expr, files, notes = microarray_maxmean("GSE104580", "GPL570")
    c = B.base_clin(meta)
    sub = meta["ch_subject subgroup"].astype(str)
    c["response"] = sub.str.replace("TACE ", "", regex=False)
    c["responder"] = sub.map({"TACE responders": 1, "TACE non-responders": 0})
    c["arm"] = "TACE"
    c["drug"] = "TACE"
    c["drugs"] = "TACE"
    c["setting"] = "locoregional (TACE)"
    c = ensure_clin(c, "hepatocellular carcinoma")
    notes += ["Deposited data: GCRMA, log2, batch-corrected on scan date in Partek (per GEO).",
              "responder from 'subject subgroup' (TACE responders = 1 / non-responders = 0). GEO and the "
              "submitter (National Cancer Centre Singapore) give no response criterion (presumably "
              "mRECIST/EASL radiological response) and no publication is linked; no survival data.",
              "Chemotherapeutic agent(s) used in TACE are not stated; drugs = 'TACE' (placeholder)."]
    return expr, c, dict(accession="GSE104580", cancer="hepatocellular carcinoma", setting="locoregional",
                         platform="GPL570", endpoint="TACE response (responder/non-responder, criterion not stated)",
                         files=files, notes=notes, urls=B.geo_urls("GSE104580"),
                         paper="unpublished (GEO submitter: Kam Hui, National Cancer Centre Singapore)",
                         short_note="TACE agent unspecified (drugs=TACE)")


EMTAB = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/267/E-MTAB-3267/Files/"


def EMTAB3267():
    sdrf = fetch(EMTAB + "E-MTAB-3267.sdrf.txt", "E-MTAB-3267/E-MTAB-3267.sdrf.txt")
    s = pd.read_csv(sdrf, sep="\t")
    s.columns = [re.sub(r"^Characteristics\[(.*)\]$", r"\1", c) for c in s.columns]
    celdir = CACHE / "E-MTAB-3267" / "cel"
    for f in s["Array Data File"]:
        fetch(EMTAB + f, f"E-MTAB-3267/cel/{f}")
    pmfile = CACHE / "E-MTAB-3267" / "pm_raw.tsv.gz"
    if not pmfile.exists():
        rs = CACHE / "E-MTAB-3267" / "pm.R"
        rs.write_text(R_PM)
        subprocess.run([RSCRIPT, str(rs), str(celdir), str(pmfile)], check=True)
    pm = pd.read_csv(pmfile, sep="\t")
    sets = pm.pop("fsetid").values
    m = pm.values.astype(float)
    m = np.column_stack([rma_bg(m[:, j]) for j in range(m.shape[1])])
    m = np.log2(quantile_normalise(m))
    ts = median_polish_sets(m, sets)
    ts.columns = pm.columns
    ts.index = ts.index.astype(str)
    mapping, gp = B.probe_to_symbol("GPL6244")
    expr = collapse_maxmean(ts, mapping)

    s = s[s["disease"] == "Tumor"].copy()
    s.index = s["Assay Name"].astype(str)
    c = pd.DataFrame({"sample": s.index}, index=s.index)
    c["patient"] = "EMTAB3267_" + s["individual"].astype(str)
    c["arm"] = "SUNITINIB"
    c["drug"] = "SUNITINIB"
    c["drugs"] = "SUNITINIB"
    c["setting"] = "metastatic first-line"
    c["response"] = s["sunitinib response"]
    c["responder"] = s["sunitinib response"].map({"PR": 1, "CR": 1, "SD": 0, "PD": 0})
    c["clinical_benefit"] = s["sunitinib response"].map({"PR": 1, "CR": 1, "SD": 1, "CLINICAL BENEFIT": 1, "PD": 0})
    c["pfs_months"] = num(s["progression free survival"])
    c["pfs_event"] = num(s["progression"])
    c["age"] = num(s["age"])
    c["sex"] = s["sex"]
    c["histology"] = s["histology type"]
    c = ensure_clin(c, "kidney renal clear cell carcinoma (metastatic)")
    notes = [f"59 CEL files (HuGene-1_0-st, A-AFFY-141; 53 tumours + 6 normals) from BioStudies/ArrayExpress.",
             "Core-target PM intensities read with Bioconductor oligo (pd.hugene.1.0.st.v1); RMA done in numpy "
             "over all 59 arrays: affy-style convolution background (KDE-mode estimates), quantile "
             "normalisation, log2, median polish per transcript cluster -> "
             f"{ts.shape[0]} transcript clusters (approximates oligo::rma(target='core')).",
             f"Transcript clusters mapped with {gp.name}; probe with highest mean kept per symbol -> "
             f"{expr.shape[0]} genes. Normals dropped (53 tumours kept).",
             "responder: RECIST best response PR = 1, SD/PD = 0; 'CLINICAL BENEFIT' (n=10, response category "
             "not given in ArrayExpress) = NaN. clinical_benefit: PR/SD/CLINICAL BENEFIT = 1, PD = 0.",
             "PFS in months; pfs_event = 'progression' (1 = progressed). No OS deposited.",
             "Beuselinck et al. derived ccrcc1-4 subtypes on these data; the paper's primary endpoints were PFS/OS."]
    return expr, c, dict(accession="E-MTAB-3267", cancer="kidney renal clear cell carcinoma (metastatic)",
                         setting="metastatic", platform="A-AFFY-141 (HuGene 1.0 ST; GPL6244 annotation)",
                         endpoint="RECIST best response (PR vs SD/PD); PFS",
                         files=[sdrf, gp], notes=notes,
                         urls=["https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3267", EMTAB],
                         paper="Beuselinck B et al. 2015 Clin Cancer Res, PMID 25593300",
                         short_note="10 'clinical benefit' pts responder=NaN; RMA re-implemented")


BRAUN_XLSX = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/MediaObjects/"
              "41591_2020_839_MOESM2_ESM.xlsx")


def BRAUN2020_CHECKMATE():
    x = fetch(BRAUN_XLSX, "Braun2020/41591_2020_839_MOESM2_ESM.xlsx")
    cl = pd.read_excel(x, sheet_name="S1_Clinical_and_Immune_Data", header=1, na_values=["NA"])
    ex = pd.read_excel(x, sheet_name="S4A_RNA_Expression", header=1, index_col=0)
    ex.index = ex.index.astype(str)
    expr = rnaseq_genes(ex.astype(float))
    cl = cl[cl["RNA_ID"].notna() & cl["RNA_ID"].astype(str).isin(expr.columns)].copy()
    cl.index = cl["RNA_ID"].astype(str)
    c = pd.DataFrame({"sample": cl.index}, index=cl.index)
    c["patient"] = cl["SUBJID"]
    c["arm"] = cl["Arm"]
    c["drug"] = cl["Arm"]
    c["drugs"] = cl["Arm"]
    c["setting"] = "advanced, previously treated (anti-angiogenic)"
    c["response"] = cl["ORR"]
    c["responder"] = cl["ORR"].map({"CR": 1, "PR": 1, "CRPR": 1, "SD": 0, "PD": 0})
    c["clinical_benefit"] = cl["Benefit"].map({"CB": 1, "NCB": 0})  # ICB (intermediate) -> NaN
    c["benefit"] = cl["Benefit"]
    c["pfs_months"] = num(cl["PFS"])
    c["pfs_event"] = num(cl["PFS_CNSR"])
    c["os_months"] = num(cl["OS"])
    c["os_event"] = num(cl["OS_CNSR"])
    c["age"] = num(cl["Age"])
    c["sex"] = cl["Sex"]
    c["trial"] = cl["Cohort"]
    c["randomised"] = cl["Cohort"].eq("CM-025")
    c["control_arm"] = cl["Arm"].eq("EVEROLIMUS")
    for k in ["MSKCC", "IMDC", "Tumor_Sample_Primary_or_Metastasis", "Number_of_Prior_Therapies",
              "Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy", "Purity", "TMB_Counts", "PBRM1"]:
        c[k.lower()] = cl[k]
    c = ensure_clin(c, "kidney renal clear cell carcinoma (advanced)")
    notes = ["Braun DA et al. 2020 Nat Med Supplementary Table S1 (clinical) and S4A (normalised RNA "
             "expression matrix, 311 samples) downloaded from the Springer static-content server.",
             "Expression kept as deposited (log2-scale normalised, batch-corrected per paper; values ~0-68, "
             "median ~24 so absolute levels are shifted vs log2(TPM+1)); gene symbols restricted to UniProt "
             f"gene names -> {expr.shape[0]} genes.",
             "responder: ORR CR/PR/CRPR = 1, SD/PD = 0, NE = NaN (RECIST 1.1 per trial).",
             "clinical_benefit per paper: CB (CR/PR or SD with tumour shrinkage and PFS >= 6 mo) = 1, "
             "NCB (PD with PFS < 3 mo) = 0, ICB (intermediate) = NaN.",
             "PFS/OS in months; *_CNSR = 1 is an event (871/1006 PFS events) as in the paper's code.",
             "CM-025 is the randomised phase III (nivolumab vs everolimus); CM-009 and CM-010 are nivolumab-only "
             "phase I/II cohorts. Everolimus arm = comparator for treatment-interaction analyses."]
    return expr, c, dict(accession="Braun2020 NatMed suppl (dbGaP phs001493 / EGA for raw)",
                         cancer="kidney renal clear cell carcinoma (advanced)", setting="metastatic",
                         platform="RNA-seq (Illumina)", endpoint="RECIST ORR; clinical benefit; PFS; OS",
                         randomised=True, has_control=True, files=[x], notes=notes,
                         urls=["https://www.nature.com/articles/s41591-020-0839-y", BRAUN_XLSX],
                         paper="Braun DA et al. 2020 Nat Med, PMID 32472114",
                         short_note="CM-025 randomised NIVO vs EVE; CM-009/010 NIVO only")


CGGA_ZENODO_EXPR = "https://zenodo.org/api/records/8193658/files/CGGA.mRNAseq_693.csv/content"
CGGA_GH_CLIN = ("https://raw.githubusercontent.com/JackWJW/LGG_Prognosis_Prediction/HEAD/CGGA_Data/"
                "CGGA.mRNAseq_693_clinical.20200506.txt")


def CGGA693():
    ep = fetch(CGGA_ZENODO_EXPR, "CGGA/CGGA.mRNAseq_693.csv")
    cp = fetch(CGGA_GH_CLIN, "CGGA/CGGA.mRNAseq_693_clinical.20200506.txt")
    ex = pd.read_csv(ep, index_col=0).T  # samples x genes in the mirror -> genes x samples
    ex.index = ex.index.astype(str)
    expr = rnaseq_genes(np.log2(ex.astype(float).clip(lower=0) + 1))
    txt = Path(cp).read_text().replace("\r\n", "\n").replace("\r", "\n")
    import io
    cl = pd.read_csv(io.StringIO(txt), sep="\t", na_values=["NA"])
    cl.columns = [c.split(" (")[0] for c in cl.columns]
    cl.index = cl["CGGA_ID"].astype(str)
    cl = cl[cl.index.isin(expr.columns)]
    c = pd.DataFrame({"sample": cl.index}, index=cl.index)
    chemo = num(cl["Chemo_status"])
    c["arm"] = chemo.map({1: "TMZ", 0: "no TMZ"})
    c["drug"] = chemo.map({1: "TEMOZOLOMIDE", 0: "NONE"})
    c["drugs"] = c["drug"]
    c["setting"] = cl["PRS_type"].str.lower()
    c["responder"] = np.nan
    c["os_months"] = num(cl["OS"]) / 30.4375
    c["os_event"] = num(cl["Censor"])
    c["radio_status"] = num(cl["Radio_status"])
    c["age"] = num(cl["Age"])
    c["sex"] = cl["Gender"]
    c["histology"] = cl["Histology"]
    c["grade"] = cl["Grade"]
    c["prs_type"] = cl["PRS_type"]
    c["idh_status"] = cl["IDH_mutation_status"]
    c["codel_1p19q"] = cl["1p19q_codeletion_status"]
    c["mgmt_status"] = cl["MGMTp_methylation_status"]
    c["randomised"] = False
    c["control_arm"] = chemo.eq(0)
    gbm = cl["Histology"].isin(["GBM", "rGBM", "sGBM"])
    c = ensure_clin(c, pd.Series(np.where(gbm, "glioblastoma", "lower-grade glioma"), index=c.index))
    notes = ["CGGA mRNAseq_693 (release 20200506). cgga.org.cn is not reachable from the build environment "
             "(egress allow-list), so mirrors were used: expression = Zenodo record 8193658 "
             "'CGGA.mRNAseq_693.csv' (samples x genes; values are the CGGA RSEM FPKM, 23987 genes, "
             "despite the record text mentioning z-normalisation), clinical = the original CGGA "
             "clinical file mirrored in GitHub JackWJW/LGG_Prognosis_Prediction.",
             f"log2(FPKM + 1); symbols restricted to UniProt gene names -> {expr.shape[0]} genes.",
             "No response endpoint: responder = NaN. arm = 'TMZ' (Chemo_status 1) vs 'no TMZ' (0); NA if "
             "unknown. Treatment was not randomised (TMZ given by clinical indication; confounded by grade/"
             "era/IDH) - use as a treatment-interaction cohort with OS, adjusting for grade/IDH/MGMT/radio.",
             "os_months = OS days / 30.4375; os_event = Censor (1 = dead). cancer = 'glioblastoma' for "
             "GBM/rGBM/sGBM histology else 'lower-grade glioma'. Primary and recurrent tumours included "
             "(prs_type); patients may contribute >1 sample."]
    return expr, c, dict(accession="CGGA mRNAseq_693", cancer="glioma (WHO II-IV; incl. GBM)",
                         setting="primary/recurrent, adjuvant TMZ vs none", platform="RNA-seq (Illumina HiSeq)",
                         endpoint="OS (treatment-interaction, TMZ vs no TMZ)", randomised=False,
                         has_control=True, files=[ep, cp], notes=notes,
                         urls=["http://www.cgga.org.cn/download.jsp", "https://zenodo.org/records/8193658",
                               "https://github.com/JackWJW/LGG_Prognosis_Prediction"],
                         paper="Zhao Z et al. 2021 Genomics Proteomics Bioinformatics (CGGA), PMID 33662628",
                         short_note="no response labels; TMZ arm not randomised; mirror sources")


def GSE7696():
    meta, expr, files, notes = microarray_maxmean("GSE7696", "GPL570")
    meta = meta[meta["ch_disease status"] == "GBM"]
    c = B.base_clin(meta)
    tr = meta["ch_treatment"]
    c["patient"] = meta["ch_patient"]
    c["arm"] = tr.map({"TMZ/radiotherapy": "RT+TMZ", "radiotherapy": "RT"})
    c["drug"] = tr.map({"TMZ/radiotherapy": "TEMOZOLOMIDE", "radiotherapy": "NONE"})
    c["drugs"] = c["drug"]
    c["setting"] = "adjuvant (newly diagnosed)"
    c["responder"] = np.nan
    c["os_months"] = num(meta["ch_survival time in months"])
    c["os_event"] = num(meta["ch_survival status"])
    c["age"] = num(meta["ch_age"])
    c["sex"] = meta["ch_gender"]
    c["mgmt_status"] = meta["ch_mgmt status"].replace({"": np.nan, "?": np.nan})
    c["randomised"] = False
    c["control_arm"] = c["arm"].eq("RT")
    c = ensure_clin(c, "glioblastoma")
    notes += ["Primary GBM only (disease status 'GBM', n=70); recurrent / re-recurrent and non-tumoral samples "
              "dropped.",
              "arm: RT+TMZ (concomitant/adjuvant temozolomide) vs RT alone. Patients came from the EORTC 26981/"
              "NCIC CE.3 trial and associated studies; the expression subset is not a randomised sample, so arm is "
              "treated as non-randomised.",
              "No response endpoint (responder = NaN); OS in months, survival status 1 = dead."]
    return expr, c, dict(accession="GSE7696", cancer="glioblastoma", setting="adjuvant",
                         platform="GPL570", endpoint="OS (treatment-interaction, RT+TMZ vs RT)",
                         has_control=True, files=files, notes=notes, urls=B.geo_urls("GSE7696"),
                         paper="Murat A et al. 2008 J Clin Oncol, PMID 18565887",
                         short_note="no response labels; RT+TMZ vs RT")


def GSE78220():
    p = B.series_matrix("GSE78220")
    meta, _ = B.read_series_matrix(p)
    x = B.series_suppl("GSE78220", "GSE78220_PatientFPKM.xlsx")
    ex = pd.read_excel(x, index_col=0)
    ex.index = ex.index.astype(str)
    # FPKM columns are '<title>.baseline' / '<title>.OnTx'; GEO titles are Pt1, ..., Pt27A, Pt27B, ...
    t2g = dict(zip(meta["title"].astype(str), meta.index))
    ex = ex.rename(columns={c: t2g.get(c.split(".")[0], c) for c in ex.columns})
    ex = ex[[c for c in ex.columns if c in meta.index]]
    expr = rnaseq_genes(np.log2(ex.astype(float).clip(lower=0) + 1))
    c = B.base_clin(meta)
    c["patient"] = meta["ch_patient id"]
    c["arm"] = "ANTI-PD-1"
    tr = meta.get("ch_treatment", pd.Series(index=meta.index, dtype=object)).fillna("Pembrolizumab")
    c["drug"] = tr.str.upper()
    c["drugs"] = c["drug"]
    c["setting"] = "metastatic"
    resp = meta["ch_anti-pd-1 response"]
    c["response"] = resp
    c["responder"] = resp.map({"Complete Response": 1, "Partial Response": 1, "Progressive Disease": 0})
    c["os_months"] = num(meta["ch_overall survival (days)"]) / 30.4375
    c["os_event"] = meta["ch_vital status"].map({"Dead": 1, "Alive": 0})
    c["age"] = num(meta["ch_age (yrs)"])
    c["sex"] = meta["ch_gender"]
    c["biopsy_time"] = meta.get("ch_biopsy time")
    c["previous_mapki"] = meta["ch_previous mapki"]
    c = ensure_clin(c, "skin cutaneous melanoma (metastatic)")
    notes = [f"Series matrix {p.name} (no expression) + supplementary {x.name} (FPKM, 28 patients); "
             f"log2(FPKM + 1), symbols restricted to UniProt gene names -> {expr.shape[0]} genes.",
             "responder: irRECIST CR/PR = 1, PD = 0 (no SD in this cohort).",
             "Treatment from the GEO 'treatment' field (pembrolizumab for all 28); one biopsy is annotated "
             "on-treatment (biopsy_time; Pt16), the rest pre-treatment. Pt27A/Pt27B are two baseline "
             "lesions of one patient (same `patient`).",
             "OS days / 30.4375; vital status Dead = event."]
    return expr, c, dict(accession="GSE78220", cancer="skin cutaneous melanoma (metastatic)",
                         setting="metastatic", platform="RNA-seq (GPL11154)",
                         endpoint="irRECIST response (CR/PR vs PD); OS", files=[p, x], notes=notes,
                         urls=B.geo_urls("GSE78220"), paper="Hugo W et al. 2016 Cell, PMID 26997480",
                         short_note="small (n=28)")


def GSE140901():
    p = B.series_matrix("GSE140901")
    meta, _ = B.read_series_matrix(p)
    f = B.series_suppl("GSE140901", "GSE140901_processed_data.txt.gz")
    ex = pd.read_csv(f, sep="\t", skiprows=4, index_col=0)
    ex = ex.loc[:, ~ex.columns.str.startswith("Unnamed")]
    ex.index = ex.index.astype(str).str.strip()
    t2g = {t.replace("HCC ", ""): g for t, g in zip(meta["title"].astype(str), meta.index)}
    ex = ex.rename(columns=t2g)[[g for g in meta.index if g in t2g.values()]].astype(float)
    ex = ex[ex.index != "nan"]
    expr = ex.groupby(level=0).max() if ex.index.duplicated().any() else ex
    expr.index.name = "gene"
    c = B.base_clin(meta)
    c["arm"] = "ICI"
    c["drug"] = "ANTI-PD-1/PD-L1"
    c["drugs"] = "ICI_ANTI_PD1_PDL1"
    c["setting"] = "advanced/metastatic"
    br = meta["ch_best_response"]
    c["response"] = br
    c["responder"] = br.map({"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    c["clinical_benefit"] = meta["ch_clinical_benefit_response"].map({"Yes": 1, "No": 0})
    wk = 7 / 30.4375
    c["pfs_months"] = num(meta["ch_pfs_time"]) * wk
    c["pfs_event"] = num(meta["ch_pfs_event"])
    c["os_months"] = num(meta["ch_os_time"]) * wk
    c["os_event"] = num(meta["ch_os_event"])
    c["age"] = num(meta["ch_age"])
    c["sex"] = meta["ch_gender"]
    c["etiology"] = meta["ch_etiology"]
    c = ensure_clin(c, "hepatocellular carcinoma (advanced)")
    notes = [f"Supplementary {f.name}: NanoString PanCancer IO 360 panel (GPL19965), TMM-normalised log2 CPM "
             f"as deposited; {expr.shape[0]} genes only (targeted immune panel, not genome-wide).",
             "Archival pre-treatment tumour; anti-PD-1/PD-L1-based ICI (agent/combination per patient not "
             "deposited) -> drugs = 'ICI_ANTI_PD1_PDL1' placeholder.",
             "responder: best response PR = 1, SD/PD = 0 (no CR); clinical_benefit as deposited.",
             "pfs_time / os_time are in weeks (shortest PFS 5.1 = first restaging) -> months = weeks * 7 / 30.4375."]
    return expr, c, dict(accession="GSE140901", cancer="hepatocellular carcinoma (advanced)",
                         setting="metastatic", platform="GPL19965 (NanoString IO360, ~770 genes)",
                         endpoint="best response (PR vs SD/PD); clinical benefit; PFS; OS", files=[p, f],
                         notes=notes, urls=B.geo_urls("GSE140901"),
                         paper="Hsu CL et al. 2021 Liver Cancer, PMID 34414122",
                         short_note="NanoString panel (~770 genes); ICI agent unspecified")


RIAZ_CLIN = "https://raw.githubusercontent.com/riazn/bms038_analysis/master/data/bms038_clinical_data.csv"
GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def GSE91061():
    f = B.series_suppl("GSE91061", "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")
    gi = fetch(GENE_INFO, "NCBI/Homo_sapiens.gene_info.gz")
    rc = fetch(RIAZ_CLIN, "GSE91061/bms038_clinical_data.csv")
    p = B.series_matrix("GSE91061")
    meta, _ = B.read_series_matrix(p)
    ex = pd.read_csv(f, index_col=0)
    g = pd.read_csv(gi, sep="\t", usecols=["GeneID", "Symbol"])
    e2s = dict(zip(g["GeneID"].astype(str), g["Symbol"]))
    ex.index = ex.index.astype(str).map(lambda i: e2s.get(i))
    ex = ex[ex.index.notna()]
    # FPKM column names equal the GEO sample titles (PtX_Pre_<id> / PtX_On_<id>); pre-treatment only
    t2g = dict(zip(meta["title"].astype(str), meta.index))
    pre = [c for c in ex.columns if "_Pre_" in c and c in t2g]
    gsm = {c.split("_")[0]: t2g[c] for c in pre}
    ex = ex[pre].rename(columns=t2g)
    expr = rnaseq_genes(np.log2(ex.astype(float).clip(lower=0) + 1))
    cl = pd.read_csv(rc)
    cl.index = cl["PatientID"].map(gsm)
    cl = cl[cl.index.notna() & cl.index.isin(expr.columns)]
    c = pd.DataFrame({"sample": cl.index}, index=cl.index)
    c["patient"] = cl["PatientID"]
    c["arm"] = cl["Cohort"]
    c["drug"] = "NIVOLUMAB"
    c["drugs"] = "NIVOLUMAB"
    c["setting"] = "metastatic"
    c["response"] = cl["BOR"]
    c["responder"] = cl["BOR"].map({"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    c["os_months"] = num(cl["OS"]) / 30.4375
    c["os_event"] = 1 - num(cl["OS_SOR"])
    c["pfs_months"] = num(cl["PFS"]) / 30.4375
    c["pfs_event"] = 1 - num(cl["PFS_SOR"])
    c["prior_ipilimumab"] = cl["Cohort"].eq("NIV3-PROG")
    c["subtype"] = cl["SubtypeEZ"]
    c = ensure_clin(c, "skin cutaneous melanoma (metastatic)")
    notes = [f"Supplementary {f.name} (FPKM, Entrez gene IDs) mapped to symbols with NCBI "
             "Homo_sapiens.gene_info; pre-treatment biopsies only; log2(FPKM + 1); symbols restricted to "
             f"UniProt gene names -> {expr.shape[0]} genes.",
             "Clinical data (BOR, OS, PFS in days) from the authors' repository riazn/bms038_analysis "
             "(data/bms038_clinical_data.csv); *_SOR = 1 is censored, so event = 1 - SOR.",
             "responder: RECIST BOR CR/PR = 1, SD/PD = 0, NE = NaN. arm: NIV3-NAIVE (ipilimumab-naive) vs "
             "NIV3-PROG (progressed on ipilimumab); all received nivolumab (CheckMate 038)."]
    return expr, c, dict(accession="GSE91061", cancer="skin cutaneous melanoma (metastatic)",
                         setting="metastatic", platform="RNA-seq (GPL9052)",
                         endpoint="RECIST BOR (CR/PR vs SD/PD); OS; PFS", files=[f, gi, rc, p],
                         notes=notes, urls=B.geo_urls("GSE91061") + [RIAZ_CLIN],
                         paper="Riaz N et al. 2017 Cell, PMID 29033130",
                         short_note="pre-treatment only; arm = ipi-naive vs ipi-progressed")


def GSE67501():
    meta, expr, files, notes = microarray_maxmean("GSE67501", "GPL14951")
    c = B.base_clin(meta)
    k = [x for x in meta.columns if x.startswith("ch_response") and "complete response" in x][0]
    c["response"] = meta[k]
    c["responder"] = meta[k].map({"CR": 1, "PR": 1, "SD": 0, "NR": 0})
    c["arm"] = "NIVOLUMAB"
    c["drug"] = "NIVOLUMAB"
    c["drugs"] = "NIVOLUMAB"
    c["setting"] = "metastatic"
    c["sex"] = meta["ch_gender"]
    c["site"] = meta["ch_primary tumor or metastasis"]
    c = ensure_clin(c, "renal cell carcinoma (metastatic)")
    notes += ["responder: CR/PR = 1, SD/NR = 0 (matches the series' response/no_response split).",
              "Archival tumour (collected 2-81 months before nivolumab). n=11 only; no survival data."]
    return expr, c, dict(accession="GSE67501", cancer="renal cell carcinoma (metastatic)", setting="metastatic",
                         platform="GPL14951", endpoint="RECIST response (CR/PR vs SD/NR)", files=files,
                         notes=notes, urls=B.geo_urls("GSE67501"),
                         paper="Ascierto ML et al. 2016 Cancer Immunol Res, PMID 27491898",
                         short_note="n=11; archival tissue")


COHORTS = OrderedDict((f.__name__, f) for f in [
    GSE104580, GSE140901, EMTAB3267, BRAUN2020_CHECKMATE, GSE67501, CGGA693, GSE7696, GSE78220, GSE91061])


def main(ids):
    cat_path = B.CATALOG
    cat = pd.read_csv(cat_path, sep="\t")
    for cid in ids:
        print(f"== {cid}", flush=True)
        expr, clin, info = COHORTS[cid]()
        row = B.write_cohort(cid, expr, clin, info)
        src = HERE / cid / "SOURCE.md"
        src.write_text(src.read_text().replace("build_geo_trials.py", "build_extra_trials.py")
                       .replace("- GEO accession:", "- Accession:")
                       .replace("abs(AUC-0.5)>0.35", "AUC<0.3 or AUC>0.7"))
        size = (HERE / cid / "expression.tsv.gz").stat().st_size / 1e6
        assert size < MAX_MB, f"{cid}: expression.tsv.gz {size:.1f} MB"
        print({k: row[k] for k in ["n", "responders", "n_labelled", "n_genes", "qc_pc1_auc", "qc_flag",
                                   "qc_arm_pc1_auc"]}, f"{size:.1f} MB", flush=True)
        new = pd.DataFrame([row])
        if cid in set(cat["cohort_id"]):
            cat = cat[cat["cohort_id"] != cid]
        cat = pd.concat([cat, new], ignore_index=True)
        for col in ["n", "responders", "n_labelled", "n_genes"]:
            cat[col] = cat[col].astype("Int64")
        cat.to_csv(cat_path, sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:] or list(COHORTS))
