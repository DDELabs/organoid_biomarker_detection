import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
"""Build data/external/trials/<cohort>/ from harvested GitHub-hosted sources (scratch build script)."""
import datetime as dt
import gzip
import hashlib
import io
import pickle
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
SCR = Path("/tmp/obd_raw")
GH = SCR / "gh"
OUT = Path(_REPO + '/data/external/trials')
OUT.mkdir(parents=True, exist_ok=True)
CATALOG_ROWS = []
ONLY = set(sys.argv[1:])


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- gene symbol repair
MON = {"Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6, "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12}


def _md_symbol(month, day):
    if month == 3:
        return f"MARCHF{day}"
    if month == 9:
        return f"SEPTIN{day}"
    if month == 12 and day == 1:
        return "DELEC1"
    return None


def fix_symbol(s):
    s = str(s).strip()
    m = re.fullmatch(r"(\d{1,2})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)", s)
    if m:
        return _md_symbol(MON[m.group(2)], int(m.group(1))) or s
    m = re.fullmatch(r"\d{4}-(\d{2})-(\d{2})( 00:00:00)?", s)  # Excel date -> yyyy-mm-dd (day = gene number)
    if m:
        return _md_symbol(int(m.group(1)), int(m.group(2))) or s
    if re.fullmatch(r"4\d{4}", s):  # Excel serial date
        d = dt.date(1899, 12, 30) + dt.timedelta(days=int(s))
        return _md_symbol(d.month, d.day) or s
    return s


def tidy_expr(expr):
    """genes x samples -> numeric, repaired symbols, junk rows dropped, duplicate symbols averaged."""
    expr = expr.apply(pd.to_numeric, errors="coerce")
    idx = pd.Index([fix_symbol(g) for g in expr.index])
    expr.index = idx
    bad = idx.isin(["---", "", "nan", "?", "NA"]) | idx.isna()
    expr = expr.loc[~bad]
    expr = expr.loc[expr.notna().any(axis=1)]
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0, sort=False).mean()
    expr.index.name = "gene"
    return expr


def to_log(expr, kind):
    """kind: 'log' (already log), 'linear' (log2(x+1)), 'counts' (log2(CPM+1))."""
    if kind == "log":
        return expr, "values as distributed (already log-scale)"
    if kind == "counts":
        lib = expr.sum(axis=0)
        return np.log2(expr.div(lib, axis=1) * 1e6 + 1), "log2(CPM+1) computed from distributed read counts"
    return np.log2(expr.clip(lower=0) + 1), "log2(x+1) of distributed linear values"


# ---------------------------------------------------------------- QC
def auc(score, y):
    y = np.asarray(y, dtype=float)
    s = np.asarray(score, dtype=float)
    ok = ~np.isnan(y) & ~np.isnan(s)
    y, s = y[ok], s[ok]
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().values
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def pc_qc(expr, y, ntop=5000):
    y = pd.Series(y).dropna()
    X = expr.loc[:, [s for s in y.index if s in expr.columns]]
    y = y.loc[X.columns]
    X = X.loc[X.notna().all(axis=1)]
    X = X.loc[X.var(axis=1) > 0]
    X = X.loc[X.var(axis=1).nlargest(min(ntop, len(X))).index]
    Z = (X.sub(X.mean(axis=1), axis=0)).T.values
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    pcs = U[:, :3] * S[:3]
    var = (S ** 2 / (S ** 2).sum())[:3]
    return [auc(pcs[:, i], y.values) for i in range(3)], var, len(y)


# ---------------------------------------------------------------- writer
def write_cohort(cid, expr, clin, meta, qc_label="responder"):
    if ONLY and cid not in ONLY:
        return
    d = OUT / cid
    d.mkdir(parents=True, exist_ok=True)
    clin = clin.copy()
    clin.index = clin.index.astype(str)
    expr.columns = expr.columns.astype(str)
    keep = [s for s in clin.index if s in expr.columns]
    clin = clin.loc[keep]
    expr = expr.loc[:, keep]
    expr = expr.loc[expr.notna().any(axis=1)]
    expr = expr.round(4)
    clin.index.name = "sample"
    with gzip.open(d / "expression.tsv.gz", "wt", compresslevel=6) as f:
        expr.to_csv(f, sep="\t", na_rep="NA")
    clin.reset_index().to_csv(d / "clinical.tsv", sep="\t", index=False, na_rep="NA")
    lab = clin[qc_label] if qc_label in clin else clin["responder"]
    try:
        aucs, var, nq = pc_qc(expr, lab)
    except Exception as e:  # noqa
        aucs, var, nq = [np.nan] * 3, [np.nan] * 3, 0
    flag = bool(abs(aucs[0] - 0.5) > 0.35) if not np.isnan(aucs[0]) else False
    within = []
    if "arm" in clin and clin["arm"].nunique() > 1 and qc_label == "responder":
        for a, sub in clin.groupby("arm"):
            if sub["responder"].nunique() == 2 and len(sub) >= 10:
                aa, _, nn = pc_qc(expr[sub.index], sub["responder"])
                within.append(f"{a}: n={nn}, PC1 AUC={aa[0]:.3f}, PC2={aa[1]:.3f}, PC3={aa[2]:.3f}")
    n = len(clin)
    nresp = int((clin["responder"] == 1).sum()) if "responder" in clin else 0
    nlab = int(clin["responder"].notna().sum()) if "responder" in clin else 0
    usable = meta.get("usable_override")
    if usable is None:
        usable = (not flag) and nlab >= 15 and min(nresp, nlab - nresp) >= 4
    qc_txt = (f"PCA on top-{min(5000, expr.shape[0])} most variable genes (gene-centred, log scale), samples with "
              f"non-missing `{qc_label}` (n={nq}). AUC of PC score vs `{qc_label}` (1 = higher score in label 1):\n"
              f"- PC1 AUC = {aucs[0]:.3f} (var {var[0]:.1%})\n- PC2 AUC = {aucs[1]:.3f} (var {var[1]:.1%})\n"
              f"- PC3 AUC = {aucs[2]:.3f} (var {var[2]:.1%})\n"
              f"- Flag |PC1 AUC - 0.5| > 0.35: **{'FAILED (batch/label aligned)' if flag else 'pass'}**\n")
    if within:
        qc_txt += "- Within-arm PC AUCs vs responder:\n" + "".join(f"  - {w}\n" for w in within)
    src_lines = "".join(f"- {u}\n  sha256 `{h}`\n" for u, h in meta["sources"])
    md = f"""# {cid} — {meta['title']}

{meta.get('citation', '')}

- Cancer: {meta['cancer']}
- Drug(s)/regimen: {meta['drugs']}
- Setting: {meta['setting']}
- Control/comparator arm: {meta.get('control', 'none')}
- Endpoint / label: {meta['endpoint']}
- Original accession: {meta['accession']}
- Platform: {meta['platform']}
- n = {n} samples in these files; responders = {nresp} of {nlab} labelled

## Provenance
NCBI GEO / EBI are blocked from the build environment; files were taken from the GitHub-hosted copies below
(downloaded {dt.date.today().isoformat()}).

{src_lines}
## Processing
- expression.tsv.gz: {expr.shape[0]} genes x {expr.shape[1]} samples, first column `gene` (HUGO symbol).
  {meta['units']}
  Duplicate symbols averaged; Excel-mangled symbols repaired (e.g. 1-Mar -> MARCHF1, 1-Sep -> SEPTIN1).
- clinical.tsv: {', '.join(clin.columns)}.
{meta.get('notes', '')}

## QC
{qc_txt}
Usable flag (QC pass, >=15 labelled samples, >=4 per class): **{usable}**
"""
    (d / "SOURCE.md").write_text(md)
    row = dict(cohort_id=cid, cancer=meta["cancer"], drugs=meta["drugs"], setting=meta["setting"], n=n,
               responders=nresp if nlab else np.nan, n_labelled=nlab,
               has_control_arm=meta.get("has_control", False), arms=";".join(map(str, sorted(clin["arm"].dropna().unique()))) if "arm" in clin else "",
               endpoint=meta["endpoint"], platform=meta["platform"], accession=meta["accession"],
               source=meta["source_short"], qc_label=qc_label,
               qc_pc1_auc=round(aucs[0], 3), qc_pc2_auc=round(aucs[1], 3), qc_pc3_auc=round(aucs[2], 3),
               qc_flag=flag, usable=bool(usable))
    CATALOG_ROWS.append(row)
    print(cid, expr.shape, n, nresp, [round(a, 3) for a in aucs], "FLAG" if flag else "", "usable" if usable else "NOT usable", flush=True)


# ======================================================================= ENLIGHT
EN = GH / "enlight-data" / "Data"
EN_RAW = "https://raw.githubusercontent.com/PangeaResearch/enlight-data/main/Data/"
# name: (drug, cancer, setting, accession, confidence, platform, kind, extra note)
ENLIGHT = {
    "Anti-PD1 +- Anti-CTLA4": ("anti-PD1 +/- anti-CTLA4", "HCC", "advanced", "GSE140901 (inferred from GSM4190059-82)", "medium", "NanoString PanCancer Immune (784 genes)", "log", ""),
    "Anti-PD1": ("nivolumab (anti-PD1)", "Melanoma", "advanced", "Riaz 2017 GSE91061 (inferred; Pt*_On = ON-TREATMENT biopsies)", "high", "RNA-seq read counts", "counts", "Samples are on-treatment biopsies (Pt*_On), not baseline."),
    "Anti-PD1_2": ("anti-PD1 (pembrolizumab/nivolumab)", "GBM", "recurrent", "Zhao et al. 2019 Nat Med, SRA PRJNA482620 (inferred from SRR8281xxx)", "medium", "RNA-seq read counts", "counts", ""),
    "Anti-PD1_3": ("nivolumab (anti-PD1)", "RCC", "metastatic", "GSE67501 (Ascierto 2016; inferred from GSM1648114-24)", "high", "Illumina microarray (probe-collapsed)", "log", ""),
    "Anti-PD1_4": ("durvalumab + olaparib + paclitaxel (I-SPY2 arm)", "Breast", "neoadjuvant", "GSE173839 (I-SPY2 durvalumab/olaparib arm; inferred from GSM5281524-628)", "high", "Agilent microarray", "log", ""),
    "Anti-PD1_5": ("anti-PD1", "Melanoma", "advanced", "unresolved (sample IDs 190511-xxF)", "low", "RNA-seq FPKM (log-like values)", "log", ""),
    "BRAFi_1": ("BRAF inhibitor (vemurafenib/dabrafenib)", "Melanoma", "advanced", "Hugo et al. 2015 Cell, GSE65185 (inferred from Pt*-baseline IDs)", "medium", "RNA-seq FPKM", "linear", ""),
    "BRAFi_2": ("BRAF inhibitor +/- MEK inhibitor", "Melanoma", "advanced", "GSE99898 (inferred from GSM2663924-60; low confidence)", "low", "microarray", "log", ""),
    "BRAFi_3": ("BRAF inhibitor (vemurafenib/dabrafenib)", "Melanoma", "advanced", "GSE50509 (Rizos 2014; inferred from GSM1220412-71)", "high", "Illumina microarray (background-subtracted log values incl. negatives)", "log", ""),
    "Bevacizumab": ("mFOLFOX6 + bevacizumab", "Colorectal", "metastatic", "GSE19860 (bevacizumab subset; inferred from GSM496016-43)", "medium", "Affymetrix (centred log ratios)", "log", ""),
    "Bevacizumab_2": ("bevacizumab-containing chemo", "Colorectal", "metastatic", "unresolved (GSM1282761-78)", "low", "microarray (linear)", "linear", ""),
    "Bevacizumab_3": ("bevacizumab + chemo", "Breast (per drug_targets.csv; README lists Breast/read counts but values are log-scale microarray)", "neoadjuvant/advanced", "unresolved (GSM2778766-86)", "low", "microarray (log)", "log", ""),
    "Bevacizumab_4": ("bevacizumab + chemo", "Colorectal", "metastatic", "unresolved (GSM1471839-86)", "low", "microarray (log, probe-summed)", "log", ""),
    "Cetuximab": ("cetuximab (+ platinum/5-FU)", "Head and neck", "recurrent/metastatic", "GSE65021 (Bossi 2016; inferred from GSM1585940-79)", "high", "Affymetrix (linear)", "linear", ""),
    "Lapatinib": ("lapatinib (+/- trastuzumab) + paclitaxel->FEC (CHER-LOB arms B/C)", "Breast (HER2+)", "neoadjuvant", "GSE66399 CHER-LOB (inferred from GSM1619387-474)", "medium", "Affymetrix (log, probe-summed)", "log", "Same series as Trastuzumab_2 (CHER-LOB arm A)."),
    "MGH_Alpelisib": ("alpelisib", "Breast", "metastatic", "MGH institutional (not public GEO)", "n/a", "RNA-seq TPM", "linear", ""),
    "MGH_Ribociclib": ("ribociclib", "Breast", "metastatic", "MGH institutional (not public GEO)", "n/a", "RNA-seq TPM", "linear", ""),
    "MK2206": ("MK-2206 + paclitaxel (I-SPY2 arm)", "Breast", "neoadjuvant", "I-SPY2 MK-2206 arm (GSM4497593-686; series unresolved)", "medium", "Agilent microarray (log, probe-summed)", "log", ""),
    "Rituximab": ("rituximab-containing (FCR)", "CLL", "first-line", "unresolved (GSM877692-753)", "low", "microarray (linear)", "linear", ""),
    "Selinexor": ("selinexor", "GBM", "recurrent", "KING trial (Lassman 2022), not in GEO", "n/a", "RNA-seq TPM (8508 genes)", "linear", ""),
    "Sorafenib": ("sorafenib (STORM adjuvant, sorafenib arm only)", "HCC", "adjuvant", "GSE109211 (STORM; 67 of 140, placebo arm absent)", "high", "Illumina WG-DASL (linear)", "linear", "Known batch/label alignment (see QC)."),
    "Sorafenib_2": ("sorafenib (BATTLE)", "NSCLC", "metastatic, pretreated", "GSE33072 (BATTLE sorafenib arm)", "high", "Affymetrix HuGene 1.0 ST (log)", "log", "Label = PFS > ~2.5 months (ENLIGHT; verified against ISLE PFS)."),
    "Tipifarnib_1": ("tipifarnib", "AML", "newly diagnosed (elderly)", "GSE8970 (GSM227xxx; ENLIGHT file Tipifarnib_1.csv holds GSM227* although its labels are filed under Tipifarnib_2)", "high", "Affymetrix U133A (MAS5 linear)", "linear", "Labels matched by sample ID; the ENLIGHT dataset names of Tipifarnib_1/2 are swapped between expression files and the label table."),
    "Tipifarnib_2": ("tipifarnib", "AML", "relapsed/refractory", "GSE5122 (GSM115xxx; labels filed under Tipifarnib_1 in ENLIGHT table)", "high", "Affymetrix U133A (MAS5 linear)", "linear", "Labels matched by sample ID (see Tipifarnib_1)."),
    "Trastuzumab": ("trastuzumab + chemo (NOAH, trastuzumab arm only)", "Breast (HER2+)", "neoadjuvant", "GSE50948 NOAH (inferred from GSM1232994-3146; chemo-only control arm NOT included)", "medium", "Affymetrix U133Plus2 (log, probe-summed)", "log", ""),
    "Trastuzumab_2": ("trastuzumab + paclitaxel->FEC (CHER-LOB arm A)", "Breast (HER2+)", "neoadjuvant", "GSE66399 CHER-LOB (inferred)", "medium", "Affymetrix (log, probe-summed)", "log", ""),
    "Trastuzumab_3": ("trastuzumab + chemo", "Breast (HER2+)", "neoadjuvant", "GSE37946 (inferred from GSM930525-74)", "high", "Affymetrix U133A (log)", "log", ""),
    "Trastuzumab_4": ("trastuzumab + chemo", "Breast (HER2+)", "neoadjuvant", "unresolved (GSM1050643-67)", "low", "Affymetrix U133A (log)", "log", ""),
    "Trastuzumab_5": ("trastuzumab-containing", "Breast (HER2+)", "unknown", "unresolved (T### IDs)", "low", "RNA-seq read counts", "counts", ""),
    "Vismodegib": ("vismodegib", "Basal cell carcinoma", "advanced", "unresolved (Resistance/Sensitive IDs; possibly Atwood 2015)", "low", "RNA-seq FPKM", "linear", ""),
}


def build_enlight():
    lab_file = EN / "drug_response_classifications.csv"
    labs = pd.read_csv(lab_file)
    labmap = labs.drop_duplicates("Sample ID").set_index("Sample ID")
    lab_sha = sha256(lab_file)
    for name, (drug, cancer, setting, acc, conf, plat, kind, note) in ENLIGHT.items():
        f = EN / f"{name}.csv"
        raw = pd.read_csv(f, index_col=0)
        expr = tidy_expr(raw)
        expr, units = to_log(expr, kind)
        samples = [s for s in expr.columns if s in labmap.index]
        lab = labmap.loc[samples]
        clin = pd.DataFrame(index=pd.Index(samples, name="sample"))
        clin["drug"] = drug
        clin["arm"] = drug
        clin["response"] = lab["Response"].map({"Responder": "R", "Non-responder": "NR"})
        clin["responder"] = (lab["Response"] == "Responder").astype(int)
        clin["enlight_dataset"] = lab["Dataset"]
        cid = "enlight_" + re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
        meta = dict(title=f"ENLIGHT dataset '{name}'", cancer=cancer, drugs=drug, setting=setting, endpoint="ENLIGHT binary R/NR (per-dataset definition, Dinstag 2023 Table S1)",
                    accession=f"{acc} [accession confidence: {conf}]", platform=plat,
                    units=f"{units}. ENLIGHT notes: raw mRNA matrices as deposited; for several microarray sets the gene value is a SUM of log2 probe values (gene max >> 16), so per-gene standardisation is advisable.",
                    citation="Dinstag G, et al. Clinically oriented prediction of patient response to targeted and immunotherapies from the tumor transcriptome. Med 2023;4:15-30.",
                    sources=[(EN_RAW + name.replace(' ', '%20').replace('+', '%2B') + ".csv", sha256(f)), (EN_RAW + "drug_response_classifications.csv", lab_sha)],
                    source_short="PangeaResearch/enlight-data", has_control=False,
                    notes=("- " + note + "\n") if note else "")
        write_cohort(cid, expr, clin, meta)


# ======================================================================= ISLE (GSE25066, GSE33072)
def build_isle():
    import rdata
    p = rdata.read_rda(str(GH / "ISLE/data/GSE25055.RData"))["prob"]
    genes = np.asarray(p["genes"]).astype(str)
    m = pd.DataFrame(np.asarray(p["mRNA"], dtype=float), index=genes, columns=np.asarray(p["samples"]).astype(str))
    expr = tidy_expr(m)
    sv = p["surv.dt"]
    clin = pd.DataFrame(index=m.columns)
    clin["drug"] = "taxane -> anthracycline (T/FAC or T/FEC)"
    clin["arm"] = "taxane-anthracycline"
    clin["response"] = np.asarray(p["response"]).astype(str)
    clin.loc[~clin["response"].isin(["pCR", "RD"]), "response"] = np.nan
    clin["responder"] = clin["response"].map({"pCR": 1, "RD": 0})
    clin["rcb_class"] = np.asarray(p["class"], dtype=object)
    for k in ["er", "pr", "her2", "subtypes"]:
        clin[k] = np.asarray(p[k], dtype=object)
    clin["age"] = pd.to_numeric(np.asarray(p["age"]), errors="coerce")
    clin["drfs_years"] = sv["time"].values
    clin["drfs_event"] = sv["status"].values
    f = GH / "ISLE/data/GSE25055.RData"
    meta = dict(title="Hatzis 2011 neoadjuvant taxane-anthracycline breast cohort (MDACC discovery+validation)",
                citation="Hatzis C, et al. A genomic predictor of response and survival following taxane-anthracycline chemotherapy for invasive breast cancer. JAMA 2011;305:1873-81.",
                cancer="Breast (HER2- mostly)", drugs="taxane -> anthracycline (T/FAC, T/FEC)", setting="neoadjuvant", endpoint="pCR vs RD (+ DRFS)",
                accession="GSE25066 (= GSE25055 + GSE25065; ISLE file is named GSE25055 but holds all 508)", platform="Affymetrix U133A (log2, gene-level as in ISLE)",
                units="values as distributed in ISLE (log2 MAS5/gene-level).", sources=[("https://github.com/jooslee/ISLE/blob/master/data/GSE25055.RData (git clone)", sha256(f))],
                source_short="jooslee/ISLE", has_control=False, notes="- No untreated arm (single regimen); DRFS in years.\n")
    write_cohort("GSE25066_hatzis", expr, clin, meta)

    # BATTLE GSE33072: erlotinib + sorafenib arms
    parts = []
    for fn, drug in [("GSE33072.RData", "erlotinib"), ("GSE33072.sorafenib.RData", "sorafenib")]:
        p = rdata.read_rda(str(GH / "ISLE/data" / fn))["prob"]
        genes = np.asarray(p["genes"]).astype(str)
        m = pd.DataFrame(np.asarray(p["mRNA"], dtype=float), index=genes, columns=np.asarray(p["samples"]).astype(str))
        c = pd.DataFrame(index=m.columns)
        c["drug"] = drug
        c["arm"] = drug
        c["pfs_months"] = p["surv.dt"]["time"].values
        c["pfs_event"] = p["surv.dt"]["status"].values
        for k in ["egfr", "kras", "sex", "race"]:
            c[k] = np.asarray(p[k], dtype=object)
        parts.append((m, c, GH / "ISLE/data" / fn))
    m = pd.concat([parts[0][0], parts[1][0]], axis=1)
    expr = tidy_expr(m)
    clin = pd.concat([parts[0][1], parts[1][1]])
    labs = pd.read_csv(EN / "drug_response_classifications.csv").set_index("Sample ID")["Response"]
    clin["enlight_label"] = labs.reindex(clin.index).values
    # ENLIGHT sorafenib labels coincide with PFS > 2.5 months (all NR <= 2.46, all R >= 2.56); apply to both arms
    clin["responder"] = (clin["pfs_months"] > 2.5).astype(int)
    clin.loc[(clin["pfs_months"] <= 2.5) & (clin["pfs_event"] == 0), "responder"] = np.nan
    clin["response"] = clin["responder"].map({1: "R", 0: "NR"})
    clin["label_source"] = np.where(clin["enlight_label"].notna(), "ENLIGHT (== PFS>2.5mo)", "derived: PFS>2.5mo")
    meta = dict(title="BATTLE trial NSCLC, erlotinib and sorafenib arms",
                citation="Kim ES, et al. The BATTLE trial. Cancer Discov 2011;1:44-53; Blumenschein GR, et al. Clin Cancer Res 2013 (GSE33072).",
                cancer="NSCLC", drugs="erlotinib; sorafenib", setting="metastatic, previously treated (adaptive randomisation)", control="two active arms (erlotinib vs sorafenib) allow drug-by-biomarker interaction, no placebo",
                endpoint="PFS > 2.5 months (binary; matches ENLIGHT sorafenib labels exactly) + PFS time", accession="GSE33072", platform="Affymetrix HuGene 1.0 ST (log2)",
                units="values as distributed in ISLE (log2 RMA, gene-level).",
                sources=[(f"https://github.com/jooslee/ISLE/blob/master/data/{x[2].name} (git clone)", sha256(x[2])) for x in parts],
                source_short="jooslee/ISLE", has_control=True,
                notes="- Erlotinib labels are DERIVED (PFS>2.5 mo, threshold inferred from ENLIGHT sorafenib labels); BATTLE primary endpoint was 8-week disease control.\n- Erlotinib arm here is n=25 (ISLE subset), not the full BATTLE erlotinib arm. Vandetanib/bexarotene arms absent.\n")
    write_cohort("GSE33072_battle", expr, clin, meta)


# ======================================================================= BrighTNess GSE164458
def build_brightness():
    base = GH / "multi_label_classification/brightness_input/GSE164458"
    files = [("GSE164458_paclitaxel.csv", "paclitaxel (+placebo)"), ("GSE164458_carboplatin_paclitaxel.csv", "carboplatin + paclitaxel"),
             ("GSE164458_veliparib_carboplatin_paclitaxel.csv", "veliparib + carboplatin + paclitaxel")]
    exprs, clins, srcs = [], [], []
    for fn, arm in files:
        d = pd.read_csv(base / fn).set_index("Sample")
        c = d[["pcr", "HR", "HER2"]].copy()
        c["arm"] = arm
        c["drug"] = arm + " -> AC"
        exprs.append(d.drop(columns=["pcr", "HR", "HER2"]).T)
        clins.append(c)
        srcs.append((f"https://github.com/moonchangin/multi_label_classification/blob/main/brightness_input/GSE164458/{fn} (git clone)", sha256(base / fn)))
    expr = pd.concat(exprs, axis=1, join="inner")
    expr.index = [g.replace(".", "-") if re.search(r"\.(AS|IT|DT|OT)\d", g) else g for g in expr.index]
    expr = tidy_expr(expr)
    clin = pd.concat(clins)
    clin["responder"] = clin["pcr"].astype(int)
    clin["response"] = clin["responder"].map({1: "pCR", 0: "RD"})
    clin = clin[["drug", "arm", "response", "responder", "HR", "HER2"]]
    meta = dict(title="BrighTNess phase 3 (TNBC), three randomised arms",
                citation="Loibl S, et al. Lancet Oncol 2018 (BrighTNess); Filho OM, et al. Clin Cancer Res 2021 (GSE164458 RNA-seq).",
                cancer="Breast (TNBC)", drugs="paclitaxel; carboplatin+paclitaxel; veliparib+carboplatin+paclitaxel (all followed by AC)", setting="neoadjuvant",
                control="paclitaxel + placebo arm (randomised 1:1:2)", endpoint="pCR vs RD", accession="GSE164458", platform="RNA-seq (log2-scale, as processed by moonchangin)",
                units="values as distributed (log2-scale RNA-seq expression; genes = intersection of the three per-arm files).", sources=srcs,
                source_short="moonchangin/multi_label_classification", has_control=True,
                notes="- Randomised: treatment-by-biomarker interaction tests possible (carboplatin or veliparib vs paclitaxel control).\n- Gene names with '.' (R make.names) for antisense genes restored to '-' where obvious.\n")
    write_cohort("GSE164458_brightness", expr, clin, meta)


# ======================================================================= I-SPY2 GSE194040 (bhklab ICB_Wolf)
def build_ispy2():
    x = pickle.load(open(SCR / "wolf.pkl", "rb"))
    cd = x.colData.listData
    df = pd.DataFrame({str(k): np.asarray(v) for k, v in cd.items()}, index=np.asarray(x.colData.rownames).astype(str))
    e = x.ExperimentList.listData["expr"]
    ex = e.assays.data.listData["expr"]
    names = np.asarray(e.rowRanges.elementMetadata.listData["gene_name"]).astype(str)
    ex = pd.DataFrame(np.asarray(ex, dtype=float), index=names, columns=np.asarray(e.colData.rownames).astype(str))
    expr = tidy_expr(ex)
    clin = pd.DataFrame(index=df.index)
    clin["drug"] = df["treatmentid"]
    arm_short = df["Arm..short.name."].astype(str)
    clin["arm"] = arm_short
    clin["is_control_arm"] = (arm_short == "Ctr").astype(int)
    clin["response"] = df["pCR"].astype(int).map({1: "pCR", 0: "RD"})
    clin["responder"] = df["pCR"].astype(int)
    for k in ["HR", "HER2", "MP2", "Receptor.Subtype", "PAM50.Subtype", "I.SPY2.Subtypes", "RPS.5"]:
        clin[k.replace(".", "_").rstrip("_")] = df[k].values
    f = GH / "ispy-readii/workflow/Source Data/ICB_Wolf.rds"
    meta = dict(title="I-SPY2 adaptive neoadjuvant breast trial, 10 arms incl. paclitaxel control",
                citation="Wolf DM, et al. Redefining breast cancer subtypes to guide treatment prioritization and maximize response: predictive biomarkers across 10 cancer therapies. Cancer Cell 2022;40:609-623.",
                cancer="Breast (stage II-III, high risk)", drugs="paclitaxel +/- {neratinib, MK-2206, ganitumab, ganetespib, trebananib(AMG386), veliparib+carboplatin, pembrolizumab, pertuzumab/trastuzumab}; T-DM1+pertuzumab; then AC",
                setting="neoadjuvant", control="Ctr = paclitaxel (HER2-) or paclitaxel+trastuzumab (HER2+), n=210", endpoint="pCR vs RD",
                accession="GSE194040", platform="Agilent 44K microarray (gene-level, log2; bhklab ORCESTRA curation)",
                units="values as distributed in bhklab ICB_Wolf MultiAssayExperiment (GEO gene-level matrix, log2 scale); Ensembl rows mapped to gene_name.",
                sources=[("https://github.com/bhklab/ispy-readii/blob/main/workflow/Source%20Data/ICB_Wolf.rds (git clone; bhklab ICB_Wolf curation)", sha256(f))],
                source_short="bhklab/ispy-readii (ICB_Wolf.rds)", has_control=True,
                notes="- Arms are adaptively randomised against a shared control: arm x biomarker interaction tests possible within HER2-/HER2+ strata.\n- Arm short names: Ctr, AMG386 (trebananib), N (neratinib), Ganitumab, MK2206, Ganetespib, VC (veliparib+carboplatin), Pembro, TDM1/P, Pertuzumab.\n")
    write_cohort("GSE194040_ispy2", expr, clin, meta)


# ======================================================================= cBioPortal iAtlas ICI cohorts
CB = SCR / "cbio"
CB_URL = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/"


def read_cbio_expr(s):
    fs = sorted((CB / s).glob("data_mrna_seq_*.txt"))
    f = [x for x in fs if "zscore" not in x.name][0]
    d = pd.read_csv(f, sep="\t", low_memory=False)
    if "Hugo_Symbol" in d:
        d = d.drop(columns=[c for c in ["Entrez_Gene_Id"] if c in d]).set_index("Hugo_Symbol")
    else:
        mp = pd.read_csv(SCR / "hgnc2.txt", sep="\t", usecols=["symbol", "entrez_id"], low_memory=False).dropna()
        mp = dict(zip(mp["entrez_id"].astype(int), mp["symbol"]))
        d["Entrez_Gene_Id"] = d["Entrez_Gene_Id"].map(lambda v: mp.get(int(v)) if pd.notna(v) else None)
        d = d.dropna(subset=["Entrez_Gene_Id"]).set_index("Entrez_Gene_Id")
    return d, f


def cbio(s, cid, title, cancer, drugs, setting, endpoint, platform, citation, arm_fn=None, sample_filter=None, has_control=False, control="none", notes="", resp_map=None, qc_label="responder", extra_src=None):
    expr, f = read_cbio_expr(s)
    pt = pd.read_csv(CB / s / "data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
    sm = pd.read_csv(CB / s / "data_clinical_sample.txt", sep="\t", comment="#").set_index("SAMPLE_ID")
    sm = sm.join(pt, on="PATIENT_ID", rsuffix="_pt")
    if sample_filter is not None:
        sm = sm.loc[sample_filter(sm)]
    sm = sm.loc[[x for x in sm.index if str(x) in set(map(str, expr.columns))]]
    expr.columns = expr.columns.astype(str)
    sm.index = sm.index.astype(str)
    expr = expr[sm.index]
    expr = tidy_expr(expr)
    expr = expr.loc[(expr > 0).mean(axis=1) >= 0.1]
    expr, units = to_log(expr, "linear")
    clin = pd.DataFrame(index=sm.index)
    clin["patient"] = sm["PATIENT_ID"]
    clin["arm"] = arm_fn(sm) if arm_fn else drugs
    clin["drug"] = clin["arm"]
    if resp_map is None:
        rmap = {"Complete Response": "CR", "Partial Response": "PR", "Stable Disease": "SD", "Progressive Disease": "PD"}
        clin["response"] = sm["RESPONSE"].map(rmap)
    else:
        clin["response"] = resp_map(sm)
    clin["responder"] = clin["response"].map({"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    if "CLINICAL_BENEFIT" in sm:
        clin["clinical_benefit"] = sm["CLINICAL_BENEFIT"].map({True: 1, False: 0, "True": 1, "False": 0})
    if "SAMPLE_TREATMENT" in sm:
        clin["sample_timing"] = sm["SAMPLE_TREATMENT"]
    for col, new in [("OS_MONTHS", "os_months"), ("OS_STATUS", "os_event"), ("PFS_MONTHS", "pfs_months"), ("PFS_STATUS", "pfs_event"), ("DFS_MONTHS", "pfs_months"), ("DFS_STATUS", "pfs_event")]:
        if col in sm:
            v = sm[col]
            if "STATUS" in col:
                v = v.astype(str).str.split(":").str[0].map({"1": 1, "0": 0})
            clin[new] = pd.to_numeric(v, errors="coerce")
    for col in ["PRIOR_ICI_RX", "NON_ICI_RX", "NEOICI_RX", "TMB_NONSYNONYMOUS", "BIOPSY_SITE", "IMMUNE_SUBTYPE"]:
        if col in sm:
            clin[col.lower()] = sm[col]
    meta = dict(title=title, citation=citation, cancer=cancer, drugs=drugs, setting=setting, endpoint=endpoint, platform=platform,
                accession=f"cBioPortal datahub `{s}`", control=control,
                units=f"{units}; genes kept if expressed (>0) in >=10% of samples.",
                sources=[(CB_URL + f"{s}/{f.name}", sha256(f)), (CB_URL + f"{s}/data_clinical_sample.txt", sha256(CB / s / "data_clinical_sample.txt")),
                         (CB_URL + f"{s}/data_clinical_patient.txt", sha256(CB / s / "data_clinical_patient.txt"))],
                source_short="cBioPortal/datahub (Git LFS)", has_control=has_control, notes=notes)
    if extra_src:
        meta["sources"].append(extra_src)
    if qc_label != "responder":
        meta["usable_override"] = None
    write_cohort(cid, expr, clin, meta, qc_label=qc_label)


def build_cbio():
    pre = lambda sm: sm["SAMPLE_TREATMENT"].eq("Pre")
    cbio("rcc_iatlas_immotion150_2018", "IMmotion150_rcc", "IMmotion150 randomised phase 2, atezolizumab +/- bevacizumab vs sunitinib", "RCC (clear cell)",
         "atezolizumab; atezolizumab+bevacizumab; sunitinib", "first-line metastatic", "RECIST best response (CR/PR vs SD/PD) + PFS", "RNA-seq TPM (iAtlas harmonised)",
         "McDermott DF, et al. Nat Med 2018;24:749-757.",
         arm_fn=lambda sm: np.where(sm["NON_ICI_RX"].eq("Sunitinib"), "sunitinib", np.where(sm["NON_ICI_RX"].eq("Bevacizumab"), "atezolizumab+bevacizumab", "atezolizumab")),
         has_control=True, control="sunitinib arm (n~89) vs atezolizumab +/- bevacizumab", notes="- Randomised 1:1:1 - arm x biomarker interaction possible.\n")
    cbio("blca_iatlas_imvigor210_2017", "IMvigor210_blca", "IMvigor210 phase 2, atezolizumab in metastatic urothelial carcinoma", "Bladder urothelial",
         "atezolizumab", "metastatic (cohort 1 cisplatin-ineligible + cohort 2 post-platinum)", "RECIST best response + OS", "RNA-seq TPM (iAtlas harmonised, FFPE)",
         "Mariathasan S, et al. Nature 2018;554:544-548; Rosenberg et al. ESMO Open 2024.")
    cbio("mel_iatlas_liu_2019", "Liu2019_mel", "Liu 2019 anti-PD1 melanoma (DFCI)", "Melanoma", "pembrolizumab or nivolumab", "metastatic (+/- prior ipilimumab)",
         "RECIST best response + PFS/OS", "RNA-seq TPM (iAtlas harmonised)", "Liu D, et al. Nat Med 2019;25:1916-1927.",
         arm_fn=lambda sm: np.where(sm["PRIOR_ICI_RX"].notna(), "anti-PD1 (ipi-experienced)", "anti-PD1 (ipi-naive)"))
    cbio("mel_iatlas_gide_2019", "Gide2019_mel", "Gide 2019 anti-PD1 +/- anti-CTLA4 melanoma (pre-treatment)", "Melanoma", "pembrolizumab/nivolumab +/- ipilimumab",
         "metastatic", "RECIST best response + PFS/OS", "RNA-seq (iAtlas harmonised)", "Gide TN, et al. Cancer Cell 2019;35:238-255.",
         arm_fn=lambda sm: np.where(sm["ICI_RX"].astype(str).str.contains("Ipilimumab"), "anti-PD1+anti-CTLA4", "anti-PD1"), sample_filter=pre,
         notes="- Two non-randomised regimens (mono vs combination).\n")
    cbio("mel_iatlas_riaz_nivolumab_2017", "Riaz2017_mel_pre", "Riaz 2017 nivolumab melanoma (pre-treatment biopsies)", "Melanoma", "nivolumab", "metastatic (+/- prior ipilimumab)",
         "RECIST best response + OS", "RNA-seq TPM (iAtlas harmonised)", "Riaz N, et al. Cell 2017;171:934-949.",
         arm_fn=lambda sm: np.where(sm["PRIOR_ICI_RX"].notna(), "nivolumab (ipi-progressed)", "nivolumab (ipi-naive)"), sample_filter=pre)
    cbio("mel_iatlas_hugo_ucla_2016", "Hugo2016_mel", "Hugo 2016 anti-PD1 melanoma", "Melanoma", "pembrolizumab", "metastatic", "RECIST best response + OS",
         "RNA-seq TPM (iAtlas harmonised)", "Hugo W, et al. Cell 2016;165:35-44.")
    cbio("skcm_dfci_2015", "VanAllen2015_mel", "Van Allen 2015 ipilimumab melanoma (RNA subset)", "Melanoma", "ipilimumab", "metastatic", "RECIST best response (X = not evaluable -> NA) + OS/PFS",
         "RNA-seq RPKM (cBioPortal; Entrez IDs mapped to HGNC symbols)", "Van Allen EM, et al. Science 2015;350:207-211.",
         resp_map=lambda sm: sm["DURABLE_CLINICAL_BENEFIT"].where(sm["DURABLE_CLINICAL_BENEFIT"].isin(["CR", "PR", "SD", "PD"])),
         extra_src=("https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/annotation/hgnc_complete_set.txt (HGNC Entrez->symbol map)", sha256(SCR / "hgnc2.txt")),
         notes="- Only 40 of 110 patients have RNA-seq.\n")
    cbio("paad_iatlas_prince_2022", "PRINCE_paad", "PRINCE randomised phase 2, gem/nab-paclitaxel + nivolumab and/or sotigalimab", "Pancreatic adenocarcinoma",
         "gemcitabine+nab-paclitaxel with nivolumab (+/- sotigalimab) or sotigalimab alone", "first-line metastatic", "RECIST best response + OS/PFS",
         "RNA-seq (iAtlas harmonised, FFPE)", "Padron LJ, et al. Nat Med 2022;28:1167-1177.",
         arm_fn=lambda sm: np.where(sm["ICI_RX"].eq("Nivolumab"), "chemo+nivolumab(+/-sotigalimab)", "chemo+sotigalimab"), sample_filter=pre, has_control=True,
         control="randomised arms differ by nivolumab (chemo backbone shared); no chemo-only arm", notes="- Arm label from iAtlas ICI_RX ('None Ici Rx' = sotigalimab+chemo arm).\n")
    cbio("gbm_iatlas_prins_2019", "Cloughesy2019_gbm", "Cloughesy 2019 neoadjuvant vs adjuvant pembrolizumab, recurrent GBM", "GBM (recurrent)", "pembrolizumab (neoadjuvant+adjuvant vs adjuvant only)",
         "recurrent, surgical", "RECIST-like best response (no responders) + OS/PFS", "RNA-seq (iAtlas harmonised)", "Cloughesy TF, et al. Nat Med 2019;25:477-486.",
         arm_fn=lambda sm: np.where(sm["NEOICI_RX"].eq("Pembrolizumab"), "neoadjuvant+adjuvant pembrolizumab", "adjuvant pembrolizumab only"), has_control=True,
         control="adjuvant-only arm (randomised); NOTE arm is confounded with sample timing (neoadjuvant arm samples are on-treatment)",
         notes="- No CR/PR: responder label is all 0; use OS/PFS and arm. QC uses clinical_benefit (CB = SD).\n", qc_label="clinical_benefit")
    cbio("ccrcc_iatlas_choueiri_2016", "Choueiri2016_rcc", "Choueiri 2016 nivolumab RCC biomarker trial", "RCC (clear cell)", "nivolumab", "metastatic", "RECIST best response + OS/PFS",
         "RNA-seq (iAtlas harmonised)", "Choueiri TK, et al. Clin Cancer Res 2016;22:5461-5471.")
    cbio("brca_iatlas_anders_2022", "Anders2022_tnbc", "Anders 2022 pembrolizumab + cyclophosphamide mTNBC", "Breast (TNBC, metastatic)", "pembrolizumab + low-dose cyclophosphamide",
         "metastatic, pretreated", "RECIST best response + OS/PFS", "RNA-seq TPM (iAtlas harmonised)", "Anders CK, et al. J Immunother Cancer 2022 (PMID 35121644).")


# ======================================================================= I-SPY2 pembro TNBC (TimiGP) -- subset of GSE194040; skipped (superseded)
# ======================================================================= lung adjuvant (osun24/nsclc-adj-chemo)
def build_lung():
    big = SCR / "lung_affymetrix.merged.csv"
    m = pd.read_csv(big, index_col=0, low_memory=False)
    clincols = ["Adjuvant Chemo", "Age", "OS_STATUS", "Stage", "Histology", "Race", "Smoked?", "IS_MALE", "OS_MONTHS", "RFS_MONTHS", "PFS_MONTHS"]
    gsm = m.index.str.extract(r"GSM(\d+)")[0].astype(int).values // 1000
    cohorts = {
        "GSE14814_jbr10": (set(pd.read_csv(SCR / "lung_GSE14814_data.csv", index_col=0).index), "JBR.10 randomised adjuvant cisplatin/vinorelbine vs observation", "GSE14814",
                           "Zhu CQ, et al. J Clin Oncol 2010;28:4417-24 (JBR.10).", "Affymetrix U133A (GPL96)", True,
                           "Randomised ACT vs OBS: the canonical predictive-interaction cohort."),
        "GSE68465_lung_adj": (set(m.index[gsm == 1672]), "Director's Challenge lung adenocarcinoma (reprocessed CELs), adjuvant chemo yes/no", "GSE68465",
                              "Shedden K, et al. Nat Med 2008;14:822-7.", "Affymetrix U133A (GPL96)", True,
                              "Adjuvant chemo NOT randomised (observational ACT vs OBS); confounded by stage."),
        "GSE37745_lung_adj": (set(pd.read_csv(SCR / "lung_GSE37745_data.csv", index_col=0).index), "Uppsala NSCLC cohort, adjuvant chemo yes/no", "GSE37745",
                              "Botling J, et al. Clin Cancer Res 2013;19:194-204.", "Affymetrix U133 Plus 2 (GPL570)", True,
                              "Adjuvant chemo NOT randomised."),
    }
    src = [("https://media.githubusercontent.com/media/osun24/nsclc-adj-chemo/main/affymetrix.merged.csv", sha256(big))]
    for cid, (ids, title, acc, cit, plat, ctrl, note) in cohorts.items():
        sub = m.loc[[i for i in m.index if i in ids]]
        expr = sub.drop(columns=clincols).T
        expr = tidy_expr(expr)
        expr.columns = sub.index.str.extract(r"(GSM\d+)")[0].values
        clin = pd.DataFrame(index=expr.columns)
        clin["arm"] = sub["Adjuvant Chemo"].map({"ACT": "adjuvant chemo", "OBS": "observation"}).values
        clin["drug"] = np.where(clin["arm"] == "adjuvant chemo", "platinum-based adjuvant chemo" if cid != "GSE14814_jbr10" else "cisplatin + vinorelbine", "none")
        clin["treated"] = (sub["Adjuvant Chemo"] == "ACT").astype(int).values
        clin["response"] = np.nan
        clin["responder"] = np.nan
        clin["os_months"] = sub["OS_MONTHS"].values
        clin["os_event"] = sub["OS_STATUS"].values
        clin["rfs_months"] = sub["RFS_MONTHS"].values
        for k in ["Age", "Stage", "Histology", "IS_MALE", "Smoked?"]:
            clin[k.lower().replace("?", "")] = sub[k].values
        clin["cel_file"] = sub.index.values
        meta = dict(title=title, citation=cit, cancer="NSCLC (resected stage I-III)", drugs=clin["drug"].iloc[0] if cid == "GSE14814_jbr10" else "platinum-based adjuvant chemo",
                    setting="adjuvant", control="observation arm", endpoint="OS (time-to-event); no response label", accession=acc, platform=plat,
                    units="values as distributed in osun24/nsclc-adj-chemo affymetrix.merged.csv (RMA log2, jointly processed across Affymetrix lung cohorts; symbols collapsed by the repo author).",
                    sources=src, source_short="osun24/nsclc-adj-chemo (Git LFS)", has_control=True, usable_override=True,
                    notes=f"- {note}\n- QC label = `treated` (ACT vs OBS) since there is no response label; a PC aligned with arm would indicate batch/arm confounding.\n- Sample->cohort assignment from the repo's per-GSE clinical files (GSE14814_data.csv, GSE37745_data.csv) or GSM prefix (GSE68465: GSM1672xxx).\n")
        write_cohort(cid, expr, clin, meta, qc_label="treated")


if __name__ == "__main__":
    build_enlight()
    build_isle()
    build_brightness()
    build_ispy2()
    build_cbio()
    build_lung()
    cat = pd.DataFrame(CATALOG_ROWS)
    cat.to_csv(SCR / ("catalog_partial.tsv" if ONLY else "catalog_built.tsv"), sep="\t", index=False)
    print(cat.to_string())
