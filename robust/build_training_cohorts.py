"""Prepare the new GEO cohorts for the combined training table.

Reads the files fetched by fetch_geo_matrices.py (data/external/geo_train/<GSE>/) and writes, per cohort,
data/external/train_cohorts/<ID>/expression.tsv.gz (gene symbols x samples) and clinical.tsv in the same
column layout as data/curated/trials. Only cohorts whose public record holds both the treatment and an
outcome are built; the others are listed in SKIPPED with the reason.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "data" / "curated" / "trials"))
import build_geo_trials as B  # noqa: E402
from obd.geo_matrix import read_expr, read_meta  # noqa: E402
from obd.reference import ensembl_to_symbol  # noqa: E402

RAW = HERE.parent / "data" / "external" / "geo_train"
OUT = HERE.parent / "data" / "external" / "train_cohorts"

SKIPPED = {
    "GSE279390": "VESPER bladder trial: gene data only; arm, response and survival are not in the public record",
    "GSE156915": "FOCUS colorectal trial: gene data only; arm and survival are not in the public record",
    "GSE319641": "NeoTRIP breast trial: gene data only; arm and response are not in the public record",
    "GSE234519": "NSABP B-41 breast trial: gene data only; arm and response are not in the public record",
    "GSE116335": "CALGB 40601 breast trial: gene data only; arm and response are not in the public record",
    "GSE87304": "bladder chemo cohort: subtype only; response and survival are not in the public record",
    "GSE235067": "bladder chemo cohort: 66-gene panel and no response in the public record",
    "GSE186102": "eribulin breast cohort: targeted gene panel, only 10 complete responses",
}


def m(gse, name=None):
    return RAW / gse / (name or f"{gse}_series_matrix.txt.gz")


def num(s):
    return pd.to_numeric(s, errors="coerce")


def array_expr(gse, gpl, name=None):
    raw = read_expr(m(gse, name))
    raw.index = raw.index.astype(str)
    raw = raw.apply(pd.to_numeric, errors="coerce")
    mapping = probe_map(gpl)
    raw, _ = B.to_log2(raw)
    return B.collapse(raw, mapping)


def probe_map(gpl):
    try:
        return B.probe_to_symbol(gpl)[0]
    except ValueError:
        import gzip, io
        p = B.gpl_file(gpl)
        txt = (gzip.open(p, "rt", errors="replace") if str(p).endswith(".gz") else open(p, errors="replace")).read().split("\n")
        a = next(i for i, l in enumerate(txt) if l.startswith("!platform_table_begin")) + 1
        z = next((i for i, l in enumerate(txt) if l.startswith("!platform_table_end")), len(txt))
        t = pd.read_csv(io.StringIO("\n".join(txt[a:z])), sep="\t", dtype=str, quoting=3)
        col = next(c for c in ("GENE_SYMBOL", "ILMN_Gene", "Symbol", "SYMBOL") if c in t.columns)
        s = pd.Series(t[col].values, index=t["ID"].astype(str)).dropna().str.strip()
        return s[(s != "") & ~s.str.contains("///|//|,| ", regex=True)]


def by_symbol(e):
    e = e[~e.index.isna()].select_dtypes("number")
    e.index = e.index.astype(str)
    e = e.groupby(level=0).mean()
    e.index.name = "gene"
    return e


def base(meta, cancer, setting, randomised=False):
    c = pd.DataFrame(index=meta.index)
    c.index.name = "sample"
    c["cancer"], c["setting"], c["randomised"], c["control_arm"] = cancer, setting, randomised, False
    return c


# ------------------------------------------------------------------------------------------- breast
ISPY_ARM = {
    "Paclitaxel": "PACLITAXEL", "Paclitaxel + Ganitumab": "PACLITAXEL;GANITUMAB;METFORMIN",
    "Paclitaxel + Ganetespib": "PACLITAXEL;GANETESPIB", "Paclitaxel + AMG 386": "PACLITAXEL;TREBANANIB",
    "Paclitaxel + AMG-386": "PACLITAXEL;TREBANANIB", "Paclitaxel + Pembrolizumab": "PACLITAXEL;PEMBROLIZUMAB",
    "T-DM1 + Pertuzumab": "TRASTUZUMAB EMTANSINE;PERTUZUMAB", "Paclitaxel + MK-2206": "PACLITAXEL;MK-2206",
    "Paclitaxel + Pertuzumab + Trastuzumab": "PACLITAXEL;PERTUZUMAB;TRASTUZUMAB",
    "Paclitaxel + MK-2206 + Trastuzumab": "PACLITAXEL;MK-2206;TRASTUZUMAB",
    "Paclitaxel + AMG 386 + Trastuzumab": "PACLITAXEL;TREBANANIB;TRASTUZUMAB",
    "Paclitaxel + Trastuzumab": "PACLITAXEL;TRASTUZUMAB", "Paclitaxel + Neratinib": "PACLITAXEL;NERATINIB",
    "Paclitaxel + ABT 888 + Carboplatin": "PACLITAXEL;VELIPARIB;CARBOPLATIN",
}


def ISPY2():
    meta = pd.concat([read_meta(m("GSE194040", f"GSE194040-{g}_series_matrix.txt.gz")) for g in ("GPL20078", "GPL30493")])
    e = pd.read_csv(RAW / "GSE194040" / "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz",
                    sep="\t", index_col=0)
    e.columns = e.columns.astype(str)
    meta = meta[meta["patient id"].isin(e.columns)]
    e = e[meta["patient id"]]
    e.columns = meta.index
    c = base(meta, "breast", "neoadjuvant", True)
    c["patient"], c["arm"] = meta["patient id"], meta["arm"].str.replace("AMG-386", "AMG 386")
    c["drugs"] = meta["arm"].map(ISPY_ARM)
    c["control_arm"] = c["arm"].isin(["Paclitaxel", "Paclitaxel + Trastuzumab"])
    c["responder"] = num(meta["pcr"])
    c["er"], c["her2"], c["platform_batch"] = num(meta["hr"]), num(meta["her2"]), meta["platform_id"]
    return by_symbol(e), c, "all arms are followed by doxorubicin/cyclophosphamide; outcome is complete response at surgery"


def BRIGHTNESS():
    meta = read_meta(m("GSE164458"))
    e = pd.read_csv(RAW / "GSE164458" / "GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz", sep="\t", index_col=0)
    ids = meta["title"].str.replace("_RNAseq", "")
    keep = ids[ids.isin(e.columns.astype(str))]
    e.columns = e.columns.astype(str)
    e = e[keep.values]
    e.columns = keep.index
    meta = meta.loc[keep.index]
    c = base(meta, "breast (triple-negative)", "neoadjuvant", True)
    c["patient"] = keep
    c["arm"] = meta["planned_arm_code"].map({"A": "PACLITAXEL+CARBOPLATIN+VELIPARIB", "B": "PACLITAXEL+CARBOPLATIN", "C": "PACLITAXEL"})
    c["drugs"] = meta["planned_arm_code"].map({"A": "PACLITAXEL;CARBOPLATIN;VELIPARIB", "B": "PACLITAXEL;CARBOPLATIN", "C": "PACLITAXEL"})
    c["control_arm"] = meta["planned_arm_code"].eq("C")
    c["responder"] = meta["pathologic_complete_response"].map({"pCR": 1, "RD": 0})
    c["er"], c["her2"] = 0, 0
    return by_symbol(e), c, "all arms are followed by doxorubicin/cyclophosphamide"


def GSE25066():
    meta = read_meta(m("GSE25066"))
    e = array_expr("GSE25066", "GPL96")
    c = base(meta, "breast (HER2-negative)", "neoadjuvant")
    tax = meta["type_taxane"].map({"Taxol": "PACLITAXEL", "Taxotere": "DOCETAXEL"}).fillna("TAXANE")
    c["arm"] = "taxane + anthracycline"
    c["drugs"] = tax + ";ANTHRACYCLINE;CYCLOPHOSPHAMIDE"
    c["responder"] = meta["pathologic_response_pcr_rd"].map({"pCR": 1, "RD": 0})
    c["rfs_months"], c["rfs_event"] = num(meta["drfs_even_time_years"]) * 12, num(meta["drfs_1_event_0_censored"])
    c["age"], c["stage"] = num(meta["age_years"]), meta["clinical_ajcc_stage"]
    c["er"] = meta["er_status_ihc"].map({"P": 1, "N": 0})
    c["her2"] = meta["her2_status"].map({"P": 1, "N": 0})
    c["site"] = meta["source"]
    return e, c, "relapse = distant relapse; may share patients with GSE20194/GSE20271/GSE22093/GSE23988 (same centre)"


def GSE163882():
    meta = read_meta(m("GSE163882"))
    e = pd.read_csv(RAW / "GSE163882" / "GSE163882_all.data.tpms_222Samples.csv.gz", index_col=0)
    e.index = e.index.map(ensembl_to_symbol())
    e = np.log2(by_symbol(e.apply(pd.to_numeric, errors="coerce")) + 1)
    ids = meta["title"].str.split(":").str[0]
    keep = ids[ids.isin(e.columns)]
    e = e[keep.values]
    e.columns = keep.index
    meta = meta.loc[keep.index]
    c = base(meta, "breast", "neoadjuvant")
    c["arm"], c["drugs"] = "neoadjuvant chemotherapy (regimen not given)", "CHEMOTHERAPY"
    c["responder"] = meta["response to nac"].map({"pCR": 1, "RD": 0})
    c["age"], c["stage"] = num(meta["age"]), meta["breast cancer stage"]
    c["er"] = meta["estrogen receptor status"].map({"P": 1, "N": 0})
    c["her2"] = meta["her2 receptor status"].map({"P": 1, "N": 0})
    return e, c, "regimen not recorded per patient"


def GSE123845():
    meta = read_meta(m("GSE123845"))
    meta = meta[meta["time"] == "T1"]
    e = pd.read_csv(RAW / "GSE123845" / "GSE123845_exp_tpm_matrix.csv.gz", index_col=0)
    keep = meta["title"][meta["title"].isin(e.columns)]
    e = np.log2(by_symbol(e[keep.values]) + 1)
    e.columns = keep.index
    meta = meta.loc[keep.index]
    c = base(meta, "breast", "neoadjuvant")
    c["patient"] = meta["title"].str.rsplit("_", n=1).str[0]
    c["arm"] = meta["nac_drug"]
    c["drugs"] = meta["nac_drug"].map({"AC+T": "DOXORUBICIN;CYCLOPHOSPHAMIDE;TAXANE",
                                       "AC+T+Trastuzumab": "DOXORUBICIN;CYCLOPHOSPHAMIDE;TAXANE;TRASTUZUMAB",
                                       "TCHP": "DOCETAXEL;CARBOPLATIN;TRASTUZUMAB;PERTUZUMAB", "Others": "CHEMOTHERAPY"})
    c["responder"] = num(meta["pcr_status"]).map({0: 1, 1: 0})   # 0 = no residual tumour (all size-0 tumours are coded 0)
    c["age"] = num(meta["age"])
    c["er"], c["her2"] = num(meta["er_status_diagnosis"]), num(meta["her2_status_diagnosis"]).where(lambda s: s < 2)
    return e, c, "pre-treatment samples only; response coding inferred from residual tumour size (0 = complete response)"


def GSE240671():
    meta = read_meta(m("GSE240671"))
    meta = meta[meta["timing biopsies"] == "pre-surgery"]
    e = pd.read_csv(RAW / "GSE240671" / "GSE240671_raw_count_all_libraries.csv.gz", sep="\t", index_col=0)
    e = e.drop(columns=["Chromosome", "Region", "Identifier"])
    keep = meta["description"][meta["description"].isin(e.columns)]
    e = by_symbol(e[keep.values])
    e = np.log2(e / e.sum() * 1e6 + 1)
    e.columns = keep.index
    meta = meta.loc[keep.index]
    c = base(meta, "breast", "neoadjuvant")
    c["patient"] = meta["patient id"]
    c["arm"] = meta["nac category"]
    her = np.where(meta["nac herceptin_binary"] == "yes", ";TRASTUZUMAB", "")
    c["drugs"] = meta["nac category"].map({"anthracyclines and taxanes": "ANTHRACYCLINE;TAXANE", "anthracyclin alone": "ANTHRACYCLINE",
                                           "taxanes alone": "TAXANE", "others": "CHEMOTHERAPY"}) + her
    c["responder"] = meta["rcb category"].eq("no residual disease").astype(int)
    c["age"] = num(meta["age diagnostic"])
    c["er"] = meta["molecular category"].isin(["Luminal A", "Luminal B", "HER2+_HR+"]).astype(int)
    c["her2"] = meta["molecular category"].str.contains("HER2").astype(int)
    return e, c, "counts converted to log2 counts-per-million; response = no residual disease"


def SCANB():
    meta = pd.concat([read_meta(m("GSE96058", f"GSE96058-{g}_series_matrix.txt.gz")) for g in ("GPL11154", "GPL18573")])
    f = RAW / "GSE96058" / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
    cols = pd.read_csv(f, nrows=0).columns
    use = [cols[0]] + [x for x in cols[1:] if not x.endswith("repl")]
    e = pd.read_csv(f, index_col=0, usecols=use, dtype={x: "float32" for x in use[1:]})
    keep = meta["title"][meta["title"].isin(e.columns)]
    e = by_symbol(e[keep.values])
    e.columns = keep.index
    meta = meta.loc[keep.index]
    c = base(meta, "breast", "adjuvant")
    ch, en = num(meta["chemo treated"]), num(meta["endocrine treated"])
    c["arm"] = np.select([(ch == 1) & (en == 1), ch == 1, en == 1, (ch == 0) & (en == 0)],
                         ["chemo + endocrine", "chemo", "endocrine", "no systemic therapy"], default="unknown")
    c["drugs"] = np.select([(ch == 1) & (en == 1), ch == 1, en == 1, (ch == 0) & (en == 0)],
                           ["CHEMOTHERAPY;ENDOCRINE", "CHEMOTHERAPY", "ENDOCRINE", ""], default="")
    c["control_arm"] = (ch == 0) & (en == 0)
    c["os_months"], c["os_event"] = num(meta["overall survival days"]) / 30.44, num(meta["overall survival event"])
    c["age"] = num(meta["age at diagnosis"])
    c["er"], c["her2"] = num(meta["er status"]), num(meta["her2 status"])
    c["grade"], c["node_positive"] = meta["nhg"], meta["lymph node status"].map({"NodePositive": 1, "NodeNegative": 0})
    c["tumor_size_mm"] = num(meta["tumor size"])
    return e, c, "treatment recorded as chemo / endocrine yes-no only, not by drug; not randomised"


# --------------------------------------------------------------------------------------- colorectal
def GSE39582():
    meta = read_meta(m("GSE39582"))
    meta = meta[meta["dataset"] != "Non Tumoral"]
    e = array_expr("GSE39582", "GPL570")[meta.index]
    c = base(meta, "colon", "adjuvant")
    y, t = meta["chemotherapy.adjuvant"], meta["chemotherapy.adjuvant.type"]
    c["arm"] = np.where(y == "Y", t.where(t != "N/A", "chemo (type not given)"), np.where(y == "N", "no chemotherapy", "unknown"))
    c["drugs"] = np.where(y == "Y", t.map({"5FU": "FLUOROURACIL", "FUFOL": "FLUOROURACIL;LEUCOVORIN",
                                           "FOLFOX": "FLUOROURACIL;LEUCOVORIN;OXALIPLATIN",
                                           "FOLFIRI": "FLUOROURACIL;LEUCOVORIN;IRINOTECAN"}).fillna("CHEMOTHERAPY"), "")
    c["control_arm"] = y.eq("N")
    c.loc[y.eq("N/A"), "drugs"] = np.nan
    c["rfs_months"], c["rfs_event"] = num(meta["rfs.delay"]), num(meta["rfs.event"])
    c["os_months"], c["os_event"] = num(meta["os.delay (months)"]), num(meta["os.event"])
    c["age"], c["sex"], c["stage"] = num(meta["age.at.diagnosis (year)"]), meta["sex"], meta["tnm.stage"].replace("N/A", np.nan)
    return e, c, "not randomised: chemotherapy was given mostly to stage III patients"


def GSE14333():
    meta = read_meta(m("GSE14333"))
    s = meta["location"].str.extract(r"DukesStage: (\w);\s*Age_Diag: ([\d.NA]+);\s*Gender: (\w);\s*DFS_Time: ([\d.NA]+);\s*DFS_Cens: (\w+);\s*AdjXRT: (\w);\s*AdjCTX: (\w)")
    e = array_expr("GSE14333", "GPL570")
    c = base(meta, "colorectal", "adjuvant")
    c["arm"] = s[6].map({"Y": "adjuvant chemotherapy (type not given)", "N": "no chemotherapy"})
    c["drugs"] = s[6].map({"Y": "CHEMOTHERAPY", "N": ""})
    c["control_arm"] = s[6].eq("N")
    c["dfs_months"] = num(s[3])
    c["dfs_event"] = num(s[4]).map({0: 1, 1: 0})   # DFS_Cens 0 = relapse (rises with stage, shorter times)
    c["age"], c["sex"] = num(s[1]), s[2].map({"M": "Male", "F": "Female"})
    c["stage"] = s[0].map({"A": "1", "B": "2", "C": "3", "D": "4"})
    c["radiotherapy"] = s[5]
    return e, c, "regimen not recorded (5-FU-based by era); event coding inferred from stage pattern; shares patients with GSE17538"


def GSE104645():
    meta = read_meta(m("GSE104645"))
    e = array_expr("GSE104645", "GPL6480")
    c = base(meta, "colorectal (metastatic)", "metastatic")
    reg = meta["1st-line chemotherapy regimens"]
    c["arm"] = reg
    fu = "FLUOROURACIL;LEUCOVORIN"
    parts = {"FOLFOX": fu + ";OXALIPLATIN", "FOLFIRI": fu + ";IRINOTECAN", "IRIS": "IRINOTECAN;TEGAFUR", "SOX": "TEGAFUR;OXALIPLATIN",
             "XELOX": "CAPECITABINE;OXALIPLATIN", "Others": "CHEMOTHERAPY"}
    add = {"Bev": "BEVACIZUMAB", "Pmab": "PANITUMUMAB", "Cmab": "CETUXIMAB"}
    c["drugs"] = [";".join([parts.get(r.split("+")[0], "CHEMOTHERAPY")] + [add[x] for x in r.split("+")[1:] if x in add]) for r in reg]
    br = meta["best response of 1st-line chemotherapy"]
    c["response"] = br.map({"Complete response": "CR", "Partial response": "PR", "Stable disease": "SD", "Progressive disease": "PD"})
    c["responder"] = c["response"].map({"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    c["sex"] = meta["gender"]
    return e, c, ""


def GSE87211():
    meta = read_meta(m("GSE87211"))
    meta = meta[meta["tissue"] == "rectal tumor"]
    e = array_expr("GSE87211", "GPL13497")[meta.index]
    c = base(meta, "rectal", "neoadjuvant chemoradiotherapy")
    r = meta["preoperative radiochemotherapy (rct)"]
    c["arm"] = r
    c["drugs"] = r.map({"5-FU + RT": "FLUOROURACIL;RADIOTHERAPY", "5-FU + Oxaliplatin + RT": "FLUOROURACIL;OXALIPLATIN;RADIOTHERAPY",
                        "5-FU + Oxaliplatin + Cetuximab + RT": "FLUOROURACIL;OXALIPLATIN;CETUXIMAB;RADIOTHERAPY"})
    c["control_arm"] = r.eq("5-FU + RT")
    t, n = num(meta["depth of invasion after rct"]), num(meta["lymph node metastasis after rct"])
    c["responder"] = ((t == 0) & (n == 0)).astype(float).where(t.notna() & n.notna())
    c["dfs_months"], c["dfs_event"] = num(meta["disease free time (month)"]), num(meta["cancer recurrance after surgery"])
    c["age"], c["sex"] = num(meta["age"]), meta["gender"].map({"M": "Male", "F": "Female"})
    return e, c, "drug given together with radiotherapy; response = no tumour left at surgery (ypT0N0)"


def GSE28702():
    meta = read_meta(m("GSE28702"))
    e = array_expr("GSE28702", "GPL570")
    c = base(meta, "colorectal (metastatic)", "metastatic")
    c["arm"], c["drugs"] = "mFOLFOX6", "FLUOROURACIL;LEUCOVORIN;OXALIPLATIN"
    c["responder"] = meta["mfolfox6"].map({"responder": 1, "non-responder": 0})
    c["sex"] = meta["gender"].map({"M": "Male", "F": "Female"})
    c["sample_site"] = meta["lesion"]
    return e, c, ""


# ------------------------------------------------------------------------------------------ bladder
MVAC = "METHOTREXATE;VINBLASTINE;DOXORUBICIN;CISPLATIN"


def GSE70691():
    meta = read_meta(m("GSE70691"))
    e = array_expr("GSE70691", "GPL14951")
    c = base(meta, "bladder (muscle-invasive)", "perioperative", True)
    c["arm"], c["drugs"] = meta["treatment"], MVAC
    c["os_months"], c["os_event"] = num(meta["survival (months)"]), meta["censor.os"].map({"Uncensored": 1, "Censored": 0})
    return e, c, "everyone received MVAC (before or after surgery); may share patients with GSE48276/GSE52219"


def GSE48276():
    meta = read_meta(m("GSE48276"))
    meta = meta[meta["nac"].isin(["no", "yes-MVAC/A"]) & meta["os censor"].notna()]
    e = array_expr("GSE48276", "GPL14951")[meta.index]
    c = base(meta, "bladder (muscle-invasive)", "neoadjuvant")
    c["arm"] = meta["nac"].map({"no": "no neoadjuvant chemotherapy", "yes-MVAC/A": "MVAC"})
    c["drugs"] = meta["nac"].map({"no": "", "yes-MVAC/A": MVAC})
    c["control_arm"] = meta["nac"].eq("no")
    c["os_months"], c["os_event"] = num(meta["survival (mo)"]), meta["os censor"].map({"uncensored": 1, "censored": 0})
    c["age"], c["sex"], c["stage"] = num(meta["age (at diagnosis)"]), meta["gender"].map({"M": "Male", "F": "Female"}), meta["cstage"]
    return e, c, "not randomised; may share patients with GSE70691/GSE52219"


def GSE31684():
    meta = read_meta(m("GSE31684"))
    e = array_expr("GSE31684", "GPL570")
    c = base(meta, "bladder", "adjuvant / salvage")
    y = meta["post rc_chemo"].str.lower()
    c["arm"] = np.where(y == "yes", "chemo after surgery (" + meta["postchemo type"].fillna("?") + ")", "no chemo after surgery")
    c["drugs"] = np.where(y == "yes", "CHEMOTHERAPY", "")
    c["control_arm"] = y.eq("no")
    c["os_months"], c["os_event"] = num(meta["survival.months"]), meta["last known status"].isin(["DOD", "DOC"]).astype(int)
    c["rfs_months"] = num(meta["recurrence free survival months (distant and local)"])
    c["rfs_event"] = meta["recurrence/dod"].map({"Yes": 1, "No": 0})
    c["age"], c["sex"], c["stage"] = num(meta["age at rc"]), meta["gender"].str.capitalize(), meta["rc_stage"]
    return e, c, "regimen not recorded; salvage chemo was given after relapse, so it is not a clean treated group"


BUILDERS = [ISPY2, BRIGHTNESS, GSE25066, GSE163882, GSE123845, GSE240671, SCANB,
            GSE39582, GSE14333, GSE104645, GSE87211, GSE28702, GSE70691, GSE48276, GSE31684]

if __name__ == "__main__":
    want = set(sys.argv[1:])
    rows = []
    for fn in BUILDERS:
        cid = fn.__name__
        if want and cid not in want:
            continue
        d = OUT / cid
        expr, clin, note = fn()
        expr = expr.loc[:, clin.index].dropna(how="all")
        clin["drugs"] = clin["drugs"].replace("", "NONE")   # NONE = known to have had no drug; empty = unknown
        d.mkdir(parents=True, exist_ok=True)
        expr.round(4).to_csv(d / "expression.tsv.gz", sep="\t")
        clin.to_csv(d / "clinical.tsv", sep="\t")
        (d / "NOTE.txt").write_text(note + "\n")
        ev = [x for x in ("os_event", "rfs_event", "dfs_event", "pfs_event") if x in clin]
        rows.append(dict(cohort=cid, samples=len(clin), genes=len(expr),
                         with_response=int(clin["responder"].notna().sum()) if "responder" in clin else 0,
                         responders=int(clin["responder"].sum()) if "responder" in clin else 0,
                         with_survival=int(clin[ev[0]].notna().sum()) if ev else 0,
                         events=int(clin[ev[0]].sum()) if ev else 0, note=note))
        print(rows[-1], flush=True)
