#!/usr/bin/env python3
"""Build curated GEO chemotherapy / treated-cohort datasets for ATLAS.

For every cohort this writes data/curated/trials/<COHORT_ID>/:
  expression.tsv.gz  genes x samples, first column 'gene' (HGNC symbol), log2 scale.
                     Microarrays: probes mapped to symbols with the GEO platform (GPL) annotation,
                     probes mapping to >1 symbol dropped, median over probes per gene.
  clinical.tsv       one row per sample (sample, arm, drug, drugs, setting, response, responder,
                     survival columns, covariates, randomised, control_arm)
  SOURCE.md          provenance, sha256 of downloaded files, derivations, PCA-vs-label QC
and data/curated/trials/CATALOG_GEO.tsv.

Usage:  python build_geo_trials.py [COHORT_ID ...]      (default: all)
Raw downloads are cached in $GEO_CACHE (default ~/.cache/geo_trials).
Requires pandas, numpy, scikit-learn, xlrd (only for GSE20194) and curl on PATH.
"""
import gzip
import hashlib
import io
import os
import re
import subprocess
import sys
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CACHE = Path(os.environ.get("GEO_CACHE", Path.home() / ".cache" / "geo_trials"))
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
QC_THRESH = 0.35

CLIN_COLS = ["sample", "patient", "arm", "drug", "drugs", "setting", "response", "responder",
             "os_months", "os_event", "pfs_months", "pfs_event", "rfs_months", "rfs_event",
             "dfs_months", "dfs_event", "age", "sex", "stage", "randomised", "control_arm"]

# generic drug-name expansions of common regimen abbreviations
FU_LV = ["FLUOROURACIL", "LEUCOVORIN"]
REGIMENS = {
    "FOLFOX": FU_LV + ["OXALIPLATIN"],
    "FOLFIRI": FU_LV + ["IRINOTECAN"],
    "FOLFIRINOX": FU_LV + ["IRINOTECAN", "OXALIPLATIN"],
    "FAC": ["FLUOROURACIL", "DOXORUBICIN", "CYCLOPHOSPHAMIDE"],
    "FEC": ["FLUOROURACIL", "EPIRUBICIN", "CYCLOPHOSPHAMIDE"],
    "AC": ["DOXORUBICIN", "CYCLOPHOSPHAMIDE"],
    "MVAC": ["METHOTREXATE", "VINBLASTINE", "DOXORUBICIN", "CISPLATIN"],
}


# ----------------------------------------------------------------------------------------- download
def stub(acc):
    return re.sub(r"\d{1,3}$", "nnn", acc)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch(url, dest, tries=6):
    """Download url to dest (cached); validates gzip magic for .gz files."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for i in range(tries):
        r = subprocess.run(["curl", "-sfL", "--retry", "3", "-o", str(tmp), url])
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
            if dest.suffix != ".gz" or tmp.open("rb").read(2) == b"\x1f\x8b":
                tmp.replace(dest)
                return dest
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"download failed: {url}")


def series_matrix(gse, name=None):
    name = name or f"{gse}_series_matrix.txt.gz"
    return fetch(f"{GEO_FTP}/series/{stub(gse)}/{gse}/matrix/{name}", CACHE / gse / name)


def series_suppl(gse, name):
    return fetch(f"{GEO_FTP}/series/{stub(gse)}/{gse}/suppl/{name}", CACHE / gse / name)


def gpl_file(gpl):
    """Curated GEO .annot.gz when it exists, otherwise the full SOFT platform table."""
    try:
        return fetch(f"{GEO_FTP}/platforms/{stub(gpl)}/{gpl}/annot/{gpl}.annot.gz",
                     CACHE / "GPL" / f"{gpl}.annot.gz", tries=2)
    except RuntimeError:
        return fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data",
                     CACHE / "GPL" / f"{gpl}.soft.txt")


# ------------------------------------------------------------------------------------------ parsing
def read_series_matrix(path):
    """Return (sample metadata DataFrame indexed by GSM, expression DataFrame probes x GSM)."""
    meta, rows, header = OrderedDict(), [], None
    chars = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t")
                v = [x.strip('"') for x in v]
                if k.startswith("!Sample_characteristics"):
                    chars.append(v)
                else:
                    kk = k[len("!Sample_"):]
                    while kk in meta:
                        kk += "_"
                    meta[kk] = v
            elif line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                header = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            elif header is not None and not line.startswith("!") and line.strip():
                rows.append(line.rstrip("\n"))
    m = pd.DataFrame(meta)
    m.index = m["geo_accession"]
    # characteristics are parsed per sample by key, since GEO rows can be misaligned across samples
    for v in chars:
        for gsm, x in zip(m.index, v):
            if ": " in x:
                k, val = x.split(": ", 1)
            elif ":" in x:
                k, val = x.split(":", 1)
            else:
                continue
            k = "ch_" + re.sub(r"\s+", " ", k.strip().lower())
            if k not in m:
                m[k] = None
            m.at[gsm, k] = val.strip()
    expr = None
    if rows:
        expr = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", header=None, index_col=0,
                           na_values=["null", "NA", "", "NaN"], low_memory=False)
        expr.columns = header[1:]
        expr.index = expr.index.astype(str).str.strip('"')
        expr = expr.apply(pd.to_numeric, errors="coerce")
    return m, expr


def _ensg_symbols():
    t = pd.read_csv(REPO / "data" / "ENSG_GENESYMBOL.txt", sep="\t").dropna()
    return dict(zip(t.iloc[:, 0], t.iloc[:, 1]))


def probe_to_symbol(gpl):
    """Series probe ID -> single HGNC symbol (probes annotated with several symbols are dropped)."""
    p = gpl_file(gpl)
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    end = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    t = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str, low_memory=False,
                    quoting=3)
    t.columns = [c.strip('"') for c in t.columns]
    col = next((c for c in ["Gene symbol", "Gene Symbol", "Gene.Symbol", "GeneSymbol", "Symbol"]
                if c in t.columns), None)
    if col is not None:
        sym = t[col]
    elif "ENSG_ID" in t.columns:  # Almac ovarian DSA (GPL8414)
        sym = t["ENSG_ID"].map(_ensg_symbols())
    else:
        raise ValueError(f"{gpl}: no symbol column in {list(t.columns)}")
    s = pd.Series(sym.values, index=t["ID"].astype(str).str.strip('"')).dropna().astype(str).str.strip()
    s = s[(s != "") & ~s.str.contains("///|//|,| ", regex=True)]
    return s, p


def collapse(expr, mapping):
    """Median over probes per gene symbol."""
    e = expr.loc[expr.index.intersection(mapping.index)]
    e = e.groupby(mapping.loc[e.index].values).median()
    e.index.name = "gene"
    return e


def to_log2(expr):
    """Return (log2 matrix, note). Data with max > 100 are treated as linear intensities."""
    v = expr.values[np.isfinite(expr.values)]
    if np.nanmax(v) > 100:
        lo = max(np.nanmin(v[v > 0]) if (v > 0).any() else 1.0, 1.0)
        return np.log2(expr.clip(lower=lo)), f"linear values (max {np.nanmax(v):.0f}) -> log2 (floor {lo:g})"
    return expr, f"already log scale (range {np.nanmin(v):.2f}..{np.nanmax(v):.2f})"


def drugs_from(*names):
    out = []
    for n in names:
        for d in REGIMENS.get(n, [n]):
            if d not in out:
                out.append(d)
    return ";".join(out)


def num(x):
    return pd.to_numeric(x, errors="coerce")


# --------------------------------------------------------------------------------------------- QC
def pca_qc(expr, clin):
    x = expr.loc[:, clin["sample"]].astype(float)
    x = x.loc[x.notna().all(axis=1)]
    x = x.loc[x.var(axis=1).sort_values(ascending=False).index[:2000]]
    z = x.sub(x.mean(axis=1), axis=0).T.values
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    pcs = u[:, :3] * s[:3]
    ve = (s ** 2 / (s ** 2).sum())[:3]
    res = {"var_explained": ve}
    for lab in ["responder", "arm"]:
        if lab not in clin:
            continue
        y = clin[lab].reset_index(drop=True)
        ok = y.notna() & (y.astype(str) != "NA")
        if lab == "arm":
            vals = y[ok].astype(str).value_counts()
            if len(vals) != 2:
                continue
            yy = (y[ok].astype(str) == vals.index[0]).astype(int)
        else:
            yy = num(y[ok]).astype(int)
        if yy.nunique() != 2:
            continue
        res[lab] = [roc_auc_score(yy, pcs[ok.values, k]) for k in range(3)]
    return res


# ---------------------------------------------------------------------------------------- writing
def write_cohort(cid, expr, clin, info):
    d = HERE / cid
    d.mkdir(parents=True, exist_ok=True)
    clin = clin.copy()
    clin["sample"] = clin["sample"].astype(str)
    expr = expr.loc[:, [s for s in clin["sample"] if s in expr.columns]]
    clin = clin[clin["sample"].isin(expr.columns)]
    expr = expr.dropna(how="all")
    expr = expr[~expr.index.isna() & (expr.index.astype(str) != "")]
    expr.index.name = "gene"
    cols = [c for c in CLIN_COLS if c in clin] + [c for c in clin if c not in CLIN_COLS]
    clin = clin[cols]
    expr.round(4).to_csv(d / "expression.tsv.gz", sep="\t", compression={"method": "gzip", "mtime": 0})
    clin.to_csv(d / "clinical.tsv", sep="\t", index=False, na_rep="NA")
    size = (d / "expression.tsv.gz").stat().st_size / 1e6
    assert size < 90, f"{cid} expression too large: {size:.1f} MB"

    qc = pca_qc(expr, clin)
    resp = num(clin["responder"]) if "responder" in clin else pd.Series(dtype=float)
    pc_auc = qc.get("responder", [np.nan] * 3)
    arm_auc = qc.get("arm", [np.nan] * 3)
    flag = any(abs(a - 0.5) > QC_THRESH for a in pc_auc if np.isfinite(a))
    flag_arm = any(abs(a - 0.5) > QC_THRESH for a in arm_auc if np.isfinite(a))

    lines = [f"# {cid}", "", f"- GEO accession: {info['accession']}",
             f"- Paper: {info['paper']}", f"- Cancer: {info['cancer']}; setting: {info['setting']}",
             f"- Platform: {info['platform']}", f"- Samples: {expr.shape[1]}; genes: {expr.shape[0]}",
             f"- Responders (responder==1): {int((resp == 1).sum())}; non-responders: {int((resp == 0).sum())}; "
             f"unlabelled: {int(resp.isna().sum()) if len(resp) else expr.shape[1]}",
             f"- Randomised: {info.get('randomised', False)}; control/comparator arm: {info.get('has_control', False)}",
             "", "## URLs", ""]
    lines += [f"- {u}" for u in info["urls"]]
    lines += ["", "## Downloaded files (sha256)", ""]
    lines += [f"- `{Path(p).name}`: {sha256(p)}" for p in info["files"]]
    lines += ["", "## Processing", ""] + [f"- {x}" for x in info["notes"]]
    lines += ["", "## Arms", ""]
    if "arm" in clin:
        for a, n in clin["arm"].astype(str).value_counts().items():
            sub = num(clin.loc[clin["arm"].astype(str) == a, "responder"]) if "responder" in clin else None
            r = "" if sub is None else f" (responders {int((sub == 1).sum())}/{int(sub.notna().sum())} labelled)"
            lines.append(f"- {a}: n={n}{r}")
    lines += ["", "## QC: PCA of top-2000-variance genes (complete genes, centred)", "",
              "PC variance explained: " + ", ".join(f"PC{k+1} {v:.3f}" for k, v in enumerate(qc["var_explained"])),
              "", "| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |", "|---|---|---|---|---|"]
    for lab, a, f in [("responder", pc_auc, flag), ("arm", arm_auc, flag_arm)]:
        if all(not np.isfinite(v) for v in a):
            lines.append(f"| {lab} | NA | NA | NA | not computable |")
        else:
            lines.append(f"| {lab} | " + " | ".join(f"{v:.3f}" for v in a) + f" | {'FLAG' if f else 'ok'} |")
    lines += ["", f"Built by `data/curated/trials/build_geo_trials.py {cid}`.", ""]
    (d / "SOURCE.md").write_text("\n".join(lines))

    drugs = sorted({x for s in clin["drugs"].dropna().astype(str) for x in s.split(";") if x and x != "NA"})
    return {
        "cohort_id": cid, "accession": info["accession"], "cancer": info["cancer"],
        "drugs": ";".join(drugs), "setting": info["setting"], "n": expr.shape[1],
        "responders": int((resp == 1).sum()), "n_labelled": int(resp.notna().sum()) if len(resp) else 0,
        "endpoint": info["endpoint"], "randomised": bool(info.get("randomised", False)),
        "has_control_arm": bool(info.get("has_control", False)),
        "arms": ";".join(clin["arm"].astype(str).value_counts().index) if "arm" in clin else "",
        "platform": info["platform"], "n_genes": expr.shape[0],
        "qc_pc1_auc": round(pc_auc[0], 3), "qc_pc2_auc": round(pc_auc[1], 3), "qc_pc3_auc": round(pc_auc[2], 3),
        "qc_flag": flag, "qc_arm_pc1_auc": round(arm_auc[0], 3), "qc_arm_flag": flag_arm,
        "usable": bool(info.get("usable", True)) and not flag, "note": info.get("short_note", ""),
    }


def microarray(gse, gpl, name=None):
    """Series matrix -> (meta, log2 gene matrix, files, notes)."""
    p = series_matrix(gse, name)
    meta, raw = read_series_matrix(p)
    mapping, gp = probe_to_symbol(gpl)
    raw, lognote = to_log2(raw)
    expr = collapse(raw, mapping)
    notes = [f"Series matrix {p.name} ({raw.shape[0]} probes); {lognote}.",
             f"Probes mapped with {gp.name}; probes with no / multiple symbols dropped; "
             f"median over probes per symbol -> {expr.shape[0]} genes."]
    return meta, expr, [p, gp], notes


def base_clin(meta):
    c = pd.DataFrame({"sample": meta.index})
    c.index = meta.index
    return c


def resp_map(series, mapping):
    return series.map(lambda v: mapping.get(str(v).strip(), np.nan) if pd.notna(v) else np.nan)


def geo_urls(gse):
    return [f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            f"{GEO_FTP}/series/{stub(gse)}/{gse}/"]


# ------------------------------------------------------------------------------------------ cohorts
def GSE72970():
    meta, expr, files, notes = microarray("GSE72970", "GPL570")
    c = base_clin(meta)
    reg = meta["ch_regimen"]
    part = {"FOLFIRI": ["FOLFIRI"], "FOLFOX": ["FOLFOX"], "FOLFIRINOX": ["FOLFIRINOX"],
            "BEVACIZUMAB": ["BEVACIZUMAB"], "ERBITUX": ["CETUXIMAB"], "XELIRI": ["CAPECITABINE", "IRINOTECAN"]}
    c["arm"] = reg
    c["drug"] = reg
    c["drugs"] = reg.map(lambda r: drugs_from(*[d for p in r.split("+") for d in part[p]]))
    c["setting"] = "metastatic"
    c["response"] = meta["ch_response category"]
    c["responder"] = resp_map(meta["ch_response status"], {"R": 1, "NR": 0})
    c["pfs_months"], c["pfs_event"] = num(meta["ch_pfs"]), num(meta["ch_pfs censored"])
    c["os_months"], c["os_event"] = num(meta["ch_os"]), num(meta["ch_os censored"])
    c["age"], c["sex"] = num(meta["ch_age"]), meta["ch_sex"]
    c["center"] = meta["ch_center"]
    c["tumor_location"] = meta["ch_tumor location"]
    c["synchronous_metastasis"] = meta["ch_synchronous metastase"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["responder from 'response status' (R = CR/PR, NR = SD/PD; RECIST best response 'response category').",
              "'pfs censored'/'os censored' = 1 taken as event (1 for 114/124 PFS, i.e. progression observed).",
              "PFS/OS units: months as deposited.",
              "Regimens expanded to generic names (FOLFIRI = FLUOROURACIL;LEUCOVORIN;IRINOTECAN, ERBITUX = CETUXIMAB, "
              "XELIRI = CAPECITABINE;IRINOTECAN)."]
    return expr, c, dict(accession="GSE72970", cancer="colorectal (metastatic)", setting="metastatic",
                         platform="GPL570", endpoint="RECIST response; PFS; OS", files=files, notes=notes,
                         urls=geo_urls("GSE72970"),
                         paper="Del Rio M et al. 2017 J Clin Oncol, PMID 28284171 (also PMID 28659146, 30863148)")


def GSE109211():
    meta, expr, files, notes = microarray("GSE109211", "GPL13938")
    c = base_clin(meta)
    c["arm"] = meta["ch_treatment"].map({"Sor": "sorafenib", "Plac": "placebo"})
    c["drug"] = c["arm"]
    c["drugs"] = c["arm"].map({"sorafenib": "SORAFENIB", "placebo": "PLACEBO"})
    c["setting"] = "adjuvant"
    c["response"] = meta["ch_outcome"].map({"responder": "R", "non-responder": "NR"})
    c["responder"] = resp_map(meta["ch_outcome"], {"responder": 1, "non-responder": 0})
    c["randomised"] = True
    c["control_arm"] = c["arm"] == "placebo"
    c["title"] = meta["title"]
    notes += ["Phase 3 STORM trial (adjuvant sorafenib vs placebo after resection/ablation of HCC); "
              "FFPE tumour, Illumina WG-DASL HumanHT-12 v4 (GPL13938).",
              "response/responder from the deposited 'outcome' (responder vs non-responder as defined by Pinyol et al.; "
              "recurrence-based - see paper). Labels exist in both arms, so the placebo arm is a prognostic control.",
              "No survival times in GEO."]
    return expr, c, dict(accession="GSE109211", cancer="hepatocellular carcinoma", setting="adjuvant",
                         platform="GPL13938", endpoint="responder (recurrence-based, per paper)", files=files,
                         notes=notes, urls=geo_urls("GSE109211"), randomised=True, has_control=True,
                         paper="Pinyol R et al. 2019 J Hepatol, PMID 30108162")


def _er_her2(c, meta, er=None, pr=None, her2=None):
    pn = {"P": "positive", "N": "negative", "positive": "positive", "negative": "negative",
          "ERpos": "positive", "ERneg": "negative"}
    for col, key in [("er", er), ("pr", pr), ("her2", her2)]:
        if key and key in meta:
            c[col] = meta[key].map(lambda v: pn.get(str(v).strip(), np.nan) if pd.notna(v) else np.nan)


def GSE20271():
    meta, expr, files, notes = microarray("GSE20271", "GPL96")
    c = base_clin(meta)
    rnd = meta["ch_randomized (1=fac, 2=t/fac)"].map({"1": "FAC", "2": "T/FAC"})
    rec = meta["ch_treatment received (1=fac, 2=t/fac)"].map({"1": "FAC", "2": "T/FAC"})
    pre = meta["ch_preoperative treatment"].fillna("")
    c["arm"] = rnd
    c["treatment_received"] = rec
    c["drug"] = pre.where(pre != "", rec)

    def drugs(row):
        s = str(row["drug"]).upper()
        s = re.sub(r"\([^)]*POST SX[^)]*\)", "", s)  # taxol given after surgery does not count
        s = re.sub(r"[^;,]*\(ADJ\)", "", s)
        if not s.strip() or s == "NAN":
            s = str(row["treatment_received"]).upper()
        anth = "EPIRUBICIN" if "FEC" in s or "EPIRUBICIN" in s else "DOXORUBICIN"
        out = ["FLUOROURACIL", anth, "CYCLOPHOSPHAMIDE"]
        if re.search(r"TAXOL|PACLITAXEL|T/FAC", s):
            out = ["PACLITAXEL"] + out
        return ";".join(out)
    c["drugs"] = c.apply(drugs, axis=1)
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pcr or rd"]
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"] = num(meta["ch_age"])
    c["sex"] = "F"
    c["stage"] = "T" + meta["ch_prechemo t"].astype(str) + "N" + meta["ch_prechemo n"].astype(str)
    _er_her2(c, meta, "ch_er status", "ch_pr status", "ch_her 2 status")
    c["grade"] = meta["ch_bmn grade"]
    c["race"] = meta["ch_race"]
    c["randomised"] = rnd.notna()
    c["control_arm"] = rnd == "FAC"
    notes += ["arm = randomised arm ('randomized (1=fac, 2=t/fac)'); treatment_received and the free-text "
              "'preoperative treatment' (drug) are also kept. drugs derived from the regimen actually given "
              "(FEC -> EPIRUBICIN instead of DOXORUBICIN; T = weekly PACLITAXEL x12). drugs include preoperative switches "
              "after non-response (14 FAC-arm patients also got paclitaxel before surgery); use 'arm' for ITT.",
              "responder = 1 for pCR, 0 for RD.",
              "Possible sample overlap with other MDACC series (GSE20194/GSE22093) - deduplicate before pooling."]
    return expr, c, dict(accession="GSE20271", cancer="breast", setting="neoadjuvant", platform="GPL96",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE20271"),
                         randomised=True, has_control=True,
                         paper="Tabchy A et al. 2010 Clin Cancer Res, PMID 20829329 (also PMID 23185353)")


def GSE41998():
    meta, expr, files, notes = microarray("GSE41998", "GPL571")
    c = base_clin(meta)
    arm = meta["ch_treatment arm"]
    c["arm"] = arm.map({"Ixabepilone": "AC->ixabepilone", "Paclitaxel": "AC->paclitaxel", "none": "AC only"})
    c["drug"] = c["arm"]
    c["drugs"] = arm.map({"Ixabepilone": "DOXORUBICIN;CYCLOPHOSPHAMIDE;IXABEPILONE",
                          "Paclitaxel": "DOXORUBICIN;CYCLOPHOSPHAMIDE;PACLITAXEL",
                          "none": "DOXORUBICIN;CYCLOPHOSPHAMIDE"})
    c["setting"] = "neoadjuvant"
    c["response"] = resp_map(meta["ch_pcr"], {"Yes": "pCR", "No": "RD"})
    c["responder"] = resp_map(meta["ch_pcr"], {"Yes": 1, "No": 0})
    c["pcr_rcb1"] = resp_map(meta["ch_pcrrcb1"], {"Yes": 1, "No": 0})
    c["ac_response"] = meta["ch_ac response"]
    c["age"], c["sex"] = num(meta["ch_age"]), "F"
    _er_her2(c, meta, "ch_er", "ch_pr", "ch_her2stat")
    c["menopause"] = meta["ch_mnssl"]
    c["tumor_size"] = meta["ch_basetms"]
    c["randomised"] = arm != "none"
    c["control_arm"] = arm == "Paclitaxel"
    notes += ["All patients received 4x AC, then were randomised to ixabepilone or paclitaxel; 'none' = not "
              "randomised / AC only (kept, randomised=False).",
              "responder from 'pcr' (Yes=1, No=0); value '0' (n=20) and blanks are ambiguous and set to NA. "
              "pcr_rcb1 (pCR or RCB-I) and clinical response to AC (ac_response) kept as extra columns.",
              "control_arm = paclitaxel (standard) arm."]
    return expr, c, dict(accession="GSE41998", cancer="breast", setting="neoadjuvant", platform="GPL571",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE41998"),
                         randomised=True, has_control=True,
                         paper="Horak CE et al. 2013 Clin Cancer Res, PMID 23340299")


def GSE32646():
    meta, expr, files, notes = microarray("GSE32646", "GPL570")
    c = base_clin(meta)
    c["arm"] = "P->FEC"
    c["drug"] = "weekly paclitaxel x12 then FEC x4"
    c["drugs"] = "PACLITAXEL;FLUOROURACIL;EPIRUBICIN;CYCLOPHOSPHAMIDE"
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pathologic response pcr ncr"].map({"pCR": "pCR", "nCR": "RD"})
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"], c["sex"] = num(meta["ch_age"]), "F"
    c["stage"] = meta["ch_clinical stage"]
    c["grade"] = meta["ch_histological grade"]
    _er_her2(c, meta, "ch_er status ihc", "ch_pr status ihc", "ch_her2 status fish")
    c["node"] = meta["ch_lymph node status"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Single-arm: paclitaxel (80 mg/m2 weekly x12) followed by FEC x4 (Miyake et al.).",
              "responder = 1 for pCR, 0 for nCR (non-pCR)."]
    return expr, c, dict(accession="GSE32646", cancer="breast", setting="neoadjuvant", platform="GPL570",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE32646"),
                         paper="Miyake T et al. 2012 Cancer Sci, PMID 22320227")


def GSE20194():
    meta, expr, files, notes = microarray("GSE20194", "GPL96")
    xls = series_suppl("GSE20194", "GSE20194_MDACC_Sample_Info.xls.gz")
    files.append(xls)
    c = base_clin(meta)
    code = meta["ch_treatment code"].replace({"NA": np.nan})
    comments = meta.get("ch_treatments comments")

    code_drugs = {
        "TFAC": ["PACLITAXEL", "FAC"], "TFEC": ["PACLITAXEL", "FEC"], "TXFAC": ["PACLITAXEL", "CAPECITABINE", "FAC"],
        "TH/FAC": ["PACLITAXEL", "TRASTUZUMAB", "FAC"], "TH/FEC": ["PACLITAXEL", "TRASTUZUMAB", "FEC"],
        "FECT": ["FEC", "PACLITAXEL"], "FACT": ["FAC", "PACLITAXEL"], "FACT+XRT/X": ["FAC", "PACLITAXEL", "CAPECITABINE"],
        "TONLY": ["PACLITAXEL"], "FAC": ["FAC"], "FEC": ["FEC"], "TFAC/HT": ["PACLITAXEL", "FAC"]}

    def drugs(cd, comment):
        if pd.isna(cd):
            return np.nan
        parts = list(code_drugs[str(cd).upper()])
        cm = str(comment).upper() if pd.notna(comment) else ""
        if "FAC" in parts and "FEC" in cm and "FAC" not in cm:  # comment says FEC although code says FAC
            parts[parts.index("FAC")] = "FEC"
        return drugs_from(*parts)
    c["arm"] = code
    c["drug"] = code if comments is None else code.astype(str) + np.where(comments.notna(), " | " + comments.astype(str), "")
    c["drugs"] = [drugs(cd, cm) for cd, cm in zip(code, comments if comments is not None else [np.nan] * len(code))]
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pcr_vs_rd"]
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"], c["sex"] = num(meta["ch_age"]), "F"
    c["stage"] = "T" + meta["ch_tbefore"].astype(str) + "N" + meta["ch_nbefore"].astype(str)
    _er_her2(c, meta, "ch_er_status", "ch_pr_status", "ch_her2 status")
    c["grade"] = meta["ch_bmngrd"]
    c["race"] = meta["ch_race"]
    c["maqc_set"] = meta["description"].str.replace("MAQC_Distribution_Status: ", "", regex=False)
    c["randomised"], c["control_arm"] = False, False
    notes += ["MAQC-II breast set (MDACC). Characteristics parsed per sample by key (GEO rows are misaligned).",
              "Treatment codes: T = paclitaxel, X = capecitabine (TXFAC: assumed paclitaxel+capecitabine then FAC), FAC/FEC, "
              "H = trastuzumab, /HT = hormone therapy (not listed); drugs derived from code, FAC->FEC when the free-text "
              "comment says FEC. "
              "Most patients received TFAC (non-randomised).",
              "responder = 1 for pCR, 0 for RD.",
              "Overlaps with GSE20271/GSE22093 MDACC samples are possible - deduplicate before pooling.",
              f"Supplementary {xls.name} downloaded for provenance (not needed: same fields as series matrix)."]
    return expr, c, dict(accession="GSE20194", cancer="breast", setting="neoadjuvant", platform="GPL96",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE20194"),
                         paper="Popovici V et al. 2010 Breast Cancer Res (MAQC-II), PMID 20064235")


def GSE22093():
    meta, expr, files, notes = microarray("GSE22093", "GPL96")
    c = base_clin(meta)
    c["arm"] = "FAC/FEC"
    c["drug"] = "5-fluorouracil, doxorubicin (or epirubicin), cyclophosphamide"
    c["drugs"] = "FLUOROURACIL;DOXORUBICIN;CYCLOPHOSPHAMIDE"
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pcr.v.rd"].replace({"NA": np.nan})
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"], c["sex"] = num(meta["ch_age"]), "F"
    c["stage"] = "T" + meta["ch_prechemo t"].astype(str) + "N" + meta["ch_prechemo n"].astype(str)
    _er_her2(c, meta, "ch_er positive vs negative by immunohistochemistry")
    c["tp53"] = meta["ch_p53 status"]
    c["grade"] = meta["ch_bmn.grade"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Anthracycline-based FAC/FEC (no taxane); per-sample anthracycline (doxorubicin vs epirubicin) "
              "not annotated, drugs lists DOXORUBICIN.",
              "responder = 1 for pCR, 0 for RD; 'NA' -> NA.",
              "Overlaps with other MDACC series possible."]
    return expr, c, dict(accession="GSE22093", cancer="breast", setting="neoadjuvant", platform="GPL96",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE22093"),
                         paper="Iwamoto T et al. 2011 J Natl Cancer Inst, PMID 21191116")


def GSE23988():
    meta, expr, files, notes = microarray("GSE23988", "GPL96")
    c = base_clin(meta)
    c["arm"] = "FAC->TX"
    c["drug"] = "FAC x4 then docetaxel + capecitabine x4"
    c["drugs"] = "FLUOROURACIL;DOXORUBICIN;CYCLOPHOSPHAMIDE;DOCETAXEL;CAPECITABINE"
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pcr.v.rd"]
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"], c["sex"] = num(meta["ch_age"]), "F"
    c["stage"] = "T" + meta["ch_prechemo t stage"].astype(str) + "N" + meta["ch_prechemo nodal status"].astype(str)
    _er_her2(c, meta, "ch_er positive vs negative")
    c["grade"] = meta["ch_grade"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Single regimen (US Oncology trial): FAC x4 followed by docetaxel + capecitabine x4.",
              "responder = 1 for pCR, 0 for RD."]
    return expr, c, dict(accession="GSE23988", cancer="breast", setting="neoadjuvant", platform="GPL96",
                         endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE23988"),
                         paper="Iwamoto T et al. 2011 J Natl Cancer Inst, PMID 21191116")


def GSE45670():
    meta, expr, files, notes = microarray("GSE45670", "GPL570")
    meta = meta[meta["ch_tissue"].str.contains("carcinoma")]
    c = base_clin(meta)
    c["arm"] = "CRT"
    c["drug"] = "vinorelbine + cisplatin with concurrent radiotherapy"
    c["drugs"] = "VINORELBINE;CISPLATIN"
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_pathological response to crp"].map(
        {"pathological complete response": "pCR", "not pathological complete response": "RD"})
    c["responder"] = resp_map(c["response"], {"pCR": 1, "RD": 0})
    c["age"], c["sex"] = num(meta["ch_age"]), meta["ch_gender"]
    c["stage"] = meta["ch_tumor stage"]
    c["radiotherapy"] = True
    c["randomised"], c["control_arm"] = False, False
    notes += ["Pre-treatment endoscopic ESCC biopsies (n=28); 10 normal-epithelium samples dropped.",
              "Preoperative chemoradiotherapy (vinorelbine + cisplatin + RT); responder = 1 for pCR."]
    return expr, c, dict(accession="GSE45670", cancer="oesophageal squamous cell carcinoma", setting="neoadjuvant",
                         platform="GPL570", endpoint="pCR", files=files, notes=notes, urls=geo_urls("GSE45670"),
                         paper="Wen J et al. 2014 Dis Esophagus / Oncotarget, PMID 24907633")


def GSE15622():
    meta, expr, files, notes = microarray("GSE15622", "GPL8414")
    meta = meta[meta["title"].str.contains("pre-treatment")]
    c = base_clin(meta)
    c["patient"] = meta["title"].str.extract(r"Patient (\d+)")[0]
    c["arm"] = meta["ch_treatment"].str.lower()
    c["drug"] = c["arm"] + " (single agent, 3 cycles)"
    c["drugs"] = meta["ch_treatment"].map({"Paclitaxel": "PACLITAXEL", "Carboplatin": "CARBOPLATIN",
                                           "Both": "PACLITAXEL;CARBOPLATIN"})
    c["setting"] = "neoadjuvant"
    c["response"] = meta["ch_response"]
    c["responder"] = resp_map(c["response"], {"sensitive": 1, "resistant": 0})
    c["ca125_coefficient"] = num(meta.get("ch_ca-125 coefficient"))
    c["randomised"], c["control_arm"] = False, False
    notes += ["CTCR-OV01: 3 cycles of single-agent paclitaxel or carboplatin before debulking; only pre-treatment "
              "biopsies kept (post-treatment samples dropped). Allocation is not described as randomised in GEO, "
              "so randomised=False; the two single-agent arms are comparators.",
              "response sensitive/resistant (CA-125 based, per GEO); responder = 1 for sensitive.",
              "Almac Ovarian Cancer DSA (GPL8414) probesets annotated with ENSG ids; mapped to HGNC symbols "
              "with data/ENSG_GENESYMBOL.txt."]
    return expr, c, dict(accession="GSE15622", cancer="ovarian (high-grade serous)", setting="neoadjuvant",
                         platform="GPL8414", endpoint="CA-125 response", files=files, notes=notes,
                         urls=geo_urls("GSE15622"), has_control=True,
                         paper="Ahmed AA et al. 2007 Cancer Cell, PMID 18068629 (also PMID 25560085)")


def GSE51373():
    meta, expr, files, notes = microarray("GSE51373", "GPL570")
    c = base_clin(meta)
    c["arm"] = "platinum-taxane"
    c["drug"] = "carboplatin + paclitaxel (first-line, post debulking)"
    c["drugs"] = "CARBOPLATIN;PACLITAXEL"
    c["setting"] = "adjuvant"
    c["response"] = meta["ch_classification"].map({"chemotherapy sensitive": "sensitive",
                                                   "chemotherapy resistant": "resistant"})
    c["responder"] = resp_map(c["response"], {"sensitive": 1, "resistant": 0})
    c["stage"] = meta["ch_stage"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["High-grade serous ovarian cancer; chemo sensitive vs resistant (platinum-free interval based, Koti et al.).",
              "Regimen per paper: carboplatin/paclitaxel after primary debulking (not per-sample annotated)."]
    return expr, c, dict(accession="GSE51373", cancer="ovarian (high-grade serous)", setting="adjuvant",
                         platform="GPL570", endpoint="platinum sensitivity", files=files, notes=notes,
                         urls=geo_urls("GSE51373"), paper="Koti M et al. 2013 BMC Cancer, PMID 24237932")


def GSE63885():
    meta, expr, files, notes = microarray("GSE63885", "GPL570")
    c = base_clin(meta)
    key = lambda s: next(k for k in meta.columns if k.startswith(s))
    chemo = meta[key("ch_adjuwant chemotherapy")].replace({"NA": np.nan})
    c["arm"] = chemo
    c["drug"] = chemo
    c["drugs"] = chemo.map({"taxane/platinum": "PACLITAXEL;PLATINUM", "platinum/cyclophosphamide": "PLATINUM;CYCLOPHOSPHAMIDE"})
    c["setting"] = "adjuvant"
    resp = meta[key("ch_clinical status post 1st line chemotherapy")].replace({"NA": np.nan, "P": "PD"})
    c["response"] = resp
    c["responder"] = resp_map(resp, {"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    c["platinum_sensitivity"] = meta[key("ch_platinium sensitivity")].replace({"NA": np.nan})
    c["dfs_months"] = num(meta[key("ch_dfs")]) / 30.4375
    c["os_months"] = num(meta[key("ch_os")]) / 30.4375
    status = meta[key("ch_clinical status at last follow-up")]
    c["os_event"] = status.map(lambda s: 1 if s == "DOD" else (0 if s in ("AWD", "NED") else np.nan))
    c["histology"] = meta[key("ch_histophatological type")]
    c["stage"] = meta[key("ch_figo_stage")].replace({"NA": np.nan})
    c["grade"] = meta[key("ch_tumor grade")]
    c["residual"] = meta[key("ch_residual tumor size")]
    c["brca1"] = meta[key("ch_brca1 mutation")]
    c["tp53"] = meta[key("ch_somatic tp53 mutation")]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Surgical tumour samples before first-line chemotherapy (taxane/platinum or platinum/cyclophosphamide); "
              "platinum agent (cisplatin vs carboplatin) not specified -> 'PLATINUM'.",
              "response = clinical status after 1st line chemo (CR/PR/SD/P->PD); responder = 1 for CR/PR.",
              "OS/DFS days converted to months (/30.4375); os_event = 1 for DOD (dead of disease), 0 for AWD/NED. "
              "DFS deposited as 0 for non-CR patients, so dfs_months is only meaningful for CR patients; no DFS event "
              "indicator deposited.",
              "platinum_sensitivity (resistant/moderately/highly sensitive) kept as extra column."]
    return expr, c, dict(accession="GSE63885", cancer="ovarian", setting="adjuvant", platform="GPL570",
                         endpoint="clinical response; OS; platinum sensitivity", files=files, notes=notes,
                         urls=geo_urls("GSE63885"),
                         paper="Lisowska KM et al. 2014 Front Oncol, PMID 24478986 (also PMID 27028324)")


def GSE52219():
    meta, expr, files, notes = microarray("GSE52219", "GPL14951")
    c = base_clin(meta)
    c["arm"] = "MVAC"
    c["drug"] = "neoadjuvant MVAC"
    c["drugs"] = drugs_from("MVAC")
    c["setting"] = "neoadjuvant"
    k = next(x for x in meta.columns if x.startswith("ch_response to mvac"))
    c["response"] = meta[k].map({"yes": "R", "no": "NR"})
    c["responder"] = resp_map(meta[k], {"yes": 1, "no": 0})
    c["stage"] = meta["ch_cstage"]
    c["pathologic_stage"] = meta["ch_pstage"]
    c["subtype"] = meta["ch_subset type"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Pre-treatment TURBT of MIBC; responder = downstaged to pT0 or pT1 after neoadjuvant MVAC."]
    return expr, c, dict(accession="GSE52219", cancer="bladder (MIBC)", setting="neoadjuvant", platform="GPL14951",
                         endpoint="downstaging (<=pT1)", files=files, notes=notes, urls=geo_urls("GSE52219"),
                         paper="Choi W et al. 2014 Cancer Cell, PMID 24525232")


def GSE169455():
    p = series_matrix("GSE169455")
    meta, _ = read_series_matrix(p)
    xp = series_suppl("GSE169455", "GSE169455_normalized_by_gene.txt.gz")
    with gzip.open(xp, "rt") as fh:
        txt = [l for l in fh if l.strip() and not l.startswith("#")]
    raw = pd.read_csv(io.StringIO("".join(txt)), sep="\t", index_col=0)
    raw = raw.loc[raw.index.notna(), [c for c in raw.columns if str(c).startswith("BLCA")]]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    raw, lognote = to_log2(raw)
    expr = raw.groupby(raw.index.astype(str)).median()
    t2g = dict(zip(meta["title"], meta.index))
    expr.columns = [t2g.get(c, c) for c in expr.columns]
    c = base_clin(meta)
    c["patient"] = meta["title"]
    c["arm"] = meta["ch_neoadjuvant vs induction chemotherapy"]
    c["drug"] = "cisplatin-based combination chemotherapy (regimen not annotated per sample)"
    c["drugs"] = "CISPLATIN"
    c["setting"] = "neoadjuvant"
    pt = meta["ch_pathologic tn stage rc specimen"].fillna("")
    c["pathologic_stage"] = pt
    pcr = pt.str.match(r"^pT0N0") | pt.str.match(r"^pT0\s*N0")
    known = pt.str.match(r"^p?T")
    c["response"] = pd.Series(np.where(pcr, "pCR", "RD"), index=c.index).where(known)
    c["responder"] = pcr.astype(float).where(known)
    ds = pt.str.match(r"^p?T(0|a|is|1)N0") | pt.str.match(r"^pTisN0") | pt.str.match(r"^pTaN0")
    c["downstaged"] = ds.astype(float).where(known)
    c["stage"] = meta["ch_clinical tnm-staging"]
    c["subtype_lundtax"] = meta["ch_lundtax rna subtype"]
    c["subtype_consensus"] = meta["ch_consensus classifier subtype"]
    c["labeling_kit"] = meta["ch_labeling kit"]
    c["labeling_batch"] = meta["ch_labeling batch"]
    c["randomised"], c["control_arm"] = False, False
    notes = [f"Expression from supplementary {xp.name} (gene-level RMA, Affymetrix Human Gene 1.0 ST, GPL6244; "
             f"rows already gene symbols; duplicated symbols median-collapsed); {lognote}. Columns renamed from "
             f"sample titles to GSM ids.",
             "arm = neoadjuvant vs induction (cN+) cisplatin-based chemotherapy (Sjodahl et al.: mostly "
             "gemcitabine/cisplatin or (dd)MVAC; per-sample regimen not deposited, drugs lists CISPLATIN only).",
             "responder = 1 for pT0N0 at cystectomy (pCR), 0 otherwise; 'downstaged' = <=pT1N0 (incl. pTa/pTis).",
             "Labeling kit / batch kept (two amplification kits) - check batch before pooling."]
    return expr, c, dict(accession="GSE169455", cancer="bladder (MIBC)", setting="neoadjuvant", platform="GPL6244",
                         endpoint="pCR (pT0N0)", files=[p, xp], notes=notes, urls=geo_urls("GSE169455"),
                         paper="Sjodahl G et al. 2022 Eur Urol, PMID 34782206")


def GSE42127():
    meta, expr, files, notes = microarray("GSE42127", "GPL6884")
    c = base_clin(meta)
    chemo = meta["ch_had_adjuvant_chemo"] == "TRUE"
    c["arm"] = np.where(chemo, "adjuvant chemo", "observation")
    c["drug"] = np.where(chemo, "adjuvant platinum-based chemotherapy", "none")
    c["drugs"] = np.where(chemo, "PLATINUM", "NONE")
    c["setting"] = "adjuvant"
    c["os_months"] = num(meta["ch_overall survival months"])
    c["os_event"] = meta["ch_survival status"].map({"D": 1, "A": 0})
    c["age"], c["sex"] = num(meta["ch_age at surgery"]), meta["ch_gender"]
    c["stage"] = meta["ch_final.pat.stage"]
    c["histology"] = meta["ch_histology"]
    c["randomised"] = False
    c["control_arm"] = ~chemo
    notes += ["Resected NSCLC (MDACC); adjuvant chemotherapy yes/no, non-randomised. Regimen not annotated "
              "(platinum doublets per Tang et al.) -> 'PLATINUM'. No response labels; OS only."]
    return expr, c, dict(accession="GSE42127", cancer="NSCLC", setting="adjuvant", platform="GPL6884",
                         endpoint="OS", files=files, notes=notes, urls=geo_urls("GSE42127"), has_control=True,
                         paper="Tang H et al. 2013 Clin Cancer Res, PMID 23357979")


def GSE103479():
    meta, expr, files, notes = microarray("GSE103479", "GPL23985")
    c = base_clin(meta)
    adj = meta["ch_adjuvanttreatment"]
    c["arm"] = adj.map({"Yes": "adjuvant chemo", "No": "surgery only"})
    c["drug"] = adj.map({"Yes": "adjuvant fluoropyrimidine-based chemotherapy", "No": "none"})
    c["drugs"] = adj.map({"Yes": "FLUOROURACIL", "No": "NONE"})
    c["setting"] = "adjuvant"
    c["os_months"] = num(meta["ch_overall survival time"])
    c["os_event"] = meta["ch_status alive.dead"].map({"Dead": 1, "Alive": 0})
    c["pfs_months"] = num(meta["ch_progression-free survival time"])
    c["pfs_event"] = meta["ch_recurrence"].map({"Yes": 1, "No": 0})
    c["age"] = num(meta["ch_age diagnosis"])
    c["sex"] = meta["ch_gender"]
    c["stage"] = meta["ch_stage ii\\iii"] if "ch_stage ii\\iii" in meta else meta.filter(like="ch_stage").iloc[:, 0]
    c["site"] = meta["ch_site"]
    c["cms"] = meta["ch_cms subgroup"]
    c["kras"] = meta["ch_kras (codons 12, 13, 61)"]
    c["braf"] = meta["ch_braf (v600e)"]
    c["randomised"] = False
    c["control_arm"] = adj == "No"
    notes += ["Stage II/III colorectal cancer, adjuvant chemotherapy yes/no (non-randomised); regimen not annotated "
              "(fluoropyrimidine +/- oxaliplatin per Allen et al.) -> 'FLUOROURACIL'.",
              "pfs_event from 'recurrence'; os_event from vital status; times as deposited (months).",
              "Almac Xcel array (GPL23985) symbols from the GEO platform table."]
    return expr, c, dict(accession="GSE103479", cancer="colorectal (stage II/III)", setting="adjuvant",
                         platform="GPL23985", endpoint="OS; PFS", files=files, notes=notes,
                         urls=geo_urls("GSE103479"), has_control=True,
                         paper="Allen WL et al. 2018 JCO Precis Oncol, PMID 30088816")


def GSE19860():
    meta, expr, files, notes = microarray("GSE19860", "GPL570")
    c = base_clin(meta)
    tr = meta["ch_treatment response"]
    c["arm"] = "mFOLFOX6"
    c["drug"] = "mFOLFOX6"
    c["drugs"] = drugs_from("FOLFOX")
    c["setting"] = "metastatic"
    fl = tr.str.extract(r"FL_(Responder|Non_responder)")[0]
    c["response"] = fl.map({"Responder": "R", "Non_responder": "NR"})
    c["responder"] = fl.map({"Responder": 1, "Non_responder": 0})
    bv = tr.str.extract(r"BV_(Responder|Non_reponder)")[0]
    c["bevacizumab_response"] = bv.map({"Responder": "R", "Non_reponder": "NR"})
    c["description"] = meta["description"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Advanced CRC treated with mFOLFOX6; responder from 'FL_Responder' / 'FL_Non_responder'.",
              "Some patients later received bevacizumab; that response is kept in bevacizumab_response."]
    return expr, c, dict(accession="GSE19860", cancer="colorectal (advanced)", setting="metastatic",
                         platform="GPL570", endpoint="response (R/NR)", files=files, notes=notes,
                         urls=geo_urls("GSE19860"),
                         paper="Watanabe T et al. 2011 Int J Cancer (FOLFOX response signature); no PMID in GEO")


def GSE19862():
    meta, expr, files, notes = microarray("GSE19862", "GPL570")
    c = base_clin(meta)
    tr = meta["ch_treatment response"]
    c["arm"] = "bevacizumab-containing"
    c["drug"] = "bevacizumab-containing chemotherapy (backbone not annotated)"
    c["drugs"] = "BEVACIZUMAB"
    c["setting"] = "metastatic"
    c["response"] = tr.str.extract(r"BV_(Responder|Non-responder)")[0].map({"Responder": "R", "Non-responder": "NR"})
    c["responder"] = c["response"].map({"R": 1, "NR": 0})
    c["randomised"], c["control_arm"] = False, False
    notes += ["Advanced CRC, response to bevacizumab-based therapy (BV_Responder / BV_Non-responder). The "
              "chemotherapy backbone is not annotated in GEO; drugs lists BEVACIZUMAB only.",
              "Deposited values are on an unusual scale (about -4.6..7.4, unlike GSE19860 from the same group); kept "
              "as deposited - standardise per gene before use. n=14 only."]
    return expr, c, dict(accession="GSE19862", cancer="colorectal (advanced)", setting="metastatic",
                         platform="GPL570", endpoint="response (R/NR)", files=files, notes=notes,
                         urls=geo_urls("GSE19862"), paper="Watanabe T et al. (Teikyo Univ.); no PMID in GEO")


def GSE62080():
    meta, expr, files, notes = microarray("GSE62080", "GPL570")
    c = base_clin(meta)
    c["arm"] = "FOLFIRI"
    c["drug"] = "FOLFIRI (LV5FU2 + irinotecan), first line"
    c["drugs"] = drugs_from("FOLFIRI")
    c["setting"] = "metastatic"
    r = meta["ch_response"]
    c["response"] = r.map(lambda s: "R" if "S (sensitive)" in s else ("NR" if "R (resistant)" in s else np.nan))
    c["responder"] = c["response"].map({"R": 1, "NR": 0})
    c["randomised"], c["control_arm"] = False, False
    notes += ["First-line FOLFIRI in metastatic CRC; sensitive (objective response) -> responder 1, resistant -> 0."]
    return expr, c, dict(accession="GSE62080", cancer="colorectal (metastatic)", setting="metastatic",
                         platform="GPL570", endpoint="response (sensitive/resistant)", files=files, notes=notes,
                         urls=geo_urls("GSE62080"),
                         paper="Del Rio M et al. 2007 J Clin Oncol, PMID 17327601 (also PMID 30863148)")


CURATED_OV = "https://bioconductor.org/packages/release/data/experiment/src/contrib/curatedOvarianData_1.50.0.tar.gz"


def GSE9891():
    meta, expr, files, notes = microarray("GSE9891", "GPL570")
    tgz = fetch(CURATED_OV, CACHE / "curatedOvarianData_1.50.0.tar.gz")
    ph = CACHE / "GSE9891" / "curatedOvarianData_GSE9891_pheno.tsv"
    if not ph.exists():  # needs R (base only): extract phenoData slot of the ExpressionSet
        subprocess.run(["tar", "-xzf", str(tgz), "-C", str(CACHE), "curatedOvarianData/data/GSE9891_eset.rda"],
                       check=True)
        rda = CACHE / "curatedOvarianData" / "data" / "GSE9891_eset.rda"
        subprocess.run(["Rscript", "-e", f'load("{rda}"); p<-GSE9891_eset@phenoData@data; '
                        f'write.table(p,"{ph}",sep="\\t",quote=F,row.names=T,col.names=NA)'], check=True)
    files += [tgz]
    p = pd.read_csv(ph, sep="\t", index_col=0).drop(columns=["uncurated_author_metadata"], errors="ignore")
    p = p[(p["sample_type"] == "tumor") & p.index.isin(meta.index)]
    c = base_clin(meta.loc[p.index])
    plat, tax = p["pltx"] == "y", p["tax"] == "y"
    c["arm"] = np.select([plat & tax, plat & ~tax, ~plat & p["pltx"].notna()],
                         ["platinum+taxane", "platinum", "no platinum"], default="unknown")
    c["drug"] = c["arm"]
    c["drugs"] = np.select([plat & tax, plat & ~tax], ["PLATINUM;PACLITAXEL", "PLATINUM"], default="NA")
    c["neoadjuvant"] = p["neo"]
    c["setting"] = np.where(p["neo"] == "y", "neoadjuvant", "adjuvant")
    c["rfs_months"] = num(p["days_to_tumor_recurrence"]) / 30.4375
    c["rfs_event"] = p["recurrence_status"].map({"recurrence": 1, "norecurrence": 0})
    c["os_months"] = num(p["days_to_death"]) / 30.4375
    c["os_event"] = p["vital_status"].map({"deceased": 1, "living": 0})
    # proxy platinum response: recurrence-free at 12 months from surgery (~>=6 months after chemotherapy)
    early = (c["rfs_event"] == 1) & (c["rfs_months"] < 12)
    late = c["rfs_months"] >= 12
    c["responder"] = np.where(early, 0.0, np.where(late, 1.0, np.nan))
    c.loc[~plat, "responder"] = np.nan
    c["response"] = c["responder"].map({1.0: "R (RFS>=12m)", 0.0: "NR (RFS<12m)"})
    c["age"], c["sex"] = num(p["age_at_initial_pathologic_diagnosis"]), "F"
    c["stage"] = p["tumorstage"].astype(str) + p["substage"].fillna("").astype(str)
    c["grade"] = p["grade"]
    c["histology"] = p["histological_type"]
    c["primary_site"] = p["primarysite"]
    c["debulking"] = p["debulking"]
    c["batch"] = p["batch"]
    c["randomised"], c["control_arm"] = False, False
    notes += ["Clinical annotation from Bioconductor curatedOvarianData 1.50.0 (GSE9891_eset phenoData; Tothill et al. "
              "clinical data curated by Ganzfried et al. 2013); borderline/LMP tumours (n=18) excluded.",
              "arm from pltx/tax flags (platinum agent not specified -> 'PLATINUM'; taxane assumed paclitaxel).",
              "rfs/os from days_to_tumor_recurrence / days_to_death (/30.4375) with recurrence_status / vital_status.",
              "responder is a DERIVED platinum-response proxy (no RECIST/CA-125 response deposited): 1 = recurrence-free "
              ">=12 months after surgery, 0 = recurrence <12 months, NA if censored <12 months or not platinum-treated.",
              "Array hybridisation batch (date) kept in 'batch'."]
    return expr, c, dict(accession="GSE9891", cancer="ovarian (serous/endometrioid)", setting="adjuvant",
                         platform="GPL570", endpoint="RFS; OS; derived RFS>=12m response", files=files, notes=notes,
                         urls=geo_urls("GSE9891") + [CURATED_OV],
                         paper="Tothill RW et al. 2008 Clin Cancer Res, PMID 18698038")


COHORTS = OrderedDict((f.__name__, f) for f in [
    GSE72970, GSE109211, GSE20271, GSE41998, GSE32646, GSE20194, GSE22093, GSE23988, GSE45670,
    GSE15622, GSE9891, GSE51373, GSE63885, GSE52219, GSE169455, GSE42127, GSE103479, GSE19860, GSE19862, GSE62080])

CATALOG = HERE / "CATALOG_GEO.tsv"


def main(ids):
    rows = []
    old = pd.read_csv(CATALOG, sep="\t") if CATALOG.exists() else pd.DataFrame(columns=["cohort_id"])
    for cid in ids:
        print(f"== {cid}", flush=True)
        expr, clin, info = COHORTS[cid]()
        row = write_cohort(cid, expr, clin, info)
        print({k: row[k] for k in ["n", "responders", "n_labelled", "n_genes", "qc_pc1_auc", "qc_flag"]}, flush=True)
        rows.append(row)
    new = pd.DataFrame(rows)
    keep = old[~old["cohort_id"].isin(new["cohort_id"])]
    cat = pd.concat([keep, new], ignore_index=True) if len(keep) else new
    for col in ["n", "responders", "n_labelled", "n_genes"]:
        cat[col] = cat[col].astype("Int64")
    order = {c: i for i, c in enumerate(COHORTS)}
    cat = cat.sort_values("cohort_id", key=lambda s: s.map(order).fillna(999))
    cat.to_csv(CATALOG, sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:] or list(COHORTS))
