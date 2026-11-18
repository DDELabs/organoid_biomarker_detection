"""Load a user's own organoid expression + drug-response files and turn them into a lab-model set.

What it does (plain words):
  1. Finds the files in one folder (expression, drug response, optional samples table).
  2. Works out by itself: genes in rows or in columns, gene symbols or Ensembl ids, raw counts / TPM-FPKM / log
     values, long or wide drug table.
  3. Scores the same Reactome pathways the rest of this repo uses (rank score, >= 5 genes per pathway; 673 when
     the file covers the full gene universe) with exactly the same function as the patient table
     (robust/build_training_table.score).
  4. Writes the same layout as data/curated/preclinical/<set>/ plus pathway_scores.tsv.gz, organoid_set.json and a
     plain-language IMPORT_REPORT.md that lists what was found and every problem.

Use from Python:   import_folder(input_folder, set_name, out_dir, cancer=None, ...)
Use from the shell: python3 nar/21_user_organoids/code/import_organoids.py <input_folder> <set_name>
"""
import datetime
import gzip
import hashlib
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
REPO = FOLDER.parents[1]
sys.path.insert(0, str(REPO / "robust"))
warnings.filterwarnings("ignore")

from obd.cohorts import to_canonical  # noqa: E402
from obd.reference import common_drug_name, ensembl_to_symbol, gene_sets, uniprot_maps  # noqa: E402
from obd.scoring import rank_score  # noqa: E402

EXTS = (".tsv", ".txt", ".csv", ".xlsx", ".xls", ".tab")
ENSG = re.compile(r"^ENSG\d{6,}(\.\d+)?$", re.I)
DATE_MANGLED = re.compile(r"^\d{1,2}-(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)$|^(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)-\d{1,2}$", re.I)
N_REACTOME_FULL = 673

SAMPLE_COLS = ("sample", "sample_id", "sampleid", "model", "model_id", "organoid", "organoid_line", "organoid_id", "line",
               "cell_line", "name", "id", "patient")
DRUG_COLS = ("drug", "drug_name", "drugname", "compound", "compound_name", "treatment", "agent")
VALUE_COLS = ("value", "response", "auc", "ic50", "ln_ic50", "gi50", "ec50", "viability", "score", "activity", "dss", "metric_value")
METRIC_COLS = ("metric", "measure", "readout", "response_type", "response_metric")
LOWER_WORDS = re.compile(r"auc|ic50|ic_50|gi50|gr50|ec50|viab|survival|growth|relative|ln_|log", re.I)
HIGHER_WORDS = re.compile(r"inhib|kill|death|cytotox|dss|sensitiv|efficacy|effect|apopto|cell_death", re.I)
IC_LIKE = re.compile(r"ic50|ic_50|gi50|gr50|ec50", re.I)
PRIMARY_ORDER = ("AUC", "IC50", "GI50", "EC50")


class ImportProblem(Exception):
    """The folder cannot be turned into a lab-model set; the message says why in plain words."""


class Report:
    def __init__(self):
        self.found, self.did, self.warn, self.problem = [], [], [], []

    def f(self, t):
        self.found.append(t)

    def d(self, t):
        self.did.append(t)

    def w(self, t):
        self.warn.append(t)

    def p(self, t):
        self.problem.append(t)


# ------------------------------------------------------------------------------------------------------ reading
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _base(path):
    """File name without .gz and without the table extension, lower case."""
    n = path.name.lower()
    if n.endswith(".gz"):
        n = n[:-3]
    return n, Path(n).suffix


def read_table(path, rep=None):
    """Read tsv/csv/txt/xlsx (optionally .gz) into a DataFrame with a header row; delimiter is detected."""
    name, ext = _base(path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=0)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        head = [next(fh, "") for _ in range(5)]
    txt = "".join(head)
    counts = {"\t": txt.count("\t"), ",": txt.count(","), ";": txt.count(";")}
    sep = max(counts, key=counts.get)
    if counts[sep] == 0:
        sep = "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)


def find_files(folder, rep):
    """Pick the expression, drug-response and (optional) samples file by name, then by content."""
    files = sorted(p for p in Path(folder).iterdir() if p.is_file() and _base(p)[1] in EXTS)
    if not files:
        raise ImportProblem(f"No table files (.tsv .csv .txt .xlsx, optionally .gz) were found in {folder}.")
    out = {}
    rules = {"samples": r"^(samples?|sample_info|metadata|meta|annotation|pheno)", "response": r"(drug|response|ic50|auc|viab|sensitiv|dose)",
             "expression": r"(expr|rna|tpm|fpkm|counts?|gene|transcript|matrix)"}
    for kind in ("samples", "response", "expression"):
        hits = [p for p in files if re.search(rules[kind], _base(p)[0]) and p not in out.values()]
        if kind == "expression":                            # 'drug_response' must not be taken as expression
            hits = [p for p in hits if not re.search(rules["response"], _base(p)[0])] or hits
        if hits:
            out[kind] = hits[0]
            if len(hits) > 1:
                rep.w(f"More than one file looks like the {kind} file ({', '.join(h.name for h in hits)}). I used {hits[0].name}.")
    left = [p for p in files if p not in out.values()]
    if "expression" not in out and left:                    # fall back on content: the widest/longest numeric table
        best = max(left, key=lambda p: p.stat().st_size)
        out["expression"] = best
        rep.w(f"No file name said 'expression'. I assumed {best.name} (the largest file) is the expression table.")
        left.remove(best)
    if "response" not in out and left:
        out["response"] = left[0]
        rep.w(f"No file name said 'drug response'. I assumed {left[0].name} is the drug-response table.")
    for k in ("expression", "response"):
        if k not in out:
            raise ImportProblem(f"I could not find the {k} file in {folder}. Name it like '{'expression' if k == 'expression' else 'drug_response'}.tsv' "
                                "(csv, xlsx and .gz also work).")
    for k, p in out.items():
        rep.f(f"{k.capitalize()} file: {p.name} ({p.stat().st_size / 1e6:.2f} MB).")
    return out


# --------------------------------------------------------------------------------------------- gene-id helpers
_G2U = None


def _known_symbol(s):
    global _G2U
    if _G2U is None:
        _G2U = uniprot_maps()[0]
    return s in _G2U


def gene_like_fraction(labels):
    labels = [str(x).strip() for x in labels]
    if not labels:
        return 0.0
    hit = sum(1 for x in labels if ENSG.match(x) or _known_symbol(x))
    return hit / len(labels)


def key(s):
    """Loose comparison key for sample names: 'P-19 TA' == 'p19ta'."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


# ----------------------------------------------------------------------------------------------- expression
def parse_expression(df, rep, value_type=None):
    """Return (genes x samples DataFrame of log-scale values with gene symbols, info dict)."""
    info = {}
    df = df.copy()
    first = df.columns[0]
    labels = df[first].astype(str).str.strip()
    body = df.drop(columns=[first])
    # annotation columns (text) are set aside; a column called symbol/gene_name can rescue Ensembl rows
    alias = None
    num = body.apply(pd.to_numeric, errors="coerce")
    text_cols = [c for c in body.columns if num[c].notna().mean() < 0.8]
    for c in text_cols:
        if re.search(r"symbol|gene_?name|hgnc", str(c), re.I):
            alias = body[c].astype(str).str.strip()
    if text_cols:
        rep.w(f"Ignored {len(text_cols)} text column(s) in the expression file that are not samples: {', '.join(map(str, text_cols[:5]))}.")
    M = num.drop(columns=text_cols)
    M.index = labels.values
    if M.shape[1] == 0 or M.shape[0] == 0:
        raise ImportProblem("The expression file has no numeric values.")
    row_g, col_g = gene_like_fraction(M.index), gene_like_fraction(M.columns)
    info["gene_like_rows"], info["gene_like_columns"] = round(row_g, 3), round(col_g, 3)
    if max(row_g, col_g) < 0.3:
        raise ImportProblem("I could not recognise gene names or Ensembl ids in the expression file (only "
                            f"{max(row_g, col_g):.0%} of the labels look like genes). If your genes are Entrez numbers or probe ids, "
                            "convert them to gene symbols or Ensembl ids first.")
    if row_g >= col_g:
        info["orientation"] = "genes in rows, samples in columns"
        if alias is not None:
            alias = pd.Series(alias.values, index=M.index)
    else:
        info["orientation"] = "samples in rows, genes in columns (flipped by me)"
        M = M.T
        alias = None
    rep.f(f"Expression table is laid out as: {info['orientation']} ({M.shape[0]} genes x {M.shape[1]} samples).")
    M.columns = [str(c).strip() for c in M.columns]
    if M.columns.duplicated().any():
        d = M.columns[M.columns.duplicated()].unique()
        rep.w(f"{len(d)} sample name(s) appear more than once in the expression file ({', '.join(d[:5])}); I took the median of the copies.")
        M = M.T.groupby(level=0).median().T
    # ---- gene id type
    g = pd.Series(M.index.astype(str), index=M.index)
    frac_ens = g.map(lambda x: bool(ENSG.match(x))).mean()
    frac_num = g.str.fullmatch(r"\d+").mean()
    if frac_num > 0.5:
        raise ImportProblem("The gene labels are plain numbers (Entrez ids). This repo has no offline Entrez map, so please convert them "
                            "to gene symbols or Ensembl ids and try again.")
    n_date = int(g.map(lambda x: bool(DATE_MANGLED.match(x))).sum())
    if n_date:
        rep.w(f"{n_date} gene label(s) look like dates (Excel turned names such as SEPT2 or MARCH1 into dates). They cannot be recovered and are dropped.")
    # ---- value type, decided on the raw numbers before genes are removed
    vals = M.values.astype(float)
    finite = vals[np.isfinite(vals)]
    vmin, vmax = float(finite.min()), float(finite.max())
    nonneg = vmin >= 0
    intlike = nonneg and np.mean(np.abs(finite[:200000] - np.round(finite[:200000])) < 1e-9) > 0.99
    if value_type:
        vt = value_type
        rep.f(f"Value type given by you: {vt}.")
    elif intlike and vmax > 100:
        vt = "counts"
    elif nonneg and vmax > 50:
        vt = "tpm"
    else:
        vt = "log"
    if not value_type:
        rep.f({"counts": f"Values look like raw read counts (whole numbers up to {vmax:,.0f}).",
               "tpm": f"Values look like TPM/FPKM (non-negative, up to {vmax:,.0f}, not on a log scale).",
               "log": f"Values look like log-scale expression (range {vmin:.2f} to {vmax:.2f})."}[vt])
    if vt == "counts":
        lib = np.nansum(vals, axis=0)
        if (lib <= 0).any():
            raise ImportProblem("Some samples have zero total counts.")
        M = pd.DataFrame(np.round(np.log2(vals / lib * 1e6 + 1), 9), index=M.index, columns=M.columns)
        rep.d("Converted raw counts to log2(counts per million + 1).")
        if lib.min() < 1e6 / 10:
            rep.w(f"Smallest library has only {lib.min():,.0f} counts in total; scores for that sample may be noisy.")
    elif vt == "tpm":
        M = pd.DataFrame(np.round(np.log2(vals + 1), 9), index=M.index, columns=M.columns)
        rep.d("Converted TPM/FPKM to log2(value + 1).")
    else:
        if vmin < 0:
            rep.w("Values include negatives, so they look already centred (z-scores or log-ratios). Kept as they are; pathway scores use only ranks within each sample.")
        rep.d("Kept values as they are (already on a log scale).")
    info["value_type"] = vt
    # ---- map genes to symbols
    labels = pd.Series(M.index.astype(str), index=M.index)
    if frac_ens > 0.5:
        info["gene_id_type"] = "Ensembl gene ids"
        e2s = ensembl_to_symbol()
        base = labels.str.replace(r"\.\d+$", "", regex=True).str.upper()
        sym = base.map(e2s)
        if alias is not None:
            sym = sym.fillna(alias.reindex(M.index).replace({"nan": np.nan}))
        n_un = int(sym.isna().sum())
        rep.f(f"Genes are Ensembl ids ({frac_ens:.0%}); {len(sym) - n_un} of {len(sym)} map to a gene symbol with the repo's protein-coding table.")
        if n_un:
            rep.w(f"{n_un} Ensembl id(s) have no symbol in the repo's protein-coding table (non-coding genes, retired ids or other species); they are dropped.")
        M = M.loc[sym.notna().values]
        M.index = sym[sym.notna()].values
    else:
        info["gene_id_type"] = "gene symbols"
        keep = labels.map(lambda x: not DATE_MANGLED.match(x) and x not in ("", "nan"))
        M = M.loc[keep.values]
        M.index = labels[keep].str.strip().values
        rep.f("Genes are gene symbols.")
    if M.index.duplicated().any():
        n = int(M.index.duplicated().sum())
        rep.w(f"{n} gene label(s) occur more than once (several ids or probes for one gene); I took the median of each.")
        M = M.groupby(level=0).median()
    # ---- missing values
    nan_frac = M.isna().mean(axis=1)
    bad = nan_frac > 0.2
    if bad.any():
        rep.w(f"Dropped {int(bad.sum())} gene(s) with more than 20% missing values.")
        M = M.loc[~bad]
    n_nan = int(M.isna().sum().sum())
    if n_nan:
        M = M.apply(lambda r: r.fillna(r.median()), axis=1)
        rep.w(f"Filled {n_nan} missing expression value(s) with that gene's median over samples.")
    # ---- how many genes the pipeline can use
    canon = to_canonical(M)
    info["n_genes_in_file"] = int(M.shape[0])
    info["n_genes_canonical"] = int(canon.shape[0])
    rep.f(f"{M.shape[0]} genes kept; {canon.shape[0]} of them are in the repo's gene-name system and are used for scoring.")
    if canon.shape[0] < 5000:
        rep.w(f"Only {canon.shape[0]} usable genes. Pathway scores are most comparable with patients when at least about 10,000 genes are present.")
    M.index.name = "gene"
    return M, info


# ---------------------------------------------------------------------------------------------- drug response
def _pick(cols, names):
    low = {str(c).strip().lower().replace(" ", "_"): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


def infer_direction(metric, hint=None):
    """('lower' | 'higher', guessed?) : which direction means more sensitive."""
    if hint:
        return hint, False
    if LOWER_WORDS.search(metric):
        return "lower", False
    if HIGHER_WORDS.search(metric):
        return "higher", False
    return "lower", True


def parse_response(df, rep, expr_samples, metric_hint=None, direction_hint=None, file_name=""):
    """Return (long DataFrame sample, drug, drug_original, response, metric, direction, ...; info dict)."""
    info = {}
    cols = list(df.columns)
    s_col, d_col = _pick(cols, SAMPLE_COLS), _pick(cols, DRUG_COLS)
    v_col = _pick(cols, VALUE_COLS)
    m_col = _pick(cols, METRIC_COLS)
    if s_col is not None and d_col is not None and s_col != d_col and (v_col is not None or m_col is not None or len(cols) >= 3):
        info["layout"] = "long"
        if v_col is None:
            numeric = [c for c in cols if c not in (s_col, d_col, m_col) and pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.8]
            if not numeric:
                raise ImportProblem("The drug-response file has sample and drug columns but no numeric value column.")
            v_col = numeric[0]
        rep.f(f"Drug-response table is long format (sample column '{s_col}', drug column '{d_col}', value column '{v_col}'"
              + (f", metric column '{m_col}')." if m_col is not None else ")."))
        others = [c for c in cols if c not in (s_col, d_col, v_col, m_col) and str(c).strip().lower() in VALUE_COLS]
        if others:
            rep.w(f"Other value columns ({', '.join(map(str, others))}) were ignored; I used '{v_col}'.")
        L = pd.DataFrame({"sample": df[s_col].astype(str).str.strip(), "drug_original": df[d_col].astype(str).str.strip(),
                          "value": pd.to_numeric(df[v_col], errors="coerce")})
        if m_col is not None:
            L["metric"] = df[m_col].astype(str).str.strip()
        else:
            vc = str(v_col).strip()
            L["metric"] = metric_hint or (vc.upper() if vc.lower() in VALUE_COLS and vc.lower() not in ("value", "response", "score", "metric_value")
                                          else (_metric_from_name(file_name) or "response"))
    else:
        info["layout"] = "wide"
        first = cols[0]
        body = df.drop(columns=[first])
        num = body.apply(pd.to_numeric, errors="coerce")
        keep = [c for c in num.columns if num[c].notna().mean() > 0.3]
        num = num[keep]
        num.index = df[first].astype(str).str.strip().values
        ek = {key(x) for x in expr_samples}
        over_rows = np.mean([key(x) in ek for x in num.index])
        over_cols = np.mean([key(x) in ek for x in num.columns])
        if over_cols > over_rows:
            num = num.T
            rep.f(f"Drug-response table is wide with drugs in rows and samples in columns (flipped by me; {over_cols:.0%} of the column names match expression samples).")
        else:
            rep.f(f"Drug-response table is wide with samples in rows and drugs in columns ({over_rows:.0%} of the row names match expression samples).")
        num.columns = [str(c).strip() for c in num.columns]
        L = num.stack().rename("value").reset_index()
        L.columns = ["sample", "drug_original", "value"]
        L["metric"] = metric_hint or _metric_from_name(file_name) or "response"
        if L["metric"].iloc[0] == "response":
            rep.w("A wide drug table does not say what the numbers are (AUC, IC50, viability, ...). Please give --metric, or put it in the file name (e.g. drug_response_AUC.tsv). Called it 'response'.")
    n0 = len(L)
    L = L.dropna(subset=["value"])
    L = L[~L["drug_original"].isin(["", "nan", "NaN"])]
    if len(L) < n0:
        rep.w(f"Dropped {n0 - len(L)} response row(s) with no number.")
    if L.empty:
        raise ImportProblem("The drug-response file has no usable numbers.")
    L["metric"] = L["metric"].map(lambda m: m.upper() if len(m) <= 6 else m)
    # ---- more than one metric: keep one so values are never mixed
    mets = L["metric"].value_counts()
    if len(mets) > 1:
        want = (metric_hint or "").upper()
        primary = want if want in mets.index else next((m for m in PRIMARY_ORDER if m in mets.index), mets.index[0])
        rep.w(f"The drug table has several metrics ({', '.join(f'{m}: {n} rows' for m, n in mets.items())}). I kept only '{primary}' in response.tsv "
              "so that different scales are never mixed; the others are in response_other_metrics.tsv.")
        info["other_metrics"] = L[L["metric"] != primary]
        L = L[L["metric"] == primary]
    metric = L["metric"].iloc[0]
    info["metric"] = metric
    # ---- IC50-like values are analysed on the ln scale (as in the paper)
    if IC_LIKE.search(metric) and not re.search(r"ln|log", metric, re.I):
        if (L["value"] > 0).all():
            L["value_raw"] = L["value"]
            L["value"] = np.log(L["value"])
            metric = "ln_" + metric
            L["metric"] = metric
            info["metric"] = metric
            rep.d(f"Took the natural log of the {L['metric'].iloc[0][3:]} values (ln scale, as in the paper).")
        else:
            rep.w("IC50-like values include zero or negatives, so I could not take the log. Used them as they are.")
    direction, guessed = infer_direction(metric, direction_hint)
    info["direction"] = direction
    if guessed:
        rep.w(f"I could not tell from the name '{metric}' whether a LOWER or a HIGHER number means more sensitive. I assumed lower = more sensitive "
              "(true for AUC, IC50, viability). If that is wrong, run again with --direction higher.")
    else:
        rep.f(f"Metric: {metric}; {direction} value = more sensitive" + (" (as you said)." if direction_hint else "."))
    # ---- drug names
    L["drug"] = L["drug_original"].map(common_drug_name)
    nd, nd0 = L["drug"].nunique(), L["drug_original"].nunique()
    if nd < nd0:
        rep.d(f"Merged {nd0} drug names into {nd} after matching synonyms with the repo's drug-name map.")
    # ---- repeats of the same sample x drug
    dup = L.duplicated(["sample", "drug"], keep=False)
    if dup.any():
        rep.w(f"{int(L.duplicated(['sample', 'drug']).sum())} repeated sample x drug measurement(s) (replicates); the pipeline takes their median.")
    return L, info


def _metric_from_name(name):
    for m in ("AUC", "IC50", "GI50", "EC50", "GR50", "viability", "inhibition"):
        if re.search(m, name, re.I):
            return m.upper() if m.upper() in PRIMARY_ORDER or m == "GR50" else m
    return None


# ---------------------------------------------------------------------------------------------- samples table
def parse_samples(df, rep):
    cols = {str(c).strip().lower().replace(" ", "_"): c for c in df.columns}
    s = cols.get("sample") or cols.get("sample_id") or df.columns[0]
    out = pd.DataFrame({"sample": df[s].astype(str).str.strip()})
    for want in ("patient", "cancer_type", "organoid_line", "replicate"):
        if want in cols and cols[want] != s:
            out[want] = df[cols[want]].values
    for c in df.columns:
        k = str(c).strip().lower().replace(" ", "_")
        if k not in out.columns and c != s:
            out[k] = df[c].values
    rep.f(f"Samples table: {len(out)} rows with columns {', '.join(out.columns)}.")
    return out.drop_duplicates("sample").set_index("sample")


# ---------------------------------------------------------------------------------------------- reconciliation
def reconcile(expr_samples, other, label, rep):
    """Map names in `other` (response / samples file) onto expression names, tolerating case, spaces, dashes."""
    ek = {}
    for s in expr_samples:
        ek.setdefault(key(s), []).append(s)
    mapping, fuzzy, missing = {}, [], []
    for s in dict.fromkeys(other):
        if s in expr_samples:
            mapping[s] = s
        elif key(s) in ek and len(ek[key(s)]) == 1:
            mapping[s] = ek[key(s)][0]
            fuzzy.append(s)
        else:
            missing.append(s)
    if fuzzy:
        rep.d(f"Matched {len(fuzzy)} {label} sample name(s) to expression names ignoring upper/lower case, spaces and dashes (e.g. '{fuzzy[0]}' -> '{mapping[fuzzy[0]]}').")
    if missing:
        rep.w(f"{len(missing)} sample(s) in the {label} file have no expression data and are left out: {', '.join(map(str, missing[:8]))}"
              + (" ..." if len(missing) > 8 else "") + ".")
    return mapping, missing


# --------------------------------------------------------------------------------------------- scoring
def score_pathways(expr):
    """Samples x pathways Reactome rank scores; same as robust/build_training_table.score()."""
    e = to_canonical(expr)
    gs = {k: [g for g in v if g in e.index] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    return rank_score(e, gs).T


# ---------------------------------------------------------------------------------------------- main entry
def import_folder(input_folder, set_name, out_dir=None, cancer=None, metric=None, direction=None, value_type=None, write=True):
    """Run the whole import. Returns a dict with the tables and the report text; writes files unless write=False."""
    input_folder = Path(input_folder)
    out = Path(out_dir) if out_dir else FOLDER / "imports" / set_name
    rep = Report()
    status = "ok"
    try:
        res = _import(input_folder, set_name, cancer, metric, direction, value_type, rep)
    except ImportProblem as e:
        rep.p(str(e))
        status = "stopped"
        res = {}
    text = report_text(set_name, input_folder, rep, res, status)
    res.update(report=text, status=status, set_name=set_name, out=out)
    if write:
        out.mkdir(parents=True, exist_ok=True)
        (out / "IMPORT_REPORT.md").write_text(text)
        if status == "ok":
            write_set(res, out)
    return res


def _import(folder, set_name, cancer, metric_hint, direction_hint, value_type, rep):
    files = find_files(folder, rep)
    expr, einfo = parse_expression(read_table(files["expression"]), rep, value_type)
    names = list(expr.columns)
    rep.f(f"{len(names)} samples in the expression file: {', '.join(names[:8])}" + (" ..." if len(names) > 8 else "") + ".")
    L, rinfo = parse_response(read_table(files["response"]), rep, names, metric_hint, direction_hint, files["response"].name)
    samples = parse_samples(read_table(files["samples"]), rep) if "samples" in files else None
    # --- sample names
    mp, miss = reconcile(names, L["sample"], "drug-response", rep)
    L["sample"] = L["sample"].map(mp)
    n_before = len(L)
    L = L.dropna(subset=["sample"])
    if L.empty:
        rep.p("No sample name in the drug-response file matches a sample in the expression file. Example expression names: "
              f"{', '.join(names[:4])}; example response names: {', '.join(map(str, pd.unique(read_table(files['response']).iloc[:, 0])[:4]))}.")
        raise ImportProblem("Sample names in the two files do not match, so no organoid has both expression and drug response.")
    if samples is not None:
        mps, _ = reconcile(names, samples.index, "samples-table", rep)
        samples.index = [mps.get(s, s) for s in samples.index]
        samples = samples[samples.index.isin(names)]
    no_resp = [s for s in names if s not in set(L["sample"])]
    if no_resp:
        rep.w(f"{len(no_resp)} expression sample(s) have no drug response: {', '.join(no_resp[:8])}" + (" ..." if len(no_resp) > 8 else "")
              + ". They stay in the expression file and are scored, but cannot be used for training.")
    # --- replicates of the same organoid line collapse to one column (only when the samples table says so)
    if samples is not None and "organoid_line" in samples:
        line = samples["organoid_line"].astype(str)
        if line.nunique() < len(line):
            m = {s: (line.get(s, s)) for s in names}
            expr = expr.T.groupby(pd.Series(m)).median().T
            L["sample"] = L["sample"].map(m)
            rep.d(f"Collapsed {len(names)} expression columns into {expr.shape[1]} organoid lines (median of replicates) using the 'organoid_line' column.")
            samples = samples.reset_index().drop_duplicates("organoid_line").assign(sample=lambda d: d["organoid_line"].astype(str)).set_index("sample")
            names = list(expr.columns)
    # --- cancer type
    if cancer is None and samples is not None and "cancer_type" in samples and samples["cancer_type"].notna().any():
        cancer = str(samples["cancer_type"].mode().iloc[0])
        rep.f(f"Cancer type taken from the samples table: {cancer}.")
    elif cancer:
        rep.f(f"Cancer type given by you: {cancer}.")
    else:
        rep.w("No cancer type was given. The data are imported, but the patient analysis needs --cancer (for example colorectal, bladder, pancreas, liver).")
    # --- scoring
    scores = score_pathways(expr)
    npw = scores.shape[1]
    rep.f(f"Scored {npw} Reactome pathways in {scores.shape[0]} samples (rank score).")
    if npw < N_REACTOME_FULL:
        rep.w(f"Only {npw} of the {N_REACTOME_FULL} Reactome pathways the repo uses could be scored (some pathways have fewer than 5 genes in your file). "
              "Scores for the missing pathways are not available.")
    if len(set(L["sample"])) < 8:
        rep.w(f"Only {len(set(L['sample']))} organoids have drug response. Models trained on so few organoids are very unstable; treat results as exploratory.")
    drugs = sorted(L["drug"].unique())
    # --- which drugs the paper's network table covers
    from obd import network as N
    try:
        prox = set(N.precomputed_proximity().columns)
    except Exception:
        prox = set()
    covered = [d for d in drugs if d in prox or d == "OXALIPLATIN"]
    rep.f(f"{len(drugs)} drugs. {len(covered)} of them have a network-proximity score in the paper's table and can go through the paper method"
          + (f": {', '.join(covered[:12])}" if covered else "") + ".")
    files_info = {k: dict(name=p.name, size_bytes=p.stat().st_size, sha256=sha256(p)) for k, p in files.items()}
    return dict(expr=expr, scores=scores, response=L, samples=samples, cancer=cancer, einfo=einfo, rinfo=rinfo, drugs=drugs,
                drugs_with_network_score=covered, files=files_info, n_pathways=npw)


def write_set(res, out):
    expr, L, samples = res["expr"], res["response"], res["samples"]
    expr.round(6).to_csv(out / "expression.tsv.gz", sep="\t")
    r = L.copy()
    r["response"] = r["value"]
    r["direction"] = "lower = more sensitive" if res["rinfo"]["direction"] == "lower" else "higher = more sensitive"
    cols = ["sample", "drug", "response", "metric", "direction", "drug_original"] + (["value_raw"] if "value_raw" in r else [])
    r[cols].to_csv(out / "response.tsv", sep="\t", index=False)
    if "other_metrics" in res["rinfo"]:
        res["rinfo"]["other_metrics"].to_csv(out / "response_other_metrics.tsv", sep="\t", index=False)
    s = samples.copy() if samples is not None else pd.DataFrame(index=pd.Index(expr.columns, name="sample"))
    s = s.reindex(expr.columns)
    s.index.name = "sample"
    if "organoid_line" not in s:
        s["organoid_line"] = s.index
    s.to_csv(out / "samples.tsv", sep="\t")
    pd.DataFrame({"geo_accession": "", "title": list(expr.columns), "sample": list(expr.columns)}).to_csv(out / "gsm_map.tsv", sep="\t", index=False)
    res["scores"].to_csv(out / "pathway_scores.tsv.gz", sep="\t")
    meta = dict(set=res["set_name"], cancer=res["cancer"], created=datetime.datetime.now().isoformat(timespec="seconds"),
                n_samples=int(expr.shape[1]), n_genes=int(expr.shape[0]), n_pathways_scored=int(res["n_pathways"]),
                n_samples_with_response=int(L["sample"].nunique()), n_drugs=len(res["drugs"]), drugs=res["drugs"],
                drugs_with_network_score=res["drugs_with_network_score"], metric=res["rinfo"]["metric"], direction=res["rinfo"]["direction"],
                expression_layout=res["einfo"]["orientation"], gene_id_type=res["einfo"]["gene_id_type"], value_type=res["einfo"]["value_type"],
                response_layout=res["rinfo"]["layout"], input_files=res["files"], scoring="rank_score (singscore-style), Reactome sets >= 5 genes, as robust/build_training_table.score")
    (out / "organoid_set.json").write_text(json.dumps(meta, indent=2))
    (out / "SOURCE.md").write_text(f"# {res['set_name']}\n\nUser-supplied organoid set imported by nar/21_user_organoids/code/load_user_organoids.py.\n"
                                   f"Input files (sha256 in organoid_set.json): {', '.join(v['name'] for v in res['files'].values())}.\n"
                                   "See IMPORT_REPORT.md for everything that was detected and changed.\n")


def report_text(set_name, folder, rep, res, status):
    L = [f"# Import report: {set_name}", "", f"Folder read: `{folder}`", ""]
    if status == "ok":
        n = res["expr"].shape[1]
        L += [f"**Result: imported.** {n} organoid samples, {len(res['drugs'])} drugs, {res['n_pathways']} pathways scored.", ""]
    else:
        L += ["**Result: stopped.** The data could not be imported. See 'Problems' below.", ""]
    L += ["## Problems (must be fixed or checked)", ""] + ([f"- {x}" for x in rep.problem] or ["- None."])
    L += ["", "## Warnings (imported, but please check)", ""] + ([f"- {x}" for x in rep.warn] or ["- None."])
    L += ["", "## What I found", ""] + ([f"- {x}" for x in rep.found] or ["- Nothing."])
    L += ["", "## What I changed in your data", ""] + ([f"- {x}" for x in rep.did] or ["- Nothing was changed."])
    if status == "ok":
        L += ["", "## Files written", "", "- `expression.tsv.gz`: genes x samples, gene symbols, log scale.",
              "- `response.tsv`: one row per sample x drug (lower or higher = more sensitive is stated in the `direction` column).",
              "- `samples.tsv`, `gsm_map.tsv`: sample information.", f"- `pathway_scores.tsv.gz`: samples x {res['n_pathways']} Reactome pathway scores.",
              "- `organoid_set.json`: all settings and the sha256 of your input files."]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit("Run import_organoids.py instead.")
