"""Read GEO series-matrix files: sample annotation table and expression matrix."""
import gzip, io
import pandas as pd

def read_meta(path):
    """Sample table: one row per GSM; 'key: value' characteristics become columns."""
    rows = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"): break
            if not line.startswith("!Sample_"): continue
            parts = line.rstrip("\n").split("\t")
            rows.setdefault(parts[0][8:], []).append([p.strip('"') for p in parts[1:]])
    gsm = rows["geo_accession"][0]
    out = pd.DataFrame(index=gsm)
    for k in ("title", "source_name_ch1", "platform_id", "description"):
        if k in rows: out[k] = rows[k][0]
    for vals in rows.get("characteristics_ch1", []):
        for g, v in zip(gsm, vals):
            if ":" in v:
                key, val = v.split(":", 1)
                out.loc[g, key.strip().lower()] = val.strip()
    return out

def read_expr(path):
    """Expression table (probes x samples); empty frame when the matrix holds no values."""
    with gzip.open(path, "rt", errors="replace") as f:
        lines, on = [], False
        for line in f:
            if line.startswith("!series_matrix_table_begin"): on = True; continue
            if line.startswith("!series_matrix_table_end"): break
            if on: lines.append(line)
    if len(lines) < 2: return pd.DataFrame()
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0, low_memory=False)
