"""Isogenic 5-FU-sensitive / -resistant colorectal cells (HCT-8 vs HCT-8/5-FU; GEO GSE81008 -> GSE81005, PrimeView).

Paper (Fig. 4b): sensitive = HCT8 0 h (n = 3), resistant = HCT8/5Fu 0 h, 24 h, 48 h (n = 9); ssGSEA NES; two-sided
two-sample Student's t-test. Two versions:
  authors : the authors' own NES (Source Data, MOESM3 sheet 'Figure 4b');
  current : recomputed here from the GEO series matrix with the same ssGSEA wrapper as everything else.
"""
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from common import DATA, RAW, REPO, C, gene_sets, scoring

MOESM3 = REPO / "nar/_cache/41467_2020_19313_MOESM3_ESM.xlsx"
GEO = RAW / "GSE81008"


def _labels(cols):
    sens = [c for c in cols if c.startswith("HCT8-0h")]
    res = [c for c in cols if c.startswith("HCT8/5Fu-") and c.split("-")[1] in ("0h", "24h", "48h")]
    return sens, res


def authors_nes():
    d = pd.read_excel(MOESM3, sheet_name="Figure 4b", header=None)
    hdr = d.iloc[4].tolist()
    nes = d.iloc[5:].copy()
    nes.columns = hdr
    nes = nes.set_index(hdr[0]).astype(float)
    nes.index.name = "pathway"
    return nes


def current_nes():
    f = DATA / "scores" / "isogenic_GSE81005_ssgsea.tsv.gz"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0)
    import gzip
    titles, rows = None, []
    with gzip.open(GEO / "GSE81005_series_matrix.txt.gz", "rt") as fh:
        inside = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"\n') for t in line.split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if inside:
                rows.append(line.rstrip("\n").split("\t"))
    m = pd.DataFrame(rows[1:], columns=[c.strip('"') for c in rows[0]]).set_index("ID_REF")
    m.index = m.index.str.strip('"')
    m = m.apply(pd.to_numeric, errors="coerce")
    m.columns = titles
    ann = pd.read_csv(GEO / "GPL15207_id_symbol.tsv.gz", sep="\t").dropna()
    ann = ann[~ann["Gene Symbol"].str.contains("///")]
    sym = ann.set_index("ID")["Gene Symbol"]
    m = m.loc[m.index.intersection(sym.index)]
    m["gene"] = sym.loc[m.index].values
    m["mean"] = m.drop(columns="gene").mean(1)
    m = m.sort_values("mean", ascending=False).drop_duplicates("gene").set_index("gene").drop(columns="mean")
    if m.max().max() > 100:
        m = np.log2(m + 1)
    e = C.to_canonical(m)
    gs = {k: [g for g in v if g in set(e.index)] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    nes = scoring.ssgsea(e, gs, threads=2)
    f.parent.mkdir(parents=True, exist_ok=True)
    nes.to_csv(f, sep="\t")
    return nes


def iso_table(nes):
    sens, res = _labels(nes.columns)
    t, p = ttest_ind(nes[sens].values, nes[res].values, axis=1)
    return pd.DataFrame({"t": t, "p": p, "mean_sensitive": nes[sens].mean(1).values,
                         "mean_resistant": nes[res].mean(1).values}, index=nes.index), sens, res


if __name__ == "__main__":
    for name, fn in (("authors", authors_nes), ("current", current_nes)):
        nes = fn()
        t, s, r = iso_table(nes)
        print(name, nes.shape, len(s), len(r), t.loc["REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS"].to_dict())
