"""RESULTS_INDEX.md: every figure (PNG) and table (TSV/JSON) under nar/ (except _data, _cache, code), with a one-line
description (figure: first heading of its caption; table: its columns), grouped by paper section and by cancer."""
from common import *
import re
skip = ("/_data/", "/_cache/", "/code/", "/_test/", "/cache/", "/data/")
files = sorted(p for ext in ("png", "tsv", "json", "tsv.gz") for p in NAR.rglob(f"*.{ext}")
               if not any(s in str(p) for s in skip) and p.name != "completeness.tsv" or p.name == "completeness.tsv")
files = [p for p in files if not any(s in str(p) for s in skip)]
def desc(p):
    if p.suffix == ".png":
        c = p.with_suffix("").with_suffix(".caption.md") if False else p.parent / (p.stem + ".caption.md")
        if c.exists():
            for l in c.read_text().splitlines():
                if l.strip():
                    return l.lstrip("# ").strip()[:150]
        return "figure (no caption file)"
    try:
        if p.suffix == ".json":
            return "JSON summary: " + ", ".join(list(json.load(open(p)).keys())[:6])[:140]
        cols = pd.read_csv(p, sep="\t", nrows=0).columns
        return "table: " + ", ".join(map(str, cols[:7]))[:140]
    except Exception:
        return "table"
SEC = [("Fig. 1 (method, schematics, data map)", r"^[Ff]ig1|workflow|data_map|S1_data|preregistration"),
       ("Fig. 2 (network pathway biomarker, survival)", r"^Fig2|_full|forest|composite|KM|cross_set|rank1|verify|summary|zcutoff|cutoff|main_|interaction|stage|km_"),
       ("Fig. 3 / Supp 6-9 (comparison feature sets, centrality, deep learning)", r"^Fig3|SuppFig[6-9]|comparison|neighbour|target|bench|central"),
       ("Fig. 4 (isogenic cells, bootstrap)", r"^Fig4|isogenic|bootstrap|S7_|permutation|random"),
       ("Fig. 5 / Supp 12, 13, 18 (mutation concordance)", r"^Fig5|ERCC2|BRAF|(?<!per)mutation|SuppFig1[238]"),
       ("Supp. Figs 1, 14 (batch / PCA)", r"SuppFig1_|PCA|alignment|SuppFig14"),
       ("Supp. Figs 2-5, 10, 11, 15-17 (models, multi-pathway, size, cell lines, proximity, inner validation)",
        r"SuppFig[2-5]|SuppFig1[015-7]|coefficients|multipathway|size_vs|GDSC|proximity|loo|split|learning"),
       ("Supp. Tables 1-4 (clinical balance, isogenic, bootstrap)", r"SuppTable|clinical"),
       ("DepMap / registered test / other", r".")]
CAN = [("PDX (Gao) / breast / melanoma", r"PDX|BRCA_|SKCM|_CM_"), ("colon / rectum", r"COAD|COLON|colorectal|CRC|5FU|FLUORO|vandewetering|isella|HCT8"), ("bladder", r"BLCA|bladder"),
       ("pancreas", r"PAAD|pancrea|PDAC|tiriac|shi"), ("liver", r"LIHC|liver|broutier|_ji|SORAF"), ("ovarian", r"OV_|ovar|vias|GSE229119"),
       ("lung", r"LUAD|lung|NSCLC|GSE165611"), ("sarcoma", r"SARC|sarcoma"),
       ("all sets / cross-cancer", r".")]
rows = []
for p in files:
    rel = str(p.relative_to(NAR)); low = rel
    order = [SEC[i] for i in (7, 5, 4, 2, 3, 6, 0, 1, 8)]
    sec = next(s for s, rx in order if re.search(rx, p.name))
    can = next(s for s, rx in CAN if re.search(rx, rel, re.I))
    rows.append((sec, can, rel, desc(p)))
T = pd.DataFrame(rows, columns=["section", "cancer", "path", "description"])
T.to_csv(TAB / "results_index.tsv", sep="\t", index=False)
out = ["# Results index: every figure and table under nar/", "",
       f"{(T.path.str.endswith('.png')).sum()} figures (each also has a PDF and a .caption.md next to it) and "
       f"{(~T.path.str.endswith('.png')).sum()} tables. Paths are relative to `nar/`. Descriptions come from the caption's first line "
       "(figures) or the column names (tables). Built by `10_completeness/code/c09_index.py`; the same list is in "
       "`10_completeness/tables/results_index.tsv`. Items being made by the B3 agent (`07_B3_models/`) will appear after a re-run.", "",
       "## Part A. By paper section", ""]
for s, _ in SEC:
    sub = T[T.section == s]
    if len(sub) == 0:
        continue
    out += [f"### {s}", ""]
    for c, _ in CAN:
        ss = sub[sub.cancer == c]
        if len(ss):
            out.append(f"**{c}**"); out.append("")
            out += [f"- `{r.path}` — {r.description}" for r in ss.itertuples()]
            out.append("")
out += ["## Part B. By cancer (counts; full list in Part A)", "", "| Cancer | Figures | Tables | Folders |", "|---|---|---|---|"]
for c, _ in CAN:
    ss = T[T.cancer == c]
    out.append(f"| {c} | {ss.path.str.endswith('.png').sum()} | {(~ss.path.str.endswith('.png')).sum()} | "
               + ", ".join(sorted({'/'.join(x.split('/')[:2]) for x in ss.path}))[:300] + " |")
(OUT / "RESULTS_INDEX.md").write_text("\n".join(out) + "\n")
print(T.groupby(["section"]).size().to_string()); print(T.groupby("cancer").size().to_string())
