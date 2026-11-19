"""Turn the user's uploaded files into the folder layout the loader (nar/21_user_organoids) reads.
Only reshaping: no value is changed.  Output: nar/22_user_data_results/_data/inputs/<set>/ (small files).
Run:  python3 -I prepare_inputs.py
"""
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
UP = REPO / "data" / "user_upload"
OUT = Path(__file__).resolve().parents[1]
INP = OUT / "_data" / "inputs"
TAB = OUT / "tables"
for p in (INP, TAB):
    p.mkdir(parents=True, exist_ok=True)
log = []

# ------------------------------------------------------------------ LUAD (GSE165611)
L = UP / "organoid_LUAD"
em = L / "expression/GSE165611/GSE165611_series_matrix.txt"
med = pd.read_csv(em / "geneID_expression_median.txt", sep="\t")
log.append(("LUAD", "geneID_expression_median.txt rows", len(med)))
log.append(("LUAD", "duplicate gene symbols in median file", int(med["geneID"].duplicated().sum())))
med = med.drop(columns=["uniprotID"]).set_index("geneID")
med.columns = [c.upper() for c in med.columns]
# cross-check against the GEO series matrix (Ensembl_symbol ids) -> same numbers?
sm = pd.read_csv(em / "GSE165611_series_matrix(LUAD).txt", sep="\t", comment="!", quotechar='"', index_col=0)
sym = sm.index.str.replace(r"^ENSG\d+\.\d+_", "", regex=True)
titles = None
for line in open(em / "GSE165611_series_matrix(LUAD).txt"):
    if line.startswith("!Sample_title"):
        titles = [t.strip('"').upper() for t in line.rstrip("\n").split("\t")[1:]]
        break
sm.columns = titles
sm.index = sym
sm = sm[~sm.index.duplicated(keep=False)]
common_g = med.index.intersection(sm.index)
a = med.loc[common_g, titles]
b = sm.loc[common_g, titles]
d = (a.values - b.values)
log.append(("LUAD", "genes compared median-file vs GEO series matrix", len(common_g)))
log.append(("LUAD", "max abs difference of values", float(np.nanmax(np.abs(d)))))
log.append(("LUAD", "share of values identical", float((np.abs(d) < 1e-9).mean())))
log.append(("LUAD", "series-matrix probes", sm.shape[0]))
d0 = INP / "luad_gse165611"
d0.mkdir(parents=True, exist_ok=True)
med.to_csv(d0 / "expression.tsv.gz", sep="\t")

s2 = pd.read_csv(L / "drug_response/supp_table_S2b.txt", sep="\t")
mr = pd.read_csv(L / "drug_response/median_response.txt", sep="\t")
log.append(("LUAD", "median_response.txt equals supp_table_S2b IC50 column (row by row)", bool((mr["IC50"].values == s2["IC50"].values).all())))
log.append(("LUAD", "rows per organoid x drug (repeats)", int(s2.groupby(["Organoid", "Drug Name"]).size().iloc[0])))
cl = s2.drop_duplicates(["Organoid"]).set_index("Organoid")["Cluster"]
# one value per cluster x drug (check)
cv = s2.groupby(["Cluster", "Drug Name"])["IC50"].agg(["nunique", "first", "size"])
log.append(("LUAD", "cluster x drug cells with more than one distinct value", int((cv["nunique"] > 1).sum())))
log.append(("LUAD", "cluster x drug cells listed", len(cv)))
cvt = cv["first"].unstack(0)
cvt.to_csv(TAB / "luad_cluster_by_drug_value_listed.tsv", sep="\t")
cl.to_frame("cluster").to_csv(TAB / "luad_cluster_labels.tsv", sep="\t")
drugs = sorted(s2["Drug Name"].unique())
rows = []
for o, c in cl.items():
    for dr in drugs:
        v = cvt.loc[dr, c] if c in cvt.columns else np.nan
        rows.append((o, dr, v))
long = pd.DataFrame(rows, columns=["sample", "drug", "value"])
long_listed = long.dropna()
long_zero = long.fillna({"value": 0.0})        # user's own Excel matrix writes 0 where the cluster has no entry
long_zero.assign(metric="log2FC").to_csv(d0 / "drug_response.tsv", sep="\t", index=False)
d1 = INP / "luad_gse165611_listed"
d1.mkdir(exist_ok=True)
long_listed.assign(metric="log2FC").to_csv(d1 / "drug_response.tsv", sep="\t", index=False)
med.to_csv(d1 / "expression.tsv.gz", sep="\t")
# does the user's xlsx matrix match? (zero fill)
x = pd.read_excel(L / "drug_response/luad drug.xlsx", sheet_name="Sheet1", header=0).iloc[:, :5].dropna(subset=["Drug"])
x = x.set_index("Drug")[["C1", "C2", "C3", "C4"]]
cmp = cvt.reindex(x.index).fillna(0)[["C1", "C2", "C3", "C4"]]
log.append(("LUAD", "luad drug.xlsx Sheet1 matrix equals zero-filled cluster x drug table", bool(np.allclose(x.values.astype(float), cmp.values))))
samp = pd.DataFrame({"sample": cl.index, "patient": cl.index, "cancer_type": "lung", "cluster": cl.values})
samp.to_csv(d0 / "samples.tsv", sep="\t", index=False)
samp.to_csv(d1 / "samples.tsv", sep="\t", index=False)
# histology from GEO matrix
hist = {}
for line in open(em / "GSE165611_series_matrix(LUAD).txt"):
    if line.startswith("!Sample_characteristics_ch1") and "histopathology" in line:
        hist = dict(zip(titles, [t.strip('"').replace("histopathology: ", "") for t in line.rstrip("\n").split("\t")[1:]]))
pd.Series(hist, name="histology").to_frame().to_csv(TAB / "luad_histology.tsv", sep="\t")
log.append(("LUAD", "organoids in expression / in response / both", f"{len(titles)} / {cl.shape[0]} / {len(set(titles) & set(o.upper() for o in cl.index))}"))

# ------------------------------------------------------------------ OVARY (GSE229119)
O = UP / "organoid_OVARY"
oe = pd.read_csv(O / "expression/GSE229119/GSE229119_series_matrix.txt/GSE229119_log2TPM.gene.expr.table.txt", sep="\t", index_col=0)
log.append(("OVARY", "log2TPM table genes x samples", f"{oe.shape[0]} x {oe.shape[1]}"))
log.append(("OVARY", "minimum value (floor for TPM = 0)", float(oe.values.min())))
log.append(("OVARY", "share of values at the floor", float((oe.values == oe.values.min()).mean())))
base = {c: c.replace("_NC_D3", "") for c in oe.columns if c.endswith("_NC_D3")}
name_map = {"Hey": "Hey", "MCF7": "MCF7", "ORG2414": "Organoid 2414", "ORG2445": "Organoid 2445", "ORG3008": "Organoid 3008", "OVCA432": "OVCA432", "SKOV3": "SKOV3"}
base_expr = oe[list(base)].rename(columns={k: name_map[v] for k, v in base.items()})
od = INP / "ovary_gse229119"
od.mkdir(parents=True, exist_ok=True)
base_expr.to_csv(od / "expression.tsv.gz", sep="\t")
ores = pd.read_csv(O / "drug_response/median_response.txt", sep="\t")
ores["line"] = ores["Organoid"].str.split(",").str[0].str.strip()
g = ores.groupby(["line", "Drug Name"])["IC50"].agg(["nunique", "first", "size"])
log.append(("OVARY", "lines with more than one distinct IC50 per drug (should be 0)", int((g["nunique"] > 1).sum())))
log.append(("OVARY", "rows collapsed into line x drug", f"{len(ores)} -> {len(g)}"))
one = g["first"].reset_index().rename(columns={"line": "sample", "Drug Name": "drug", "first": "value"})
one["metric"] = "IC50"
one.to_csv(od / "drug_response.tsv", sep="\t", index=False)
one.pivot(index="sample", columns="drug", values="value").to_csv(TAB / "ovary_ic50_per_line.tsv", sep="\t")
ls = pd.DataFrame({"sample": list(name_map.values()), "patient": list(name_map.values()), "cancer_type": "ovarian",
                   "line_type": ["cell line", "cell line (breast)", "organoid", "organoid", "organoid", "cell line", "cell line"]})
ls.to_csv(od / "samples.tsv", sep="\t", index=False)
# ovarian-only version (drop MCF7, a breast line)
od2 = INP / "ovary_gse229119_noMCF7"
od2.mkdir(exist_ok=True)
base_expr.drop(columns=["MCF7"]).to_csv(od2 / "expression.tsv.gz", sep="\t")
one[one["sample"] != "MCF7"].to_csv(od2 / "drug_response.tsv", sep="\t", index=False)
ls[ls["sample"] != "MCF7"].to_csv(od2 / "samples.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ LIHC (GSE188541)
H = UP / "organoid_LIHC"
hr = pd.read_csv(H / "drug_response/median_response.txt", sep="\t")
hr.to_csv(TAB / "lihc_ic50_long.tsv", sep="\t", index=False)
hr.pivot(index="Drug Name", columns="Organoid", values="median_IC50").to_csv(TAB / "lihc_ic50_per_organoid.tsv", sep="\t")
txt = open(H / "expression/GSE188541/GSE188541_series_matrix.txt/GSE188541_series_matrix.txt").read()
n_rows = txt.split("!series_matrix_table_begin")[1].split("!series_matrix_table_end")[0].strip().count("\n")
log.append(("LIHC", "data rows in the series matrix (only the header line = 0)", n_rows))
pd.DataFrame(log, columns=["set", "check", "value"]).to_csv(TAB / "input_checks.tsv", sep="\t", index=False)
print(pd.DataFrame(log, columns=["set", "check", "value"]).to_string())
