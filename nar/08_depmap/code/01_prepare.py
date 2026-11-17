"""Step 1: parse the downloaded (untrusted) GDSC / Cell Model Passports / Project Score files
into clean pickles under nar/_cache/depmap/derived/ and record file versions.

python -I nar/08_depmap/code/01_prepare.py
"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # -I drops the script dir; add our own code folder only

import numpy as np
import pandas as pd

from common import CACHE, CUR, DERIVED, TAB

from obd import scoring  # noqa: E402  (robust/ package, added to sys.path by common)
from obd.cohorts import to_canonical  # noqa: E402
from obd.reference import canonical_symbol, common_drug_name, gene_sets  # noqa: E402


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def versions():
    rows = []
    for p in sorted(CACHE.rglob("*")):
        if p.is_file() and "derived" not in p.parts:
            rows.append({"file": str(p.relative_to(CACHE)), "bytes": p.stat().st_size, "sha256": sha(p)})
    for p in ("ctrpv2_ccle/response.tsv", "ctrpv2_ccle/expression.tsv.gz", "ctrpv2_ccle/samples.tsv",
              "prism_repurposing/response.tsv", "prism_repurposing/expression.tsv.gz", "prism_repurposing/samples.tsv",
              "bladder_lee2018/response.tsv", "bladder_lee2018/expression.tsv.gz"):
        f = CUR / p
        rows.append({"file": f"data/curated/preclinical/{p}", "bytes": f.stat().st_size, "sha256": sha(f)})
    pd.DataFrame(rows).to_csv(TAB / "data_versions.tsv", sep="\t", index=False)


def gdsc():
    parts = []
    for ds in ("gdsc1", "gdsc2"):
        f = next((CACHE / ds).glob("*.xlsx"))
        d = pd.read_excel(f, engine="openpyxl")
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    d.columns = [c.lower() for c in d.columns]
    alias = {"OBATOCLAX MESYLATE": "OBATOCLAX", "ABT737": "ABT-737", "VENOTOCLAX": "VENETOCLAX"}
    d["drug"] = d["drug_name"].map(lambda s: common_drug_name(alias.get(str(s).upper(), str(s).upper())))
    keep = ["dataset", "cosmic_id", "cell_line_name", "sanger_model_id", "tcga_desc", "drug_id", "drug_name", "drug",
            "putative_target", "pathway_name", "min_conc", "max_conc", "ln_ic50", "auc", "rmse", "z_score"]
    d = d[keep]
    d.to_pickle(DERIVED / "gdsc_long.pkl")
    print("GDSC rows", len(d), d.groupby("dataset").size().to_dict())


def models():
    m = pd.read_csv(next((CACHE / "cmp_models").glob("*.csv.gz")), low_memory=False)
    m = m.set_index("model_id")
    m[["model_name", "model_type", "msi_status", "tissue", "cancer_type", "cancer_type_detail", "COSMIC_ID",
       "BROAD_ID", "mutation_data", "rnaseq_data", "crispr_ko_data", "drug_data"]].to_pickle(DERIVED / "models.pkl")


def expression():
    f = CACHE / "cmp_rnaseq_x" / "rnaseq_tpm_20220624.csv"
    head = pd.read_csv(f, nrows=4, header=None, low_memory=False)
    model_ids = head.iloc[0, 2:].tolist()
    x = pd.read_csv(f, skiprows=5, header=None, low_memory=False)
    sym = x.iloc[:, 1].astype(str)
    vals = x.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").astype(np.float32)
    vals.columns = model_ids
    vals.index = sym
    vals = np.log2(vals.fillna(0) + 1)
    vals = vals.T.groupby(level=0).mean().T          # Sanger and Broad RNA-seq of the same model -> mean
    vals = to_canonical(vals).astype(np.float32)
    vals.to_pickle(DERIVED / "expr_cmp.pkl")
    print("CMP expression", vals.shape)
    return vals


def score(expr, tag):
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    sc = scoring.rank_score(expr, gs).astype(np.float32)
    sc.to_pickle(DERIVED / f"reactome_{tag}.pkl")
    print("scores", tag, sc.shape)


def ccle():
    for tag, folder in (("ctrp", "ctrpv2_ccle"), ("prism", "prism_repurposing")):
        e = pd.read_csv(CUR / folder / "expression.tsv.gz", sep="\t", index_col=0).astype(np.float32)
        e = to_canonical(e).astype(np.float32)
        e.to_pickle(DERIVED / f"expr_{tag}.pkl")
        score(e, tag)
        r = pd.read_csv(CUR / folder / "response.tsv", sep="\t", low_memory=False)
        r["drug"] = r["drug"].map(lambda s: common_drug_name(str(s).upper()))
        r.to_pickle(DERIVED / f"{tag}_long.pkl")
        pd.read_csv(CUR / folder / "samples.tsv", sep="\t").to_pickle(DERIVED / f"{tag}_samples.pkl")


def crispr():
    f = next((CACHE / "project_score_x").glob("*fold_change_values*.tsv"))
    head = pd.read_csv(f, sep="\t", nrows=4, header=None, index_col=0)
    ids, qc = head.loc["model_id"].values, head.loc["qc_pass"].astype(str).str.upper().values == "TRUE"
    x = pd.read_csv(f, sep="\t", skiprows=5, header=None, index_col=0, low_memory=False)
    x = x.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    x.columns = ids
    x = x.loc[:, qc]
    x = x.T.groupby(level=0).mean().T
    sym = x.index.map(lambda s: canonical_symbol(str(s)) or str(s))
    x.index = sym
    x = x.groupby(level=0).mean().astype(np.float32)
    x.to_pickle(DERIVED / "crispr_lfc.pkl")
    print("CRISPR", x.shape, "qc-pass screens", int(qc.sum()))


def mutations():
    m = pd.read_csv(next((CACHE / "cmp_mutsummary").glob("*.csv.gz")), low_memory=False)
    m = m[(m["coding"] == "t")]
    models_ = pd.read_pickle(DERIVED / "models.pkl")
    have = models_.index[models_["mutation_data"].astype(str).str.upper() == "TRUE"].union(pd.Index(m["model_id"].unique()).intersection(models_.index))
    out = pd.DataFrame(0, index=have, columns=["TP53", "KRAS", "BRAF", "APC", "PIK3CA"], dtype=float)
    for g in out.columns:
        hit = set(m.loc[m["gene_symbol"] == g, "model_id"])
        out[g] = [1.0 if i in hit else 0.0 for i in out.index]
    out["BRAF_V600E"] = [1.0 if i in set(m.loc[(m.gene_symbol == "BRAF") & (m.protein_mutation == "p.V600E"), "model_id"])
                         else 0.0 for i in out.index]
    out["MSI"] = models_.loc[out.index, "msi_status"].map({"MSI": 1.0, "MSS": 0.0})
    out.to_pickle(DERIVED / "mutations.pkl")
    print("mutations", out.shape, out.sum().to_dict())


if __name__ == "__main__":
    import sys as _s
    if "--rest" not in _s.argv:
        versions()
        gdsc()
        models()
        e = expression()
        score(e, "cmp")
        del e
        ccle()
    if "--mut" not in _s.argv:
        crispr()
    mutations()
