import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
"""Build data/external/gi_networks/{name}.tsv.gz + SOURCE.md (scratch build script)."""
import hashlib
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
SCR = Path("/tmp/obd_raw")
GH = SCR / "gh"
OUT = Path(_REPO + '/data/external/gi_networks')
OUT.mkdir(parents=True, exist_ok=True)
ROWS = []


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(name, df, desc, srcs):
    df = df.copy()
    df["gene_a"] = df["gene_a"].astype(str).str.strip()
    df["gene_b"] = df["gene_b"].astype(str).str.strip()
    df = df[(df.gene_a != "") & (df.gene_b != "") & (df.gene_a != "nan") & (df.gene_b != "nan")]
    cols = ["gene_a", "gene_b", "type", "score", "source"] + [c for c in df.columns if c not in ("gene_a", "gene_b", "type", "score", "source")]
    df = df[cols]
    df.to_csv(OUT / f"{name}.tsv.gz", sep="\t", index=False, na_rep="NA", compression="gzip")
    genes = len(set(df.gene_a) | set(df.gene_b))
    ROWS.append(dict(name=name, n_pairs=len(df), n_genes=genes, types=";".join(f"{k}:{v}" for k, v in df["type"].value_counts().items()), description=desc,
                     sources=" | ".join(f"{u} (sha256 {h})" for u, h in srcs)))
    print(name, len(df), genes, df["type"].value_counts().to_dict())


SL = GH / "synthLethal/data"
SL_URL = "https://raw.githubusercontent.com/bhklab/synthLethal/main/data/"

# 1. SynLethDB 2.0 human SL (as mirrored in bhklab/synthLethal Human_SL.csv)
d = pd.read_csv(SL / "Human_SL.csv", encoding="utf-8-sig", low_memory=False)
df = pd.DataFrame({"gene_a": d["n1.name"], "gene_b": d["n2.name"], "type": "SL",
                   "score": pd.to_numeric(d["r.statistic_score"], errors="coerce"),
                   "source": "SynLethDB2.0:" + d["r.source"].astype(str), "pubmed_id": d["r.pubmed_id"], "cell_line": d["r.cell_line"]})
save("synlethdb_human_sl", df, "SynLethDB 2.0 human SL pairs (all evidence types: CRISPR/RNAi screens, computational prediction, text mining, ...). score = SynLethDB statistic_score (0-1).",
     [(SL_URL + "Human_SL.csv", sha(SL / "Human_SL.csv"))])

# 2. SLKB (CRISPR combinatorial screens)
d = pd.read_csv(SL / "SLKB_pairs.csv", low_memory=False)
df = pd.DataFrame({"gene_a": d["gene1"], "gene_b": d["gene2"], "type": np.where(d["SL_or_not"].astype(str).str.upper().str.startswith("SL"), "SL", "notSL"),
                   "score": d["SL_score"], "source": "SLKB:PMID" + d["study_origin"].astype(str) + ":" + d["cell_line_origin"].astype(str),
                   "statistical_score": d["statistical_score"], "sl_score_cutoff": d["SL_score_cutoff"]})
save("slkb_crispr_pairs", df, "SLKB combinatorial-CRISPR SL calls (as mirrored in bhklab/synthLethal SLKB_pairs.csv). score = SL_score (more negative = stronger SL).",
     [(SL_URL + "SLKB_pairs.csv", sha(SL / "SLKB_pairs.csv"))])


# 3-5. ISLE (Lee et al. 2018) Cytoscape sessions
def cys_edges(path):
    txt = Path(path).read_text()
    out = []
    for lab in re.findall(r'<edge id="\d+" label="([^"]*)"', txt):
        m = re.match(r"(.+?) \((.+?)\) (.+)", lab)
        if m:
            out.append((m.group(1), m.group(3)))
    return pd.DataFrame(out, columns=["gene_a", "gene_b"])


ISLE_URL = "https://github.com/jooslee/ISLE/blob/master/networks/"
for name, sub, cys, desc in [
    ("isle_clinical_sl_fdr0.1", "ISLE_clinical_SL_network_FDR_0.1", "ISLE_clinical_SL_network_FDR_0.1.cys", "ISLE genome-wide clinically relevant SL network, FDR 0.1 (Lee et al. Nat Commun 2018; 2,326 pairs)."),
    ("isle_clinical_sl_fdr0.2", "ISLE_clinical_SL_network_FDR_0.2", "ISLE_clinical_SL_network_FDR_0.2.cys", "ISLE genome-wide clinically relevant SL network, FDR 0.2 (Lee et al. 2018)."),
    ("isle_drug_target_csl", "ISLE_drug_cSL_network", "ISLE_drug_cSL_network.cys", "ISLE drug-target clinical SL (cSL) network used for TCGA drug-response prediction (gene_b = drug target)."),
]:
    x = next((SCR / "cys" / sub).glob("*/networks/*.xgmml"))
    df = cys_edges(x)
    df["type"] = "SL"
    df["score"] = np.nan
    df["source"] = "ISLE:" + name
    save(name, df, desc + " Edges parsed from the Cytoscape session (.cys) XGMML; no per-edge scores in the session.",
         [(ISLE_URL + cys, sha(GH / "ISLE/networks" / cys))])

# 6. ISLE gold-standard experimental SL set (sl.golden.set.RData)
import rdata  # noqa: E402

g = rdata.read_rda(str(GH / "ISLE/data/sl.golden.set.RData"))["gd"]
names = np.asarray(g["sr0"]).astype(str).reshape(2, -1).T if np.asarray(g["sr0"]).size == 2 * len(g["flag"]) else None
dat = np.asarray(g["dat"]).astype(str).reshape(4, -1).T
# sr0 is the column-major flattening of the 2-column gene-name matrix
df = pd.DataFrame({"gene_a": names[:, 0], "gene_b": names[:, 1], "type": np.where(np.asarray(g["flag"]) == 1, "SL", "notSL"), "score": np.asarray(g["flag"], dtype=float),
                   "source": ["ISLE_gold:" + ";".join(sorted(set(r) - {"", "NA"})) for r in dat]})
save("isle_gold_standard_sl", df, "ISLE experimentally reported SL gold-standard set (literature screens: shRNA/sgRNA/drug/mutation screens), sl.golden.set.RData. source = screen(s).",
     [("https://github.com/jooslee/ISLE/blob/master/data/sl.golden.set.RData", sha(GH / "ISLE/data/sl.golden.set.RData"))])

# 7. bhklab ISLE re-implementation (pan-cancer) incl. SR (DD) statistics
f = SL / "SL_pairs_pan-cancer_binarize_expression_cox.csv"
d = pd.read_csv(f)
df = pd.DataFrame({"gene_a": d["gene1"], "gene_b": d["gene2"], "type": "SL", "score": d["q_value"], "source": "bhklab_synthLethal_ISLE-like_pancancer"})
for c in ["p_value", "depletion_q_value", "SR_DD_q_value", "survival_coef", "survival_p_value", "phylo_coefficient"]:
    df[c] = d[c]
save("bhklab_pancancer_sl", df, "bhklab/synthLethal ISLE-style pan-cancer SL inference (TCGA co-inactivation depletion + CRISPR + survival + phylogeny); score = q_value. Columns SR_DD_* are the synthetic-rescue (down-down) test.",
     [(SL_URL + f.name, sha(f))])

# 8. bhklab TCGA drug cSL (identical to ISLE drug cSL) - skipped as duplicate; CSL_pairs.csv identical to ISLE FDR0.2 (checked below)
d = pd.read_csv(SL / "CSL_networks/CSL_pairs.csv")
print("bhklab CSL_pairs rows", len(d))

pd.DataFrame(ROWS).to_csv(SCR / "gi_catalog.tsv", sep="\t", index=False)
