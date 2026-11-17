"""User's lung attempt: GSE165611 lung cancer patient-derived cells (PDC), expression only.

The per-PDC drug table ('supp_table_S2b.txt') is not public: GEO has expression only (sample characteristics list tissue,
disease, histology), and the paper's supplement (Yu et al. 2023 J Exp Clin Cancer Res 42:37, Additional files 1-3 on
static-content.springer.com) lists the 48 screened drugs (Table S2, incl. cisplatin, pemetrexed, gemcitabine, paclitaxel,
etoposide, olaparib) but no per-sample response values. So the expression is downloaded, assembled and compared with
TCGA-LUAD (where the PDCs sit relative to patients), and the lung analysis itself uses GDSC/CTRP cell lines (run_sets.py LUAD).
"""
import io
import re
import tarfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import common as K  # noqa: E402
import engine as E  # noqa: E402

OUT = K.PARENT / "LUAD_GSE165611_substitute"
FIG, TAB = OUT / "figures", OUT / "tables"


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    f = K.CACHE / "GSE165611_log2FPKM.tsv.gz"
    if not f.exists():
        cols = {}
        with tarfile.open(K.RAW / "geo" / "GSE165611_RAW.tar") as t:
            for m in t.getmembers():
                if not m.name.endswith(".txt.gz"):
                    continue
                name = re.search(r"(NCCLu-\d+)", m.name).group(1)
                d = pd.read_csv(io.BytesIO(t.extractfile(m).read()), sep="\t", compression="gzip")
                gid = d.iloc[:, 0].astype(str)
                sym = gid.str.split("_").str[-1] if gid.str.contains("_").all() else gid
                cols[name] = pd.Series(d["FPKM"].values, index=sym.values).groupby(level=0).sum()
        X = np.log2(pd.DataFrame(cols) + 1)
        X.to_csv(f, sep="\t")
    X = pd.read_csv(f, sep="\t", index_col=0)
    brief = open(K.RAW / "geo" / "GSE165611_gsm_brief.txt").read().split("^SAMPLE")[1:]
    meta = {}
    for b in brief:
        title = re.search(r"Sample_title = (\S+)", b).group(1)
        meta[title] = {k: re.search(rf"{k}: ([^\n]+)", b).group(1).strip() for k in ("tissue", "histopathology")}
    meta = pd.DataFrame(meta).T
    meta.to_csv(TAB / "GSE165611_sample_annotation.tsv", sep="\t")
    Xc = K.C.to_canonical(X)
    texpr = K.tcga_expr("LUAD")
    genes = Xc.index.intersection(texpr.index)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    rows = []
    for ax, scaled, L in zip(axes, (False, True), "ab"):
        A, Bm = Xc.loc[genes].T.values, texpr.loc[genes].T.values
        if scaled:
            A, Bm = StandardScaler().fit_transform(A), StandardScaler().fit_transform(Bm)
        M = np.vstack([Bm, A])
        M = M[:, M.std(axis=0) > 0]
        pc = PCA(2).fit(M)
        Z = pc.transform(M)
        ax.scatter(Z[: len(Bm), 0], Z[: len(Bm), 1], s=5, color=E.OK["grey"], alpha=0.5, label=f"TCGA-LUAD (n={len(Bm)})")
        h = meta.reindex(Xc.columns)["histopathology"].fillna("NA")
        for hh, col in zip(sorted(h.unique()), [E.OK["verm"], E.OK["blue"], E.OK["green"], E.OK["purple"]]):
            m = (h == hh).values
            ax.scatter(Z[len(Bm):][m, 0], Z[len(Bm):][m, 1], s=18, color=col, label=f"PDC: {hh} (n={m.sum()})")
        ev = pc.explained_variance_ratio_ * 100
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.set_title("genes, " + ("z-scored within each dataset" if scaled else "not scaled"))
        ax.legend(frameon=False, fontsize=7)
        E.panel_letter(ax, L)
        rows.append(dict(scaled=scaled, pc1=ev[0], pc2=ev[1]))
    fig.tight_layout()
    fig.savefig(FIG / "GSE165611_PDC_vs_TCGA_PCA.png", dpi=300)
    fig.savefig(FIG / "GSE165611_PDC_vs_TCGA_PCA.pdf")
    pd.DataFrame(rows).to_csv(TAB / "GSE165611_PCA_variance.tsv", sep="\t", index=False)
    open(FIG / "GSE165611_PDC_vs_TCGA_PCA.caption.md", "w").write(
        f"# Lung patient-derived cells (GSE165611) next to TCGA-LUAD tumours\n\n{X.shape[1]} PDC RNA-seq profiles (RSEM FPKM from GEO, "
        f"log2(FPKM+1)) and TCGA-LUAD primary tumours, on {len(genes)} shared genes. a: raw; b: each gene z-scored within its own dataset "
        "(the paper's batch correction). Colours: PDC histology from GEO. The drug-response values for these PDCs are not public, "
        "so they cannot be used to train the model; the lung analysis uses GDSC/CTRP lung adenocarcinoma cell lines instead.\n")
    print(X.shape, meta.histopathology.value_counts().to_dict())


if __name__ == "__main__":
    main()
