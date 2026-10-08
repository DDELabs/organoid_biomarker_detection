"""Heatmap of recurrent top pathways across cancer x drug pairs (standard + relaxed tiers).

Rows: the 20 pathways most often in a pair's top-5 response list and the 20 most often in a top-5
resistance list. Columns: cancer x drug pairs, standard tier (n >= 20 RECIST) then relaxed tiers,
grouped by drug. Colour: combined z (blue = higher pathway -> response/longer PFI, red = resistance),
clipped at +/-4. Dot = pathway is in that pair's top 10.

Also writes the same layout coloured by local z (that cancer's patients only): *_local.png.

Usage: python robust/plot_cancer_drug_heatmap.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from obd import RESULTS  # noqa: E402
from obd.atlas import DRUG_CLASS  # noqa: E402

D = RESULTS / "atlas" / "cancer_specific"
CMAP = LinearSegmentedColormap.from_list("div", ["#e34948", "#f0efec", "#2a78d6"])


def load():
    s = pd.read_csv(D / "all_pairs_local.tsv.gz", sep="\t").assign(tier="A20")
    r = pd.read_csv(D / "relaxed_all_pairs.tsv.gz", sep="\t")
    a = pd.concat([s, r], ignore_index=True)
    st = pd.read_csv(D / "top10_by_cancer_drug.tsv", sep="\t").assign(tier="A20")
    rt = pd.read_csv(D / "relaxed_top10.tsv", sep="\t")
    return a, pd.concat([st, rt], ignore_index=True)


def short(p):
    p = p.replace("REACTOME_", "").replace("_", " ").lower()
    return p if len(p) <= 48 else p[:46] + "…"


def main(value="combined_z"):
    a, top = load()
    a["pair"] = a["cancer"] + " · " + a["drug"].str.title().str.replace("Interferon Alfa-2B, Recombinant", "IFN-α2b")
    top["pair"] = top["cancer"] + " · " + top["drug"].str.title().str.replace("Interferon Alfa-2B, Recombinant", "IFN-α2b")
    t5 = top[top["rank"] <= 5]
    resp = t5[t5.direction == "response"].pathway.value_counts().head(20).index
    res = [p for p in t5[t5.direction == "resistance"].pathway.value_counts().index if p not in resp][:20]
    rows = list(resp) + res
    M = a.pivot_table(index="pathway", columns="pair", values=value).reindex(rows)
    meta = a.drop_duplicates("pair").set_index("pair")
    meta["cls"] = meta["drug"].map(DRUG_CLASS).fillna("other / targeted / hormonal")
    meta["std"] = meta["tier"].eq("A20")
    order = meta.sort_values(["std", "cls", "drug", "cancer"], ascending=[False, True, True, True]).index
    M = M[order]
    M = M.loc[M.mean(axis=1).sort_values(ascending=False).index]
    inset = set(zip(top.pathway, top.pair))
    nr, nc = M.shape
    fig, ax = plt.subplots(figsize=(0.2 * nc + 5.5, 0.24 * nr + 3.2))
    im = ax.imshow(M.clip(-4, 4).values, cmap=CMAP, vmin=-4, vmax=4, aspect="auto", interpolation="nearest")
    for i, p in enumerate(M.index):
        for j, c in enumerate(M.columns):
            if (p, c) in inset:
                ax.plot(j, i, "o", ms=2.6, mfc="#1a1a19", mec="none")
    ax.set_yticks(range(nr))
    ax.set_yticklabels([short(p) for p in M.index], fontsize=7)
    tiers = meta.loc[M.columns, "tier"]
    ax.set_xticks(range(nc))
    ax.set_xticklabels([f"{c} [{t}]" if t != "A20" else c for c, t in zip(M.columns, tiers)], rotation=90, fontsize=6)
    nstd = int(meta.loc[M.columns, "std"].sum())
    ax.axvline(nstd - 0.5, color="#1a1a19", lw=1.2)
    ax.text(nstd / 2, -0.9, "standard tier (RECIST, n ≥ 20)", ha="center", fontsize=8)
    ax.text(nstd + (nc - nstd) / 2, -0.9, "relaxed tiers: A RECIST n ≥ 10 · B first-course outcome · C PFI",
            ha="center", fontsize=8)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
    cb.set_label(f"{value.replace('_', ' ')}  (blue = response / longer PFI, red = resistance)", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    what = "combined z (local + pan-cancer ATLAS prior)" if value == "combined_z" else "local z (that cancer's patients only)"
    ax.set_title(f"Recurrent pathway biomarkers across cancer × drug pairs, TCGA — {what}; dot = in that pair's top 10",
                 fontsize=10, loc="left", pad=30)
    fig.tight_layout()
    tag = "" if value == "combined_z" else "_local"
    fig.savefig(D / f"heatmap_cancer_drug_pathways{tag}.png", dpi=170)
    fig.savefig(D / f"heatmap_cancer_drug_pathways{tag}.pdf")
    M.round(3).to_csv(D / f"heatmap_cancer_drug_pathways{tag}.tsv", sep="\t")
    print(M.shape)


if __name__ == "__main__":
    main("combined_z")
    main("local_z")
