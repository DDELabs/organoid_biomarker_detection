"""G13: gene composition of the 'odd' frozen pathways (APOBEC3G-mediated resistance to HIV-1, HIV infection) and of all
35 frozen pathways: overlap with proteasome, nuclear-pore, ribosome, RNA polymerase II and immune gene sets, the
most-overlapping canonical pathways, and correlation of the pathway score with tumour purity (ABSOLUTE) and with an
immune pathway score in TCGA (pan-cancer, within each cancer type).
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

from gcommon import NAR, OI, REPO, TAB, panel, savefig, style, version_rows

ODD = ["REACTOME_APOBEC3G_MEDIATED_RESISTANCE_TO_HIV1_INFECTION", "REACTOME_HIV_INFECTION"]
FROZEN = NAR / "08_depmap/tables/frozen_hypotheses.tsv"
PUR = REPO / "robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
SCORES = REPO / "robust/results/training_table/pathway_scores_tcga.tsv.gz"
IMMUNE = "REACTOME_ADAPTIVE_IMMUNE_SYSTEM"


def gmt():
    out = {}
    for line in open(REPO / "data/msigdb.v6.1.symbols.gmt.txt"):
        p = line.rstrip("\n").split("\t")
        out[p[0]] = set(p[2:])
    return out


def families(genes):
    fam = {"proteasome (PSM*)": r"^PSM[ABCDEFG]\d", "nuclear pore (NUP*, RAN*, KPNB1, XPO1, SEH1L, SEC13, TPR, AAAS, RAE1)":
           r"^(NUP\d|NUPL|RANBP2|RANGAP1|RAN$|KPNB1|XPO1|SEH1L|SEC13|TPR$|AAAS|RAE1|POM121|NDC1)",
           "ribosomal (RPL*/RPS*)": r"^RP[LS]\d", "RNA polymerase II (POLR2*)": r"^POLR2", "ubiquitin (UB*)": r"^(UBA52|UBB|UBC|RPS27A)$",
           "HLA / immune receptors": r"^(HLA-|CD\d|TRA|TRB|IG[HKL])", "HIV accessory host factors (PPIA, BANF1, HMGA1, PSIP1, XRCC*, LIG4)":
           r"^(PPIA|BANF1|HMGA1|PSIP1|XRCC[456]|LIG4|APOBEC3G)$"}
    return {k: sorted(g for g in genes if re.match(v, g)) for k, v in fam.items()}


def main():
    G = gmt()
    universe = set().union(*G.values())
    N = len(universe)
    fz = pd.read_csv(FROZEN, sep="\t")
    paths = list(dict.fromkeys(fz.pathway))
    refs = {"KEGG_PROTEASOME": "proteasome", "REACTOME_TRANSPORT_OF_MATURE_TRANSCRIPT_TO_CYTOPLASM": "nuclear export / pore",
            "REACTOME_NEP_NS2_INTERACTS_WITH_THE_CELLULAR_EXPORT_MACHINERY": "nuclear pore (NEP/NS2)",
            "KEGG_RIBOSOME": "ribosome", "REACTOME_IMMUNE_SYSTEM": "immune system", IMMUNE: "adaptive immune",
            "REACTOME_INNATE_IMMUNE_SYSTEM": "innate immune", "REACTOME_INTERFERON_SIGNALING": "interferon",
            "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION": "RNA pol II transcription", "REACTOME_CELL_CYCLE": "cell cycle"}
    refs = {k: v for k, v in refs.items() if k in G}
    comp, top = [], []
    for p in paths:
        g = G[p]
        fam = families(g)
        row = dict(pathway=p, n_genes=len(g), odd=p in ODD)
        for k, v in fam.items():
            row[f"n_{k}"] = len(v)
            row[f"share_{k}"] = len(v) / len(g)
        for k, lab in refs.items():
            ov = len(g & G[k])
            row[f"overlap_{lab}"] = ov
            row[f"share_in_{lab}"] = ov / len(g)
            row[f"hyper_p_{lab}"] = float(hypergeom.sf(ov - 1, N, len(G[k]), len(g)))
        comp.append(row)
        # most-overlapping canonical pathways (KEGG/Reactome/Biocarta/PID), by Jaccard
        cand = [(k, len(g & s) / len(g | s), len(g & s)) for k, s in G.items()
                if k != p and k.split("_")[0] in ("KEGG", "REACTOME", "BIOCARTA", "PID") and len(g & s) > 0]
        for k, j, ov in sorted(cand, key=lambda x: -x[1])[:5]:
            top.append(dict(pathway=p, other=k, jaccard=j, overlap=ov, n_other=len(G[k])))
        if p in ODD:
            (TAB / f"G13_genes_{p.replace('REACTOME_', '')}.txt").write_text(
                "\n".join(f"{x}\t" + ";".join(k for k, v in fam.items() if x in v) for x in sorted(g)) + "\n")
    C = pd.DataFrame(comp)
    C.to_csv(TAB / "G13_frozen_pathway_composition.tsv", sep="\t", index=False)
    pd.DataFrame(top).to_csv(TAB / "G13_most_overlapping_pathways.tsv", sep="\t", index=False)

    # TCGA: purity and immune-score correlation within each cancer type
    S = pd.read_csv(SCORES, sep="\t", index_col=0, usecols=lambda c: c in (["sample_id", IMMUNE] + paths))
    proj = S.index.str.split(":").str[0]
    pat = S.index.str.split(":").str[1].str[:12]
    pur = pd.read_csv(PUR, sep="\t")
    pur = pur.assign(patient=pur["array"].str[:12]).drop_duplicates("patient").set_index("patient")["purity"]
    S["purity"] = pur.reindex(pat).values
    S["project"] = proj
    rows = []
    for p in paths:
        if p not in S:
            continue
        for who, col in (("purity", "purity"), ("adaptive immune score", IMMUNE)):
            rr = []
            for pr, d in S.groupby("project"):
                d = d[[p, col]].dropna()
                if len(d) >= 30:
                    rr.append(spearmanr(d[p], d[col])[0])
            rr = np.array(rr)
            rows.append(dict(pathway=p, odd=p in ODD, vs=who, n_cancer_types=len(rr), median_rho=float(np.median(rr)),
                             q10=float(np.quantile(rr, 0.1)), q90=float(np.quantile(rr, 0.9))))
    P = pd.DataFrame(rows)
    P.to_csv(TAB / "G13_purity_immune_correlation.tsv", sep="\t", index=False)
    pd.DataFrame(version_rows([FROZEN, PUR, SCORES, REPO / "data/msigdb.v6.1.symbols.gmt.txt"])).to_csv(
        TAB / "G13_data_versions.tsv", sep="\t", index=False)
    figure(C, P, pd.DataFrame(top))


def figure(C, P, top):
    plt = style()
    fig, axs = plt.subplots(1, 3, figsize=(14, 5.2), gridspec_kw=dict(width_ratios=[1.1, 1.2, 1]))
    ax = axs[0]
    fams = [c for c in C.columns if c.startswith("share_") and not c.startswith("share_in_")]
    cols = [OI["blue"], OI["orange"], OI["green"], OI["purple"], OI["yellow"], OI["red"], OI["sky"]]
    for i, p in enumerate(C[C.odd].pathway):
        r = C[C.pathway == p].iloc[0]
        left = 0
        for f, c in zip(fams, cols):
            ax.barh(i, r[f], left=left, color=c, label=f.replace("share_", "") if i == 0 else None)
            left += r[f]
        ax.barh(i, 1 - left, left=left, color=OI["grey"], alpha=0.4, label="other" if i == 0 else None)
        ax.text(1.01, i, f"{int(r.n_genes)} genes", va="center", fontsize=7)
    ax.set_yticks(range(int(C.odd.sum())))
    ax.set_yticklabels([p.replace("REACTOME_", "").replace("_", " ").lower()[:40] for p in C[C.odd].pathway], fontsize=7.5)
    ax.set_xlim(0, 1.2)
    ax.set_xlabel("share of genes")
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    ax.set_title("Gene families in the two 'HIV' pathways")
    panel(ax, "a", -0.3)
    ax = axs[1]
    t = top[top.pathway.isin(C[C.odd].pathway)].copy()
    t["lab"] = t.pathway.str.replace("REACTOME_", "").str[:14] + " ~ " + t.other.str.replace("REACTOME_", "R_").str.lower().str[:44]
    ax.barh(range(len(t)), t.jaccard, color=[OI["blue"] if "APOBEC" in p else OI["orange"] for p in t.pathway])
    for i, r in enumerate(t.itertuples()):
        ax.text(r.jaccard + 0.01, i, f"{r.overlap} shared", va="center", fontsize=6.5)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(t.lab, fontsize=6.3)
    ax.invert_yaxis()
    ax.set_xlabel("Jaccard overlap (shared genes / all genes of both)")
    ax.set_title("Most-overlapping canonical pathways")
    panel(ax, "b", -0.9)
    ax = axs[2]
    for k, (who, c) in enumerate((("purity", OI["red"]), ("adaptive immune score", OI["green"]))):
        q = P[P.vs == who].reset_index(drop=True)
        y = np.arange(len(q))
        ax.scatter(q.median_rho, y + (k - 0.5) * 0.3, color=c, s=[30 if o else 10 for o in q.odd], label=who,
                   edgecolor=["k" if o else "none" for o in q.odd], lw=0.6)
    q = P[P.vs == "purity"].reset_index(drop=True)
    ax.set_yticks(range(len(q)))
    ax.set_yticklabels([p.replace("REACTOME_", "").lower()[:32] for p in q.pathway], fontsize=5.5)
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("median within-cancer Spearman rho (33 TCGA types)")
    ax.set_title("All 35 frozen pathways (big dots = 'HIV' ones)")
    ax.legend(frameon=False, fontsize=7)
    panel(ax, "c", -0.7)
    fig.tight_layout()
    a = C[C.odd]
    lines = []
    for r in a.itertuples():
        d = r._asdict()
        fam_txt = ", ".join(f"{k.replace('share_', '')}: {d[k]*1:.0%}" for k in d if isinstance(k, str) and k.startswith("share_")
                            and not k.startswith("share_in_") and d[k] > 0)
        lines.append(f"- {r.pathway}: {r.n_genes} genes; {fam_txt}")
    pp = P[P.odd]
    lines += [f"- {r.pathway} vs {r.vs}: median rho {r.median_rho:.2f} (10th-90th percentile {r.q10:.2f} to {r.q90:.2f})"
              for r in pp.itertuples()]
    savefig(fig, "G13_odd_pathway_composition", """
# What is inside the "HIV" pathways among the frozen hypotheses?

Two of the 35 frozen cell-line hypotheses carry virus names (Reactome "APOBEC3G-mediated resistance to HIV-1
infection" and "HIV infection"). Reactome describes how the virus uses human proteins, so these sets are made of human
genes. Here we show which human genes they contain.

- **a** Share of genes from common protein families (by gene-name pattern).
- **b** The five canonical pathways (KEGG, Reactome, BioCarta, PID; MSigDB v6.1) that overlap each set most.
- **c** For all 35 frozen pathways: how the pathway score tracks tumour purity (ABSOLUTE) and an adaptive-immune pathway
  score, as the median correlation within each TCGA cancer type. Large outlined dots = the two "HIV" sets.

""" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
