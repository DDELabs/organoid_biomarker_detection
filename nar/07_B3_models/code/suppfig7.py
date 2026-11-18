"""Supp Fig 6/7 replica (centrality-based pathway selection) and a network-centrality profile of the selected
pathways and drug targets, on the paper's STRING>700 network.

Paper Supp Fig 6 (COAD 5-FU) / 7 (BLCA cisplatin): pick the |proximal| most central pathways (mean degree,
betweenness or closeness of their genes) instead of the network-proximal ones, rank them by ridge weight on the
organoids, and draw KM curves of the top 3 in drug-treated TCGA patients.
"""
import numpy as np
import pandas as pd

import b3common as B
from bench import Ctx, ridge_paper, zy
from figures import km, letter

plt = B.plot_style()
PAIRS = [("COAD_paper", "FLUOROURACIL", "SuppFig6"), ("BLCA_paper", "CISPLATIN", "SuppFig7")]


def replica():
    rows = []
    for b, d, stem in PAIRS:
        C = Ctx(b, d)
        fig, axs = plt.subplots(3, 4, figsize=(14, 9.5))
        for i, cen in enumerate(("degree", "betweenness", "closeness")):
            cols = C.central[cen]
            X = B.zcols(C.Xo_pw[:, cols])[0]
            c = ridge_paper(X, zy(C.y))
            o = np.argsort(-np.abs(c))
            ax = axs[i, 0]
            top = o[:10]
            ax.barh(range(len(top))[::-1], c[top], color=[B.OK["verm"] if v > 0 else B.OK["blue"] for v in c[top]])
            ax.set_yticks(range(len(top))[::-1])
            ax.set_yticklabels([B.short(C.U[cols[j]], 38) for j in top], fontsize=6.5)
            ax.axvline(0, color="black", lw=0.6)
            ax.set_xlabel("ridge weight (organoids)")
            ax.set_title(f"{cen}: {len(cols)} most central pathways", fontsize=8.5)
            letter(ax, "abc"[i], x=-0.9)
            for j in range(3):
                jj = o[j]
                risk = c[jj] * C.Xp_pw[C.ti, cols[jj]]
                km(axs[i, j + 1], C.T[C.ti], C.E[C.ti], risk, f"rank {j + 1}: {B.short(C.U[cols[jj]], 40)}")
                ms, _ = B.median_split(C.T[C.ti], C.E[C.ti], risk)
                rows.append(dict(pair=f"{b} {d}", centrality=cen, n_selected=len(cols), rank=j + 1, pathway=C.U[cols[jj]],
                                 coef=c[jj], n=len(C.ti), events=int(C.E[C.ti].sum()), **ms,
                                 C_index=B.cindex(C.T[C.ti], C.E[C.ti], risk)))
        fig.tight_layout()
        paper = {"SuppFig6": "paper: degree 0.997 / 0.71 / 0.802; betweenness 0.79 / 0.503 / 0.91; closeness 0.17 / 0.71 / 0.32",
                 "SuppFig7": "paper: degree 0.094 / 0.96 / 0.51; betweenness 0.23 / 0.31 / 0.091; closeness 0.56 / 0.29 / 0.74"}[stem]
        B.savefig(fig, f"{stem}_replica_{b}_{d}", f"""
**Paper {stem.replace('SuppFig', 'Supp Fig ')} replica ({b.replace('_paper', '')} organoids, {d.lower()}; paper source data).** Instead of the
network-proximal pathways, the {len(C.central['degree'])} pathways whose genes are most central in the STRING>700 network are used
(rows: mean degree, betweenness, closeness of the pathway genes). Left: ridge weights of the 10 top-ranked pathways on the organoids
(red = higher score, more resistant). Right: Kaplan-Meier curves of drug-treated TCGA patients split at the median of each of the top 3
pathways (log-rank p, HR of predicted non-responders). Log-rank p of the top 3 in the {paper}.""")
        plt.close("all")
    pd.DataFrame(rows).to_csv(B.TABLES / "SuppFig6_7_centrality_selection.tsv", sep="\t", index=False)


def profile():
    G = pd.read_csv(B.DATA / "centrality_genes.tsv", sep="\t", index_col=0)
    Pw = pd.read_csv(B.DATA / "centrality_pathways.tsv", sep="\t", index_col=0)
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
    rows = []
    for i, (b, d, _) in enumerate(PAIRS):
        C = Ctx(b, d)
        tg = [t for t in B.targets(d) if t in G.index]
        Fp = [C.U[j] for j in C.F]
        genes_F = sorted({g for p in Fp for g in B.reactome().get(p, []) if g in G.index})
        allpw = sorted({g for p in C.U for g in B.reactome().get(p, []) if g in G.index})
        for j, cen in enumerate(("degree", "betweenness")):
            ax = axs[i, j]
            groups = [("all network genes", G[cen]), ("genes in any pathway", G.loc[allpw, cen]),
                      ("genes in proximal pathways", G.loc[genes_F, cen]), ("drug targets", G.loc[tg, cen])]
            data = [np.log10(v.values + (1 if cen == "degree" else 1e-7)) for _, v in groups]
            ax.boxplot(data, showfliers=False, widths=0.6)
            for k, v in enumerate(data):
                if len(v) <= 30:
                    ax.scatter(np.full(len(v), k + 1) + np.random.default_rng(0).uniform(-0.1, 0.1, len(v)), v, s=14,
                               color=B.OK["verm"], zorder=3)
            ax.set_xticks(range(1, 5))
            ax.set_xticklabels([g for g, _ in groups], rotation=20, fontsize=7)
            ax.set_ylabel(f"log10 {cen}")
            ax.set_title(f"{b.replace('_paper', '')} {d.lower()}: {cen}", fontsize=8.5)
            letter(ax, "abcdef"[i * 3 + j])
            from scipy.stats import mannwhitneyu
            rows.append(dict(pair=f"{b} {d}", centrality=cen, n_targets=len(tg), targets=",".join(tg),
                             median_targets=G.loc[tg, cen].median(), median_proximal_genes=G.loc[genes_F, cen].median(),
                             median_all_pathway_genes=G.loc[allpw, cen].median(),
                             mwu_p_proximal_vs_all=mannwhitneyu(G.loc[genes_F, cen], G.loc[allpw, cen]).pvalue,
                             n_proximal_genes=len(genes_F)))
        ax = axs[i, 2]
        z = C.z.dropna()
        pw = Pw.reindex(z.index)
        prox = z <= B.CUTOFF
        ax.scatter(pw.loc[~prox, "degree"], z[~prox], s=6, color=B.OK["grey"], alpha=0.5, label="other pathways")
        ax.scatter(pw.loc[prox, "degree"], z[prox], s=12, color=B.OK["blue"], label=f"proximal (z <= {B.CUTOFF})")
        ax.axhline(B.CUTOFF, color="black", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("mean degree of pathway genes")
        ax.set_ylabel("network proximity z-score")
        from scipy.stats import spearmanr
        r, p = spearmanr(pw["degree"], z, nan_policy="omit")
        ax.set_title(f"pathway degree vs proximity: Spearman {r:.2f} (p = {p:.2g})", fontsize=8)
        ax.legend(fontsize=6.5, frameon=False)
        letter(ax, "abcdef"[i * 3 + 2])
        rows.append(dict(pair=f"{b} {d}", centrality="pathway degree vs z", spearman=r, p=p))
    fig.tight_layout()
    pd.DataFrame(rows).to_csv(B.TABLES / "SuppFig7_centrality_profile.tsv", sep="\t", index=False)
    B.savefig(fig, "SuppFig7_network_centrality_profile", """
**Network centrality of the drug targets and of the selected pathways (STRING>700, the paper's network).** Rows: colorectal
5-FU (top) and bladder cisplatin (bottom). (a, b, d, e) Degree (number of direct partners) and betweenness (how often a protein lies
on shortest paths between others; estimated from 1,000 sampled sources), log scale, for all network genes, genes in any Reactome
pathway, genes in the network-proximal pathways, and the drug targets (red dots). (c, f) Each Reactome pathway's mean gene degree
against its network-proximity z-score; blue = proximal pathways used by the paper. A strong correlation would mean the proximity
score mostly picks hub-rich pathways; the Spearman correlation is printed in the title. Numbers: tables/SuppFig7_centrality_profile.tsv.""")
    plt.close("all")


if __name__ == "__main__":
    replica()
    profile()
