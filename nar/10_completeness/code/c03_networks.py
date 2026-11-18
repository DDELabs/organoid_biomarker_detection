"""Fig 2c / 2g equivalents: the drug's target proteins and the genes of the rank-1 pathway on the paper's STRING>700
network, for the paper's two cases and every other set's headline drug."""
from common import *
import networkx as nx

G = paper_network()
R = reactome()
keys = [k for k in REG if k != "LIHC_broutier_SOR"] + ["LIHC_broutier_SOR"]
extra = {"BLCA_paper_AMINO": ("Bladder, paper's reported biomarker | cisplatin", "CISPLATIN",
                              "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION")}
jobs = [(k, REG[k][0], REG[k][2], rank1(k)[0]) for k in keys]
jobs.insert(2, ("BLCA_paper_AMINO",) + extra["BLCA_paper_AMINO"])
rows = []
n = len(jobs); nc = 4; nr = int(np.ceil(n / nc))
fig, axs = plt.subplots(nr, nc, figsize=(11, 2.9 * nr))
for i, (k, lab, drug, pw) in enumerate(jobs):
    ax = axs.flat[i]
    T = [t for t in targets(drug) if t in G]
    P = [g for g in R.get(pw, []) if g in G]
    d = {}
    for g in P:
        dl = [nx.shortest_path_length(G, g, t) for t in T] if T else []
        d[g] = min(dl) if dl else np.nan
    dist = pd.Series(d, dtype=float)
    show = list(dist.sort_values().index[:30])
    nodes = set(T) | set(show)
    # connectors: one shortest path from each shown gene to its nearest target (max 2 intermediate genes shown)
    for g in show:
        if T and np.isfinite(dist[g]) and dist[g] <= 3:
            tt = min(T, key=lambda t: nx.shortest_path_length(G, g, t))
            nodes |= set(nx.shortest_path(G, g, tt))
    H = G.subgraph(nodes)
    pos = nx.spring_layout(H, seed=1, k=0.6)
    col = [OI[3] if v in T else (OI[0] if v in P else "#CCCCCC") for v in H.nodes]
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.35, width=0.5)
    nx.draw_networkx_nodes(H, pos, ax=ax, node_color=col, node_size=28, edgecolors="black", linewidths=0.25)
    lab_nodes = {v: v for v in H.nodes if v in T or (v in P and dist.get(v, 9) <= 1)}
    nx.draw_networkx_labels(H, pos, labels=lab_nodes, ax=ax, font_size=5.5)
    direct = int((dist == 1).sum()); md = float(dist.mean()) if len(dist) else np.nan
    ax.set_title(f"{lab}\n{short(pw, 52)}\n{len(T)} targets, {len(P)} pathway genes, {direct} touch a target, mean d = {md:.2f}",
                 fontsize=6.5)
    ax.axis("off"); letter(ax, "abcdefghijklmnopqrstuvwxyz"[i])
    rows.append(dict(key=k, label=lab, drug=drug, pathway=pw, targets=",".join(T), n_targets_in_network=len(T),
                     n_pathway_genes=len(R.get(pw, [])), n_pathway_genes_in_network=len(P),
                     n_genes_directly_linked_to_target=direct, mean_closest_distance=md, genes_shown=len(show)))
for j in range(i + 1, nr * nc):
    axs.flat[j].axis("off")
fig.tight_layout()
t = pd.DataFrame(rows); t.to_csv(TAB / "Fig2c_2g_network_neighbourhoods.tsv", sep="\t", index=False)
save(fig, "Fig2c_Fig2g_network_neighbourhoods_all_sets", f"""
# Fig. 2c and 2g (and the same picture for every other set): where the biomarker pathway sits in the protein network

Each panel shows part of the paper's STRING protein network (interaction score > 700). Red dots = the drug's target
proteins (DrugBank; carboplatin and oxaliplatin use the cisplatin targets, as in the analyses). Blue dots = genes of
the rank-1 pathway (the biomarker chosen by the organoid model); up to 30 pathway genes closest to the targets are
drawn. Grey dots = connector proteins on the shortest path between a pathway gene and its nearest target.
Labels are shown for targets and for pathway genes that touch a target directly.
Panel title: number of targets, pathway genes in the network, pathway genes touching a target, and mean closest
distance d (the quantity behind the proximity z-score).
Panels a-c are the paper's own cases: colon 5-FU (TYMS and the BH3-only pathway = paper Fig. 2c), bladder cisplatin
with our rank-1 pathway (b) and with the paper's reported amino acid pathway (c, = paper Fig. 2g).
Numbers: tables/Fig2c_2g_network_neighbourhoods.tsv. Broutier uses the same drug as Ji (sorafenib) but its own rank-1 pathway.
""")
print(t[["key", "n_targets_in_network", "n_pathway_genes_in_network", "n_genes_directly_linked_to_target", "mean_closest_distance"]].to_string())
