"""Gene and pathway network centrality on the paper's STRING>700 network (MOESM3 'Figure 2a-f').

Degree: exact. Betweenness: networkx sampled estimate (k = 1000 source nodes, seed 0).
Closeness: harmonic-free classic closeness estimated from 1000 sampled BFS sources (mean shortest-path distance
from the sampled sources to each node; closeness = 1 / mean distance).
Pathway centrality = mean centrality of the pathway's genes that are in the network (as Supp Fig 6/7 of the paper).
Writes _data/centrality_genes.tsv and _data/centrality_pathways.tsv.
"""
import numpy as np
import pandas as pd
import networkx as nx

from b3common import DATA, paper_network, reactome

g = paper_network()
nodes = list(g.nodes)
print("network", g.number_of_nodes(), g.number_of_edges(), flush=True)
deg = pd.Series(dict(g.degree()), dtype=float)
btw = pd.Series(nx.betweenness_centrality(g, k=1000, seed=0, normalized=True))
rng = np.random.default_rng(0)
src = rng.choice(nodes, 1000, replace=False)
tot = pd.Series(0.0, index=nodes)
for s in src:
    d = nx.single_source_shortest_path_length(g, s)
    tot += pd.Series(d).reindex(nodes).fillna(0)
clo = 1.0 / (tot / len(src)).replace(0, np.nan)
G = pd.DataFrame({"degree": deg, "betweenness": btw, "closeness": clo})
G.to_csv(DATA / "centrality_genes.tsv", sep="\t")
rows = {}
for p, genes in reactome().items():
    gg = [x for x in genes if x in G.index]
    if gg:
        rows[p] = G.loc[gg].mean().to_dict() | {"n_genes_in_network": len(gg)}
pd.DataFrame(rows).T.to_csv(DATA / "centrality_pathways.tsv", sep="\t")
print("done", flush=True)
