"""Drug-target to pathway network proximity.

Two complementary measures on a protein interaction network:
  * closest-distance proximity z-score (Guney et al. 2016), the measure used
    by the original pipeline, with degree-matched random target/pathway sets;
  * random-walk-with-restart (RWR) propagation from the drug targets, scored as
    the mean propagated signal in the pathway against degree-matched random
    seeds. RWR uses all paths rather than only the shortest one, so it is less
    sensitive to single missing or spurious edges.
A pathway is called proximal when it passes the cut-off on every measure
available (consensus).
"""
import random

import networkx as nx
import numpy as np
import pandas as pd

from . import DATA


def load_network(path, min_score=700, sep=None):
    """Edge list with columns (gene_a, gene_b[, score]); keeps the largest component."""
    df = pd.read_csv(path, sep=sep, engine="python", header=None, comment="#")
    if df.shape[1] >= 3 and np.issubdtype(df[2].dtype, np.number):
        df = df[df[2] >= min_score]
    g = nx.Graph()
    g.add_edges_from(zip(df[0].astype(str), df[1].astype(str)))
    g.remove_edges_from(nx.selfloop_edges(g))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def degree_bins(g, min_bin_size=100):
    by_deg = {}
    for n, d in g.degree():
        by_deg.setdefault(d, []).append(n)
    bins, cur = [], []
    for d in sorted(by_deg):
        cur += by_deg[d]
        if len(cur) >= min_bin_size:
            bins.append(cur)
            cur = []
    if cur:
        if bins:
            bins[-1] += cur
        else:
            bins.append(cur)
    lookup = {n: b for b in bins for n in b}
    return lookup


def _random_like(nodes, lookup, rng):
    out = set()
    for n in nodes:
        pool = lookup[n]
        for _ in range(20):
            c = rng.choice(pool)
            if c not in out:
                break
        out.add(c)
    return out


def proximity_z(g, targets, pathways, n_random=1000, seed=452456, min_bin_size=100, chunk=400):
    """Closest-distance proximity z-score of every pathway to one drug's targets.

    d = mean over targets of the shortest-path distance to the nearest pathway
    gene; the null draws degree-matched random target and pathway sets.
    Vectorised: BFS distances from all (real + random) target nodes are computed
    once with scipy's C implementation and reused for every pathway.
    """
    from scipy.sparse.csgraph import shortest_path

    order = list(g)
    idx = {n: i for i, n in enumerate(order)}
    targets = [t for t in targets if t in idx]
    if not targets:
        return pd.Series(dtype=float)
    lookup = degree_bins(g, min_bin_size)
    bins = {id(b): np.array([idx[n] for n in b]) for b in {id(b): b for b in lookup.values()}.values()}
    node_bin = {idx[n]: bins[id(b)] for n, b in lookup.items()}
    rng = np.random.default_rng(seed)

    def random_sets(members, n):
        return np.column_stack([rng.choice(node_bin[m], n) for m in members])

    t_idx = np.array([idx[t] for t in targets])
    rand_t = random_sets(t_idx, n_random)                       # n_random x |targets|
    sources = np.unique(np.r_[t_idx, rand_t.ravel()])
    row = {s: i for i, s in enumerate(sources)}
    A = nx.to_scipy_sparse_array(g, nodelist=order, format="csr")
    D = np.empty((len(sources), len(order)), dtype=np.uint8)
    for i in range(0, len(sources), chunk):
        d = shortest_path(A, unweighted=True, directed=False, indices=sources[i:i + chunk])
        D[i:i + chunk] = np.where(np.isinf(d), 255, d).astype(np.uint8)
    t_rows = np.array([row[s] for s in t_idx])
    rt_rows = np.vectorize(row.get)(rand_t)
    out = {}
    for name, genes in pathways.items():
        p_idx = np.array([idx[x] for x in genes if x in idx])
        if len(p_idx) == 0:
            continue
        obs = D[np.ix_(t_rows, p_idx)].min(axis=1).mean()
        rp = random_sets(p_idx, n_random)                       # n_random x |pathway|
        null = np.array([D[np.ix_(rt_rows[r], rp[r])].min(axis=1).mean() for r in range(n_random)])
        sd = null.std()
        out[name] = 0.0 if sd == 0 else (obs - null.mean()) / sd
    return pd.Series(out)


def rwr_z(g, targets, pathways, n_random=200, restart=0.5, seed=1, min_bin_size=100):
    """RWR propagation z-score (higher = closer; sign flipped to match proximity)."""
    nodes = set(g)
    targets = [t for t in targets if t in nodes]
    if not targets:
        return pd.Series(dtype=float)
    lookup = degree_bins(g, min_bin_size)
    rng = random.Random(seed)
    order = list(g)
    idx = {n: i for i, n in enumerate(order)}
    pw = {k: [idx[x] for x in v if x in idx] for k, v in pathways.items()}
    pw = {k: v for k, v in pw.items() if v}

    def propagate(seeds):
        p = nx.pagerank(g, alpha=1 - restart, personalization={s: 1.0 for s in seeds}, tol=1e-8)
        return np.array([p[n] for n in order])

    obs = propagate(targets)
    null = np.vstack([propagate(_random_like(targets, lookup, rng)) for _ in range(n_random)])
    out = {}
    for k, ii in pw.items():
        o, nl = obs[ii].mean(), null[:, ii].mean(axis=1)
        out[k] = 0.0 if nl.std() == 0 else -(o - nl.mean()) / nl.std()
    return pd.Series(out)


def load_string(min_score=700):
    """STRING v12 human network (combined score >= min_score) on canonical gene symbols."""
    from . import EXTERNAL
    from .reference import canonical_symbol

    info = pd.read_csv(EXTERNAL / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
    ensp2sym = {e: canonical_symbol(s) for e, s in zip(info.iloc[:, 0], info.iloc[:, 1])}
    links = pd.read_csv(EXTERNAL / "9606.protein.links.v12.0.txt.gz", sep=" ")
    links = links[links.combined_score >= min_score]
    a, b = links.protein1.map(ensp2sym), links.protein2.map(ensp2sym)
    ok = a.notna() & b.notna() & (a != b)
    g = nx.Graph()
    g.add_edges_from(zip(a[ok], b[ok]))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def precomputed_proximity():
    """The original STRING>700 Reactome proximity z-scores (pathways x drugs)."""
    df = pd.read_csv(DATA / "coad_blca_organoid_drugs_zscore_result_reactome.txt", sep="\t", index_col=0)
    cols = {}
    for c in df.columns:
        for d in c.split("_"):
            cols[d] = df[c]
    return pd.DataFrame(cols)


def proximal(zs, cutoff=-1.2816):
    """Consensus: pathways at or below the cut-off in every supplied z-score Series."""
    zs = [z.dropna() for z in zs if z is not None and len(z)]
    common = set.intersection(*(set(z.index) for z in zs))
    return sorted(p for p in common if all(z[p] <= cutoff for z in zs))
