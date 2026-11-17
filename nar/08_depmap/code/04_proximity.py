"""Network proximity (STRING v12, combined score >= 700) of Reactome pathways to each Stage-0 drug's targets.

Reuses robust/obd/network.proximity_z (closest distance, degree-matched null, 1000 randomisations).
Targets: DrugBank (robust/obd/reference.drugbank_targets); platinum drugs have no protein target in
DrugBank, so carboplatin and oxaliplatin use the cisplatin target set, as in robust/run_survlearn.py (PROXY).
The paper's own STRING>700 table (data/coad_blca_organoid_drugs_zscore_result_reactome.txt) is used as well
where it has the drug.

python -I nar/08_depmap/code/04_proximity.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import networkx as nx  # noqa: E402
import pandas as pd  # noqa: E402

from common import CACHE, TAB  # noqa: E402
from obd import network as N  # noqa: E402
from obd.reference import canonical_symbol, drugbank_targets, gene_sets  # noqa: E402

DRUGS = ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CISPLATIN", "CARBOPLATIN", "GEMCITABINE", "PACLITAXEL",
         "DOXORUBICIN", "PEMETREXED", "SORAFENIB"]
PROXY = {"CARBOPLATIN": "CISPLATIN", "OXALIPLATIN": "CISPLATIN"}


def string_graph(min_score=700):
    d = CACHE / "string"
    info = pd.read_csv(d / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
    ensp2sym = {e: canonical_symbol(s) for e, s in zip(info.iloc[:, 0], info.iloc[:, 1])}
    links = pd.read_csv(d / "9606.protein.links.v12.0.txt.gz", sep=" ")
    links = links[links.combined_score >= min_score]
    a, b = links.protein1.map(ensp2sym), links.protein2.map(ensp2sym)
    ok = a.notna() & b.notna() & (a != b)
    g = nx.Graph()
    g.add_edges_from(zip(a[ok], b[ok]))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def main():
    g = string_graph()
    print("STRING v12 nodes/edges", g.number_of_nodes(), g.number_of_edges(), flush=True)
    gs = gene_sets(("REACTOME",))
    db = drugbank_targets()
    out, tg = {}, []
    done = {}
    for d in DRUGS:
        src = PROXY.get(d, d)
        targets = db.get(src, [])
        tg.append({"drug": d, "target_source": f"DrugBank ({src})", "targets": ",".join(targets)})
        if src not in done:
            done[src] = N.proximity_z(g, targets, gs, n_random=1000)
            print(d, "done", flush=True)
        out[d] = done[src]
    z = pd.DataFrame(out)
    paper = N.precomputed_proximity()
    for d in DRUGS:
        src = PROXY.get(d, d)
        z[f"{d}__paper_table"] = paper[src] if src in paper.columns else float("nan")
    z.to_csv(TAB / "B_network_proximity_z.tsv", sep="\t")
    pd.DataFrame(tg).to_csv(TAB / "B_drug_targets_used.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
