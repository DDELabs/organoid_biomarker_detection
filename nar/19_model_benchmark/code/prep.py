"""Pre-build things that do not use outcomes: shortlist, pathway graph, proliferation list, prior."""
import sys; sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import run_survlearn as RS
from obd import network as N
from obd.reference import gene_sets
d, Z, P = breast_data()
cols = list(Z.columns)
DRUGS = ["PACLITAXEL","DOCETAXEL","DOXORUBICIN","EPIRUBICIN","CYCLOPHOSPHAMIDE","FLUOROURACIL","CARBOPLATIN","TRASTUZUMAB"]
paper = N.precomputed_proximity()
zs, zp = {}, {}
for dr in DRUGS:
    z = RS.proximity(dr)
    if z is not None: zs[dr] = z.reindex(cols)
    if dr in paper.columns: zp[dr] = paper[dr].reindex(cols)
ZS, ZP = pd.DataFrame(zs), pd.DataFrame(zp)
CUT = -1.2816
out = pd.DataFrame(index=cols)
out["min_z_stringv12"] = ZS.min(axis=1); out["min_z_paper"] = ZP.min(axis=1)
out["shortlist_stringv12"] = out.min_z_stringv12 <= CUT
out["shortlist_paper"] = out.min_z_paper <= CUT
# proliferation list (same keywords as W4)
kw = ["CELL_CYCLE","MITOTIC","G1_S","G2_M","E2F","DNA_REPLICATION","PRE_REPLICATIVE","CHECKPOINT"]
out["proliferation"] = [any(k in c for k in kw) for c in cols]
# pathway graph from gene overlap (Jaccard, 10 nearest neighbours)
gs = gene_sets(("REACTOME",))
miss = [c for c in cols if c not in gs]
print("pathways without gene set:", len(miss))
sets = [set(gs.get(c, [])) for c in cols]
n = len(cols); J = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        u = len(sets[i] | sets[j]); J[i, j] = J[j, i] = len(sets[i] & sets[j]) / u if u else 0
np.save(HERE/"tables"/"pathway_jaccard.npy", J)
out.to_csv(TAB/"pathway_flags.tsv", sep="\t")
print(ZS.shape, ZP.shape, out[["shortlist_stringv12","shortlist_paper","proliferation"]].sum().to_dict(), "union both:", (out.shortlist_stringv12|out.shortlist_paper).sum())
print("drugs string:", list(zs), "paper:", list(zp))
