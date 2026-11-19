"""Shared set-up for the nar/22 analyses.  Reuses nar/02_organoid_only/code unchanged and redirects every write to
nar/22_user_data_results/_data.  One CPU only.  The TCGA drug tables are the files the user uploaded
(data/user_upload/TCGA_LUAD, TCGA_OVARY); TCGA LUSC uses the repo's copy (not uploaded)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import glob, sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
NAR = REPO / "nar"
DATA = OUT / "_data"
RAW = DATA / "raw"
CACHE = DATA / "cache"
for p in (RAW, CACHE, OUT / "tables", OUT / "figures", OUT / "logs"):
    p.mkdir(parents=True, exist_ok=True)


def link(src, dst):
    dst = Path(dst)
    if dst.is_symlink() or dst.exists():
        return
    os.symlink(str(src), str(dst))


def make_raw():
    for d in (NAR / "02_organoid_only" / "_data" / "raw", NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "_data" / "raw",
              NAR / "09_registered_test" / "_data" / "tcga_response" / "raw"):
        for p in d.iterdir():
            if p.name.startswith("nationwidechildrens.org_clinical_drug_") and p.name.endswith(("luad.txt", "ov.txt")):
                continue                                  # replaced by the user's files below
            link(p, RAW / p.name)
    for tag, sub in (("luad", "TCGA_LUAD"), ("ov", "TCGA_OVARY")):
        f = glob.glob(str(REPO / "data" / "user_upload" / sub / "clinical_drug" / "*" / f"nationwidechildrens.org_clinical_drug_{tag}.txt"))[0]
        link(f, RAW / f"nationwidechildrens.org_clinical_drug_{tag}.txt")
    # repo's own copies kept under another name for the comparison table
    for tag in ("luad", "ov"):
        link(NAR / "09_registered_test" / "_data" / "tcga_response" / "raw" / f"nationwidechildrens.org_clinical_drug_{tag}.txt", RAW / f"repo_copy_clinical_drug_{tag}.txt")


make_raw()
sys.path.insert(0, str(NAR / "02_organoid_only" / "code"))
sys.path.insert(0, str(REPO / "robust"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common as C  # noqa: E402   (nar/02 code)
from obd import scoring, network as N  # noqa: E402
from obd.reference import drug_targets, drugbank_targets, gene_sets, canonical_symbol  # noqa: E402

C.RAW = RAW
C.C.EXTERNAL = RAW
C.SCORES = DATA / "scores"
_ss = scoring.ssgsea
scoring.ssgsea = lambda expr, gs, min_size=2, threads=1: _ss(expr, gs, min_size=min_size, threads=1)

MOESM3 = NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "_data" / "moesm3"
MANUAL_TARGETS = {"ETOPOSIDE": ["TOP2A", "TOP2B"], "OLAPARIB": ["PARP1", "PARP2"], "NIRAPARIB": ["PARP1", "PARP2"],
                  "ADAVOSERTIB": ["WEE1"], "AZD7762": ["CHEK1", "CHEK2"], "BERZOSERTIB": ["ATR"], "BARASERTIB": ["AURKB"],
                  "CAPIVARSERTIB": ["AKT1", "AKT2", "AKT3"], "EVEROLIMUS": ["MTOR"], "VEMURAFENIB": ["BRAF"],
                  "DINACICLIB": ["CDK1", "CDK2", "CDK5", "CDK9"], "AT7519": ["CDK1", "CDK2", "CDK4", "CDK5", "CDK9"],
                  "AZD4547": ["FGFR1", "FGFR2", "FGFR3"], "ALECTINIB": ["ALK", "RET"], "PD0166285": ["WEE1", "CHEK1"],
                  "REPOTRECTINIB": ["ROS1", "NTRK1", "NTRK2", "NTRK3", "ALK"], "CABOZANTINIB": ["MET", "VEGFR2", "KDR", "AXL", "RET", "KIT", "FLT3"],
                  "DASATINIB": ["ABL1", "SRC", "LCK", "KIT", "PDGFRB", "EPHA2"], "AFATINIB": ["EGFR", "ERBB2", "ERBB4"],
                  "PAZOPANIB": ["KDR", "FLT1", "FLT4", "PDGFRA", "PDGFRB", "KIT"]}
_PAPER_NET = None


def paper_network():
    global _PAPER_NET
    if _PAPER_NET is None:
        import networkx as nx
        e = pd.read_csv(MOESM3 / "Figure_2a_f.tsv", sep="\t")
        a = e["Gene_Symbol_1"].astype(str).map(lambda s: canonical_symbol(s) or s)
        b = e["Gene_Symbol_2"].astype(str).map(lambda s: canonical_symbol(s) or s)
        g = nx.Graph()
        g.add_edges_from(zip(a, b))
        g.remove_edges_from(nx.selfloop_edges(g))
        _PAPER_NET = g.subgraph(max(nx.connected_components(g), key=len)).copy()
    return _PAPER_NET


PROX_LOG = []


def proximity(drug):
    """(z Series, column name).  Paper table first (with the platinum proxy of nar/02); otherwise recomputed on the
    paper's own STRING network (1000 degree-matched random target sets), cached."""
    paper = N.precomputed_proximity()
    d = C.PROX_PROXY.get(drug, drug)
    if d in paper.columns:
        PROX_LOG.append((drug, d, "paper table"))
        return paper[d].dropna(), d
    f = CACHE / f"prox_{drug}.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0], drug + " (recomputed)"
    t = MANUAL_TARGETS.get(drug) or drug_targets().get(drug) or drugbank_targets().get(drug)
    if not t:
        return None, d
    pw = {k: v for k, v in gene_sets(("REACTOME",)).items() if k in paper.index}
    z = N.proximity_z(paper_network(), list(t), pw, n_random=1000)
    z.to_frame("z").to_csv(f, sep="\t")
    open(str(f) + ".targets", "w").write(",".join(t))
    return z, drug + " (recomputed)"


C.proximity = proximity
