"""Shared paths, registry of headline biomarkers and small helpers for the completeness agent (nar/10_completeness).

Reads (never writes) other agents' outputs:
  * B3 bundles (nar/07_B3_models/_data/raw/bundles/*.pkl): lab-model + TCGA pathway scores, gene matrices, survival,
    treated / untreated patients, proximity z per drug (built from the organoid-only agents' caches);
  * rank-1 Ridge pathway and coefficient of each set x drug (nar/02_organoid_only/**/coefficients*.tsv);
  * the paper's STRING>700 network (MOESM3 'Figure 2a-f', extracted by the organoid-only agent).
"""
import os, sys, pickle, json
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
NAR = REPO / "nar"
FIG = OUT / "figures"; TAB = OUT / "tables"; DATA = OUT / "_data"
for d in (FIG, TAB, DATA):
    d.mkdir(exist_ok=True)
O2 = NAR / "02_organoid_only"
BUNDLES = NAR / "07_B3_models" / "_data" / "raw" / "bundles"
MOESM3 = O2 / "BLCA_OV_LUAD_OTHER" / "_data" / "moesm3"
PROXCACHE = O2 / "BLCA_OV_LUAD_OTHER" / "_data" / "cache"
sys.path.insert(0, str(REPO / "robust"))
from obd import survival as SV  # noqa: E402
from obd.reference import gene_sets, canonical_symbol, drug_targets, drugbank_targets  # noqa: E402
SEED = 20261117
CUTOFF = -1.2816
PROXY = {"CARBOPLATIN": "CISPLATIN", "OXALIPLATIN": "CISPLATIN"}

# headline biomarker per set: key -> (label, bundle, drug, coefficient file (relative to nar/02_organoid_only), group)
REG = {
    "COAD_paper_5FU": ("Colon organoids, paper data | 5-FU", "COAD_paper", "FLUOROURACIL",
                       "COAD_vandewetering/tables/coefficients_COAD_authors_FLUOROURACIL.tsv", "colon"),
    "BLCA_paper_CIS": ("Bladder organoids, paper data | cisplatin", "BLCA_paper", "CISPLATIN",
                       "BLCA_lee2018/paper/tables/BLCA_CISPLATIN_coefficients.tsv", "bladder"),
    "COADREAD_5FU": ("Colorectal organoids, current data | 5-FU", "COADREAD_vdW", "FLUOROURACIL",
                     "COAD_vandewetering/tables/coefficients_COADREAD_FLUOROURACIL.tsv", "colon"),
    "BLCA_cur_CIS": ("Bladder organoids, current data | cisplatin", "BLCA_current", "CISPLATIN",
                     "BLCA_lee2018/current/tables/BLCA_CISPLATIN_coefficients.tsv", "bladder"),
    "PAAD_tiriac_GEM": ("Pancreas organoids (Tiriac) | gemcitabine", "PAAD_tiriac", "GEMCITABINE",
                        "PAAD_tiriac/tables/coefficients_PAAD_GEMCITABINE.tsv", "pancreas"),
    "PAAD_shi_GEM": ("Pancreas organoids (Shi) | gemcitabine", "PAAD_shi", "GEMCITABINE",
                     "PAAD_shi/tables/coefficients_PAAD_GEMCITABINE.tsv", "pancreas"),
    "LIHC_ji_SOR": ("Liver organoids (Ji) | sorafenib", "LIHC_ji", "SORAFENIB",
                    "LIHC_ji/tables/coefficients_LIHC_SORAFENIB.tsv", "liver"),
    "LIHC_broutier_SOR": ("Liver organoids (Broutier) | sorafenib", "LIHC_ji", "SORAFENIB",
                          "LIHC_broutier/tables/coefficients_LIHC_SORAFENIB.tsv", "liver"),
    "OV_vias_CARBO": ("Ovarian organoids (Vias) | carboplatin", "OV_vias", "CARBOPLATIN",
                      "OV_vias2023/tables/OV_CARBOPLATIN_coefficients.tsv", "ovarian"),
    "LUAD_cells_CIS": ("Lung adenocarcinoma cell lines (substitute) | cisplatin", "LUAD_cells", "CISPLATIN",
                       "LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_coefficients.tsv", "lung"),
    "SARC_GEM": ("Sarcoma organoids | gemcitabine", "SARC_alshihabi", "GEMCITABINE",
                 "SARC_alshihabi2024/tables/SARC_GEMCITABINE_coefficients.tsv", "sarcoma"),
    "COAD_isella_CET": ("Colorectal PDX (Isella) | cetuximab", "COAD_isella", "CETUXIMAB",
                        "COAD_isella_pdx/tables/coefficients_COAD_CETUXIMAB.tsv", "PDX"),
    "PDX_CM_DTIC": ("Melanoma PDX (Gao) | dacarbazine", "PDX_CM", "DACARBAZINE",
                    "PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_coefficients.tsv", "PDX"),
    "PDX_NSCLC_ERL": ("Lung PDX (Gao) | erlotinib", "PDX_NSCLC", "ERLOTINIB",
                      "PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_coefficients.tsv", "PDX"),
    "PDX_PDAC_GEM": ("Pancreas PDX (Gao) | gemcitabine", "PDX_PDAC", "GEMCITABINE",
                     "PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_coefficients.tsv", "PDX"),
    "PDX_BRCA_PTX": ("Breast PDX (Gao) | paclitaxel", "PDX_BRCA", "PACLITAXEL",
                     "PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_coefficients.tsv", "PDX"),
    "PDX_CRC_5FU": ("Colorectal PDX (Gao) | 5-FU", "PDX_CRC", "FLUOROURACIL",
                    "PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_coefficients.tsv", "PDX"),
}

_B = {}
def bundle(name):
    if name not in _B:
        _B[name] = pickle.load(open(BUNDLES / f"{name}.pkl", "rb"))
    return _B[name]


def coefs(key):
    f = O2 / REG[key][3]
    c = pd.read_csv(f, sep="\t", index_col=0)
    col = "Ridge" if "Ridge" in c.columns else "coef_Ridge"
    s = c[col].astype(float)
    return s.reindex(s.abs().sort_values(ascending=False).index)


def rank1(key):
    s = coefs(key)
    return s.index[0], float(s.iloc[0])


def km_split(T, E, score):
    """Median split; high score = predicted resistant ('non-responder'). Returns dict with log-rank p and plain HR."""
    T = np.asarray(T, float); E = np.asarray(E, int); s = np.asarray(score, float)
    hi = s > np.median(s)
    if hi.sum() in (0, len(hi)):
        return dict(n=len(s), events=int(E.sum()), p=np.nan, HR=np.nan, lo=np.nan, hi=np.nan), hi
    p = SV.logrank(T[hi], E[hi], T[~hi], E[~hi])
    HR = lo = up = np.nan
    if E[hi].sum() > 0 and E[~hi].sum() > 0:
        from statsmodels.duration.hazard_regression import PHReg
        try:
            r = PHReg(np.clip(T, 1e-6, None), hi[:, None].astype(float), status=E.astype(float), ties="efron").fit(disp=False)
            b, se = float(r.params[0]), float(r.bse[0]); HR, lo, up = np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
        except Exception:
            pass
    return dict(n=len(s), events=int(E.sum()), ev_resp=int(E[~hi].sum()), ev_nonresp=int(E[hi].sum()), p=p, HR=HR, lo=lo, hi=up), hi


def treated_split(b, drug, score_series):
    pats = [p for p in b["drugs"][drug]["treated"] if p in score_series.index and p in b["surv"].index]
    sv = b["surv"].loc[pats]
    return km_split(sv["T"].values, sv["E"].values, score_series.loc[pats].values)


def untreated_split(b, score_series):
    pats = [p for p in b["untreated"] if p in score_series.index and p in b["surv"].index]
    sv = b["surv"].loc[pats]
    return km_split(sv["T"].values, sv["E"].values, score_series.loc[pats].values)


def reactome():
    return dict(gene_sets(("REACTOME",)))


def paper_network():
    import networkx as nx
    e = pd.read_csv(MOESM3 / "Figure_2a_f.tsv", sep="\t")
    a = e["Gene_Symbol_1"].astype(str).map(lambda s: canonical_symbol(s) or s)
    c = e["Gene_Symbol_2"].astype(str).map(lambda s: canonical_symbol(s) or s)
    g = nx.Graph(); g.add_edges_from(zip(a, c)); g.remove_edges_from(nx.selfloop_edges(g))
    return g.subgraph(max(nx.connected_components(g), key=len)).copy()


def targets(drug):
    d = PROXY.get(drug, drug)
    f = PROXCACHE / f"proximity_papernet_{d}.tsv.targets"
    if f.exists():
        return open(f).read().strip().split(",")
    return list(drug_targets().get(d) or drugbank_targets().get(d) or [])


# ------------------------------------------------------------------ plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
# Okabe-Ito colour-blind-safe palette
OI = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]


def save(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    (FIG / f"{name}.caption.md").write_text(caption.strip() + "\n")
    plt.close(fig)


def letter(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom", ha="left")


def short(p, n=48):
    s = p.replace("REACTOME_", "").replace("_", " ").lower()
    s = s[0].upper() + s[1:]
    return s if len(s) <= n else s[:n - 3] + "..."
