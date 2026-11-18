"""Schematic figures of the paper (Fig 1a, 1b, 3a, 4a, 4d) redrawn for our study, plus a map of all data sets
(Fig 1c). Pure drawings, except Fig 1c, whose numbers are read from the organoid-only summary tables."""
from common import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import networkx as nx

C_LAB, C_NET, C_ML, C_PAT, C_CHK = "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#999999"


def box(ax, x, y, w, h, text, col, fs=8, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=col, ec="black", lw=0.8, alpha=0.9))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", wrap=True)


def arrow(ax, a, b, txt=None, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=12, lw=1.1, color="black", linestyle=ls))
    if txt:
        ax.text((a[0] + b[0]) / 2 + 0.01, (a[1] + b[1]) / 2, txt, fontsize=7.5, ha="left", va="center", style="italic")


def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


# ---------------------------------------------------------------- Fig 1a workflow
fig, ax = plt.subplots(figsize=(7.6, 5.4)); blank(ax)
box(ax, 0.19, 0.90, 0.36, 0.11, "Lab models\n(organoids, PDX, cell lines)\nRNA + drug response", C_LAB, bold=True)
box(ax, 0.81, 0.90, 0.36, 0.11, "STRING protein network\n+ drug target proteins\n(DrugBank)", C_NET, bold=True)
box(ax, 0.81, 0.70, 0.38, 0.11, "Network proximity of each\nReactome pathway to the targets\n(z <= -1.2816 = 'proximal')", C_NET)
box(ax, 0.19, 0.70, 0.36, 0.11, "Pathway activity score per\nsample (ssGSEA NES)", C_LAB)
box(ax, 0.50, 0.50, 0.52, 0.10, "Train Ridge / SVR / LinearRegression:\nproximal-pathway scores -> drug response", C_ML, bold=True)
box(ax, 0.50, 0.33, 0.46, 0.08, "Rank-1 pathway (largest Ridge weight)\n= candidate biomarker", C_ML)
box(ax, 0.22, 0.13, 0.42, 0.13, "TCGA patients GIVEN the drug:\nscore, split at median,\nKaplan-Meier + log-rank", C_PAT, bold=True)
box(ax, 0.79, 0.13, 0.42, 0.13, "Checks: untreated patients, random\npathway sets (bootstrap), isogenic\nresistant cells, mutation concordance", C_CHK)
arrow(ax, (0.19, 0.845), (0.19, 0.755)); arrow(ax, (0.81, 0.845), (0.81, 0.755))
arrow(ax, (0.17, 0.65), (0.40, 0.555)); arrow(ax, (0.83, 0.65), (0.60, 0.555), "selects features")
arrow(ax, (0.50, 0.45), (0.50, 0.37)); arrow(ax, (0.42, 0.29), (0.27, 0.19)); arrow(ax, (0.58, 0.29), (0.73, 0.19))
ax.text(0.5, 0.995, "Kong et al. 2020 method, as re-run here for every cancer set", ha="center", va="top", fontsize=10, fontweight="bold")
save(fig, "Fig1a_workflow_schematic", """
# Fig. 1a (redrawn): the whole method in one picture

Read from top to bottom. Blue = lab models. Orange = the protein network step. Green = the machine-learning step.
Pink = the patient test. Grey = extra checks.
1. Lab models (organoids, PDX or cell lines) give two things: gene activity (RNA) and how well the drug kills them.
2. Each gene set from Reactome (a "pathway") gets one activity score per sample (ssGSEA NES).
3. In the STRING protein network we measure how close every pathway is to the drug's target proteins.
   Pathways much closer than random (z <= -1.2816) are kept.
4. Three simple models learn drug response from the kept pathways. The pathway with the largest Ridge weight is the
   candidate biomarker.
5. TCGA patients who received the drug are scored with that pathway, split at the median, and their survival compared.
6. Checks: the same split in untreated patients, random pathway sets, drug-resistant cell lines and known mutations.
This is a drawing; it contains no data.
""")

# ---------------------------------------------------------------- Fig 1b network proximity
rng = np.random.default_rng(3)
G = nx.connected_watts_strogatz_graph(40, 4, 0.25, seed=3)
pos = nx.spring_layout(G, seed=4)
targ = [0, 1]
nbr = sorted(set(G[0]) | set(G[1]) - set(targ))
path_genes = nbr[:3]
_dd = {n: min(nx.shortest_path_length(G, n, t) for t in targ) for n in G.nodes}
far = sorted(G.nodes, key=lambda n: -_dd[n])[:3]
fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.8), gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
for ax, (genes, title) in zip(axs[:2], ((path_genes, "Proximal pathway\n(short distance to targets)"),
                                        (far, "Distant pathway\n(long distance)"))):
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.25, width=0.6)
    col = [OI[3] if n in targ else (OI[0] if n in genes else "#DDDDDD") for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=col, node_size=45, edgecolors="black", linewidths=0.3)
    d = np.mean([min(nx.shortest_path_length(G, g, t) for t in targ) for g in genes])
    ax.set_title(f"{title}\nmean closest distance d = {d:.1f}", fontsize=8.5); ax.axis("off")
letter(axs[0], "a"); letter(axs[1], "b")
ax = axs[2]
null = rng.normal(2.6, 0.35, 1000)
ax.hist(null, bins=30, color="#BBBBBB", edgecolor="white")
ax.axvline(1.7, color=OI[3], lw=2); ax.text(1.72, ax.get_ylim()[1] * 0.85, "observed d", color=OI[3], fontsize=8)
ax.set_xlabel("closest distance d of 1,000 random\ngene sets with matched degree"); ax.set_ylabel("count")
ax.set_title("z = (d - mean random) / SD random\nproximal if z <= -1.2816", fontsize=8.5); letter(ax, "c")
fig.tight_layout()
save(fig, "Fig1b_network_proximity_schematic", """
# Fig. 1b (redrawn): how "network proximity" is measured

a, b. A toy protein network (grey dots = proteins, lines = known interactions). Red dots = the drug's target
proteins. Blue dots = the genes of one pathway. For each pathway gene we count the steps to the nearest target and
average them (closest distance d). In (a) the pathway sits next to the targets; in (b) it sits far away.
c. To know whether d is small, the same distance is computed for 1,000 random gene sets whose genes have the same
number of network links (degree-matched). The z-score compares the real d with these random ones. A pathway is
called "proximal" when z <= -1.2816 (the lowest 10% of the random sets), the paper's cut-off.
Illustration only: the network and the histogram are simulated drawings, not data.
""")

# ---------------------------------------------------------------- Fig 3a comparison feature sets
fig, ax = plt.subplots(figsize=(7.4, 4.0)); blank(ax)
items = [("Network-proximal\npathways (paper)", C_NET), ("Whole\ntranscriptome", C_LAB), ("All Reactome\npathways", C_LAB),
         ("Network neighbours\nof targets\n(1, 2, 3 steps)", C_NET), ("Drug target\ngenes alone", C_NET),
         ("Genes correlated\nwith response\n(Bolis)", C_LAB), ("Deep learning\n(MOLI-style)", C_CHK),
         ("Most central\npathways (degree,\nbetweenness, closeness)", C_NET)]
for i, (t, c) in enumerate(items):
    x = 0.125 + (i % 4) * 0.25; y = 0.86 if i < 4 else 0.56
    box(ax, x, y, 0.225, 0.22, t, c, fs=8)
ax.add_patch(FancyBboxPatch((0.0, 0.42), 1.0, 0.56, boxstyle="round,pad=0.005", fc="none", ec="#666666", lw=0.8, ls="--"))
arrow(ax, (0.5, 0.41), (0.5, 0.27), "each feature set in turn")
box(ax, 0.5, 0.15, 0.94, 0.18, "Same task for every feature set: train on lab models -> score TCGA patients\ngiven the drug -> split at the median -> log-rank test of survival", C_ML, bold=True)
save(fig, "Fig3a_comparison_feature_sets_schematic", """
# Fig. 3a (redrawn): the other ways to build a biomarker that the network method is compared with

Every box is one way to choose the input features. Each is trained on the same lab models and tested on the same
treated TCGA patients, so the results can be compared directly (paper Fig. 3b-f and Supp. Figs. 6-9).
Orange = uses the protein network; blue = uses only expression; grey = deep learning.
Where the results are: bladder in nar/02_organoid_only/BLCA_lee2018/extras (S_BLCA_comparison_feature_sets);
every set in nar/07_B3_models (benchmark, produced by the B3 agent); target genes and 2- and 3-step neighbours for
every set in nar/10_completeness (SuppFig8_SuppFig9_*). Drawing only.
""")

# ---------------------------------------------------------------- Fig 4a isogenic
fig, ax = plt.subplots(figsize=(7.0, 3.0)); blank(ax)
box(ax, 0.14, 0.70, 0.24, 0.16, "Parental cancer cell line\n(drug-sensitive)", C_LAB, bold=True)
box(ax, 0.14, 0.25, 0.24, 0.16, "Same line made resistant\nby long drug exposure", C_PAT, bold=True)
arrow(ax, (0.14, 0.61), (0.14, 0.34), "months of drug")
box(ax, 0.50, 0.48, 0.26, 0.20, "RNA of both lines\n-> score of the biomarker\npathway (ssGSEA)", C_ML)
box(ax, 0.85, 0.48, 0.26, 0.26, "t-test sensitive vs resistant.\nPass = difference in the\ndirection the organoid\nmodel predicts", C_CHK)
arrow(ax, (0.27, 0.68), (0.37, 0.53)); arrow(ax, (0.27, 0.27), (0.37, 0.43)); arrow(ax, (0.63, 0.48), (0.72, 0.48))
save(fig, "Fig4a_isogenic_validation_schematic", """
# Fig. 4a (redrawn): the isogenic cell-line check

"Isogenic" = two versions of the same cell line that differ only in drug resistance: the parental line and a copy
made resistant by growing it in the drug for months. Because the genetic background is the same, a pathway that truly
tracks resistance should differ between the two. The biomarker pathway is scored in both and compared with a t-test.
It passes when the difference goes the way the organoid model says (for example, a pathway whose high score means
"sensitive" in organoids should be higher in the parental line). Paper: HCT-8 vs HCT-8/5-FU (colon) and T24 vs
T24-cisplatin-resistant (bladder). Our results: nar/02_organoid_only (colon, bladder), nar/08_depmap S7 (more pairs),
nar/10_completeness (Fig4bc_isogenic_other_cancers, other cancers). Drawing only.
""")

# ---------------------------------------------------------------- Fig 4d bootstrap
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
ax = axs[0]; blank(ax)
box(ax, 0.2, 0.75, 0.34, 0.2, "Real proximal\npathways (n = k)", C_NET, bold=True)
box(ax, 0.2, 0.25, 0.34, 0.2, "10,000 random sets of\nk pathways", C_CHK, bold=True)
box(ax, 0.75, 0.5, 0.42, 0.3, "same training ->\nrank-1 pathway ->\ntreated-patient log-rank p", C_ML)
arrow(ax, (0.37, 0.75), (0.54, 0.58)); arrow(ax, (0.37, 0.25), (0.54, 0.42)); letter(ax, "a")
ax = axs[1]
null = np.random.default_rng(1).uniform(0, 1, 10000) ** 0.8
ax.hist(-np.log10(null), bins=40, color="#BBBBBB", edgecolor="white")
ax.axvline(-np.log10(0.014), color=OI[3], lw=2)
ax.text(-np.log10(0.014) + 0.05, ax.get_ylim()[1] * 0.8, "real set", color=OI[3], fontsize=8)
ax.set_xlabel("-log10 p (random sets)"); ax.set_ylabel("count"); letter(ax, "b")
ax.set_title("share of random sets with p as small\n= bootstrap p", fontsize=8.5)
fig.tight_layout()
save(fig, "Fig4d_bootstrap_schematic", """
# Fig. 4d (redrawn): the random-pathway bootstrap

a. The whole pipeline is repeated with random sets of pathways of the same size as the real proximal set.
b. Each random set gives a treated-patient log-rank p. The share of random sets whose p is as small as (or smaller
than) the real one is the bootstrap p. A small share means the network-chosen pathways are special.
The histogram in (b) is simulated for illustration. Real results: Fig2a_Fig2d_Fig4e panels (colon, pancreas, liver,
PDX) and *_full panels h-i (bladder, ovary, lung, sarcoma, PDX) in nar/02_organoid_only.
""")

# ---------------------------------------------------------------- Fig 1c map of all sets (real numbers)
g = pd.read_csv(O2 / "tables/cross_set_summary.tsv", sep="\t")
o = pd.read_csv(O2 / "BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv", sep="\t")
rows = [(r.set, r.cohort, r.patient_drug, r.n_models, r.treated_n) for r in g.itertuples()]
rows += [(r.set_label, r.project, r.drug, r.n_models, r.treated_n) for r in o.itertuples()]
m = pd.DataFrame(rows, columns=["set", "cohort", "drug", "n_models", "treated_n"])
m.to_csv(TAB / "Fig1c_data_map.tsv", sep="\t", index=False)
m["row"] = m.set.str.replace(r",? ?(paper source-data NES|current GDC data)", "", regex=True).str.slice(0, 45) + " -> " + m.cohort
piv = m.pivot_table(index="row", columns="drug", values="treated_n", aggfunc="max")
fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(piv) + 1.6))
im = ax.imshow(np.log10(piv.values.astype(float)), cmap="viridis", aspect="auto")
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, int(v), ha="center", va="center", fontsize=6.5, color="white" if np.log10(v) < 2.0 else "black")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels([c.lower() for c in piv.columns], rotation=60, ha="right")
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=7)
cb = fig.colorbar(im, ax=ax, shrink=0.5); cb.set_label("log10 treated patients")
ax.set_title("Every lab-model set x patient cohort x drug analysed (number = drug-treated TCGA patients)")
save(fig, "Fig1c_data_map_all_sets", f"""
# Fig. 1c (new): map of every analysis

Rows: lab-model set and the TCGA cohort it was tested in. Columns: drugs. Each number is how many TCGA patients
received that drug and had RNA and survival data (the size of the survival test). Empty = not analysed (no screen of
that drug in the lab models, or fewer than 12 treated patients). {len(m)} set x cohort x drug analyses in total.
Source: nar/02_organoid_only/tables/cross_set_summary.tsv and
nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv (copied to tables/Fig1c_data_map.tsv).
""")
print("done", len(m))
