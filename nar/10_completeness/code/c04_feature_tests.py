"""Supp Fig 4 (cross-test matrix), Supp Fig 8 (drug-target genes alone), Supp Fig 9 (2- and 3-step network
neighbourhoods), Supp Fig 10 (pathway size vs survival p) and Supp Fig 16 (proximity of all drugs) for the paper's
two cases and every other set."""
from common import *
import networkx as nx
from scipy.stats import spearmanr, linregress
from sklearn.linear_model import RidgeCV

KEYS = [k for k in REG if k != "LIHC_broutier_SOR"]
AMINO = "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION"
bio = {k: rank1(k) for k in KEYS}
bio["BLCA_paper_AMINO"] = (AMINO, -0.14713652465773913)   # paper's reported biomarker, coef from isogenic table (paper data)
BLAB = {k: REG[k][0] for k in KEYS}; BLAB["BLCA_paper_AMINO"] = "Bladder, paper's amino acid pathway | cisplatin"

def zs(df):  # z-score each row (feature) across samples
    return df.sub(df.mean(1), axis=0).div(df.std(1).replace(0, np.nan), axis=0).fillna(0)

# ------------------------------------------------------------------ Supp Fig 4: cross matrix
rows = []
for kb, (pw, c) in bio.items():
    for kc in KEYS:
        lab, bn, drug, _, _ = REG[kc]
        b = bundle(bn)
        if pw not in b["tsc"].index:
            rows.append(dict(biomarker=kb, cohort=kc, pathway=pw, p=np.nan, note="pathway not scored in cohort")); continue
        r, _ = treated_split(b, drug, b["tsc"].loc[pw] * np.sign(c))
        rows.append(dict(biomarker=kb, cohort=kc, pathway=pw, own=(kb == kc) or (kb == "BLCA_paper_AMINO" and kc == "BLCA_paper_CIS"),
                         n=r["n"], events=r["events"], p=r["p"], HR=r["HR"], HR_lo=r["lo"], HR_hi=r["hi"]))
X = pd.DataFrame(rows); X.to_csv(TAB / "SuppFig4_cross_test_matrix.tsv", sep="\t", index=False)
P = X.pivot(index="biomarker", columns="cohort", values="p").reindex(index=list(bio), columns=KEYS)
H = X.pivot(index="biomarker", columns="cohort", values="HR").reindex(index=list(bio), columns=KEYS)
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(-np.log10(P.values.astype(float)), cmap="Blues", vmin=0, vmax=3)
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        p, h = P.values[i, j], H.values[i, j]
        if np.isfinite(p):
            own = (P.index[i] == P.columns[j]) or (P.index[i] == "BLCA_paper_AMINO" and P.columns[j] == "BLCA_paper_CIS")
            ax.text(j, i, f"{p:.2f}" if p >= 0.01 else f"{p:.0e}", ha="center", va="center", fontsize=5.5,
                    color="white" if p < 0.01 else "black", fontweight="bold" if own else "normal")
            if own:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=OI[3], lw=1.5))
ax.set_xticks(range(len(KEYS))); ax.set_xticklabels([REG[k][0] for k in KEYS], rotation=70, ha="right", fontsize=6.5)
ax.set_yticks(range(len(bio))); ax.set_yticklabels([f"{BLAB[k]}: {short(bio[k][0], 34)}" for k in bio], fontsize=6.5)
ax.set_xlabel("treated TCGA patients (cohort | drug)"); ax.set_ylabel("biomarker (lab-model set | drug: rank-1 pathway)")
cb = fig.colorbar(im, ax=ax, shrink=0.5); cb.set_label("-log10 log-rank p")
save(fig, "SuppFig4_cross_test_matrix_all_sets", """
# Supp. Fig. 4 (extended to every set): is each biomarker specific to its own cancer and drug?

Each row is one biomarker: the rank-1 pathway of one lab-model set and drug, with the sign of its organoid Ridge
weight. Each column is one group of drug-treated TCGA patients. Every cell applies the row's biomarker to the column's
patients (median split, log-rank p printed). Red frames = the biomarker in its own patients (the main result).
Darker blue = smaller p. The paper's cross-tests are: colon 5-FU patients scored with the bladder biomarker
(row "Bladder, paper's amino acid pathway", column "Colon organoids, paper data | 5-FU") and bladder cisplatin patients
scored with the colon biomarker (row "Colon organoids, paper data", column "Bladder organoids, paper data").
Paper: p = 0.16 and 0.86. Blank = pathway not scored in that cohort. Numbers incl. HR and 95% CI:
tables/SuppFig4_cross_test_matrix.tsv.
""")
spec = X.dropna(subset=["p"])
spec = spec.assign(own=spec["own"].fillna(False).astype(bool)); own = spec[spec.own]; off = spec[~spec.own]
summ = dict(n_own=len(own), own_p_below_005=int((own.p < 0.05).sum()), n_off=len(off),
            off_p_below_005=int((off.p < 0.05).sum()), off_share=float((off.p < 0.05).mean()))
pc = X[(X.biomarker == "BLCA_paper_AMINO") & (X.cohort == "COAD_paper_5FU")].p.values
pb = X[(X.biomarker == "COAD_paper_5FU") & (X.cohort == "BLCA_paper_CIS")].p.values
summ.update(paper_SuppFig4a_COAD_with_BLCA_amino=float(pc[0]) if len(pc) and np.isfinite(pc[0]) else None,
            paper_SuppFig4b_BLCA_with_COAD_BH3=float(pb[0]) if len(pb) else None)
json.dump(summ, open(TAB / "SuppFig4_cross_test_summary.json", "w"), indent=1)
print("cross", summ)

# ------------------------------------------------------------------ Supp Fig 8 + 9: target genes and k-step neighbours
G = paper_network()
rows8, rows9 = [], []
for k in KEYS:
    lab, bn, drug, _, _ = REG[k]
    b = bundle(bn); y = b["drugs"][drug]["y"]
    S = [s for s in y.index if s in b["gm"].columns]
    gm, yy = b["gm"][S], y.loc[S].values.astype(float)
    tp = [p for p in b["drugs"][drug]["treated"] if p in b["gp"].columns and p in b["surv"].index]
    gp = b["gp"][tp]
    T = targets(drug)
    for t in T:
        if t in gm.index and t in gp.index and gm.loc[t].std() > 0:
            rho, rp = spearmanr(gm.loc[t].values, yy)
            r, _ = treated_split(b, drug, gp.loc[t] * np.sign(rho if np.isfinite(rho) and rho != 0 else 1))
            rows8.append(dict(key=k, label=lab, drug=drug, gene=t, lab_rho=rho, lab_rho_p=rp, n=r["n"], events=r["events"],
                              p=r["p"], HR=r["HR"], HR_lo=r["lo"], HR_hi=r["hi"]))
        else:
            rows8.append(dict(key=k, label=lab, drug=drug, gene=t, note="not measured in lab models or patients"))
    Tn = [t for t in T if t in G]
    if not Tn:
        continue
    dist = nx.multi_source_dijkstra_path_length(G, Tn, cutoff=3)
    for hop in (1, 2, 3):
        genes = [g for g, d in dist.items() if 1 <= d <= hop and g in gm.index and g in gp.index]
        if len(genes) < 2:
            continue
        Xl = zs(gm.loc[genes]).T.values; Xp = zs(gp.loc[genes]).T.values
        ys = (yy - yy.mean()) / yy.std()
        m = RidgeCV(alphas=np.arange(0.1, 1.0, 0.1), cv=min(3, len(S))).fit(Xl, ys)
        risk = pd.Series(Xp @ m.coef_, index=tp)
        r, _ = treated_split(b, drug, risk)
        rows9.append(dict(key=k, label=lab, drug=drug, hops=hop, n_genes=len(genes), n=r["n"], events=r["events"],
                          p=r["p"], HR=r["HR"], HR_lo=r["lo"], HR_hi=r["hi"]))
T8 = pd.DataFrame(rows8); T8.to_csv(TAB / "SuppFig8_target_genes_alone.tsv", sep="\t", index=False)
T9 = pd.DataFrame(rows9); T9.to_csv(TAB / "SuppFig9_network_neighbourhood_1_2_3_steps.tsv", sep="\t", index=False)
ver = pd.read_csv(TAB / "verify_headline_biomarkers.tsv", sep="\t").set_index("key")
fig, ax = plt.subplots(figsize=(9, 4.6))
xs = np.arange(len(KEYS))
ax.scatter(xs, -np.log10(ver.loc[KEYS, "treated_p"]), marker="*", s=90, color=OI[3], label="rank-1 network pathway (paper method)", zorder=3)
for h, mk, cc in ((1, "o", OI[0]), (2, "s", OI[2]), (3, "^", OI[4])):
    t = T9[T9.hops == h].set_index("key").reindex(KEYS)
    ax.scatter(xs + (h - 2) * 0.15, -np.log10(t.p), marker=mk, color=cc, s=22, label=f"genes within {h} step(s) of targets (ridge)")
t8 = T8.dropna(subset=["p"])
for k, sub in t8.groupby("key"):
    ax.scatter([KEYS.index(k) + 0.3] * len(sub), -np.log10(sub.p), marker="x", color="#555555", s=18)
ax.scatter([], [], marker="x", color="#555555", label="single drug-target gene")
ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.8); ax.text(len(KEYS) - 0.5, -np.log10(0.05) + 0.03, "p = 0.05", fontsize=7, ha="right")
ax.set_xticks(xs); ax.set_xticklabels([REG[k][0] for k in KEYS], rotation=65, ha="right", fontsize=6.5)
ax.set_ylabel("-log10 log-rank p, treated patients"); ax.legend(fontsize=6.5, loc="upper right")
save(fig, "SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets", """
# Supp. Figs. 8 and 9 (extended to every set): simpler network features instead of the proximal pathway

For each lab-model set and drug (x axis) the treated TCGA patients are split at the median of different scores and
their survival compared (y = -log10 log-rank p; higher = stronger split; dashed line = p 0.05).
Red star: the paper's method (rank-1 network-proximal pathway).
Crosses (Supp. Fig. 8): each drug-target gene alone. Its direction comes from its correlation with drug response in the
lab models.
Circles, squares, triangles (Supp. Fig. 9): ridge regression trained on lab models using all genes within 1, 2 or 3
steps of the drug targets in the STRING network (first neighbours = paper Fig. 3d; 2 and 3 steps = Supp. Fig. 9).
Numbers: tables/SuppFig8_target_genes_alone.tsv and tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv.
Paper (bladder): neighbours 0.70 / 0.94 / 0.99 type values; colon: target TYMS alone did not split patients.
""")
print(T9.pivot(index="key", columns="hops", values="p").round(3).to_string())

# ------------------------------------------------------------------ Supp Fig 10: size vs p
R = reactome(); rows10 = []
for k in KEYS:
    lab, bn, drug, _, _ = REG[k]; b = bundle(bn)
    tp = [p for p in b["drugs"][drug]["treated"] if p in b["tsc"].columns and p in b["surv"].index]
    sv = b["surv"].loc[tp]
    for pw in b["tsc"].index:
        s = b["tsc"].loc[pw, tp].values.astype(float); hi = s > np.median(s)
        if 0 < hi.sum() < len(hi):
            rows10.append(dict(key=k, pathway=pw, size=len(R.get(pw, [])), p=SV.logrank(sv["T"].values[hi], sv["E"].values[hi],
                                                                                      sv["T"].values[~hi], sv["E"].values[~hi])))
T10 = pd.DataFrame(rows10); T10.to_csv(DATA / "SuppFig10_size_vs_p_all.tsv.gz", sep="\t", index=False)
s10 = []
nc = 4; nr = int(np.ceil(len(KEYS) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(10, 2.3 * nr))
for i, k in enumerate(KEYS):
    t = T10[(T10.key == k) & (T10["size"] > 0)]
    lr = linregress(t["size"], -np.log10(t.p))
    pw = bio[k][0]; rp = t[t.pathway == pw]
    s10.append(dict(key=k, label=REG[k][0], n_pathways=len(t), R2=lr.rvalue ** 2, slope_p=lr.pvalue,
                    share_p_below_005=float((t.p < 0.05).mean()), rank1_pathway_p_any_direction=float(rp.p.iloc[0]) if len(rp) else np.nan))
    ax = axs.flat[i]; ax.scatter(t["size"], -np.log10(t.p), s=4, color="#888888", alpha=0.6)
    if len(rp):
        ax.scatter(rp["size"], -np.log10(rp.p), s=30, color=OI[3], marker="*")
    ax.set_xscale("log"); ax.set_title(f"{REG[k][0]}\nR2 = {lr.rvalue**2:.3f}", fontsize=6.5); ax.tick_params(labelsize=6)
    letter(ax, "abcdefghijklmnopqrstuvwxyz"[i])
for j in range(i + 1, nr * nc):
    axs.flat[j].axis("off")
fig.supxlabel("pathway size (genes, log scale)", fontsize=8); fig.supylabel("-log10 log-rank p (treated)", fontsize=8)
fig.tight_layout()
pd.DataFrame(s10).to_csv(TAB / "SuppFig10_size_vs_p_summary_all_sets.tsv", sep="\t", index=False)
save(fig, "SuppFig10_size_vs_p_all_sets", """
# Supp. Fig. 10 (extended to every set): does pathway size drive the survival split?

Each grey dot is one Reactome pathway. x = number of genes in the pathway; y = how well that pathway alone splits the
drug-treated TCGA patients (-log10 log-rank p, median split). Red star = the set's rank-1 (biomarker) pathway.
R2 = share of the spread in y explained by pathway size (paper: R2 = 0.0026 for colon and 0.0008 for bladder).
Small R2 = size does not drive the result. Numbers: tables/SuppFig10_size_vs_p_summary_all_sets.tsv.
""")
print(pd.DataFrame(s10)[["key", "R2", "share_p_below_005"]].round(3).to_string())

# ------------------------------------------------------------------ Supp Fig 16: proximity of all drugs
sys.path.insert(0, str(REPO / "robust"))
from obd import network as N
Pz = N.precomputed_proximity()
Pz.to_csv(TAB / "SuppFig16_paper_proximity_44_columns.tsv.gz", sep="\t")
cols = dict()
for k in KEYS:
    b = bundle(REG[k][1])
    for d, v in b["drugs"].items():
        if v["z"] is not None and d not in Pz.columns and d not in cols:
            cols[d] = v["z"]
A = pd.concat([Pz, pd.DataFrame(cols)], axis=1)
A.to_csv(TAB / "SuppFig16_proximity_all_drugs.tsv.gz", sep="\t")
keep = A.index[(A <= CUTOFF).sum(1) >= 3]
M = A.loc[keep]
import scipy.cluster.hierarchy as hc
o1 = hc.leaves_list(hc.linkage(M.fillna(0).values, "average")); o2 = hc.leaves_list(hc.linkage(M.fillna(0).values.T, "average"))
M = M.iloc[o1, o2]
fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(M.values, cmap="PuOr", vmin=-4, vmax=4, aspect="auto")
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels([c.lower()[:22] + (" *" if c in cols else "") for c in M.columns], rotation=80, fontsize=6)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels([short(p, 45) for p in M.index], fontsize=4.5)
cb = fig.colorbar(im, ax=ax, shrink=0.4); cb.set_label("network proximity z (negative = closer)")
save(fig, "SuppFig16_proximity_all_drugs", f"""
# Supp. Fig. 16 (paper's 44 drug columns + drugs added for the other sets)

Network proximity z-score between each drug's targets and each Reactome pathway in the STRING network.
Purple / negative = the pathway is closer to the targets than random gene sets; a pathway with z <= -1.2816 is
"proximal". Columns = {Pz.shape[1]} columns of the paper's own table (MOESM3 'Figure 2a,d'; some merge drugs with the same
targets, e.g. CETUXIMAB_GEFITINIB) plus {len(cols)} drugs recomputed on the paper's network for the other sets (marked *).
Rows = the {M.shape[0]} pathways that are proximal for at least 3 drugs (all {A.shape[0]} pathways are in
tables/SuppFig16_proximity_all_drugs.tsv.gz). Rows and columns are clustered, so drugs with similar target
neighbourhoods sit next to each other (e.g. the platinum drugs, the taxanes).
""")
print("prox", Pz.shape, len(cols), M.shape)
