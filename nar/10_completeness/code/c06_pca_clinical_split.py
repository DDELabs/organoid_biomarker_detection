"""Supp Fig 1 (PCA explained variance, lab models vs TCGA, every set), Supp Fig 14 (GDSC/CMP cell lines vs TCGA,
with explained-variance numbers), Supp Tables 1-2 (clinical balance of predicted groups, every headline set) and
Supp Fig 17 (repeated 60/10/30 split inside the lab models, every headline set)."""
from common import *
from scipy.stats import chi2_contingency, spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV

def pca2(A, B, scale):
    """A, B: features x samples. Returns PC1/PC2 % variance and coordinates."""
    common = A.index.intersection(B.index)
    A, B = A.loc[common], B.loc[common]
    if scale:
        f = lambda d: d.sub(d.mean(1), axis=0).div(d.std(1).replace(0, np.nan), axis=0)
        A, B = f(A), f(B)
    X = pd.concat([A, B], axis=1).T.fillna(0)
    X = X.loc[:, X.std() > 0]
    p = PCA(2, random_state=0).fit(X.values)
    return 100 * p.explained_variance_ratio_, p.transform(X.values), A.shape[1], len(common)

# ---------------------------------------- Supp Fig 1 variance, every bundle
rows = []
for bn in sorted({REG[k][1] for k in REG}):
    b = bundle(bn)
    for lvl, A, B in (("pathways", b["osc"], b["tsc"]), ("genes", np.log2(b["gm"].clip(lower=0) + 1) if b["gm"].max().max() > 50 else b["gm"], b["gp"])):
        for sc in (False, True):
            v, _, na, nf = pca2(A, B, sc)
            rows.append(dict(bundle=bn, label=b["meta"]["label"], level=lvl, scaled=sc, n_lab=na, n_patients=B.shape[1], n_features=nf, pc1=v[0], pc2=v[1]))
P1 = pd.DataFrame(rows); P1.to_csv(TAB / "SuppFig1_PCA_variance_all_sets.tsv", sep="\t", index=False)
print(P1[P1.level == "pathways"][["bundle", "scaled", "pc1", "pc2"]].round(1).to_string())

# ---------------------------------------- Supp Fig 14: cell lines vs TCGA
E = pd.read_pickle(NAR / "_cache/depmap/derived/expr_cmp.pkl"); Mo = pd.read_pickle(NAR / "_cache/depmap/derived/models.pkl")
Mo = Mo[Mo.model_type == "Cell Line"] if "model_type" in Mo else Mo
CT = {"COAD_paper": "Colorectal Carcinoma", "BLCA_paper": "Bladder Carcinoma", "PAAD_tiriac": "Pancreatic Carcinoma",
      "LIHC_ji": "Hepatocellular Carcinoma", "OV_vias": "Ovarian Carcinoma", "LUAD_cells": "Non-Small Cell Lung Carcinoma",
      "PDX_CM": "Melanoma", "PDX_BRCA": "Breast Carcinoma"}
R = reactome()
def path_scores(G):  # mean of per-gene z-scores within the dataset (simple pathway score, same for both sides)
    Z = G.sub(G.mean(1), axis=0).div(G.std(1).replace(0, np.nan), axis=0)
    out = {p: Z.loc[[g for g in gs if g in Z.index]].mean() for p, gs in R.items() if sum(g in Z.index for g in gs) >= 5}
    return pd.DataFrame(out).T
rows = []; fig, axs = plt.subplots(len(CT), 4, figsize=(10, 2.2 * len(CT)))
for i, (bn, ct) in enumerate(CT.items()):
    b = bundle(bn); ids = [m for m in Mo.index[Mo.cancer_type == ct] if m in E.columns]
    cl = np.log2(E[ids].clip(lower=0) + 1) if E[ids].max().max() > 50 else E[ids]
    pat = b["gp"]
    for j, (lvl, sc) in enumerate((("genes", False), ("genes", True), ("pathways", False), ("pathways", True))):
        A, B = (cl, pat) if lvl == "genes" else (None, None)
        if lvl == "pathways":
            common = cl.index.intersection(pat.index)
            # raw pathway level = mean expression of the set's genes (keeps the batch difference); scaled = per-dataset z
            mean_expr = lambda G: pd.DataFrame({p: G.loc[[g for g in gs if g in G.index]].mean() for p, gs in R.items()
                                                if sum(g in G.index for g in gs) >= 5}).T
            A, B = mean_expr(cl.loc[common]), mean_expr(pat.loc[common])
        v, Xc, na, nf = pca2(A, B, sc)
        rows.append(dict(cancer=ct, tcga=b["meta"]["project"], level=lvl, scaled=sc, n_cell_lines=na, n_patients=B.shape[1], n_features=nf, pc1=v[0], pc2=v[1]))
        ax = axs[i, j]
        ax.scatter(Xc[na:, 0], Xc[na:, 1], s=3, color=OI[0], alpha=0.5, label="TCGA tumours")
        ax.scatter(Xc[:na, 0], Xc[:na, 1], s=6, color=OI[3], label="cell lines")
        ax.set_xlabel(f"PC1 {v[0]:.1f}%", fontsize=7); ax.set_ylabel(f"PC2 {v[1]:.1f}%", fontsize=7); ax.tick_params(labelsize=6)
        ax.set_title(f"{b['meta']['project']}: {lvl}, {'scaled' if sc else 'raw'}", fontsize=7.5)
        if i == 0 and j == 0:
            ax.legend(fontsize=6, markerscale=2)
        letter(ax, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"[4 * i + j]) if i < 1 else None
fig.tight_layout()
P14 = pd.DataFrame(rows); P14.to_csv(TAB / "SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv", sep="\t", index=False)
save(fig, "SuppFig14_cell_lines_vs_TCGA_PCA", """
# Supp. Fig. 14 (with the explained-variance numbers; extended to 8 cancers)

Each row is one cancer: Cell Model Passports cancer cell lines of that type (red; RNA-seq, log2 TPM) and the TCGA
tumours of the matching cohort (blue). Columns: gene level raw, gene level after scaling each gene within each data
set (z-score), pathway level raw (mean expression of the Reactome set's genes), pathway level scaled. Axis labels give
the % of variance explained by PC1 and PC2 (the numbers the paper reports). Raw: cell lines and tumours form two
separate clouds (a technical and biological batch). Scaled: the clouds overlap, which is why the method scales each
data set before training. Numbers: tables/SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv. The paper used ssGSEA NES for
the pathway level; here a simpler mean-expression score is used for the cell lines because no ssGSEA run of all
cell lines exists in the project (stated plainly; the gene-level panels are identical in kind to the paper's).
""")
print(P14.round(1).to_string())

# ---------------------------------------- Supp Tables 1-2, every headline set
cdr = pd.read_csv(NAR / "09_registered_test/_data/tcga_response/survival_cdr.tsv.gz", sep="\t").set_index("patient")
def stage(s):
    s = str(s)
    return np.nan if "Stage" not in s else ("I-II" if s.replace("Stage ", "").startswith(("I ", "IA", "IB", "IC", "II")) and not s.replace("Stage ", "").startswith(("III", "IV")) else "III-IV")
rows = []
for k, (lab, bn, drug, _, _) in REG.items():
    if k == "LIHC_broutier_SOR":
        continue
    b = bundle(bn); pw, c = rank1(k)
    r, hi = treated_split(b, drug, b["tsc"].loc[pw] * np.sign(c))
    pats = [p for p in b["drugs"][drug]["treated"] if p in b["tsc"].columns and p in b["surv"].index]
    grp = pd.Series(np.where(hi, "non_responder", "responder"), index=pats)
    cl = cdr.reindex(pats)
    feats = {"age": np.where(cl.age >= 60, ">=60", "<60"), "sex": cl.gender.values, "stage": cl.stage.map(stage).values}
    for f, v in feats.items():
        d = pd.DataFrame({"g": grp.values, "v": v}).dropna()
        d = d[d.v.astype(str) != "nan"]
        tab = pd.crosstab(d.v, d.g)
        p = chi2_contingency(tab)[1] if tab.shape[0] > 1 and tab.shape[1] > 1 else np.nan
        for lev in tab.index:
            rows.append(dict(key=k, label=lab, feature=f, level=lev, responder=int(tab.loc[lev].get("responder", 0)),
                             non_responder=int(tab.loc[lev].get("non_responder", 0)), chi2_p=p))
T1 = pd.DataFrame(rows); T1.to_csv(TAB / "SuppTable1_2_clinical_balance_all_sets.tsv", sep="\t", index=False)
print("clinical p<0.05:", T1.drop_duplicates(["key", "feature"]).query("chi2_p < 0.05")[["key", "feature", "chi2_p"]].to_string())

# ---------------------------------------- Supp Fig 17: 60/10/30 split inside lab models (validation 10% unused, as paper)
rng = np.random.default_rng(SEED); rows = []
for k, (lab, bn, drug, _, _) in REG.items():
    if k == "LIHC_broutier_SOR":
        continue
    b = bundle(bn); y = b["drugs"][drug]["y"]; S = [s for s in y.index if s in b["osc"].columns]
    z = b["drugs"][drug]["z"]; prox = [p for p in z.index[z <= CUTOFF] if p in b["osc"].index]
    X = b["osc"].loc[prox, S].T.values; yy = y.loc[S].values.astype(float); n = len(S)
    if n < 8:
        continue
    for rep in range(100):
        o = rng.permutation(n); ntr = max(3, int(round(0.6 * n))); nte = max(2, int(round(0.3 * n)))
        tr, te = o[:ntr], o[n - nte:]
        m = RidgeCV(alphas=np.arange(0.1, 1.0, 0.1), cv=min(3, ntr)).fit(X[tr], yy[tr])
        rho = spearmanr(m.predict(X[te]), yy[te])[0] if np.std(yy[te]) > 0 else np.nan
        rows.append(dict(key=k, label=lab, n_models=n, n_train=ntr, n_test=nte, rep=rep, spearman=rho))
T17 = pd.DataFrame(rows); T17.to_csv(DATA / "SuppFig17_split_603010_reps.tsv.gz", sep="\t", index=False)
s17 = T17.groupby(["key", "label", "n_models", "n_test"]).spearman.agg(["median", lambda x: float((x > 0).mean())]).reset_index()
s17.columns = ["key", "label", "n_models", "n_test", "median_spearman", "share_positive"]
s17.to_csv(TAB / "SuppFig17_split_603010_summary.tsv", sep="\t", index=False); print(s17.round(2).to_string())
fig, ax = plt.subplots(figsize=(9, 3.8))
ks = list(s17.key)
ax.boxplot([T17[T17.key == k].spearman.dropna() for k in ks], tick_labels=[REG[k][0] for k in ks])
ax.axhline(0, color="grey", ls="--", lw=0.8); ax.set_ylabel("Spearman rho, predicted vs observed\n(30% held-out lab models)")
plt.setp(ax.get_xticklabels(), rotation=65, ha="right", fontsize=6.5)
save(fig, "SuppFig17_split_603010_all_sets", """
# Supp. Fig. 17 (paper's 60/10/30 split, every headline set)

Inside each set of lab models, the models were split at random into 60% training, 10% validation (not used, as in
the paper's ridge set-up with built-in cross-validation) and 30% test, 100 times. A ridge model on the network-proximal
pathways was trained on the 60% and its predictions on the 30% compared with the measured drug response (Spearman rho).
Boxes show the 100 repeats; above 0 = the model ranks unseen lab models in the right order. With 9-19 models the test
part holds only 3-6 models, so single repeats vary a lot. Numbers: tables/SuppFig17_split_603010_summary.tsv
(leave-one-out results of the organoid-only agents give the same kind of information with all models).
""")
