"""Fig 5 equivalents (mutation concordance) for other drugs/cancers, plus Supp Fig 12 (comparison feature sets for
cetuximab vs BRAF V600E) and Supp Fig 18 (leave-k-out organoids) for the paper's colon cetuximab case.
Test: one-sided Mann-Whitney on the predicted-resistance score (coef x pathway score) in ALL TCGA tumours with
mutation calls (not only treated patients), as in the paper."""
from common import *
from scipy.stats import mannwhitneyu
from sklearn.linear_model import RidgeCV
import networkx as nx

XR = DATA / "raw" / "xena_mut"; ORAW = O2 / "BLCA_OV_LUAD_OTHER" / "_data" / "raw"
SILENT = {"synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "upstream_gene_variant",
          "downstream_gene_variant", "intergenic_region", "non_coding_transcript_exon_variant", "splice_region_variant"}

def mut(projects):
    m = pd.concat([pd.read_csv(f, sep="\t", usecols=["sample", "gene", "Amino_Acid_Change", "effect"]) for f in projects])
    m = m[m["sample"].str[13:15] == "01"].copy(); m["patient"] = m["sample"].str[:12]
    return m

def carriers(m, gene, rule):
    s = m[(m.gene == gene) & ~m.effect.isin(SILENT)]
    return set(s.loc[s.Amino_Acid_Change.fillna("").map(rule), "patient"]), set(m.patient)

def test(score, mutset, prof, expect):  # expect 'higher' or 'lower' score in mutants (score = predicted resistance)
    pats = [p for p in score.index if p in prof]
    a = score.loc[[p for p in pats if p in mutset]].values; b = score.loc[[p for p in pats if p not in mutset]].values
    if len(a) < 3:
        return dict(n_mut=len(a), n_wt=len(b), p_one_sided=np.nan)
    p = mannwhitneyu(a, b, alternative="greater" if expect == "higher" else "less").pvalue
    return dict(n_mut=len(a), n_wt=len(b), median_mut=float(np.median(a)), median_wt=float(np.median(b)), p_one_sided=p)

crc = mut([XR / "TCGA-COAD.somaticmutation_wxs.tsv.gz", XR / "TCGA-READ.somaticmutation_wxs.tsv.gz"])
lua = mut([XR / "TCGA-LUAD.somaticmutation_wxs.tsv.gz"])
ov = mut([ORAW / "TCGA-OV.somaticmutation_wxs.tsv.gz"])
KRAS = carriers(crc, "KRAS", lambda a: any(a.startswith(f"p.{x}") for x in ("G12", "G13", "Q61", "A146")))
BRAF = carriers(crc, "BRAF", lambda a: a in ("p.V600E", "p.V640E"))  # Xena uses a transcript where V600E reads V640E
EGFR = carriers(lua, "EGFR", lambda a: a in ("p.L858R", "p.L861Q") or a.startswith("p.G719") or
                ("del" in a and any(f"{x}" in a for x in ("E746", "L747", "E745"))))
brca = set().union(*[carriers(ov, g, lambda a: True)[0] for g in ("BRCA1", "BRCA2")]); BRCA = (brca, set(ov.patient))

def score_of(key, bn=None):
    pw, c = rank1(key); b = bundle(bn or REG[key][1]); return b["tsc"].loc[pw] * c, pw

cet_vdw = pd.read_csv(O2 / "COAD_vandewetering/tables/coefficients_COAD_CETUXIMAB.tsv", sep="\t", index_col=0)["Ridge"]
cet_vdw = cet_vdw.reindex(cet_vdw.abs().sort_values(ascending=False).index)
jobs = [
    ("Colorectal organoids (vdW), current | cetuximab", bundle("COADREAD_vdW")["tsc"].loc[cet_vdw.index[0]] * cet_vdw.iloc[0], cet_vdw.index[0], "BRAF V600E", BRAF, "higher"),
    ("Colorectal organoids (vdW), current | cetuximab", bundle("COADREAD_vdW")["tsc"].loc[cet_vdw.index[0]] * cet_vdw.iloc[0], cet_vdw.index[0], "KRAS G12/G13/Q61/A146", KRAS, "higher"),
    ("Colorectal PDX (Isella) | cetuximab",) + score_of("COAD_isella_CET") + ("KRAS G12/G13/Q61/A146", KRAS, "higher"),
    ("Colorectal PDX (Isella) | cetuximab",) + score_of("COAD_isella_CET") + ("BRAF V600E", BRAF, "higher"),
    ("Lung PDX (Gao) | erlotinib",) + score_of("PDX_NSCLC_ERL") + ("EGFR activating", EGFR, "lower"),
]
pw, c = pd.read_csv(O2 / "LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_coefficients.tsv", sep="\t", index_col=0).pipe(
    lambda d: (d["coef_Ridge"].abs().idxmax(), d.loc[d["coef_Ridge"].abs().idxmax(), "coef_Ridge"]))
jobs.append(("Lung adenocarcinoma cell lines | erlotinib", bundle("LUAD_cells")["tsc"].loc[pw] * c, pw, "EGFR activating", EGFR, "lower"))
for k, lab in (("OV_vias_CARBO", "Ovarian organoids (Vias) | carboplatin"),):
    jobs.append((lab,) + score_of(k) + ("BRCA1/2 somatic", BRCA, "lower"))
pw, c = pd.read_csv(O2 / "OV_vias2023/tables/OV_CISPLATIN_coefficients.tsv", sep="\t", index_col=0).pipe(
    lambda d: (d["coef_Ridge"].abs().idxmax(), d.loc[d["coef_Ridge"].abs().idxmax(), "coef_Ridge"]))
jobs.append(("Ovarian organoids (Vias) | cisplatin", bundle("OV_vias")["tsc"].loc[pw] * c, pw, "BRCA1/2 somatic", BRCA, "lower"))
rows = []
for lab, sc, pw, mname, (ms, prof), exp in jobs:
    rows.append(dict(model=lab, pathway=pw, mutation=mname, expected_in_mutants=("more resistant" if exp == "higher" else "more sensitive"), **test(sc, ms, prof, exp)))
M = pd.DataFrame(rows); M.to_csv(TAB / "Fig5_mutation_concordance_other_cancers.tsv", sep="\t", index=False)
print(M[["model", "mutation", "n_mut", "n_wt", "p_one_sided"]].to_string())

# ---------------------------------------------- Supp Fig 12 + 18: colon cetuximab vs BRAF V600E, feature sets and leave-k-out
b = bundle("COADREAD_vdW"); y = b["drugs"]["CETUXIMAB"]["y"]; S = list(y.index); yy = y.values.astype(float)
z = b["drugs"]["CETUXIMAB"]["z"]; prox = [p for p in z.index[z <= CUTOFF] if p in b["osc"].index and p in b["tsc"].index]
def ridge_rank1(cols, X, Y):
    m = RidgeCV(alphas=np.arange(0.1, 1.0, 0.1), cv=3).fit(X, Y); i = np.argmax(np.abs(m.coef_)); return cols[i], m.coef_[i]
p1, c1 = ridge_rank1(prox, b["osc"].loc[prox, S].T.values, yy)
zs = lambda d: d.sub(d.mean(1), axis=0).div(d.std(1).replace(0, np.nan), axis=0).fillna(0)
def ridge_score(Fl, Fp):
    ys = (yy - yy.mean()) / yy.std()
    m = RidgeCV(alphas=np.arange(0.1, 1.0, 0.1), cv=3).fit(zs(Fl).T.values, ys); return pd.Series(zs(Fp).T.values @ m.coef_, index=Fp.columns)
G = paper_network(); dist = nx.multi_source_dijkstra_path_length(G, [t for t in targets("CETUXIMAB") if t in G], cutoff=3)
gm, gp = b["gm"][S], b["gp"]
feat = {"network rank-1 pathway (refit here)": b["tsc"].loc[p1] * c1,
        "whole transcriptome": ridge_score(gm, gp.loc[gm.index]),
        "all pathways": ridge_score(b["osc"][S], b["tsc"].loc[b["osc"].index])}
for h in (1, 2, 3):
    gg = [g for g, d in dist.items() if 1 <= d <= h and g in gm.index and g in gp.index]
    feat[f"neighbours up to {h} step(s) ({len(gg)} genes)"] = ridge_score(gm.loc[gg], gp.loc[gg])
r12 = [dict(feature_set=k, **test(v, BRAF[0], BRAF[1], "higher")) for k, v in feat.items()]
T12 = pd.DataFrame(r12); T12.to_csv(TAB / "SuppFig12_cetuximab_BRAF_feature_sets.tsv", sep="\t", index=False)
print(p1, c1); print(T12.to_string())
rng = np.random.default_rng(SEED); r18 = []
for k in (1, 2, 3, 4):
    for rep in range(100):
        keep = sorted(rng.choice(len(S), len(S) - k, replace=False))
        pw, cc = ridge_rank1(prox, b["osc"].loc[prox, [S[i] for i in keep]].T.values, yy[keep])
        r18.append(dict(k=k, rep=rep, pathway=pw, **test(b["tsc"].loc[pw] * cc, BRAF[0], BRAF[1], "higher")))
T18 = pd.DataFrame(r18); T18.to_csv(TAB / "SuppFig18_leave_k_out_cetuximab_BRAF.tsv", sep="\t", index=False)
s18 = T18.groupby("k").agg(share_p_below_005=("p_one_sided", lambda x: float((x < 0.05).mean())), median_p=("p_one_sided", "median"),
                           share_same_rank1=("pathway", lambda x: float((x == p1).mean()))).reset_index()
s18.to_csv(TAB / "SuppFig18_leave_k_out_summary.tsv", sep="\t", index=False); print(s18.to_string())

fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), gridspec_kw=dict(width_ratios=[1.3, 1, 1], wspace=0.9))
ax = axs[0]; lbl = [f"{r.model}\n{r.mutation}" for r in M.itertuples()]
ax.barh(range(len(M)), -np.log10(M.p_one_sided), color=[OI[0] if "cetux" in m else OI[2] if "erlot" in m else OI[4] for m in M.model])
for i, r in enumerate(M.itertuples()):
    ax.text(-np.log10(r.p_one_sided) + 0.05, i, f"p={r.p_one_sided:.2g} ({r.n_mut} vs {r.n_wt})", va="center", fontsize=6)
ax.set_yticks(range(len(M))); ax.set_yticklabels(lbl, fontsize=6); ax.invert_yaxis(); ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=0.8)
ax.set_xlabel("-log10 one-sided p (mutants vs wild type)"); ax.set_title("Mutation concordance, other drugs", fontsize=8.5); letter(ax, "a")
ax = axs[1]; ax.barh(range(len(T12)), -np.log10(T12.p_one_sided), color=OI[1])
ax.set_yticks(range(len(T12))); ax.set_yticklabels(T12.feature_set, fontsize=6.5); ax.invert_yaxis()
ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=0.8); ax.set_xlabel("-log10 one-sided p"); ax.set_title("Colon cetuximab vs BRAF V600E:\nfeature sets (Supp. Fig. 12)", fontsize=8.5); letter(ax, "b")
ax = axs[2]; ax.boxplot([-np.log10(T18[T18.k == k].p_one_sided) for k in (1, 2, 3, 4)], tick_labels=["1", "2", "3", "4"])
ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.8); ax.set_xlabel("organoids left out (k), 100 draws each")
ax.set_ylabel("-log10 one-sided p"); ax.set_title("Leave-k-out organoids (Supp. Fig. 18)", fontsize=8.5); letter(ax, "c")
fig.tight_layout()
save(fig, "Fig5_SuppFig12_SuppFig18_mutation_concordance", f"""
# Fig. 5 equivalents for other cancers, Supp. Fig. 12 and Supp. Fig. 18

The idea (paper Fig. 5): a drug model should call tumours with a known resistance (or sensitivity) mutation more
resistant (or more sensitive), even without survival data. Score = rank-1 pathway x its Ridge weight (higher =
predicted more resistant), in all TCGA tumours with open somatic mutation calls (Xena GDC hub). One-sided
Mann-Whitney test in the expected direction. Dashed line = p 0.05.
a. Mutation pairs with an established clinical meaning: cetuximab vs KRAS (codons 12, 13, 61, 146) and BRAF V600E in
colon + rectum (mutants expected more resistant); erlotinib vs activating EGFR mutations (L858R, L861Q, G719X, exon 19
deletions) in lung adenocarcinoma (mutants expected more sensitive); platinum vs somatic BRCA1/2 mutations in ovarian
cancer (expected more sensitive; germline BRCA status is controlled access, so only somatic mutations are used).
No widely accepted mutation marker exists for gemcitabine (pancreas, sarcoma), sorafenib (liver) or dacarbazine
(melanoma), so no such test is made there. Bladder ERCC2: nar/02_organoid_only/BLCA_lee2018/extras.
b. Supp. Fig. 12: colon cetuximab model trained on the 19 van de Wetering organoids (current data) with other feature
sets, tested on BRAF V600E. The network rank-1 pathway refit here is {short(p1)} (same as the organoid-only agent).
c. Supp. Fig. 18: the network model refit after leaving 1-4 organoids out (100 random draws each); boxes show the
spread of the BRAF V600E p. Share with p < 0.05 and share keeping the same rank-1 pathway:
tables/SuppFig18_leave_k_out_summary.tsv. Numbers: tables/Fig5_mutation_concordance_other_cancers.tsv,
tables/SuppFig12_cetuximab_BRAF_feature_sets.tsv.
""")
