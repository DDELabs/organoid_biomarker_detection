"""Fig 4b/c equivalents for other cancers: biomarker pathways scored (ssGSEA, as the paper) in public isogenic
drug-sensitive vs drug-resistant cell-line pairs. GSE15709: ovarian A2780 parental vs cisplatin-resistant Round5
(Affymetrix U133 Plus 2, 5 vs 5). GSE322742: liver HepG2 and Huh7 parental vs sorafenib-resistant (RNA-seq FPKM, 2 vs 2
per line; resistant cells are kept in sorafenib)."""
from common import *
from scipy.stats import ttest_ind
from obd.scoring import ssgsea
from obd.reference import ensembl_to_symbol
from obd.geo_matrix import read_expr

RAW = DATA / "raw"
# ovarian A2780
ex = read_expr(RAW / "GSE15709/GSE15709-GPL570_series_matrix.txt.gz")
_hdr = next(i for i, l in enumerate(open(RAW / "GPL570/GPL570_full.txt")) if l.startswith("ID\t"))
ann = pd.read_csv(RAW / "GPL570/GPL570_full.txt", sep="\t", skiprows=_hdr, usecols=["ID", "Gene Symbol"], dtype=str, on_bad_lines="skip").dropna()
ann = ann[~ann["Gene Symbol"].str.contains("///")].set_index("ID")["Gene Symbol"]
ex = ex.loc[ex.index.intersection(ann.index)]
ex = ex.groupby(ann.loc[ex.index].values).mean()
if ex.max().max() > 50:
    ex = np.log2(ex.clip(lower=0) + 1)
import gzip
tl = [l for l in gzip.open(RAW / "GSE15709/GSE15709-GPL570_series_matrix.txt.gz", "rt") if l.startswith("!Sample_title")][0]
tl = [t.strip('"\n') for t in tl.split("\t")[1:]]
grp_ov = pd.Series(["resistant" if "Round5" in t else "sensitive" for t in tl], index=ex.columns)
# liver sorafenib
e2s = ensembl_to_symbol(); fr = []
for f in sorted((RAW / "GSE322742").glob("GSM*_FPKM_values.csv.gz")):
    d = pd.read_csv(f, index_col=0).iloc[:, 0]; d.name = f.name.split("_")[0]; fr.append(d)
lv = pd.concat(fr, axis=1); lv = lv.groupby(lv.index.map(lambda x: e2s.get(x, np.nan))).sum()
lv = np.log2(lv + 1)
tl2 = [l for l in gzip.open(RAW / "GSE322742/GSE322742_series_matrix.txt.gz", "rt") if l.startswith("!Sample_title")][0]
tl2 = [t.strip('"\n') for t in tl2.split("\t")[1:]]
gsm2 = [l for l in gzip.open(RAW / "GSE322742/GSE322742_series_matrix.txt.gz", "rt") if l.startswith("!Sample_geo_accession")][0]
gsm2 = [t.strip('"\n') for t in gsm2.split("\t")[1:]]
lab2 = dict(zip(gsm2, tl2)); lv = lv[[g for g in gsm2 if g in lv.columns]]

R = reactome()
def coef_of(f, pw=None):
    c = pd.read_csv(O2 / f, sep="\t", index_col=0); col = "Ridge" if "Ridge" in c.columns else "coef_Ridge"
    s = c[col]; pw = pw or s.abs().idxmax(); return pw, float(s.loc[pw])
BIOS = {"cisplatin": [("Bladder organoids, paper data (rank 1)",) + coef_of("BLCA_lee2018/paper/tables/BLCA_CISPLATIN_coefficients.tsv"),
                      ("Bladder, paper's amino acid pathway", "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION", -0.14713652465773913),
                      ("Bladder organoids, current data (rank 1)",) + coef_of("BLCA_lee2018/current/tables/BLCA_CISPLATIN_coefficients.tsv"),
                      ("Ovarian organoids, carboplatin (rank 1)",) + coef_of("OV_vias2023/tables/OV_CARBOPLATIN_coefficients.tsv"),
                      ("Ovarian organoids, cisplatin (rank 1)",) + coef_of("OV_vias2023/tables/OV_CISPLATIN_coefficients.tsv"),
                      ("Lung cell lines, cisplatin (rank 1)",) + coef_of("LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_coefficients.tsv")],
        "sorafenib": [("Liver organoids (Ji) (rank 1)",) + coef_of("LIHC_ji/tables/coefficients_LIHC_SORAFENIB.tsv"),
                      ("Liver organoids (Broutier) (rank 1)",) + coef_of("LIHC_broutier/tables/coefficients_LIHC_SORAFENIB.tsv")]}
rows = []
def run(expr, groups, dataset, drug):
    sets = {pw: R[pw] for _, pw, _ in BIOS[drug] if pw in R}
    sc = ssgsea(expr, sets, min_size=2, threads=1)
    for name, pw, c in BIOS[drug]:
        if pw not in sc.index:
            rows.append(dict(dataset=dataset, biomarker=name, pathway=pw, note="too few genes measured")); continue
        s, r = sc.loc[pw, groups.index[groups == "sensitive"]], sc.loc[pw, groups.index[groups == "resistant"]]
        t, p = ttest_ind(s, r, equal_var=False)
        exp = "resistant higher" if c > 0 else "resistant lower"
        ok = (r.mean() > s.mean()) == (c > 0)
        rows.append(dict(dataset=dataset, drug=drug, biomarker=name, pathway=pw, organoid_coef=c, expected=exp, n_sensitive=len(s),
                         n_resistant=len(r), mean_sensitive=s.mean(), mean_resistant=r.mean(), welch_t=t, p=p, direction_as_predicted=ok))
run(ex, grp_ov, "GSE15709 ovarian A2780 vs cisplatin-resistant Round5", "cisplatin")
for line in ("HepG2", "Huh7"):
    cols = [g for g in lv.columns if line in lab2[g]]
    grp = pd.Series(["resistant" if "SR" in lab2[g] else "sensitive" for g in cols], index=cols)
    run(lv[cols], grp, f"GSE322742 liver {line} parental vs sorafenib-resistant", "sorafenib")
T = pd.DataFrame(rows); T.to_csv(TAB / "Fig4bc_isogenic_other_cancers.tsv", sep="\t", index=False)
print(grp_ov.value_counts().to_dict()); print(T[["dataset", "biomarker", "n_sensitive", "n_resistant", "p", "direction_as_predicted"]].to_string())
fig, ax = plt.subplots(figsize=(8, 3.8))
T2 = T.dropna(subset=["p"]).reset_index(drop=True)
cols = [OI[2] if d else OI[3] for d in T2.direction_as_predicted]
ax.barh(range(len(T2)), -np.log10(T2.p), color=cols)
for i, r in T2.iterrows():
    ax.text(-np.log10(r.p) + 0.05, i, f"p={r.p:.2g}", va="center", fontsize=6.5)
ax.set_yticks(range(len(T2))); ax.set_yticklabels([f"{r.dataset.split(' vs')[0].replace('GSE15709 ', '').replace('GSE322742 ', '')}: {r.biomarker}" for r in T2.itertuples()], fontsize=6.5)
ax.invert_yaxis(); ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=0.8)
ax.set_xlabel("-log10 Welch t-test p, sensitive vs resistant cells")
ax.bar(0, 0, color=OI[2], label="difference in the predicted direction"); ax.bar(0, 0, color=OI[3], label="opposite direction"); ax.legend(fontsize=7, loc="lower right")
save(fig, "Fig4bc_isogenic_other_cancers", """
# Fig. 4b/c equivalents for other cancers: biomarker pathways in drug-resistant cell lines

Paper Fig. 4b/c compared the biomarker pathway between a drug-sensitive parental cell line and a drug-resistant copy of
the same line. We found two public pairs that match drugs of our other sets:
- GSE15709: ovarian cancer line A2780, parental (cisplatin-sensitive, 5 replicates) vs Round5 cisplatin-resistant
  (5 replicates), Affymetrix arrays. Tested: the platinum biomarkers of bladder, ovarian and lung sets
  (carboplatin uses the cisplatin network, as in the main analyses).
- GSE322742: liver cancer lines HepG2 and Huh7, parental vs sorafenib-resistant (2 vs 2 each, RNA-seq). Tested: the
  sorafenib biomarkers of the Ji and Broutier liver organoids. Note: the resistant cells were kept in sorafenib, so drug
  exposure and resistance are mixed.
Pathways were scored with ssGSEA (as in the paper). Bars: -log10 Welch t-test p; green = the difference goes the way
the organoid Ridge weight predicts (positive weight -> higher in resistant cells), red = opposite. Dashed line = p 0.05.
With 2 vs 2 samples a t-test has little power. No public isogenic pair with expression data was found for gemcitabine
in pancreatic organoid drugs within this work, nor for sarcoma, melanoma (dacarbazine) or the PDX drugs.
Numbers: tables/Fig4bc_isogenic_other_cancers.tsv.
""")
