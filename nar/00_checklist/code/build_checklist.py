#!/usr/bin/env python3
"""Build the paper checklist for Kong et al. 2020 Nat Commun (10.1038/s41467-020-19313-8).

Every number below was read by hand from the cached sources in nar/_cache/:
  - main article PDF (PMC OA S3 copy, PMC7599252.1.pdf -> 41467_2020_Article_19313.pdf), figure pages rendered at 150 dpi
  - JATS XML from NCBI efetch (PMC7599252_jats.xml): figure legends + text with figure references
  - supplementary PDF MOESM1 (18 Supp Figs + 4 Supp Tables, images extracted with pdfimages)
  - source data MOESM3 xlsx (13 sheets)
Numbers tagged [derived] were computed by derive_source_data_facts.py from MOESM3, not printed in the paper.
Numbers tagged [read from plot, approx] were read off a plotted axis by eye.

Run:  python3 nar/00_checklist/code/build_checklist.py
Writes: nar/00_checklist/checklist.tsv and nar/00_checklist/PAPER_CHECKLIST.md
"""
import csv
import re
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)

# Short codes used in the generalisation column (defined in the markdown, section 3).
G_KM = ("Same KM split for each GEN cancer-drug pair: LUAD cisplatin/carboplatin/pemetrexed; OV carboplatin/paclitaxel; "
        "BLCA gemcitabine/methotrexate; PAAD gemcitabine; LIHC sorafenib; COAD-external GSE39582 5FU.")

ROWS = [
# ------------------------------------------------------------------ MAIN FIGURE 1
dict(ID="Fig1a", kind="main-figure", cancer_drug="all (schematic)",
     shows="Cartoon of the whole method: pathways near drug targets in a protein network -> train ML on organoid drug response -> pick best pathway -> split patients into predicted responders/non-responders.",
     numbers="No numbers (schematic).",
     data="None (drawing).",
     generalise="One schematic for the NAR paper showing the same 4 steps with a box listing all cancers/models used (organoid, PDC, cell line, PDX)."),
dict(ID="Fig1b", kind="main-figure", cancer_drug="all (schematic)",
     shows="Four-step detail: (1) network distance drug targets->pathway genes; (2-3) pathway ssGSEA scores + IC50 train regression, coefficient = predictive performance; (4) biomarker expression splits patients, compare overall survival.",
     numbers="No numbers (schematic).",
     data="None (drawing).",
     generalise="Same; add step for external validation cohorts (GEO trials) and for cross-model transfer."),
# ------------------------------------------------------------------ MAIN FIGURE 2
dict(ID="Fig2a", kind="main-figure", cancer_drug="COAD / 5-fluorouracil",
     shows="Each dot is a Reactome pathway: x = average shortest network distance from the 5FU target (TYMS) to the pathway, y = proximity z-score. Pathways with z <= -1.2816 (shaded) are 'proximal'.",
     numbers="674 Reactome pathways tested; 37 proximal pathways (z <= -1.2816, i.e. alpha=0.10). [derived] 38 pathways pass the z cutoff in the source matrix; 37 of them have organoid ssGSEA scores. Biomarker 'activation of BH3-only proteins' z = -1.2857 [derived] (just past the cutoff).",
     data="STRING edges (version not stated) score>700, largest component (13,824 proteins, 323,774 links; MOESM3 sheet 'Figure 2a-f'); DrugBank targets (5FU -> TYMS); MSigDB C2 Reactome; proximity z-scores (MOESM3 sheet 'Figure 2a,d', column FLUOROURACIL). 1000 degree-matched random draws.",
     generalise="Proximity scatter for each GEN drug: cisplatin/carboplatin (LUAD, OV), pemetrexed (LUAD), paclitaxel (OV), gemcitabine (PAAD, BLCA), sorafenib (LIHC), olaparib (OV). Platinum targets are DrugBank 'targets' that are not the real DNA mechanism - note this."),
dict(ID="Fig2b", kind="main-figure", cancer_drug="COAD / 5-fluorouracil",
     shows="Bar chart of ridge-regression |coefficient| ('predictive ability') for the 37 proximal pathways, trained on 19 colorectal organoids' 5FU IC50. The chosen biomarker (orange) is rank 1.",
     numbers="37 pathways; rank 1 = 'Activation of BH3-only proteins', |coef| ~0.8 [read from plot, approx]; next pathways ~0.6. RidgeCV alpha 0.1-1.0 step 0.1, 3-fold CV.",
     data="COAD organoid expression GSE64392 (RMA) -> ssGSEA NES (MOESM3 'Figure 2b', 612 pathways x 22 columns [derived]); 5FU IC50 from van de Wetering 2015 supplement (data/curated/preclinical colorectal_vandewetering2015).",
     generalise="Ranked-coefficient bar chart per model/drug: OV organoids (Vias 2023) carboplatin/paclitaxel; PAAD organoids (Tiriac 2018) gemcitabine; LIHC organoids (Broutier 2017/Ji 2023) sorafenib; LUAD from GDSC/CTRP lung lines or GSE165611 PDCs; BLCA other drugs from Lee 2018."),
dict(ID="Fig2c", kind="main-figure", cancer_drug="COAD / 5-fluorouracil",
     shows="Network drawing: 5FU -> TYMS (red triangle) linked via TFDP1/TP53 to the genes of the BH3-only pathway (E2F1, TFDP1, TP53, AKT1, PMAIP1, DYNLL1, YWHAB, MAPK8, BBC3, BID, BCL2L11, BCL2, BAD, PPP3R1, BMF).",
     numbers="No statistics; gene list as drawn.",
     data="STRING edges (score>700) among TYMS + pathway genes; Reactome gene set REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS.",
     generalise="Subnetwork plot (drug target + chosen biomarker pathway genes + shortest-path links) for each new cancer-drug biomarker."),
dict(ID="Fig2d-left", kind="main-figure", cancer_drug="COAD / 5-fluorouracil (TCGA treated)",
     shows="Kaplan-Meier overall survival of TCGA-COAD patients who got 5FU, split at the median predicted resistance score (biomarker NES x organoid coefficient).",
     numbers="n=114 (57 predicted responders / 57 non-responders); two-sided log-rank P = 0.014. Responders: no deaths; non-responders 5-year OS ~50%.",
     data="TCGA-COAD FPKM-UQ log2(x+1) -> ssGSEA NES (MOESM3 'Figure 2g', 674 pathways x 449 patients [derived]); TCGA drug table (5FU-treated, DrugBank-normalised names); overall survival.",
     generalise=G_KM + " Brief says only ~5 TCGA-COAD got 5FU alone; most got FOLFOX - say which definition is used."),
dict(ID="Fig2d-right", kind="main-figure", cancer_drug="COAD / untreated control (TCGA)",
     shows="Same split applied to TCGA-COAD patients with no recorded drug treatment (negative control: is the biomarker just prognostic?).",
     numbers="n=298 (149/149); log-rank P = 0.16.",
     data="Same as Fig2d-left, untreated patients.",
     generalise="Untreated-patient control KM for every GEN cancer (TCGA untreated of same cancer; GSE39582 has 326 untreated). Better: formal treatment x biomarker interaction Cox test."),
dict(ID="Fig2e", kind="main-figure", cancer_drug="BLCA / cisplatin",
     shows="Proximity scatter for cisplatin targets vs Reactome pathways (as Fig2a).",
     numbers="30 proximal pathways. [derived] 47 pathways pass z <= -1.2816 in the source matrix; 30 of them have BLCA organoid ssGSEA scores. Biomarker 'amino acid synthesis and interconversion (transamination)' z = -2.115 [derived].",
     data="As Fig2a, column CISPLATIN; cisplatin DrugBank targets MPG, A2M, ATOX1, TF.",
     generalise="As Fig2a for each GEN drug."),
dict(ID="Fig2f", kind="main-figure", cancer_drug="BLCA / cisplatin",
     shows="Ridge |coefficient| bar chart for 30 proximal pathways trained on bladder organoid cisplatin IC50. Chosen biomarker (orange) is rank 2 (rank 1 is a grey bar).",
     numbers="30 pathways; biomarker |coef| ~0.2, rank-1 pathway ~0.3 [read from plot, approx]. Biomarker is rank 2, not rank 1 (same in Supp Fig 2c,d).",
     data="BLCA organoid expression + cisplatin IC50 from Lee et al. 2018 (MOESM3 'Figure 2e': 488 pathways x 11 organoid columns [derived]; text says 9 organoids; Supp Fig 17 shows 9 with IC50).",
     generalise="As Fig2b for each GEN model/drug."),
dict(ID="Fig2g", kind="main-figure", cancer_drug="BLCA / cisplatin",
     shows="Network drawing: cisplatin targets MPG, A2M, ATOX1, TF linked to amino-acid pathway genes (KYAT1, ASNS, OAT, GLUD1, GOT2, PSAT1, PSPH, PYCR1, PHGDH, GLUL, GLS, GPT, GOT1, GPT2).",
     numbers="No statistics.",
     data="STRING edges among targets + pathway genes. NOTE: MOESM3 sheet 'Figure 2g' actually holds TCGA-COAD ssGSEA NES, not this network.",
     generalise="As Fig2c."),
dict(ID="Fig2h-left", kind="main-figure", cancer_drug="BLCA / cisplatin (TCGA treated)",
     shows="KM overall survival of cisplatin-treated TCGA-BLCA patients split at median predicted resistance.",
     numbers="n=77 (39 responders / 38 non-responders); log-rank P = 0.01. 5-year OS ~75% vs ~50%.",
     data="TCGA-BLCA FPKM-UQ -> ssGSEA (MOESM3 'Figure 2h', 490 pathways x 407 patients [derived]); cisplatin-treated patients; OS.",
     generalise=G_KM),
dict(ID="Fig2h-right", kind="main-figure", cancer_drug="BLCA / untreated control (TCGA)",
     shows="Same split in TCGA-BLCA patients with no drug treatment (negative control).",
     numbers="n=294 (147/147); log-rank P = 0.066.",
     data="As Fig2h-left, untreated patients.",
     generalise="Untreated control + interaction test for each GEN cancer."),
# ------------------------------------------------------------------ MAIN FIGURE 3
dict(ID="Fig3a", kind="main-figure", cancer_drug="COAD 5FU + BLCA cisplatin (schematic)",
     shows="Scheme of comparison methods: (1) no feature selection, (2) other feature selection, both -> ML -> predict TCGA responders.",
     numbers="No numbers.",
     data="None.",
     generalise="Same scheme, add modern baselines (elastic net on all genes, random pathways, published signatures)."),
dict(ID="Fig3b-COAD", kind="main-figure", cancer_drug="COAD / 5FU", shows="KM: ridge model on whole transcriptome (all genes) of organoids, applied to 5FU-treated TCGA-COAD.",
     numbers="n=57/57; log-rank P = 0.69.", data="Organoid gene-level expression + IC50; TCGA gene expression; OS.",
     generalise="Whole-transcriptome baseline for every GEN cancer-drug pair."),
dict(ID="Fig3b-BLCA", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="KM: whole-transcriptome ridge model, cisplatin-treated TCGA-BLCA.",
     numbers="n=39/38; log-rank P = 0.79.", data="As above for BLCA.", generalise="As Fig3b-COAD."),
dict(ID="Fig3c-COAD", kind="main-figure", cancer_drug="COAD / 5FU", shows="KM: ridge model on all pathways (whole-pathway ssGSEA).",
     numbers="n=57/57; log-rank P = 0.81.", data="Organoid + TCGA ssGSEA of all Reactome pathways.", generalise="Whole-pathway baseline for every GEN pair."),
dict(ID="Fig3c-BLCA", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="KM: whole-pathway ridge model.",
     numbers="n=39/38; log-rank P = 0.82.", data="As above.", generalise="As Fig3c-COAD."),
dict(ID="Fig3d-COAD", kind="main-figure", cancer_drug="COAD / 5FU", shows="KM: ridge model on genes that are direct (first) network neighbours of the drug target.",
     numbers="n=57/57; log-rank P = 0.96.", data="STRING first neighbours of TYMS; organoid + TCGA gene expression.", generalise="First-neighbour baseline for each GEN drug."),
dict(ID="Fig3d-BLCA", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="KM: first-neighbour gene model.",
     numbers="n=39/38; log-rank P = 0.37.", data="STRING first neighbours of MPG/A2M/ATOX1/TF.", generalise="As Fig3d-COAD."),
dict(ID="Fig3e-COAD", kind="main-figure", cancer_drug="COAD / 5FU", shows="KM: genes pre-selected by Bolis et al. method (Spearman P<0.05 with IC50 over 10 leave-half-out repeats), then ridge.",
     numbers="n=57/57; log-rank P = 0.15. 15 genes selected [derived from MOESM3 'Figure 3e'].", data="Organoid gene expression + IC50; MOESM3 'Figure 3e' gene lists.", generalise="Correlation-filter baseline for each GEN pair."),
dict(ID="Fig3e-BLCA", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="KM: Bolis feature selection + ridge.",
     numbers="n=39/38; log-rank P = 0.96. 8 genes selected [derived].", data="As above.", generalise="As Fig3e-COAD."),
dict(ID="Fig3f-COAD", kind="main-figure", cancer_drug="COAD / 5FU", shows="KM: deep neural network of Sharifi-Noghabi et al. (MOLI-style, triplet + cross-entropy loss), organoids split at median IC50.",
     numbers="n=57/57; log-rank P = 0.91. Hyper-parameters + CV AUC in MOESM3 'Figure 3f'.", data="Organoid + TCGA expression; MOESM3 'Figure 3f'.", generalise="Optional DL baseline (CPU only; torch from PyPI) for GEN pairs with >=30 training models."),
dict(ID="Fig3f-BLCA", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="KM: deep-learning model.",
     numbers="n=39/38; log-rank P = 0.82.", data="As above.", generalise="As Fig3f-COAD."),
# ------------------------------------------------------------------ MAIN FIGURE 4
dict(ID="Fig4a", kind="main-figure", cancer_drug="isogenic lines (schematic)", shows="Scheme: drug-sensitive parent cells made resistant by rising drug doses; compare transcriptomes.",
     numbers="No numbers.", data="None.", generalise="Same scheme; list isogenic pairs per GEN cancer."),
dict(ID="Fig4b", kind="main-figure", cancer_drug="COAD / 5FU (HCT-8 vs HCT-8/5FU)", shows="Boxplot of BH3-only pathway ssGSEA score in 5FU-sensitive vs 5FU-resistant colorectal cells.",
     numbers="Sensitive n=3, resistant n=9 (resistant = 0h, 24h, 48h triplicates); two-sided two-sample t-test P = 0.0085; sensitive higher.",
     data="GEO GSE81008 microarray -> ssGSEA (MOESM3 'Figure 4b', 674 pathways x 12 samples).",
     generalise="Isogenic sensitive/resistant pair per GEN drug from GEO: e.g. A549 vs A549/DDP (LUAD cisplatin), A2780 vs A2780cis or PEO1 vs PEO4 (OV platinum), gemcitabine-resistant PANC-1/MIA PaCa-2 (PAAD), sorafenib-resistant Huh7/HepG2 (LIHC). Accessions must be checked, not assumed."),
dict(ID="Fig4c", kind="main-figure", cancer_drug="BLCA / cisplatin (T24S vs T24R)", shows="Boxplot of amino-acid pathway score in cisplatin-sensitive vs -resistant T24 bladder cells.",
     numbers="Sensitive n=3, resistant n=3; t-test P = 0.00022; resistant lower.",
     data="Yeon et al. data by personal request (not public); MOESM3 'Figure 4c' gives ssGSEA NES for 6 samples.",
     generalise="As Fig4b. For BLCA a public cisplatin-resistant bladder line set from GEO is the substitute for raw data."),
dict(ID="Fig4d", kind="main-figure", cancer_drug="bootstrap (schematic)", shows="Scheme of 10,000-iteration random-pathway bootstrap.",
     numbers="10,000 iterations; random pathway sets matched to number of proximal pathways; 'top rank' = best ranked pathway that is both survival-predictive (P<0.05) and isogenic-predictive (P<0.05).",
     data="None.", generalise="Same null for each GEN pair."),
dict(ID="Fig4e", kind="main-figure", cancer_drug="COAD / 5FU", shows="Histogram of the top predictive rank reached by random pathway sets; observed rank from proximal pathways in green.",
     numbers="Empirical P = 0.0012. [derived] observed rank = 1 (12/10,000 random sets reach rank 1); 8,927/10,000 random sets give NA (no pathway passes both tests).",
     data="MOESM3 'Figure 4e' (10,000 rows: top_rank, feature, coefficient); requires full pipeline + TCGA + GSE81008.",
     generalise="Random-pathway empirical P for each GEN pair (needs an isogenic or second validation set to define 'predictive'; otherwise survival-only null)."),
dict(ID="Fig4f", kind="main-figure", cancer_drug="BLCA / cisplatin", shows="Same bootstrap histogram for cisplatin/BLCA.",
     numbers="Empirical P = 0.014. [derived] observed rank = 2 (66 sets rank 1 + 74 rank 2 = 140/10,000); 6,366/10,000 NA.",
     data="MOESM3 'Figure 4f'.", generalise="As Fig4e."),
# ------------------------------------------------------------------ MAIN FIGURE 5
dict(ID="Fig5a", kind="main-figure", cancer_drug="COAD / cetuximab vs BRAF V600E (TCGA)",
     shows="Boxplot of predicted cetuximab resistance (Gastrin-CREB via PKC and MAPK pathway, trained on organoids) in BRAF wild-type vs BRAF V600E TCGA-COAD tumours.",
     numbers="WT n=328, V600E n=47; one-sided Mann-Whitney U P = 0.037 (V600E higher). Trained after removing 3 hypermutated organoids (P19a, P19b, P24a).",
     data="Organoid cetuximab IC50 (proximity column CETUXIMAB_GEFITINIB); TCGA-COAD expression + MAF mutations (BRAF).",
     generalise="Mutation-concordance test per GEN cancer: OV BRCA1/2 mut vs platinum (expect lower resistance); LUAD EGFR mut vs EGFR TKI (erlotinib/gefitinib) or KEAP1/STK11 vs platinum; BLCA ERCC2 (already); PAAD BRCA/ATM vs platinum; COAD KRAS mut vs cetuximab."),
dict(ID="Fig5b", kind="main-figure", cancer_drug="BLCA / cisplatin vs ERCC2 (TCGA)",
     shows="Boxplot of predicted cisplatin resistance (amino-acid pathway) in ERCC2-mutant vs wild-type TCGA-BLCA tumours.",
     numbers="ERCC2 mut n=37, WT n=367; one-sided Mann-Whitney U P = 0.002 (mutants lower = more sensitive).",
     data="TCGA-BLCA expression + mutations (all patients, not only treated).",
     generalise="As Fig5a."),
# ------------------------------------------------------------------ SUPPLEMENTARY FIGURES
dict(ID="SuppFig1a-d", kind="supp-figure", cancer_drug="COAD + BLCA (organoid vs TCGA)",
     shows="PCA of organoids (red) and TCGA patients (grey) before scaling: genes and pathways. Organoids sit apart from patients (batch effect).",
     numbers="PC1/PC2 explained variance: (a) COAD genes 39.8%/5.05%; (b) COAD pathways 47.3%/13.6%; (c) BLCA genes 21.6%/6.37%; (d) BLCA pathways 27.0%/17.8%.",
     data="Organoid + TCGA gene expression and ssGSEA NES.",
     generalise="Batch-effect PCA for every GEN model-vs-patient pairing (organoid/PDC/cell line vs TCGA or GEO cohort)."),
dict(ID="SuppFig1e-h", kind="supp-figure", cancer_drug="COAD + BLCA",
     shows="Same PCA after z-scoring each gene/pathway within each dataset: organoids now mix with patients.",
     numbers="(e) COAD scaled genes 12.2%/8.93%; (f) COAD scaled pathways 27.1%/10.7%; (g) BLCA scaled genes 9.50%/6.87%; (h) BLCA scaled pathways 20.7%/18.8%.",
     data="As above, z-scored per dataset.", generalise="As SuppFig1a-d."),
dict(ID="SuppFig2a", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Predictive performance ranking with linear regression (37 pathways).",
     numbers="BH3-only pathway rank 1, ~0.88 [read from plot, approx].", data="As Fig2b.", generalise="Algorithm-robustness ranking per GEN pair."),
dict(ID="SuppFig2b", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Same with support vector regression (linear kernel).",
     numbers="BH3-only rank 1, ~0.82 [approx].", data="As Fig2b.", generalise="As SuppFig2a."),
dict(ID="SuppFig2c", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Linear regression ranking (30 pathways).",
     numbers="Amino-acid pathway rank 2, ~0.23; rank 1 ~0.34 [approx].", data="As Fig2f.", generalise="As SuppFig2a."),
dict(ID="SuppFig2d", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="SVR ranking.",
     numbers="Amino-acid pathway rank 2, ~0.20; rank 1 ~0.29 [approx].", data="As Fig2f.", generalise="As SuppFig2a."),
dict(ID="SuppFig3a-c", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Scatter of regression coefficients between algorithms (37 pathways).",
     numbers="R2: ridge vs linear 0.9948; SVR vs linear 0.9979; SVR vs ridge 0.9973.", data="Coefficients from 3 models on COAD organoids.", generalise="Coefficient agreement plot per GEN pair."),
dict(ID="SuppFig3d-f", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Same for bladder (30 pathways). Legend says (d)-(e) but figure has d, e, f.",
     numbers="R2: ridge vs linear 0.9978; SVR vs linear 0.9932; SVR vs ridge 0.9977.", data="As above, BLCA.", generalise="As SuppFig3a-c."),
dict(ID="SuppFig4a", kind="supp-figure", cancer_drug="COAD / 5FU with BLCA biomarker", shows="Cross-test: the bladder cisplatin biomarker (amino-acid pathway) used on 5FU-treated COAD patients.",
     numbers="n=57/57; log-rank P = 0.16.", data="TCGA-COAD NES; BLCA organoid coefficient.", generalise="Cross-cancer/cross-drug specificity matrix: every biomarker x every treated cohort (heatmap of log-rank P / HR)."),
dict(ID="SuppFig4b", kind="supp-figure", cancer_drug="BLCA / cisplatin with COAD biomarker", shows="Cross-test: BH3-only pathway used on cisplatin-treated BLCA patients.",
     numbers="n=39/38; log-rank P = 0.86.", data="TCGA-BLCA NES; COAD organoid coefficient.", generalise="As SuppFig4a."),
dict(ID="SuppFig5a", kind="supp-figure", cancer_drug="COAD / 5FU", shows="KM using weighted sum of top-2 pathways.",
     numbers="n=57/57; log-rank P = 0.88. (Main text calls this 5b; legend/figure say 5a.)", data="Top-2 pathways x ridge coefficients; TCGA-COAD.", generalise="Top-k (1..10) multi-pathway curve per GEN pair."),
dict(ID="SuppFig5b", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="KM using top-2 pathways.",
     numbers="n=39/38; log-rank P = 0.022. (Main text calls this 5a.)", data="As above, BLCA.", generalise="As SuppFig5a."),
dict(ID="SourceData-multipathway", kind="source-data", cancer_drug="COAD 5FU + BLCA cisplatin",
     shows="Table of top-k pathway models (k=2..10) for Ridge, SVR and LinearRegression: 5-year survival in each group and log-rank P (MOESM3 sheet 'FigureS5').",
     numbers="Ridge COAD P for k=2..10: 0.882, 0.373, 0.314, 0.130, 0.253, 0.076, 0.216, 0.314, 0.314. Ridge BLCA: 0.022, 0.496, 0.161, 0.480, 0.709, 0.342, 0.639, 0.126, 0.256. Smallest non-k2 values: COAD 0.076 (k=7), BLCA LinearRegression k=10 0.0625.",
     data="MOESM3 'FigureS5'; original code python/multiple_pathway_prediction.py produces it.", generalise="Same k-curve per GEN pair."),
dict(ID="SuppFig6a", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Degree-centrality feature selection (37 most central pathways): ranking + KM of top 3 pathways.",
     numbers="Top 3: NFkB activated and signals survival P=0.997; SPRY regulation of FGF signaling P=0.71; IRAK2-mediated TAK1 activation P=0.802 (n=57/57).",
     data="STRING degree centrality averaged over pathway genes; organoid + TCGA ssGSEA.", generalise="Centrality-selected baseline for each GEN pair."),
dict(ID="SuppFig6b", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Betweenness-centrality selection.",
     numbers="Top 3: Membrane binding and targeting of GAG proteins P=0.79; Regulation of Fanconi anemia pathway P=0.503; Hormone-sensitive lipase P=0.91.",
     data="As above (betweenness).", generalise="As SuppFig6a."),
dict(ID="SuppFig6c", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Closeness-centrality selection.",
     numbers="Top 3: Signaling by constitutively active EGFR P=0.17; SPRY regulation of FGF signaling P=0.71; CTLA4 inhibitory signaling P=0.32.",
     data="As above (closeness).", generalise="As SuppFig6a."),
dict(ID="SuppFig7a", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Degree-centrality selection (30 pathways) + KM of top 3.",
     numbers="G alpha z signaling events P=0.094; mRNA capping P=0.96; GABA B receptor activation P=0.51 (n=39/38).",
     data="As SuppFig6.", generalise="As SuppFig6a."),
dict(ID="SuppFig7b", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Betweenness selection.",
     numbers="Pre-NOTCH transcription and translation P=0.23; G1 phase P=0.31; TGF-beta receptor signaling activates SMADs P=0.091.",
     data="As SuppFig6.", generalise="As SuppFig6a."),
dict(ID="SuppFig7c", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Closeness selection.",
     numbers="Signaling by WNT P=0.56; EGFR downregulation P=0.29; CTNNB1 phosphorylation cascade P=0.74.",
     data="As SuppFig6.", generalise="As SuppFig6a."),
dict(ID="SuppFig8a", kind="supp-figure", cancer_drug="COAD / 5FU", shows="KM split by expression of the drug target TYMS alone (median).",
     numbers="high n=57, low n=57; log-rank P = 0.73.", data="TCGA-COAD gene expression.", generalise="Drug-target-expression KM per GEN drug (e.g. TUBB for paclitaxel, TYMS/RRM1 for pemetrexed/gemcitabine, PARP1 for olaparib)."),
dict(ID="SuppFig8b", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="KM split by each cisplatin target gene.",
     numbers="high n=38 / low n=39: A2M P=0.92; ATOX1 P=0.18; MPG P=0.83; TF P=0.48.", data="TCGA-BLCA gene expression.", generalise="As SuppFig8a."),
dict(ID="SuppFig9a", kind="supp-figure", cancer_drug="COAD 5FU + BLCA cisplatin", shows="KM with ridge model on genes up to 2nd network neighbours of drug targets.",
     numbers="COAD P = 0.24 (57/57); BLCA P = 0.79 (39/38).", data="STRING 2-hop neighbourhood; organoid + TCGA gene expression.", generalise="k-hop neighbour baselines per GEN pair."),
dict(ID="SuppFig9b", kind="supp-figure", cancer_drug="COAD 5FU + BLCA cisplatin", shows="Same with up to 3rd neighbours.",
     numbers="COAD P = 0.36; BLCA P = 0.79.", data="STRING 3-hop.", generalise="As SuppFig9a."),
dict(ID="SuppFig10a", kind="supp-figure", cancer_drug="COAD / 5FU", shows="Scatter: pathway size (number of genes) vs log-rank P when each pathway alone splits 5FU-treated patients (median).",
     numbers="R2 = 0.0026.", data="All pathways' TCGA-COAD NES + gene-set sizes.", generalise="Pathway-size confounding check per GEN pair."),
dict(ID="SuppFig10b", kind="supp-figure", cancer_drug="BLCA / cisplatin", shows="Same for cisplatin-treated BLCA.",
     numbers="R2 = 0.00034.", data="As above, BLCA.", generalise="As SuppFig10a."),
dict(ID="SuppFig11", kind="supp-figure", cancer_drug="COAD / cetuximab (organoids)", shows="Heatmap of 37 cetuximab-proximal pathways across 19 colorectal organoids sorted by cetuximab IC50, with ridge weights. Top weight: Gastrin-CREB signalling via PKC and MAPK.",
     numbers="19 organoids (P18 most sensitive ... P31 most resistant; IC50 ~2.4 to ~8); 37 pathways; Gastrin-CREB weight ~0.75 [approx]; most negative 'role of DCC in regulating apoptosis' ~-0.4 [approx].",
     data="COAD organoid NES + cetuximab IC50; proximity column CETUXIMAB_GEFITINIB.", generalise="Organoid heatmap + weights for each GEN model/drug (e.g. OV organoid carboplatin, PAAD organoid gemcitabine)."),
dict(ID="SuppFig12a-e", kind="supp-figure", cancer_drug="COAD / cetuximab vs BRAF V600E", shows="Same BRAF test as Fig5a using other feature sets: (a) whole transcriptome, (b) whole pathway, (c) direct neighbours, (d) 2nd neighbours, (e) 3rd neighbours.",
     numbers="WT n=328, V600E n=47; one-sided MWU P: (a) 0.98, (b) 0.63, (c) 0.99, (d) 0.99, (e) 0.99.", data="As Fig5a.", generalise="Baseline comparison for each GEN mutation-concordance test."),
dict(ID="SuppFig13a-e", kind="supp-figure", cancer_drug="BLCA / cisplatin vs ERCC2", shows="Same ERCC2 test as Fig5b with other feature sets (a-e as above).",
     numbers="ERCC2 mut n=37, WT n=367; one-sided MWU P: (a) 0.04, (b) 0.037, (c) 0.017, (d) 0.021, (e) 0.028. (All < 0.05; paper's own model P = 0.002.)", data="As Fig5b.", generalise="As SuppFig12a-e."),
dict(ID="SuppFig14a-h", kind="supp-figure", cancer_drug="COAD + BLCA (GDSC cell lines vs TCGA)", shows="PCA of GDSC cancer cell lines (red) vs TCGA patients (grey), gene and pathway level, unscaled (all lines and tissue-matched lines).",
     numbers="Explained-variance values printed on axes (not transcribed).", data="GDSC/Cell Model Passports expression; TCGA.", generalise="Batch PCA for every cell-line model used (GDSC/CTRP lung, ovary, pancreas, liver lines)."),
dict(ID="SuppFig14i-p", kind="supp-figure", cancer_drug="COAD + BLCA (GDSC vs TCGA)", shows="Same PCA after per-dataset z-scoring.",
     numbers="Values on axes (not transcribed).", data="As above.", generalise="As SuppFig14a-h."),
dict(ID="SuppFig15a", kind="supp-figure", cancer_drug="COAD / 5FU from all GDSC lines", shows="Pathway ranking trained on all 893 GDSC lines' 5FU response + KM of top 3 pathways in 5FU-treated TCGA-COAD.",
     numbers="n=893 lines; Cell cycle P=0.17; CDC6 association with ORC origin complex P=0.58; Association of licensing factors with pre-replicative complex P=0.59 (n=57/57).",
     data="GDSC 5FU IC50 + expression; TCGA-COAD.", generalise="Cell-line-trained comparator for each GEN pair (all lines vs tissue-matched lines)."),
dict(ID="SuppFig15b", kind="supp-figure", cancer_drug="COAD / 5FU from colorectal GDSC lines", shows="Same with 46 colorectal lines.",
     numbers="n=46; Metabolism of nucleotides P=0.86; Cell cycle mitotic P=0.17; E2F-enabled inhibition of pre-replication complex formation P=0.51.",
     data="GDSC colorectal lines.", generalise="As SuppFig15a."),
dict(ID="SuppFig15c", kind="supp-figure", cancer_drug="BLCA / cisplatin from all GDSC lines", shows="Ranking from 831 lines' cisplatin response + KM of top 3 in cisplatin-treated TCGA-BLCA.",
     numbers="n=831; Transport of glucose and other sugars P=0.102; Transmembrane transport of small molecules P=0.0092; Platelet activation signaling and aggregation P=0.84 (n=39/38).",
     data="GDSC cisplatin; TCGA-BLCA.", generalise="As SuppFig15a."),
dict(ID="SuppFig15d", kind="supp-figure", cancer_drug="BLCA / cisplatin from bladder GDSC lines", shows="Same with 14 bladder lines.",
     numbers="n=14; BMAL1:CLOCK/NPAS2 activates circadian expression P=0.13; Circadian clock P=0.23; Transmembrane transport of small molecules P=0.0092.",
     data="GDSC bladder lines.", generalise="As SuppFig15a."),
dict(ID="SuppFig16", kind="supp-figure", cancer_drug="44 organoid drugs (COAD + BLCA screens)", shows="For each Reactome pathway, how many of the 44 drugs it is proximal to (checks that no pathway is proximal to everything).",
     numbers="674 pathways, 44 drugs; max 31 drugs [derived]; only 3 pathways (0.45%) proximal to >=22 drugs; mean 4.9 (derived 4.86), median 4 drugs per pathway; 111 pathways proximal to no drug.",
     data="MOESM3 'Figure 2a,d' (674 x 44 z-scores; some columns merge drugs with identical targets, e.g. CETUXIMAB_GEFITINIB).",
     generalise="Same hub-bias check for the full GEN drug list (all drugs used across cancers)."),
dict(ID="SuppFig17a", kind="supp-figure", cancer_drug="COAD / 5FU (organoids)", shows="Within-organoid prediction: 60/10/30 train/validation/test split; observed vs predicted IC50 in test organoids over all validation choices.",
     numbers="R2 = 0.98 (test organoids P11, P26, P8, P18, P14, P16 in green).", data="COAD organoid NES of proximal pathways + 5FU IC50.", generalise="Proper repeated nested CV on each GEN preclinical set (report Spearman/Pearson + permutation P, not one split)."),
dict(ID="SuppFig17b", kind="supp-figure", cancer_drug="BLCA / cisplatin (organoids)", shows="Same for 9 bladder organoids.",
     numbers="R2 = 0.89 (test organoids SCBO-1, SCBO-11.2, SCBO-4).", data="BLCA organoid NES + cisplatin IC50.", generalise="As SuppFig17a."),
dict(ID="SuppFig18", kind="supp-figure", cancer_drug="COAD / cetuximab vs BRAF V600E", shows="Robustness: fraction of organoid subsets (drop 0-3 organoids) that still recover the BRAF V600E-cetuximab association (one-sided MWU P<0.05).",
     numbers="Combinations tested: 969 (16 organoids), 171 (17), 19 (18), 1 (19). Network method ~0.73/0.80/0.89/1.0; whole-transcriptome ~0.64/0.71/0.79/1.0; whole-pathway ~0.40/0.38/0.31/0.0 [read from plot, approx]. Text: >70% recover with 3 organoids removed.",
     data="As Fig5a, leave-k-out over organoids.", generalise="Leave-k-models-out robustness curve for each GEN concordance/survival result."),
# ------------------------------------------------------------------ SUPPLEMENTARY TABLES
dict(ID="SuppTable1", kind="supp-table", cancer_drug="COAD / 5FU (TCGA)", shows="Clinical features of predicted responders vs non-responders (chi-squared tests).",
     numbers="n=114. Age <60: 53 (30 R/23 NR), >=60: 61 (27/34), P=0.26. Male 62 (27/35), female 52 (30/22), P=0.19. T1-T2 7 (2/5), T3-T4 107 (55/52), P=0.44. N0 27 (16/11), N1 52 (26/26), N2 35 (15/20), P=0.44.",
     data="TCGA-COAD clinical (age, gender, T, N) + group labels.", generalise="Clinical-balance table per GEN cohort; better: multivariable Cox with stage/age (and MMR/KRAS/BRAF for CRC)."),
dict(ID="SuppTable2", kind="supp-table", cancer_drug="BLCA / cisplatin (TCGA)", shows="Clinical features of predicted groups.",
     numbers="n=77. Age <60: 21 (11/10), >=60: 56 (28/28), P=0.94. Male 53 (27/26), female 24 (12/12), P=0.87. T1-T2 12 (6/6), T3-T4 65 (33/32), P=0.79. N0 39 (22/17), N1 13 (3/10), N2 25 (14/11), P=0.093.",
     data="TCGA-BLCA clinical.", generalise="As SuppTable1."),
dict(ID="SuppTable3", kind="supp-table", cancer_drug="COAD 5FU + BLCA cisplatin (isogenic lines)", shows="Top-3 most proximal pathways: proximity z and mean ssGSEA in sensitive vs resistant isogenic lines.",
     numbers="COAD: MITOTIC_G1_G1_S_PHASES z=-4.26, 0.34 vs 0.31, P=0.15; G1_S_TRANSITION z=-3.9, 0.38 vs 0.34, P=0.15; PYRIMIDINE_METABOLISM z=-3.42, 0.17 vs 0.16, P=0.11. BLCA: METAL_ION_SLC_TRANSPORTERS z=-4.03, 0.078 vs 0.06, P=0.014; RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 z=-3.75, 0.14 vs 0.12, P=0.00054; RECYCLING_OF_BILE_ACIDS_AND_SALTS z=-3.73, -0.28 vs -0.22, P=0.01.",
     data="Proximity z (MOESM3 'Figure 2a,d') + isogenic NES (MOESM3 'Figure 4b','Figure 4c').", generalise="Same table for each GEN isogenic pair."),
dict(ID="SuppTable4", kind="supp-table", cancer_drug="COAD 5FU + BLCA cisplatin (bootstrap)", shows="Random-bootstrap pathways that ranked equal to or better than the real biomarker, with their proximity z.",
     numbers="COAD: PI3K_EVENTS_IN_ERBB2_SIGNALING z=1.4. BLCA: REGULATION_OF_KIT_SIGNALING -1.51; NOTCH1_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION -0.47; SIGNALING_BY_NOTCH1 -0.16; ADAPTIVE_IMMUNE_SYSTEM 1.63; IL_2_SIGNALING 2.09; CA_DEPENDENT_EVENTS 2.94.",
     data="MOESM3 'Figure 4e','Figure 4f' + proximity matrix.", generalise="Same table per GEN pair."),
# ------------------------------------------------------------------ SOURCE-DATA-ONLY / TEXT-ONLY RESULTS
dict(ID="SourceData-STRING", kind="source-data", cancer_drug="all", shows="STRING network edge list used for proximity (MOESM3 'Figure 2a-f').",
     numbers="Text: 13,824 proteins, 323,774 interactions (score>700, largest component). [derived] sheet has 647,548 rows = each link listed in both directions; 13,825 distinct IDs across the two UniProt columns.",
     data="MOESM3 'Figure 2a-f'.", generalise="Re-use same network for every GEN cancer; sensitivity check with a second network (e.g. STRING 900, or HuRI/BioGRID)."),
dict(ID="SourceData-DLhyper", kind="source-data", cancer_drug="COAD 5FU + BLCA cisplatin", shows="Deep-learning hyper-parameter grid and mean CV AUC (MOESM3 'Figure 3f').",
     numbers="Grid columns: batch_size, #nodes, learning rates, dropouts, weight_decay, gamma, #epoch, #Folds, margin, average_AUC (highlighted rows used).",
     data="MOESM3 'Figure 3f'.", generalise="Only if DL baseline is repeated."),
dict(ID="Text-method-params", kind="text", cancer_drug="all", shows="Key method settings stated in Methods.",
     numbers="Reactome from MSigDB C2 (674 pathways usable); proximity = closest distance d_c, z from 1000 degree-matched random sets, proximal if z <= -1.2816; ridge RidgeCV alpha 0.1..1.0 step 0.1, 3-fold CV; SVR linear kernel; score = pathway NES x organoid coefficient; median split; lifelines log-rank; TCGA FPKM-UQ log2(x+1); per-dataset z-scoring; gseapy 0.10.1 ssGSEA NES.",
     data="Methods text.", generalise="Keep identical settings across GEN cancers; report any change."),
dict(ID="Text-claim-specificity", kind="text", cancer_drug="COAD + BLCA", shows="Claim that biomarkers are cancer/drug specific (cross-test) and independent of clinical factors.",
     numbers="Cross-test P = 0.16 (COAD) and 0.86 (BLCA) (= SuppFig4); no clinical factor P<0.05 (= SuppTables 1-2).",
     data="As SuppFig4, SuppTable1-2.", generalise="Full cross matrix of all GEN biomarkers x cohorts."),
]

COLS = ["ID", "kind", "cancer_drug", "shows", "numbers", "data", "generalise"]
HEAD = {"ID": "ID", "kind": "type", "cancer_drug": "cancer / drug", "shows": "what it shows (plain words)",
        "numbers": "exact numbers reported", "data": "data needed to reproduce",
        "generalise": "how to generalise to other cancers"}


def token(i):
    """Search token for the completeness checker (match it followed by a non-digit)."""
    if not re.match(r"(Supp)?Fig", i):
        return i                      # SourceData-*, Text-*, SuppTableN
    base = i.split("-")[0]
    if "-" in i and re.fullmatch(r"[a-z]", i.split("-")[1]):
        base = re.sub(r"[a-z]$", "", base)   # ranges like SuppFig1a-d -> SuppFig1
    return base


def write_tsv(path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["ID", "type", "cancer_drug", "what_it_shows", "reported_numbers", "data_needed",
                    "how_to_generalise", "search_token", "status"])
        for r in ROWS:
            w.writerow([r[c].replace("\t", " ").replace("\n", " ") for c in COLS] + [token(r["ID"]), "todo"])


def md_escape(s):
    return s.replace("|", "/").replace("\n", " ")


def write_md(path, header, footer):
    lines = [header, "", "| " + " | ".join(HEAD[c] for c in COLS) + " |", "|" + "---|" * len(COLS)]
    for r in ROWS:
        lines.append("| " + " | ".join(md_escape(r[c]) for c in COLS) + " |")
    lines += ["", footer]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    ids = [r["ID"] for r in ROWS]
    assert len(ids) == len(set(ids)), "duplicate IDs"
    write_tsv(os.path.join(OUT, "checklist.tsv"))
    header = open(os.path.join(HERE, "md_header.md")).read().rstrip()
    footer = open(os.path.join(HERE, "md_footer.md")).read().rstrip()
    write_md(os.path.join(OUT, "PAPER_CHECKLIST.md"), header, footer)
    from collections import Counter
    print(len(ROWS), "items", dict(Counter(r["kind"] for r in ROWS)))
