# Paper checklist: Kong et al. 2020, Nature Communications 11:5485

"Network-based machine learning in colorectal and bladder organoid models predicts anti-cancer drug efficacy in patients".
DOI [10.1038/s41467-020-19313-8](https://doi.org/10.1038/s41467-020-19313-8), PMID 33127883, PMC7599252 (open access, CC BY).
Full text was retrieved via PubMed (PMC) and the NCBI PMC open-access copy.

This file lists **every** main figure panel, supplementary figure panel, supplementary table, source-data sheet and
text-only number in the paper. It was built by `code/build_checklist.py`. The same rows are in `checklist.tsv`
(machine-readable, tab-separated).

## 0. How to use this checklist (for the agents that come later)

- **Name every file you make with the checklist ID at the start**, e.g. `Fig2d-left_KM_COAD_5FU.png`,
  `SuppFig8b_targets_BLCA.pdf`, `SuppTable3_isogenic.tsv`. For other cancers add the cancer after the ID,
  e.g. `Fig2d-left_LUAD_cisplatin_KM.png`. The completeness checker searches `nar/` for the `search_token`
  column of `checklist.tsv` (the ID without the `-left`/`-COAD` suffix; ranges like `SuppFig1a-d` become `SuppFig1`). Match the token followed by a non-digit, so `SuppFig1` does not match `SuppFig10`.
- `status` column starts as `todo`. The checker should set `done` / `substitute` / `blocked` with a one-line reason.
- Numbers tagged **[derived]** are not printed in the paper; they were computed from the paper's own source data
  by `code/derive_source_data_facts.py` (output: `tables/source_data_facts.tsv`). Numbers tagged
  **[read from plot, approx]** were read by eye from a plotted axis. Every other number is printed in the paper.

## 1. Sources used (all cached under `nar/_cache/`)

| file | what | sha256 (first 16) |
|---|---|---|
| `41467_2020_Article_19313.pdf` | main article PDF (PMC open-access S3 copy, `PMC7599252.1.pdf`) | 1d4d34cad1b88adc |
| `PMC7599252_jats.xml` | JATS full text from NCBI efetch (figure legends, text with figure refs) | f031c8ab4074db96 |
| `main_text_PMC7599252.txt` | PubMed MCP full text (figure numbers stripped; used only as a cross-check) | - |
| `41467_2020_19313_MOESM1_ESM.pdf` | Supplementary Information: 18 Supp Figs + 4 Supp Tables (28 pages) | baf1c7f4d603eb2c |
| `41467_2020_19313_MOESM2_ESM.pdf` | Reporting summary (no results) | bfc626963bae8f30 |
| `41467_2020_19313_MOESM3_ESM.xlsx` | Source data, 13 sheets | c1e065aa9490694a |

Note: nature.com and media.springernature.com are blocked here. Main-figure numbers were read from the article PDF
pages (rendered to `nar/_cache/main_pages/`). Supplementary figures were extracted as images into `nar/_cache/supp_imgs/`
(img-000 = Supp Fig 1 ... img-017 = Supp Fig 18, img-018..021 = Supp Tables 1-4).

## 2. The paper's key claims in plain language (with the statistics it reports)

1. **Network filter.** Of 674 Reactome pathways, 37 are "close" to the 5-fluorouracil (5FU) target and 30 are close
   to cisplatin targets in the STRING protein network (13,824 proteins, 323,774 links; close = z-score <= -1.2816,
   from 1,000 random gene sets).
2. **Organoid biomarker.** Ridge regression on 19 colorectal organoids picks "activation of BH3-only proteins" as the
   5FU biomarker (rank 1). On bladder organoids (text: 9) it picks "amino acid synthesis and interconversion"
   for cisplatin (it is rank 2 in the bar chart, chosen because it was consistent across 3 algorithms).
   Coefficients agree across ridge, linear and SVR (R2 0.993-0.998).
3. **Patient survival.** In TCGA, predicted 5FU responders live longer than predicted non-responders
   (colon, n = 114, log-rank P = 0.014; responders had no deaths). Predicted cisplatin responders live longer
   (bladder, n = 77, P = 0.01; 5-year survival about 75% vs 50%).
4. **Not just prognosis.** In untreated patients the same split gives P = 0.16 (colon, n = 298) and P = 0.066
   (bladder, n = 294).
5. **Specific to cancer and drug.** Swapping the two biomarkers gives P = 0.16 (colon) and P = 0.86 (bladder).
   Predicted groups do not differ in age, sex, T stage or N stage (all chi-squared P >= 0.093).
6. **Other methods do worse.** Whole transcriptome (P = 0.69 / 0.79), whole pathway (0.81 / 0.82), direct target
   neighbours (0.96 / 0.37), 2nd neighbours (0.24 / 0.79), 3rd neighbours (0.36 / 0.79), Bolis gene filter (0.15 / 0.96),
   deep learning (0.91 / 0.82), network-centrality pathways (all 18 P between 0.091 and 0.997), drug-target
   genes alone (0.18-0.92), and GDSC cell-line-trained models (top pathway P = 0.17, 0.86, 0.102, 0.13)
   (colon / bladder).
7. **Multi-pathway.** Using the top 2 pathways keeps the bladder result (P = 0.022) but not colon (P = 0.88).
8. **Isogenic cells agree.** The BH3 pathway is higher in 5FU-sensitive than resistant HCT-8 cells
   (n = 3 vs 9, t-test P = 0.0085). The amino-acid pathway is higher in cisplatin-sensitive than resistant T24 cells
   (n = 3 vs 3, P = 0.00022).
9. **Random pathways rarely do as well.** 10,000 random pathway sets: empirical P = 0.0012 (colon) and 0.014 (bladder).
   Pathway size does not explain survival P (R2 = 0.0026 and 0.00034).
10. **Agreement with known mutation biomarkers.** A cetuximab model from colon organoids gives higher predicted
    resistance in BRAF V600E tumours (n = 47 vs 328 wild-type, one-sided Mann-Whitney P = 0.037; other feature sets
    P = 0.63-0.99). The cisplatin model gives lower predicted resistance in ERCC2-mutant bladder tumours
    (n = 37 vs 367, P = 0.002; other feature sets P = 0.017-0.04). With 1-3 organoids removed, >70% of subsets still
    recover the BRAF result.
11. **Within-organoid fit.** A 60/10/30 split of the organoids gives observed-vs-predicted R2 = 0.98 (colon) and
    0.89 (bladder).

## 3. Generalisation targets (used in the last column)

"GEN" in the table means the set of new cancer-drug pairs below. Patient counts are from the project brief
(TCGA drug tables) unless marked "to count".

| cancer | TCGA project, treated n | lab model with drug response (repo: `data/curated/preclinical/`) | isogenic pair idea (check GEO) | known mutation biomarker for a Fig5-style test |
|---|---|---|---|---|
| colorectal (external) | GSE39582: 5FU 82, FUFOL 54, FOLFOX 23, untreated 326 | van de Wetering 2015 organoids; Isella 2017 PDX | HCT-8/5FU (GSE81008, as paper) | BRAF V600E / KRAS vs cetuximab |
| bladder (other drugs) | BLCA cisplatin 77, methotrexate 12, gemcitabine to count | Lee 2018 organoids | cisplatin-resistant bladder lines | ERCC2 mut vs cisplatin |
| lung adenocarcinoma | LUAD cisplatin 87, carboplatin 77, pemetrexed 60, gemcitabine 25 | GDSC/CTRP lung lines; GSE165611 patient-derived cells (drug table not in repo) | A549 vs A549/DDP | EGFR mut vs EGFR TKI; KEAP1/STK11 vs platinum |
| ovarian | OV carboplatin 464, paclitaxel 448, cisplatin 156, olaparib 0 | Vias 2023 organoids; GDSC ovary lines | PEO1/PEO4, A2780/A2780cis | BRCA1/2 mut vs platinum |
| pancreas | PAAD gemcitabine (to count) | Tiriac 2018, Shi 2022 organoids | gemcitabine-resistant PANC-1 / MIA PaCa-2 | BRCA/ATM (HRD) vs platinum |
| liver | LIHC sorafenib (to count, likely small) | Broutier 2017, Ji 2023 organoids | sorafenib-resistant Huh7 / HepG2 | none established; use external sorafenib-response cohort if one exists in GEO |
| pan-cancer | all of the above | CTRPv2/CCLE, PRISM, Gao 2015 PDX | - | - |

## 4. Checklist table

| ID | type | cancer / drug | what it shows (plain words) | exact numbers reported | data needed to reproduce | how to generalise to other cancers |
|---|---|---|---|---|---|---|
| Fig1a | main-figure | all (schematic) | Cartoon of the whole method: pathways near drug targets in a protein network -> train ML on organoid drug response -> pick best pathway -> split patients into predicted responders/non-responders. | No numbers (schematic). | None (drawing). | One schematic for the NAR paper showing the same 4 steps with a box listing all cancers/models used (organoid, PDC, cell line, PDX). |
| Fig1b | main-figure | all (schematic) | Four-step detail: (1) network distance drug targets->pathway genes; (2-3) pathway ssGSEA scores + IC50 train regression, coefficient = predictive performance; (4) biomarker expression splits patients, compare overall survival. | No numbers (schematic). | None (drawing). | Same; add step for external validation cohorts (GEO trials) and for cross-model transfer. |
| Fig2a | main-figure | COAD / 5-fluorouracil | Each dot is a Reactome pathway: x = average shortest network distance from the 5FU target (TYMS) to the pathway, y = proximity z-score. Pathways with z <= -1.2816 (shaded) are 'proximal'. | 674 Reactome pathways tested; 37 proximal pathways (z <= -1.2816, i.e. alpha=0.10). [derived] 38 pathways pass the z cutoff in the source matrix; 37 of them have organoid ssGSEA scores. Biomarker 'activation of BH3-only proteins' z = -1.2857 [derived] (just past the cutoff). | STRING edges (version not stated) score>700, largest component (13,824 proteins, 323,774 links; MOESM3 sheet 'Figure 2a-f'); DrugBank targets (5FU -> TYMS); MSigDB C2 Reactome; proximity z-scores (MOESM3 sheet 'Figure 2a,d', column FLUOROURACIL). 1000 degree-matched random draws. | Proximity scatter for each GEN drug: cisplatin/carboplatin (LUAD, OV), pemetrexed (LUAD), paclitaxel (OV), gemcitabine (PAAD, BLCA), sorafenib (LIHC), olaparib (OV). Platinum targets are DrugBank 'targets' that are not the real DNA mechanism - note this. |
| Fig2b | main-figure | COAD / 5-fluorouracil | Bar chart of ridge-regression /coefficient/ ('predictive ability') for the 37 proximal pathways, trained on 19 colorectal organoids' 5FU IC50. The chosen biomarker (orange) is rank 1. | 37 pathways; rank 1 = 'Activation of BH3-only proteins', /coef/ ~0.8 [read from plot, approx]; next pathways ~0.6. RidgeCV alpha 0.1-1.0 step 0.1, 3-fold CV. | COAD organoid expression GSE64392 (RMA) -> ssGSEA NES (MOESM3 'Figure 2b', 612 pathways x 22 columns [derived]); 5FU IC50 from van de Wetering 2015 supplement (data/curated/preclinical colorectal_vandewetering2015). | Ranked-coefficient bar chart per model/drug: OV organoids (Vias 2023) carboplatin/paclitaxel; PAAD organoids (Tiriac 2018) gemcitabine; LIHC organoids (Broutier 2017/Ji 2023) sorafenib; LUAD from GDSC/CTRP lung lines or GSE165611 PDCs; BLCA other drugs from Lee 2018. |
| Fig2c | main-figure | COAD / 5-fluorouracil | Network drawing: 5FU -> TYMS (red triangle) linked via TFDP1/TP53 to the genes of the BH3-only pathway (E2F1, TFDP1, TP53, AKT1, PMAIP1, DYNLL1, YWHAB, MAPK8, BBC3, BID, BCL2L11, BCL2, BAD, PPP3R1, BMF). | No statistics; gene list as drawn. | STRING edges (score>700) among TYMS + pathway genes; Reactome gene set REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS. | Subnetwork plot (drug target + chosen biomarker pathway genes + shortest-path links) for each new cancer-drug biomarker. |
| Fig2d-left | main-figure | COAD / 5-fluorouracil (TCGA treated) | Kaplan-Meier overall survival of TCGA-COAD patients who got 5FU, split at the median predicted resistance score (biomarker NES x organoid coefficient). | n=114 (57 predicted responders / 57 non-responders); two-sided log-rank P = 0.014. Responders: no deaths; non-responders 5-year OS ~50%. | TCGA-COAD FPKM-UQ log2(x+1) -> ssGSEA NES (MOESM3 'Figure 2g', 674 pathways x 449 patients [derived]); TCGA drug table (5FU-treated, DrugBank-normalised names); overall survival. | Same KM split for each GEN cancer-drug pair: LUAD cisplatin/carboplatin/pemetrexed; OV carboplatin/paclitaxel; BLCA gemcitabine/methotrexate; PAAD gemcitabine; LIHC sorafenib; COAD-external GSE39582 5FU. Brief says only ~5 TCGA-COAD got 5FU alone; most got FOLFOX - say which definition is used. |
| Fig2d-right | main-figure | COAD / untreated control (TCGA) | Same split applied to TCGA-COAD patients with no recorded drug treatment (negative control: is the biomarker just prognostic?). | n=298 (149/149); log-rank P = 0.16. | Same as Fig2d-left, untreated patients. | Untreated-patient control KM for every GEN cancer (TCGA untreated of same cancer; GSE39582 has 326 untreated). Better: formal treatment x biomarker interaction Cox test. |
| Fig2e | main-figure | BLCA / cisplatin | Proximity scatter for cisplatin targets vs Reactome pathways (as Fig2a). | 30 proximal pathways. [derived] 47 pathways pass z <= -1.2816 in the source matrix; 30 of them have BLCA organoid ssGSEA scores. Biomarker 'amino acid synthesis and interconversion (transamination)' z = -2.115 [derived]. | As Fig2a, column CISPLATIN; cisplatin DrugBank targets MPG, A2M, ATOX1, TF. | As Fig2a for each GEN drug. |
| Fig2f | main-figure | BLCA / cisplatin | Ridge /coefficient/ bar chart for 30 proximal pathways trained on bladder organoid cisplatin IC50. Chosen biomarker (orange) is rank 2 (rank 1 is a grey bar). | 30 pathways; biomarker /coef/ ~0.2, rank-1 pathway ~0.3 [read from plot, approx]. Biomarker is rank 2, not rank 1 (same in Supp Fig 2c,d). | BLCA organoid expression + cisplatin IC50 from Lee et al. 2018 (MOESM3 'Figure 2e': 488 pathways x 11 organoid columns [derived]; text says 9 organoids; Supp Fig 17 shows 9 with IC50). | As Fig2b for each GEN model/drug. |
| Fig2g | main-figure | BLCA / cisplatin | Network drawing: cisplatin targets MPG, A2M, ATOX1, TF linked to amino-acid pathway genes (KYAT1, ASNS, OAT, GLUD1, GOT2, PSAT1, PSPH, PYCR1, PHGDH, GLUL, GLS, GPT, GOT1, GPT2). | No statistics. | STRING edges among targets + pathway genes. NOTE: MOESM3 sheet 'Figure 2g' actually holds TCGA-COAD ssGSEA NES, not this network. | As Fig2c. |
| Fig2h-left | main-figure | BLCA / cisplatin (TCGA treated) | KM overall survival of cisplatin-treated TCGA-BLCA patients split at median predicted resistance. | n=77 (39 responders / 38 non-responders); log-rank P = 0.01. 5-year OS ~75% vs ~50%. | TCGA-BLCA FPKM-UQ -> ssGSEA (MOESM3 'Figure 2h', 490 pathways x 407 patients [derived]); cisplatin-treated patients; OS. | Same KM split for each GEN cancer-drug pair: LUAD cisplatin/carboplatin/pemetrexed; OV carboplatin/paclitaxel; BLCA gemcitabine/methotrexate; PAAD gemcitabine; LIHC sorafenib; COAD-external GSE39582 5FU. |
| Fig2h-right | main-figure | BLCA / untreated control (TCGA) | Same split in TCGA-BLCA patients with no drug treatment (negative control). | n=294 (147/147); log-rank P = 0.066. | As Fig2h-left, untreated patients. | Untreated control + interaction test for each GEN cancer. |
| Fig3a | main-figure | COAD 5FU + BLCA cisplatin (schematic) | Scheme of comparison methods: (1) no feature selection, (2) other feature selection, both -> ML -> predict TCGA responders. | No numbers. | None. | Same scheme, add modern baselines (elastic net on all genes, random pathways, published signatures). |
| Fig3b-COAD | main-figure | COAD / 5FU | KM: ridge model on whole transcriptome (all genes) of organoids, applied to 5FU-treated TCGA-COAD. | n=57/57; log-rank P = 0.69. | Organoid gene-level expression + IC50; TCGA gene expression; OS. | Whole-transcriptome baseline for every GEN cancer-drug pair. |
| Fig3b-BLCA | main-figure | BLCA / cisplatin | KM: whole-transcriptome ridge model, cisplatin-treated TCGA-BLCA. | n=39/38; log-rank P = 0.79. | As above for BLCA. | As Fig3b-COAD. |
| Fig3c-COAD | main-figure | COAD / 5FU | KM: ridge model on all pathways (whole-pathway ssGSEA). | n=57/57; log-rank P = 0.81. | Organoid + TCGA ssGSEA of all Reactome pathways. | Whole-pathway baseline for every GEN pair. |
| Fig3c-BLCA | main-figure | BLCA / cisplatin | KM: whole-pathway ridge model. | n=39/38; log-rank P = 0.82. | As above. | As Fig3c-COAD. |
| Fig3d-COAD | main-figure | COAD / 5FU | KM: ridge model on genes that are direct (first) network neighbours of the drug target. | n=57/57; log-rank P = 0.96. | STRING first neighbours of TYMS; organoid + TCGA gene expression. | First-neighbour baseline for each GEN drug. |
| Fig3d-BLCA | main-figure | BLCA / cisplatin | KM: first-neighbour gene model. | n=39/38; log-rank P = 0.37. | STRING first neighbours of MPG/A2M/ATOX1/TF. | As Fig3d-COAD. |
| Fig3e-COAD | main-figure | COAD / 5FU | KM: genes pre-selected by Bolis et al. method (Spearman P<0.05 with IC50 over 10 leave-half-out repeats), then ridge. | n=57/57; log-rank P = 0.15. 15 genes selected [derived from MOESM3 'Figure 3e']. | Organoid gene expression + IC50; MOESM3 'Figure 3e' gene lists. | Correlation-filter baseline for each GEN pair. |
| Fig3e-BLCA | main-figure | BLCA / cisplatin | KM: Bolis feature selection + ridge. | n=39/38; log-rank P = 0.96. 8 genes selected [derived]. | As above. | As Fig3e-COAD. |
| Fig3f-COAD | main-figure | COAD / 5FU | KM: deep neural network of Sharifi-Noghabi et al. (MOLI-style, triplet + cross-entropy loss), organoids split at median IC50. | n=57/57; log-rank P = 0.91. Hyper-parameters + CV AUC in MOESM3 'Figure 3f'. | Organoid + TCGA expression; MOESM3 'Figure 3f'. | Optional DL baseline (CPU only; torch from PyPI) for GEN pairs with >=30 training models. |
| Fig3f-BLCA | main-figure | BLCA / cisplatin | KM: deep-learning model. | n=39/38; log-rank P = 0.82. | As above. | As Fig3f-COAD. |
| Fig4a | main-figure | isogenic lines (schematic) | Scheme: drug-sensitive parent cells made resistant by rising drug doses; compare transcriptomes. | No numbers. | None. | Same scheme; list isogenic pairs per GEN cancer. |
| Fig4b | main-figure | COAD / 5FU (HCT-8 vs HCT-8/5FU) | Boxplot of BH3-only pathway ssGSEA score in 5FU-sensitive vs 5FU-resistant colorectal cells. | Sensitive n=3, resistant n=9 (resistant = 0h, 24h, 48h triplicates); two-sided two-sample t-test P = 0.0085; sensitive higher. | GEO GSE81008 microarray -> ssGSEA (MOESM3 'Figure 4b', 674 pathways x 12 samples). | Isogenic sensitive/resistant pair per GEN drug from GEO: e.g. A549 vs A549/DDP (LUAD cisplatin), A2780 vs A2780cis or PEO1 vs PEO4 (OV platinum), gemcitabine-resistant PANC-1/MIA PaCa-2 (PAAD), sorafenib-resistant Huh7/HepG2 (LIHC). Accessions must be checked, not assumed. |
| Fig4c | main-figure | BLCA / cisplatin (T24S vs T24R) | Boxplot of amino-acid pathway score in cisplatin-sensitive vs -resistant T24 bladder cells. | Sensitive n=3, resistant n=3; t-test P = 0.00022; resistant lower. | Yeon et al. data by personal request (not public); MOESM3 'Figure 4c' gives ssGSEA NES for 6 samples. | As Fig4b. For BLCA a public cisplatin-resistant bladder line set from GEO is the substitute for raw data. |
| Fig4d | main-figure | bootstrap (schematic) | Scheme of 10,000-iteration random-pathway bootstrap. | 10,000 iterations; random pathway sets matched to number of proximal pathways; 'top rank' = best ranked pathway that is both survival-predictive (P<0.05) and isogenic-predictive (P<0.05). | None. | Same null for each GEN pair. |
| Fig4e | main-figure | COAD / 5FU | Histogram of the top predictive rank reached by random pathway sets; observed rank from proximal pathways in green. | Empirical P = 0.0012. [derived] observed rank = 1 (12/10,000 random sets reach rank 1); 8,927/10,000 random sets give NA (no pathway passes both tests). | MOESM3 'Figure 4e' (10,000 rows: top_rank, feature, coefficient); requires full pipeline + TCGA + GSE81008. | Random-pathway empirical P for each GEN pair (needs an isogenic or second validation set to define 'predictive'; otherwise survival-only null). |
| Fig4f | main-figure | BLCA / cisplatin | Same bootstrap histogram for cisplatin/BLCA. | Empirical P = 0.014. [derived] observed rank = 2 (66 sets rank 1 + 74 rank 2 = 140/10,000); 6,366/10,000 NA. | MOESM3 'Figure 4f'. | As Fig4e. |
| Fig5a | main-figure | COAD / cetuximab vs BRAF V600E (TCGA) | Boxplot of predicted cetuximab resistance (Gastrin-CREB via PKC and MAPK pathway, trained on organoids) in BRAF wild-type vs BRAF V600E TCGA-COAD tumours. | WT n=328, V600E n=47; one-sided Mann-Whitney U P = 0.037 (V600E higher). Trained after removing 3 hypermutated organoids (P19a, P19b, P24a). | Organoid cetuximab IC50 (proximity column CETUXIMAB_GEFITINIB); TCGA-COAD expression + MAF mutations (BRAF). | Mutation-concordance test per GEN cancer: OV BRCA1/2 mut vs platinum (expect lower resistance); LUAD EGFR mut vs EGFR TKI (erlotinib/gefitinib) or KEAP1/STK11 vs platinum; BLCA ERCC2 (already); PAAD BRCA/ATM vs platinum; COAD KRAS mut vs cetuximab. |
| Fig5b | main-figure | BLCA / cisplatin vs ERCC2 (TCGA) | Boxplot of predicted cisplatin resistance (amino-acid pathway) in ERCC2-mutant vs wild-type TCGA-BLCA tumours. | ERCC2 mut n=37, WT n=367; one-sided Mann-Whitney U P = 0.002 (mutants lower = more sensitive). | TCGA-BLCA expression + mutations (all patients, not only treated). | As Fig5a. |
| SuppFig1a-d | supp-figure | COAD + BLCA (organoid vs TCGA) | PCA of organoids (red) and TCGA patients (grey) before scaling: genes and pathways. Organoids sit apart from patients (batch effect). | PC1/PC2 explained variance: (a) COAD genes 39.8%/5.05%; (b) COAD pathways 47.3%/13.6%; (c) BLCA genes 21.6%/6.37%; (d) BLCA pathways 27.0%/17.8%. | Organoid + TCGA gene expression and ssGSEA NES. | Batch-effect PCA for every GEN model-vs-patient pairing (organoid/PDC/cell line vs TCGA or GEO cohort). |
| SuppFig1e-h | supp-figure | COAD + BLCA | Same PCA after z-scoring each gene/pathway within each dataset: organoids now mix with patients. | (e) COAD scaled genes 12.2%/8.93%; (f) COAD scaled pathways 27.1%/10.7%; (g) BLCA scaled genes 9.50%/6.87%; (h) BLCA scaled pathways 20.7%/18.8%. | As above, z-scored per dataset. | As SuppFig1a-d. |
| SuppFig2a | supp-figure | COAD / 5FU | Predictive performance ranking with linear regression (37 pathways). | BH3-only pathway rank 1, ~0.88 [read from plot, approx]. | As Fig2b. | Algorithm-robustness ranking per GEN pair. |
| SuppFig2b | supp-figure | COAD / 5FU | Same with support vector regression (linear kernel). | BH3-only rank 1, ~0.82 [approx]. | As Fig2b. | As SuppFig2a. |
| SuppFig2c | supp-figure | BLCA / cisplatin | Linear regression ranking (30 pathways). | Amino-acid pathway rank 2, ~0.23; rank 1 ~0.34 [approx]. | As Fig2f. | As SuppFig2a. |
| SuppFig2d | supp-figure | BLCA / cisplatin | SVR ranking. | Amino-acid pathway rank 2, ~0.20; rank 1 ~0.29 [approx]. | As Fig2f. | As SuppFig2a. |
| SuppFig3a-c | supp-figure | COAD / 5FU | Scatter of regression coefficients between algorithms (37 pathways). | R2: ridge vs linear 0.9948; SVR vs linear 0.9979; SVR vs ridge 0.9973. | Coefficients from 3 models on COAD organoids. | Coefficient agreement plot per GEN pair. |
| SuppFig3d-f | supp-figure | BLCA / cisplatin | Same for bladder (30 pathways). Legend says (d)-(e) but figure has d, e, f. | R2: ridge vs linear 0.9978; SVR vs linear 0.9932; SVR vs ridge 0.9977. | As above, BLCA. | As SuppFig3a-c. |
| SuppFig4a | supp-figure | COAD / 5FU with BLCA biomarker | Cross-test: the bladder cisplatin biomarker (amino-acid pathway) used on 5FU-treated COAD patients. | n=57/57; log-rank P = 0.16. | TCGA-COAD NES; BLCA organoid coefficient. | Cross-cancer/cross-drug specificity matrix: every biomarker x every treated cohort (heatmap of log-rank P / HR). |
| SuppFig4b | supp-figure | BLCA / cisplatin with COAD biomarker | Cross-test: BH3-only pathway used on cisplatin-treated BLCA patients. | n=39/38; log-rank P = 0.86. | TCGA-BLCA NES; COAD organoid coefficient. | As SuppFig4a. |
| SuppFig5a | supp-figure | COAD / 5FU | KM using weighted sum of top-2 pathways. | n=57/57; log-rank P = 0.88. (Main text calls this 5b; legend/figure say 5a.) | Top-2 pathways x ridge coefficients; TCGA-COAD. | Top-k (1..10) multi-pathway curve per GEN pair. |
| SuppFig5b | supp-figure | BLCA / cisplatin | KM using top-2 pathways. | n=39/38; log-rank P = 0.022. (Main text calls this 5a.) | As above, BLCA. | As SuppFig5a. |
| SourceData-multipathway | source-data | COAD 5FU + BLCA cisplatin | Table of top-k pathway models (k=2..10) for Ridge, SVR and LinearRegression: 5-year survival in each group and log-rank P (MOESM3 sheet 'FigureS5'). | Ridge COAD P for k=2..10: 0.882, 0.373, 0.314, 0.130, 0.253, 0.076, 0.216, 0.314, 0.314. Ridge BLCA: 0.022, 0.496, 0.161, 0.480, 0.709, 0.342, 0.639, 0.126, 0.256. Smallest non-k2 values: COAD 0.076 (k=7), BLCA LinearRegression k=10 0.0625. | MOESM3 'FigureS5'; original code python/multiple_pathway_prediction.py produces it. | Same k-curve per GEN pair. |
| SuppFig6a | supp-figure | COAD / 5FU | Degree-centrality feature selection (37 most central pathways): ranking + KM of top 3 pathways. | Top 3: NFkB activated and signals survival P=0.997; SPRY regulation of FGF signaling P=0.71; IRAK2-mediated TAK1 activation P=0.802 (n=57/57). | STRING degree centrality averaged over pathway genes; organoid + TCGA ssGSEA. | Centrality-selected baseline for each GEN pair. |
| SuppFig6b | supp-figure | COAD / 5FU | Betweenness-centrality selection. | Top 3: Membrane binding and targeting of GAG proteins P=0.79; Regulation of Fanconi anemia pathway P=0.503; Hormone-sensitive lipase P=0.91. | As above (betweenness). | As SuppFig6a. |
| SuppFig6c | supp-figure | COAD / 5FU | Closeness-centrality selection. | Top 3: Signaling by constitutively active EGFR P=0.17; SPRY regulation of FGF signaling P=0.71; CTLA4 inhibitory signaling P=0.32. | As above (closeness). | As SuppFig6a. |
| SuppFig7a | supp-figure | BLCA / cisplatin | Degree-centrality selection (30 pathways) + KM of top 3. | G alpha z signaling events P=0.094; mRNA capping P=0.96; GABA B receptor activation P=0.51 (n=39/38). | As SuppFig6. | As SuppFig6a. |
| SuppFig7b | supp-figure | BLCA / cisplatin | Betweenness selection. | Pre-NOTCH transcription and translation P=0.23; G1 phase P=0.31; TGF-beta receptor signaling activates SMADs P=0.091. | As SuppFig6. | As SuppFig6a. |
| SuppFig7c | supp-figure | BLCA / cisplatin | Closeness selection. | Signaling by WNT P=0.56; EGFR downregulation P=0.29; CTNNB1 phosphorylation cascade P=0.74. | As SuppFig6. | As SuppFig6a. |
| SuppFig8a | supp-figure | COAD / 5FU | KM split by expression of the drug target TYMS alone (median). | high n=57, low n=57; log-rank P = 0.73. | TCGA-COAD gene expression. | Drug-target-expression KM per GEN drug (e.g. TUBB for paclitaxel, TYMS/RRM1 for pemetrexed/gemcitabine, PARP1 for olaparib). |
| SuppFig8b | supp-figure | BLCA / cisplatin | KM split by each cisplatin target gene. | high n=38 / low n=39: A2M P=0.92; ATOX1 P=0.18; MPG P=0.83; TF P=0.48. | TCGA-BLCA gene expression. | As SuppFig8a. |
| SuppFig9a | supp-figure | COAD 5FU + BLCA cisplatin | KM with ridge model on genes up to 2nd network neighbours of drug targets. | COAD P = 0.24 (57/57); BLCA P = 0.79 (39/38). | STRING 2-hop neighbourhood; organoid + TCGA gene expression. | k-hop neighbour baselines per GEN pair. |
| SuppFig9b | supp-figure | COAD 5FU + BLCA cisplatin | Same with up to 3rd neighbours. | COAD P = 0.36; BLCA P = 0.79. | STRING 3-hop. | As SuppFig9a. |
| SuppFig10a | supp-figure | COAD / 5FU | Scatter: pathway size (number of genes) vs log-rank P when each pathway alone splits 5FU-treated patients (median). | R2 = 0.0026. | All pathways' TCGA-COAD NES + gene-set sizes. | Pathway-size confounding check per GEN pair. |
| SuppFig10b | supp-figure | BLCA / cisplatin | Same for cisplatin-treated BLCA. | R2 = 0.00034. | As above, BLCA. | As SuppFig10a. |
| SuppFig11 | supp-figure | COAD / cetuximab (organoids) | Heatmap of 37 cetuximab-proximal pathways across 19 colorectal organoids sorted by cetuximab IC50, with ridge weights. Top weight: Gastrin-CREB signalling via PKC and MAPK. | 19 organoids (P18 most sensitive ... P31 most resistant; IC50 ~2.4 to ~8); 37 pathways; Gastrin-CREB weight ~0.75 [approx]; most negative 'role of DCC in regulating apoptosis' ~-0.4 [approx]. | COAD organoid NES + cetuximab IC50; proximity column CETUXIMAB_GEFITINIB. | Organoid heatmap + weights for each GEN model/drug (e.g. OV organoid carboplatin, PAAD organoid gemcitabine). |
| SuppFig12a-e | supp-figure | COAD / cetuximab vs BRAF V600E | Same BRAF test as Fig5a using other feature sets: (a) whole transcriptome, (b) whole pathway, (c) direct neighbours, (d) 2nd neighbours, (e) 3rd neighbours. | WT n=328, V600E n=47; one-sided MWU P: (a) 0.98, (b) 0.63, (c) 0.99, (d) 0.99, (e) 0.99. | As Fig5a. | Baseline comparison for each GEN mutation-concordance test. |
| SuppFig13a-e | supp-figure | BLCA / cisplatin vs ERCC2 | Same ERCC2 test as Fig5b with other feature sets (a-e as above). | ERCC2 mut n=37, WT n=367; one-sided MWU P: (a) 0.04, (b) 0.037, (c) 0.017, (d) 0.021, (e) 0.028. (All < 0.05; paper's own model P = 0.002.) | As Fig5b. | As SuppFig12a-e. |
| SuppFig14a-h | supp-figure | COAD + BLCA (GDSC cell lines vs TCGA) | PCA of GDSC cancer cell lines (red) vs TCGA patients (grey), gene and pathway level, unscaled (all lines and tissue-matched lines). | Explained-variance values printed on axes (not transcribed). | GDSC/Cell Model Passports expression; TCGA. | Batch PCA for every cell-line model used (GDSC/CTRP lung, ovary, pancreas, liver lines). |
| SuppFig14i-p | supp-figure | COAD + BLCA (GDSC vs TCGA) | Same PCA after per-dataset z-scoring. | Values on axes (not transcribed). | As above. | As SuppFig14a-h. |
| SuppFig15a | supp-figure | COAD / 5FU from all GDSC lines | Pathway ranking trained on all 893 GDSC lines' 5FU response + KM of top 3 pathways in 5FU-treated TCGA-COAD. | n=893 lines; Cell cycle P=0.17; CDC6 association with ORC origin complex P=0.58; Association of licensing factors with pre-replicative complex P=0.59 (n=57/57). | GDSC 5FU IC50 + expression; TCGA-COAD. | Cell-line-trained comparator for each GEN pair (all lines vs tissue-matched lines). |
| SuppFig15b | supp-figure | COAD / 5FU from colorectal GDSC lines | Same with 46 colorectal lines. | n=46; Metabolism of nucleotides P=0.86; Cell cycle mitotic P=0.17; E2F-enabled inhibition of pre-replication complex formation P=0.51. | GDSC colorectal lines. | As SuppFig15a. |
| SuppFig15c | supp-figure | BLCA / cisplatin from all GDSC lines | Ranking from 831 lines' cisplatin response + KM of top 3 in cisplatin-treated TCGA-BLCA. | n=831; Transport of glucose and other sugars P=0.102; Transmembrane transport of small molecules P=0.0092; Platelet activation signaling and aggregation P=0.84 (n=39/38). | GDSC cisplatin; TCGA-BLCA. | As SuppFig15a. |
| SuppFig15d | supp-figure | BLCA / cisplatin from bladder GDSC lines | Same with 14 bladder lines. | n=14; BMAL1:CLOCK/NPAS2 activates circadian expression P=0.13; Circadian clock P=0.23; Transmembrane transport of small molecules P=0.0092. | GDSC bladder lines. | As SuppFig15a. |
| SuppFig16 | supp-figure | 44 organoid drugs (COAD + BLCA screens) | For each Reactome pathway, how many of the 44 drugs it is proximal to (checks that no pathway is proximal to everything). | 674 pathways, 44 drugs; max 31 drugs [derived]; only 3 pathways (0.45%) proximal to >=22 drugs; mean 4.9 (derived 4.86), median 4 drugs per pathway; 111 pathways proximal to no drug. | MOESM3 'Figure 2a,d' (674 x 44 z-scores; some columns merge drugs with identical targets, e.g. CETUXIMAB_GEFITINIB). | Same hub-bias check for the full GEN drug list (all drugs used across cancers). |
| SuppFig17a | supp-figure | COAD / 5FU (organoids) | Within-organoid prediction: 60/10/30 train/validation/test split; observed vs predicted IC50 in test organoids over all validation choices. | R2 = 0.98 (test organoids P11, P26, P8, P18, P14, P16 in green). | COAD organoid NES of proximal pathways + 5FU IC50. | Proper repeated nested CV on each GEN preclinical set (report Spearman/Pearson + permutation P, not one split). |
| SuppFig17b | supp-figure | BLCA / cisplatin (organoids) | Same for 9 bladder organoids. | R2 = 0.89 (test organoids SCBO-1, SCBO-11.2, SCBO-4). | BLCA organoid NES + cisplatin IC50. | As SuppFig17a. |
| SuppFig18 | supp-figure | COAD / cetuximab vs BRAF V600E | Robustness: fraction of organoid subsets (drop 0-3 organoids) that still recover the BRAF V600E-cetuximab association (one-sided MWU P<0.05). | Combinations tested: 969 (16 organoids), 171 (17), 19 (18), 1 (19). Network method ~0.73/0.80/0.89/1.0; whole-transcriptome ~0.64/0.71/0.79/1.0; whole-pathway ~0.40/0.38/0.31/0.0 [read from plot, approx]. Text: >70% recover with 3 organoids removed. | As Fig5a, leave-k-out over organoids. | Leave-k-models-out robustness curve for each GEN concordance/survival result. |
| SuppTable1 | supp-table | COAD / 5FU (TCGA) | Clinical features of predicted responders vs non-responders (chi-squared tests). | n=114. Age <60: 53 (30 R/23 NR), >=60: 61 (27/34), P=0.26. Male 62 (27/35), female 52 (30/22), P=0.19. T1-T2 7 (2/5), T3-T4 107 (55/52), P=0.44. N0 27 (16/11), N1 52 (26/26), N2 35 (15/20), P=0.44. | TCGA-COAD clinical (age, gender, T, N) + group labels. | Clinical-balance table per GEN cohort; better: multivariable Cox with stage/age (and MMR/KRAS/BRAF for CRC). |
| SuppTable2 | supp-table | BLCA / cisplatin (TCGA) | Clinical features of predicted groups. | n=77. Age <60: 21 (11/10), >=60: 56 (28/28), P=0.94. Male 53 (27/26), female 24 (12/12), P=0.87. T1-T2 12 (6/6), T3-T4 65 (33/32), P=0.79. N0 39 (22/17), N1 13 (3/10), N2 25 (14/11), P=0.093. | TCGA-BLCA clinical. | As SuppTable1. |
| SuppTable3 | supp-table | COAD 5FU + BLCA cisplatin (isogenic lines) | Top-3 most proximal pathways: proximity z and mean ssGSEA in sensitive vs resistant isogenic lines. | COAD: MITOTIC_G1_G1_S_PHASES z=-4.26, 0.34 vs 0.31, P=0.15; G1_S_TRANSITION z=-3.9, 0.38 vs 0.34, P=0.15; PYRIMIDINE_METABOLISM z=-3.42, 0.17 vs 0.16, P=0.11. BLCA: METAL_ION_SLC_TRANSPORTERS z=-4.03, 0.078 vs 0.06, P=0.014; RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 z=-3.75, 0.14 vs 0.12, P=0.00054; RECYCLING_OF_BILE_ACIDS_AND_SALTS z=-3.73, -0.28 vs -0.22, P=0.01. | Proximity z (MOESM3 'Figure 2a,d') + isogenic NES (MOESM3 'Figure 4b','Figure 4c'). | Same table for each GEN isogenic pair. |
| SuppTable4 | supp-table | COAD 5FU + BLCA cisplatin (bootstrap) | Random-bootstrap pathways that ranked equal to or better than the real biomarker, with their proximity z. | COAD: PI3K_EVENTS_IN_ERBB2_SIGNALING z=1.4. BLCA: REGULATION_OF_KIT_SIGNALING -1.51; NOTCH1_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION -0.47; SIGNALING_BY_NOTCH1 -0.16; ADAPTIVE_IMMUNE_SYSTEM 1.63; IL_2_SIGNALING 2.09; CA_DEPENDENT_EVENTS 2.94. | MOESM3 'Figure 4e','Figure 4f' + proximity matrix. | Same table per GEN pair. |
| SourceData-STRING | source-data | all | STRING network edge list used for proximity (MOESM3 'Figure 2a-f'). | Text: 13,824 proteins, 323,774 interactions (score>700, largest component). [derived] sheet has 647,548 rows = each link listed in both directions; 13,825 distinct IDs across the two UniProt columns. | MOESM3 'Figure 2a-f'. | Re-use same network for every GEN cancer; sensitivity check with a second network (e.g. STRING 900, or HuRI/BioGRID). |
| SourceData-DLhyper | source-data | COAD 5FU + BLCA cisplatin | Deep-learning hyper-parameter grid and mean CV AUC (MOESM3 'Figure 3f'). | Grid columns: batch_size, #nodes, learning rates, dropouts, weight_decay, gamma, #epoch, #Folds, margin, average_AUC (highlighted rows used). | MOESM3 'Figure 3f'. | Only if DL baseline is repeated. |
| Text-method-params | text | all | Key method settings stated in Methods. | Reactome from MSigDB C2 (674 pathways usable); proximity = closest distance d_c, z from 1000 degree-matched random sets, proximal if z <= -1.2816; ridge RidgeCV alpha 0.1..1.0 step 0.1, 3-fold CV; SVR linear kernel; score = pathway NES x organoid coefficient; median split; lifelines log-rank; TCGA FPKM-UQ log2(x+1); per-dataset z-scoring; gseapy 0.10.1 ssGSEA NES. | Methods text. | Keep identical settings across GEN cancers; report any change. |
| Text-claim-specificity | text | COAD + BLCA | Claim that biomarkers are cancer/drug specific (cross-test) and independent of clinical factors. | Cross-test P = 0.16 (COAD) and 0.86 (BLCA) (= SuppFig4); no clinical factor P<0.05 (= SuppTables 1-2). | As SuppFig4, SuppTable1-2. | Full cross matrix of all GEN biomarkers x cohorts. |

## 5. Things in the paper that do not match each other (worth knowing before reproducing)

1. **Supp Fig 5 labels are swapped in the main text.** Text: "bladder P = 0.022, Supp Fig 5a; colon P = 0.88, Supp Fig 5b".
   Legend and figure: a = colon (0.88), b = bladder (0.022). Source data agrees with the figure.
2. **Supp Fig 3 legend** says bladder panels are "(d)-(e)" but the figure has d, e and f.
3. **The bladder biomarker is rank 2, not rank 1**, in ridge (Fig 2f), linear and SVR (Supp Fig 2c,d). The bootstrap
   (Fig 4f) also places the observed result at rank 2 [derived: 140 of 10,000 random sets reach rank <= 2 -> P = 0.014].
4. **The colon biomarker only just passes the proximity cut-off**: z = -1.2857 vs cut-off -1.2816 [derived].
5. **Proximal-pathway counts depend on a hidden filter.** The source z-score matrix gives 38 (5FU) and 47 (cisplatin)
   pathways with z <= -1.2816; only 37 and 30 of these have organoid ssGSEA scores [derived]. So "37" and "30" are
   after removing pathways missing from the organoid data (the bladder organoid sheet has only 488 of 674 pathways).
6. **Bladder organoid count.** Text says 9 bladder organoids; source sheet 'Figure 2e' has 11 sample columns
   (some patients have two lines, e.g. SCBO-3 and SCBO-3.2); Supp Fig 17b shows 9 with cisplatin IC50.
   Colon sheet 'Figure 2b' has 22 columns for 19 organoids.
7. **Source-data sheet names do not match panels.** 'Figure 2g' and 'Figure 2h' hold TCGA COAD/BLCA ssGSEA scores
   (inputs to Fig 2d/2h KM), not the network drawings. Sheet 'FigureS5' is titled "Figure 17" inside.
8. **Fig 5a** model was trained on 16 organoids (3 hypermutated removed: P19a, P19b, P24a), while Supp Fig 11 shows 19.
9. **Supp Fig 13**: all five alternative feature sets also give P < 0.05 for ERCC2 (0.017-0.04); the paper
   describes them as "reduced statistical power" compared with P = 0.002.
10. **Supp Fig 15c,d**: cell-line-trained models give P = 0.0092 for their rank-2 (all lines) or rank-3 (bladder lines)
    pathway; the paper reports these models as not strongly predictive because the top-ranked pathway was not significant.
11. **TCGA COAD "5FU-treated" n = 114** in the paper, but the project brief says only ~5 TCGA-COAD patients got 5FU alone
    and most got FOLFOX. The paper's 114 almost certainly includes 5FU-containing combinations.
12. **Isogenic BLCA data (T24S/T24R, Yeon et al.)** were obtained by personal request; only ssGSEA scores are public
    (MOESM3 'Figure 4c').
13. Text says pathway mean proximity count 4.9; the source matrix gives 4.86 [derived] (rounding, consistent).
14. STRING sheet lists each link twice (647,548 rows = 2 x 323,774) and has 13,825 distinct UniProt IDs vs 13,824 in the text [derived].

## 6. What the original code in `python/` already produces

| script | produces | checklist items it covers |
|---|---|---|
| `python/run_ssGSEA.py` (+ `utilities/`) | Reactome ssGSEA NES for organoids and TCGA (`python/results/COAD/...`) | inputs for Fig2b, Fig2d, Fig2f, Fig2h, Fig3c, SuppFig1 |
| `python/single_pathway_prediction.py` | per-ML (Ridge/SVR/Linear) coefficients + KM plot of the rank-k pathway; `single_pathway_predictions.txt` (COAD 5FU only configured) | Fig2b, Fig2d-left, Fig2f/Fig2h-left (if BLCA enabled), SuppFig2, SuppFig3 (coefficients) |
| `python/multiple_pathway_prediction.py` | top-k (2-10) multi-pathway KM plots + `multi_pathway_predictions.txt` (present in `python/results/multi_pathway_predictions/`) | SuppFig5, SourceData-multipathway |
| `modified_codes_py3/proximity_code.py` | network proximity z-scores (user's Python 3 port) | Fig2a, Fig2e, SuppFig16 |

Not in the original code (must be written new): untreated-control KM (Fig2d/h right), all Fig3 baselines,
isogenic t-tests (Fig4b,c, SuppTable3), bootstrap (Fig4e,f, SuppTable4), mutation concordance (Fig5, SuppFig12, 13, 18),
cross-biomarker test (SuppFig4), centrality (SuppFig6, 7), target-gene KM (SuppFig8), k-hop neighbours (SuppFig9),
pathway-size check (SuppFig10), cetuximab heatmap (SuppFig11), GDSC models + PCA (SuppFig14, 15), split-organoid test
(SuppFig17), clinical tables (SuppTable1, 2), network drawings (Fig2c, g), PCA (SuppFig1).

## 7. Counts

Main-figure items: 31. Supplementary-figure items: 37. Supplementary tables: 4. Source-data-only items: 3. Text-only items: 2.
Total: 77 checklist rows (see `checklist.tsv`).
