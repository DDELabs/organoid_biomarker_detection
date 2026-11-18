# Organoid-only strategy (Kong et al. 2020) for bladder, ovary, lung, sarcoma and PDX

Agent: "organoid-only, bladder / ovary / lung / other". Date: 2026-11-17.
Every number below was written by the code in `BLCA_OV_LUAD_OTHER/code/` into the tables named next to it.

## 1. What was done

I ran the paper's method step by step. (1) Find Reactome pathways close to the drug's target proteins in the paper's own
STRING network (proximity z-score <= -1.2816). (2) Score these pathways in the lab models (ssGSEA). (3) Train Ridge, linear
SVR and LinearRegression on the lab models' drug response. (4) Take the pathway with the largest Ridge coefficient
(rank 1), score TCGA patients with it, split them at the median and compare overall survival (Kaplan-Meier, log-rank test),
in patients given the drug and in patients with no drug record. Added: hazard ratio (HR) with 95% confidence interval (CI),
leave-one-out (LOO) prediction inside the lab models, multi-pathway models (k = 1..10), the paper's random-pathway bootstrap
(Fig. 4e/f), a permutation test (1000 shuffles of the lab-model response) and a cut-off check (paper -1.2816; user's
-2.1925 for bladder and -1.1824 for lung; also -1.0 and -1.645).

How to read the numbers: HR > 1 means "predicted non-responders died sooner", the direction the method expects.
The permutation p says how often shuffled (meaningless) responses give a patient p at least as small as the real one. It
includes the step of picking one pathway out of many, so it is the most honest single number per test.

## 2. Data (and substitutions)

| Set | Lab models | Response | Patients | Substitution / caveat |
|---|---|---|---|---|
| Bladder, Lee 2018, current data | 9-11 organoids (`data/curated/preclinical/bladder_lee2018`) | CoderData refit AUC | TCGA-BLCA, current GDC | The paper's IC50 table is not public; AUC used. |
| Bladder, paper source data | paper's organoid and TCGA NES (MOESM3 'Figure 2e', 'Figure 2h') | same AUC | TCGA-BLCA, paper-era biotab survival | Exact paper inputs except the response values. |
| Ovary, Vias 2023 | 14 HGSOC organoids (rebuilt from the raw repo files; sha256 match the provenance file) | AUC | TCGA-OV | No carboplatin/cisplatin screen: oxaliplatin = platinum proxy; docetaxel uses paclitaxel. 3 patients give 2 organoids each. |
| Ovary, GSE229119 (user's) | 3 organoids + 4 cell lines, olaparib time course | none | TCGA-OV | Not a response panel (section 5). |
| Lung, GSE165611 (user's) | 102 patient-derived cells (expression only) | not public | TCGA-LUAD | Substitute: 40-46 GDSC lung adenocarcinoma lines (ln IC50) with CCLE RNA-seq; CTRPv2 AUC for carboplatin (section 6). |
| Sarcoma, Al Shihabi 2024 | 14 organoids | % viability | TCGA-SARC | Doxorubicin not screened. |
| PDX, Gao 2015 | 24-41 PDX per tumour type | BestAvgResponse | TCGA BRCA, PAAD, COAD, LUAD, SKCM | One mouse per model. |

Network: the paper's STRING>700 edges (MOESM3 'Figure 2a-f'). For drugs not in the paper's z table (pemetrexed, etoposide,
vinorelbine, tamoxifen, trastuzumab, dacarbazine, olaparib) proximity was recomputed on that network with
`robust/obd/network.proximity_z`. Check: recomputed cisplatin vs the paper's table, Pearson r = 0.9988, 47 of 47 proximal
pathways shared (`BLCA_OV_LUAD_OTHER/SUMMARY/tables/proximity_recomputation_check_cisplatin.tsv`). Carboplatin and
oxaliplatin have no protein targets in DrugBank, so the paper's cisplatin proximity is used for them.
Treated patients: BCR biotab drug files, parsed with `robust/obd/tcga_response.py` functions (BLCA cisplatin 80 records,
79 with RNA-seq and survival; methotrexate 12; OV carboplatin 482; LUAD cisplatin 87; as in the brief).
File sizes and sha256 of every input: `BLCA_OV_LUAD_OTHER/SUMMARY/tables/data_versions.tsv` (MOESM3: 23,499,402 bytes).
Technical note: the response is z-scored before fitting. This leaves the Ridge and LinearRegression ranking unchanged but
stops linear SVR (default margin 0.1) from returning all-zero coefficients on 0-1 AUC values.

## 3. Bladder cisplatin (paper headline #2)

### 3.1 Exact reproduction on the paper's own source data (`BLCA_lee2018/paper/`)

- 47 pathways pass z <= -1.2816; **30** remain after dropping those missing from the organoid data, exactly the paper's
  number. 9 organoids have a cisplatin response (the paper's 9). These 9 lines come from **7 patients**: SCBO-3_2 and
  SCBO-11_2 are recurrence lines of the patients of SCBO-3 and SCBO-11 (audit I07). The SCBO-3 response in the curated
  table also pools SMBO-3 screens (`data/curated/preclinical/bladder_lee2018/response.tsv`, column `screens`).
  Treated patients: ours 79 (20 deaths) and 288 untreated, the paper 77 and 294 (the paper's patient list is not
  public, so the 2 extra / 6 missing patients cannot be named); organoid response = CoderData AUC, not the paper's IC50
  (audit I17). One-line-per-patient sensitivity: see section 3.1b below.
- Rank 1 in all three models: "transport of glucose and other sugars, bile salts, organic acids, metal ions and amine
  compounds". This is also first in the paper's own multi-pathway source data. The paper's reported biomarker, "amino acid
  synthesis and interconversion (transamination)", is **rank 2 in the paper** (its code tests rank 1). In our AUC refit it is rank 6.

| Test | n | deaths | p | HR (95% CI) | What I think |
|---|---|---|---|---|---|
| Rank-1 pathway, cisplatin-treated | 79 | 20 | 0.044 | 2.63 (0.99-6.96) | Same direction as the paper; 5-year survival 0.75 vs 0.44. Close to the paper's 0.01 given AUC instead of IC50. |
| Rank-1 pathway, no drug record | 288 | 79 | 0.93 | 1.02 (0.66-1.59) | No split in untreated patients, as the paper wants. |
| Paper's pathway (amino acid), treated | 79 | 20 | 0.031 | 2.66 (1.06-6.69) | The paper's biomarker reproduces on its own data (paper p = 0.01). |
| Paper's pathway, untreated | 288 | 79 | 0.087 | 0.67 (0.42-1.06) | Paper reports 0.066; similar. |
| Permutation, 1000 shuffles | | | 0.30 | | 29.6% of shuffles also reach p < 0.05: with 9 organoids and 30 pathways a p of 0.04 is easy to reach by chance. |
| Our bootstrap, 10,000 random sets (survival p<0.05 and isogenic p<0.05 in the predicted direction) | | | 0.044 | | First proximal pathway passing both criteria is at rank 6. |
| Paper's bootstrap table recounted (MOESM3 'Figure 4f') | | | 0.014 | | 66 + 74 of 10,000 random sets reach rank 1 or 2 (6,366 none). Matches the paper's 0.014, but only counting rank 2. |
| LOO inside organoids (Ridge) | 9 | | Spearman r = -0.07 (p = 0.86) | | The 30 pathways do not predict a held-out organoid's AUC. The paper's 0.89 (Supp Fig 17b) used a different split and IC50. |
| Ridge k = 2 (paper Supp Fig 5b: 0.022) | 79 | | 0.060 | | Close. Smallest of 27 multi-pathway tests: 0.017 (Ridge k = 10). |
| Cut-off -2.1925 (user's), 11 pathways | 79 | | 0.91 | | The user's stricter cut-off loses the signal. |

Clinical balance of the predicted groups (Supp Table 2): age, sex, T and N stage, chi-square p 0.11-1.0
(`BLCA_lee2018/extras/tables/SuppTable2_clinical_balance.tsv`). What I think: the split is not a hidden stage split.

### 3.1b One line per patient (audit I07; `BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv`)
The two recurrence lines (SCBO-3_2, SCBO-11_2) were dropped, leaving 7 lines from 7 patients (paper source data, same
code, 500 permutations and 2,000 bootstraps to save CPU, so the empirical p values are only comparable in size).

| Analysis | Lines | Rank-1 pathway | Treated n / deaths | Treated log-rank p | HR (95% CI) | Untreated p | Shuffled-response p | Bootstrap p (rank statistic) |
|---|---|---|---|---|---|---|---|---|
| As main run | 9 (7 patients) | glucose and other sugars transport | 79 / 20 | 0.044 | 2.63 (0.99-6.96) | 0.93 | 0.295 | 0.078 |
| One line per patient | 7 (7 patients) | same pathway | 79 / 20 | 0.044 | 2.63 (0.99-6.96) | 0.93 | 0.263 | 0.084 |

What I think: the rank-1 pathway and its treated-patient statistics do not change when the recurrence lines are
removed (the 30 candidate pathways and the rank-1 choice stay the same). The shuffled-response and bootstrap checks
move a little (0.295 -> 0.263; 0.078 -> 0.084; rank-1 p bootstrap 0.064 -> 0.083). n = 7 patients is the honest
sample size for the organoid side.

### 3.2 Current GDC data (`BLCA_lee2018/current/`)

47 proximal pathways; rank 1 "transport of organic anions". Treated n = 79, deaths = 28, **p = 0.067, HR 2.03 (0.94-4.38)**;
untreated n = 283, p = 0.11; **permutation p = 0.036**; bootstrap (rank-1 p) 0.15; Ridge k = 2 p = 0.071; smallest
multi-pathway p = 0.008. What I think: the patient p is slightly larger than on the paper's data, but this is the only
bladder test where the shuffled-response check is below 0.05. The amino acid pathway is rank 30 of 47 here; used directly it
gives treated p = 0.063, HR 2.03 (0.95-4.34), and untreated p = 0.021, HR 0.66 (`extras/tables/paper_biomarker_amino_acid_KM.tsv`).

### 3.3 Isogenic cells (Fig 4c, Supp Table 3) (`BLCA_lee2018/extras/`)

Only the paper's NES for 3 T24-sensitive and 3 T24-resistant samples exist publicly. Two-sided t-test: amino acid pathway
lower in resistant cells, p = 0.00022 (the paper's number). Our rank-1 pathway (paper data): p = 0.017, but higher in
sensitive cells while its organoid coefficient says higher = resistant, so not in the predicted direction. Top-3 most
proximal pathways: p = 0.014, 0.00054, 0.010. What I think: the isogenic result supports the amino acid pathway, not rank 1.

### 3.4 ERCC2 mutation concordance (Fig 5b, Supp Fig 13) (`extras/tables/ERCC2_mutation_concordance.tsv`)

Open GDC mutations: 36 ERCC2-mutant vs 366-367 wild-type tumours (paper 37 / 367). One-sided Mann-Whitney, mutants
predicted more sensitive.

| Score | p | What I think |
|---|---|---|
| Amino acid pathway, paper NES | 0.0030 | Reproduces the paper's 0.002. |
| Amino acid pathway, current data | 0.00055 | Stronger on current data. |
| Our rank-1 pathway, paper / current data | 0.88 / 0.33 | Rank 1 does not carry this signal. |
| Whole transcriptome | 0.065 | Paper 0.04. |
| Whole pathway, paper / current data | 0.010 / 0.022 | Paper 0.037. |
| Network neighbours, 1 / 2 / 3 steps | 0.081 / 0.035 / 0.052 | Paper 0.017 / 0.021 / 0.028. |

### 3.5 Other feature sets in cisplatin-treated patients (Fig 3b-d, Supp Fig 4b, 8b, 9)

Whole transcriptome p = 0.71; whole pathway 0.59 (paper data) / 0.86 (current); neighbours 1/2/3 steps 0.70 / 0.94 / 0.99;
colorectal BH3 biomarker cross-test 0.98 / 0.45; target genes A2M 0.74, ATOX1 0.37, MPG 0.68, TF 0.65 (paper: 0.79, 0.82,
0.37, 0.79, 0.79, 0.86, 0.92, 0.18, 0.83, 0.48). What I think: as in the paper, only the network-proximal pathway gives a small p.

### 3.6 Cell-line-trained models (Supp Fig 15c,d) and batch PCA (Supp Fig 1)

GDSC2 cisplatin, all lines with CCLE RNA-seq: top-3 pathways p = 0.82, 0.012 (HR 0.36), 0.011 (HR 2.76); bladder lines only:
0.21, 0.69, 0.24 (paper 0.102, 0.0092, 0.84 and 0.13, 0.23, 0.0092). What I think: as in the paper, cell lines give scattered
p-values below rank 1. PCA of pathway scores: PC1/PC2 27.1%/17.8% raw and 20.7%/18.8% scaled (paper 27.0/17.8 and 20.7/18.8).

### 3.7 Other bladder drugs (current data)

| Drug | organoids | treated n / deaths | p | HR (95% CI) | untreated p | perm p | What I think |
|---|---|---|---|---|---|---|---|
| Gemcitabine | 11 | 88 / 34 | 0.51 | 0.80 (0.41-1.56) | 1.00 | 0.45 | Largest bladder cohort, no split. |
| Methotrexate (user's) | 9 | 12 / 5 | 0.61 | 0.63 (0.11-3.81) | 0.041 | 0.80 | 12 patients are too few; the untreated p points to prognosis. User's cut-off -2.1925: p = 0.60. |
| Doxorubicin | 9 | 16 / 8 | 0.37 | 0.50 (0.11-2.32) | 0.90 | 0.41 | Too few patients. |
| Vinblastine | 9 | 12 / 5 | 0.60 | 1.61 (0.27-9.67) | 0.93 | 0.50 | Too few patients. |
| Paclitaxel | 9 | 10 / 6 | 0.36 | 0.46 (0.08-2.53) | 0.20 | 0.62 | Too few patients. |
| Carboplatin (cisplatin proxy) | 9 | 22 / 11 | 0.26 | 2.08 (0.57-7.54) | 0.11 | 0.32 | Same pathway and direction as cisplatin, small group. |

## 4. Ovarian organoids (Vias 2023) -> TCGA-OV (`OV_vias2023/`)

Only 37 TCGA-OV patients with RNA-seq have no drug record, so the untreated comparison is small (untreated p 0.04-0.88).

| Drug (organoid drug) | proximal | treated n / deaths | p | HR (95% CI) | perm p | LOO r (p) | What I think |
|---|---|---|---|---|---|---|---|
| Carboplatin (oxaliplatin) | 47 | 347 / 218 | 0.54 | 1.09 (0.83-1.42) | 0.63 | 0.55 (0.041) | The model predicts held-out organoids, but the polyamine pathway does not split patients. Oxaliplatin AUC is near 1 in most organoids: a weak proxy. |
| Cisplatin (oxaliplatin) | 47 | 123 / 77 | 0.26 | 0.77 (0.49-1.21) | 0.69 | 0.55 (0.041) | Same. |
| Paclitaxel | 150 | 348 / 216 | 0.38 | 1.13 (0.86-1.48) | 0.47 | -0.01 (0.97) | No split; 150 pathways for 14 organoids. |
| Docetaxel (paclitaxel) | 150 | 82 / 48 | 0.47 | 0.81 (0.45-1.44) | 0.49 | -0.01 (0.97) | No split. |
| Doxorubicin | 14 | 117 / 98 | 0.98 | 0.99 (0.66-1.50) | 0.90 | 0.02 (0.93) | No split. |
| Gemcitabine | 45 | 102 / 84 | 0.53 | 1.15 (0.74-1.78) | 0.57 | -0.04 (0.90) | No split. |

Olaparib: no TCGA-OV patient got it, so no survival test is possible.

## 5. User's ovarian attempt: GSE229119 olaparib time course (`OV_GSE229119_olaparib/`)

In plain words: every organoid got olaparib and nobody measured how much each responded. The paper's method needs one
response number per model, so it cannot be trained here. What the data show is how olaparib changes pathways over time.
Each sample was scored, the model's own DMSO control subtracted (in units of each pathway's SD), and pathways that moved the
same way in all 3 organoids at day 7 were kept: 15 up (e.g. unfolded protein response, interferon gamma signalling) and
15 down (e.g. glucose transport, NOTCH pathways, RNA export). Most changes shrink after the drug is removed. The 3 organoids
differ from each other (`tables/day7_change_similarity_spearman.tsv`).
Test in TCGA-OV (score = up pathways minus down pathways):

| Test | n | p | Effect | What I think |
|---|---|---|---|---|
| Overall survival, platinum-treated | 375 (229 deaths) | 0.32 | HR 1.14 (0.88-1.48) | No link with survival. |
| Progression-free interval | 374 (283 events) | 0.13 | HR 1.20 (0.95-1.51) | Weak trend: tumours that already look "olaparib-treated" progress a little sooner. |
| Progression within 6 months vs later | 23 vs 326 | 0.52 | medians 0.11 vs 0.01 | No difference. |
| BRCA1/2 somatic mutation vs wild type | 17 vs 271 | 0.13 | medians -0.01 vs 0.06 | Mutants score slightly lower; unclear. Germline BRCA status is controlled access. |

Can answer: which pathways olaparib moves, how fast, and whether they recover. Cannot answer: who responds to olaparib.

## 6. User's lung attempt: GSE165611 (`LUAD_GSE165611_substitute/`)

The drug table was not found. GEO holds expression only (fields: tissue, disease, histology). The paper (Yu et al. 2023,
J Exp Clin Cancer Res 42:37) has three supplementary files on static-content.springer.com; Table S2 lists the 48 screened
drugs (cisplatin, pemetrexed, gemcitabine, paclitaxel, etoposide, olaparib, ...) but no per-sample values. I assembled the
expression (102 PDCs: 84 adenocarcinoma, 6 SCLC, 12 other; RSEM FPKM) and placed it next to TCGA-LUAD (PCA figure).
The lung analysis uses **GDSC lung adenocarcinoma cell lines** instead (CTRPv2 for carboplatin).

| Drug | lines | treated n / deaths | p | HR (95% CI) | untreated p | perm p | What I think |
|---|---|---|---|---|---|---|---|
| Cisplatin | 40 | 84 / 28 | 0.58 | 1.24 (0.58-2.64) | 0.92 | 0.52 | No split. |
| Carboplatin (CTRP) | 43 | 77 / 33 | 0.58 | 0.82 (0.41-1.65) | 0.54 | 0.51 | No split. |
| Pemetrexed | 45 | 59 / 15 | 0.0048 | 0.19 (0.05-0.69) | 0.54 | 0.12 | Strong split in the **reverse** direction (predicted non-responders live longer); shuffles reach it 12% of the time. User's cut-off -1.1824 gives the same result. |
| Gemcitabine | 46 | 25 / 16 | 0.030 | 3.17 (1.06-9.44) | 0.025 | 0.14 | Right direction, but the untreated split is as strong: the G1/S pathway looks prognostic, not drug-specific. |
| Paclitaxel | 46 | 44 / 21 | 0.71 | 0.84 (0.35-2.05) | 0.080 | 0.70 | No split. |
| Vinorelbine | 46 | 24 / 13 | 0.34 | 1.74 (0.56-5.42) | 0.013 | 0.26 | Signal only in untreated patients (prognostic). |
| Docetaxel | 46 | 22 / 13 | 0.34 | 1.80 (0.53-6.09) | 0.13 | 0.63 | Small group. |
| Etoposide | 46 | 15 / 4 | 0.80 | 1.30 (0.17-9.74) | 0.51 | 0.62 | 4 deaths; cannot judge. |
| Erlotinib | 46 | 21 / 7 | 0.51 | 0.57 (0.10-3.17) | 0.35 | 0.64 | Small group. |

## 7. Sarcoma organoids (`SARC_alshihabi2024/`)

Gemcitabine: 14 organoids, 42 patients, 20 deaths, p = 0.084, HR 0.45 (0.18-1.14, reverse direction), untreated p = 0.34,
permutation p = 0.37, Ridge k = 2 p = 0.019. What I think: a reverse-direction trend in a small group of mixed subtypes.
Docetaxel: 37 patients, 17 deaths, p = 0.75, HR 1.17 (0.45-3.04), permutation p = 0.83. What I think: no split.

## 8. Novartis PDX (`PDX_novartis_gao2015/`)

| Tumour -> TCGA, drug | PDX | treated n / deaths | p | HR (95% CI) | perm p | What I think |
|---|---|---|---|---|---|---|
| CM -> SKCM, dacarbazine | 31 | 32 / 21 | 0.0011 | 0.21 (0.07-0.58) | 0.065 | Smallest p of all sets, but **reverse** direction. Bootstrap p 0.013. Worth a look; not a biomarker as it stands. |
| NSCLC -> LUAD, erlotinib | 25 | 21 / 7 | 0.040 | 7.61 (0.81-71) | 0.045 | Right direction ("EGFR downregulation"), permutation p < 0.05, but only 7 deaths. |
| BRCA, trastuzumab | 36 | 74 / 4 | 0.092 | not estimable (all 4 deaths in one group) | 0.13 | 4 deaths; cannot judge. |
| PDAC -> PAAD, gemcitabine | 33 | 104 / 51 | 0.14 | 1.54 (0.87-2.74) | 0.59 | Right direction, weak. |
| BRCA, paclitaxel | 37 | 231 / 22 | 0.64 | 1.23 (0.52-2.92) | 0.64 | Few deaths in breast cancer. |
| BRCA, tamoxifen | 37 | 253 / 22 | 0.79 | 1.13 (0.48-2.65) | 0.87 | No split. |
| CRC -> COAD, 5-FU | 41 | 129 / 25 | 0.96 | 1.02 (0.46-2.24) | 0.97 | No split. |
| CRC -> COAD, cetuximab | 41 | 14 / 6 | 0.93 | 0.93 (0.18-4.67) | 0.89 | Too few patients. |
| NSCLC -> LUAD, paclitaxel | 24 | 44 / 21 | 0.91 | 1.05 (0.42-2.63) | 0.95 | No split (untreated p = 0.001: prognostic). |

## 9. Across all sets (`BLCA_OV_LUAD_OTHER/SUMMARY/`)

34 set x drug tests. Treated log-rank p < 0.05: 5 (bladder cisplatin on paper data, LUAD gemcitabine, LUAD pemetrexed, PDX
erlotinib, PDX dacarbazine), 3 of them in the expected direction. Permutation p < 0.05: 2 (bladder cisplatin current data
0.036; PDX erlotinib 0.045). About 1.7 of 34 would reach p < 0.05 by chance.
What I think: the clearest support for the paper's idea is still bladder cisplatin, the case the paper was built on. The
paper's amino acid pathway reproduces on its own data in survival (p = 0.031), isogenic cells (p = 0.00022) and ERCC2
mutants (p = 0.003; 0.0005 on current data). That is a consistent story. Against it: the pathway was rank 2 in the paper's
own model, and the permutation test on paper data (p = 0.30) says a 9-organoid model makes small p-values easy to find.
Ovarian and sarcoma organoids give no drug-specific split; the lung and PDX hits are small groups or reverse direction.

## 10. Figures and tables

Per set x drug, in `<set>/figures/` (each with `.caption.md`, PNG 300 dpi + PDF): `<PROJECT>_<DRUG>_full`, 12 panels:
a proximity, b pathway heatmap, c Ridge/SVR/LinearRegression coefficients, d LOO predicted vs observed, e KM treated,
f KM untreated, g p vs k (1..10), h-i bootstrap (Fig 4e/f style), j permutation null, k cut-off sensitivity, l HR of the
three models. Supplement `<PROJECT>_<DRUG>_supp_zcutoff_KM` (KM at the 5 cut-offs). Tables in `<set>/tables/`:
coefficients, rank-1 tests, multi-pathway k, LOO predictions, permutation and bootstrap draws, cut-off sensitivity, summary JSON.
Bladder has `paper/` and `current/` sub-folders.

Main figure per set = headline drug `_full` figure + per-set forest `<set>/MAIN_*_forest.png`: bladder
`BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png`; ovary `OV_vias2023/figures/OV_CARBOPLATIN_full.png`; lung
`LUAD_GSE165611_substitute/figures/LUAD_CISPLATIN_full.png`; sarcoma `SARC_alshihabi2024/figures/SARC_GEMCITABINE_full.png`;
PDX `PDX_novartis_gao2015/figures/SKCM_DACARBAZINE_full.png`.
Bladder extras (`BLCA_lee2018/extras/figures/`): `S_BLCA_paper_biomarker_KM`, `S_BLCA_paper_bootstrap_recount`,
`S_BLCA_isogenic_T24`, `S_BLCA_ERCC2_concordance`, `S_BLCA_comparison_feature_sets`, `S_BLCA_GDSC_trained`,
`S_BLCA_multipathway_vs_paper`, `S_BLCA_PCA_batch`.
Other: `OV_GSE229119_olaparib/figures/OV_GSE229119_olaparib_timecourse.png`,
`LUAD_GSE165611_substitute/figures/GSE165611_PDC_vs_TCGA_PCA.png`, `BLCA_OV_LUAD_OTHER/SUMMARY/figures/cross_set_forest.png`,
`BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv`.

Checklist IDs covered (search_token): Fig2e, Fig2f, Fig2h, Fig3b, Fig3c, Fig3d (bladder), Fig4c, Fig4d, Fig4f, Fig5b,
SuppFig1 (bladder panels), SuppFig2c, SuppFig2d, SuppFig4b, SuppFig5b, SourceData-multipathway, SuppFig8b, SuppFig9a,
SuppFig9b, SuppFig13, SuppFig15c, SuppFig15d, SuppFig17b (as LOO), SuppTable2, SuppTable3, SourceData-STRING, Text-method-params.
Not done for bladder: Fig3e (Bolis selection), Fig3f (deep learning), SuppFig7 (centrality), SuppTable4 as a formatted table
(the random pathways are in `paper/tables/BLCA_CISPLATIN_bootstrap_random_pathways.tsv.gz`).

## 11. Could not be done (one line each)

- Paper's bladder IC50 values: not public; CoderData AUC used.
- T24 isogenic raw expression: not public; the paper's NES used.
- GSE165611 per-PDC drug response: not public anywhere reachable; GDSC/CTRP lung lines used.
- Olaparib survival test: no TCGA-OV patient got olaparib.
- Carboplatin/cisplatin in ovarian organoids and doxorubicin in sarcoma organoids: not screened.

## 12. Reproduce (folder `BLCA_OV_LUAD_OTHER/code/`, seeds 0)

`download_raw.sh` (1.5 GB into `_data/raw`; extra hand-fetched URLs listed at its end) -> `extract_moesm3.py` ->
`run_proximity.py` -> `run_blca.py` -> `run_sets.py OV|SARC|LUAD|PDX` -> `blca_extras.py` -> `ov_olaparib.py` ->
`luad_pdc_expression.py` -> `summary.py`. Shared engine: `engine.py`; data access (reuses `robust/obd`): `common.py`.
