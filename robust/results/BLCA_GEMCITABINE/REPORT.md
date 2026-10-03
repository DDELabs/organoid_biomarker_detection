# BLCA / GEMCITABINE: robust biomarker report

Pre-clinical model: **cell-line proxy (GDSC 2012)**. Drug targets: CMPK1, RRM1, TYMS. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.609, direction: reversed (resistant score -> better survival), treatment-interaction p = 1, events in treated patients = 34).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 17 | 35 | 88 | 34 | 0.0289 (LinearRegression, k=5) | 0.873 | 0.959 (0.674-1.37) | 0.817 | 0.532 | 1 | 1.29 |
| current_rank | 17 | 35 | 88 | 34 | 0.0083 (LinearRegression, k=4) | 0.197 | 0.752 (0.52-1.09) | 0.129 | 0.402 | 0.184 | 1.18 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.609 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0289, empirical p = 0.22

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.75 (p = 0.000526)
- Original Ridge top-7: Spearman rho = -0.446 (p = 0.0727)

## Fragility

- Leave-one-organoid-out: HR > 1 in 11.8% of refits (range 0.714-1.11)
- Patient bootstrap: adjusted HR 95% interval 0.717-1.33, HR > 1 in 42.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 50 pathways, adj HR 1.06 (p 0.737), signature: GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS, SHC1_EVENTS_IN_EGFR_SIGNALING
  - z <= -1.2816: 35 pathways, adj HR 0.959 (p 0.817), signature: SHC1_EVENTS_IN_EGFR_SIGNALING, RESOLUTION_OF_AP_SITES_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY
  - z <= -1.645: 22 pathways, adj HR 0.783 (p 0.155), signature: BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY, CYTOSOLIC_TRNA_AMINOACYLATION, PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS, PURINE_METABOLISM, EARLY_PHASE_OF_HIV_LIFE_CYCLE

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| RESOLUTION_OF_AP_SITES_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY | supported | 3/5 | -1.52 | 0.562 | 0.987 | 0.647 | 0.239 | 0.788 |
| SHC1_EVENTS_IN_EGFR_SIGNALING | supported | 3/5 | -1.34 | 0.501 | 0.925 | 0.588 | 0.216 | 1.26 |
| PURINE_METABOLISM | supported | 3/5 | -2.64 | 0.435 | 1 | 0.118 | -0.228 | 0.912 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | fragile | 2/5 | -1.86 | 0.459 | 1 | 0.471 | 0.214 | 0.903 |
| PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS | fragile | 2/5 | -2.51 | 0.412 | 0.976 | 0.353 | -0.217 | 1.15 |
| PYRIMIDINE_METABOLISM | fragile | 2/5 | -4.53 | 0.396 | 1 | 0.118 | -0.185 | 0.984 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.58 | 0.285 | 1 | 0 | -0.158 | 1.05 |
| METABOLISM_OF_NUCLEOTIDES | fragile | 2/5 | -5.21 | 0.139 | 1 | 0 | -0.121 | 0.956 |
| GLUCONEOGENESIS | fragile | 2/5 | -1.42 | 0.135 | 1 | 0 | -0.0924 | 0.729 |
| GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS | fragile | 1/5 | -1.33 | 0.424 | 1 | 0.412 | 0.237 | 0.949 |
| CYTOSOLIC_TRNA_AMINOACYLATION | fragile | 1/5 | -2.36 | 0.383 | 0.993 | 0.0588 | 0.195 | 0.944 |
| BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY | fragile | 1/5 | -1.73 | 0.364 | 0.911 | 0.0588 | 0.171 | 0.709 |
| P2Y_RECEPTORS | fragile | 1/5 | -1.54 | 0.282 | 0.708 | 0 | -0.1 | 0.967 |
| VITAMIN_B5_PANTOTHENATE_METABOLISM | fragile | 1/5 | -2.27 | 0.21 | 1 | 0.0588 | -0.133 | 1.07 |
| PURINE_SALVAGE | fragile | 1/5 | -3.82 | 0.175 | 0.671 | 0 | -0.0669 | 0.851 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
