# BRCA / DOXORUBICIN: robust biomarker report

Pre-clinical model: **cell-line proxy (GDSC 2012)**. Drug targets: NOLC1, TOP2A. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.507, direction: inconsistent, treatment-interaction p = 0.254, events in treated patients = 29 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 37 | 8 | 334 | 29 | 0.0923 (LinearRegression, k=5) | 0.853 | 1.01 (0.686-1.49) | 0.959 | 0.582 | 0.254 | 1.01 |
| current_rank | 37 | 8 | 334 | 29 | 0.166 (LinearRegression, k=2) | 0.993 | 0.841 (0.58-1.22) | 0.361 | 0.524 | 0.824 | 0.841 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.507 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0923, empirical p = 0.489

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.11 (p = 0.517)
- Original Ridge top-7: Spearman rho = -0.167 (p = 0.323)

## Fragility

- Leave-one-organoid-out: HR > 1 in 86.5% of refits (range 0.875-1.18)
- Patient bootstrap: adjusted HR 95% interval 0.665-1.58, HR > 1 in 58%
- Proximity cut-off sensitivity:
  - z <= -1.0: 16 pathways, adj HR 1.03 (p 0.894), signature: METABOLISM_OF_NON_CODING_RNA, MITOTIC_G2_G2_M_PHASES, CHROMOSOME_MAINTENANCE, SIGNAL_TRANSDUCTION_BY_L1, SYNTHESIS_SECRETION_AND_INACTIVATION_OF_GLP1
  - z <= -1.2816: 8 pathways, adj HR 1.01 (p 0.959), signature: SIGNAL_TRANSDUCTION_BY_L1
  - z <= -1.645: 0 pathways, adj HR nan (p nan), signature: 

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| SIGNAL_TRANSDUCTION_BY_L1 | robust | 5/5 | -1.59 | 0.98 | 0.985 | 1 | 0.272 | 1.01 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | supported | 3/5 | -1.36 | 0.921 | 0.776 | 0.0541 | 0.16 | 1.19 |
| EXTENSION_OF_TELOMERES | supported | 3/5 | -1.51 | 0.863 | 0.00145 | 0 | 0.0195 | 1.3 |
| LOSS_OF_NLP_FROM_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.5 | 0.94 | 0.713 | 0.0811 | -0.133 | 1.16 |
| METABOLISM_OF_NON_CODING_RNA | fragile | 2/5 | -1.52 | 0.929 | 0.779 | 0.108 | -0.164 | 1.17 |
| RECRUITMENT_OF_MITOTIC_CENTROSOME_PROTEINS_AND_COMPLEXES | fragile | 2/5 | -1.41 | 0.767 | 0.443 | 0 | -0.0542 | 1.15 |
| TELOMERE_MAINTENANCE | fragile | 2/5 | -1.54 | 0.73 | 0.13 | 0 | 0.0234 | 1.24 |
| INTEGRATION_OF_PROVIRUS | fragile | 1/5 | -1.51 | 0.87 | 0.251 | 0 | -0.0167 | 1.56 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
