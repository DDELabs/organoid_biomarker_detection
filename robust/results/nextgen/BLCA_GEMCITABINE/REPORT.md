# BLCA / GEMCITABINE: robust biomarker report

Pre-clinical model: **gdsc**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.966, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.383, events in treated patients = 34).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 17 | 181 | 86 | 34 | 0.00408 (Ridge, k=10) | 0.258 | 0.687 (0.46-1.03) | 0.0676 | 0.41 | 0.383 | 0.902 |
| nit_nodeconf | 17 | 181 | 86 | 34 | 0.0234 (Ridge, k=9) | 0.919 | 0.98 (0.637-1.51) | 0.927 | 0.476 | 0.83 | 1.01 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.966 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.00408, empirical p = 0.0679

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.463 (p = 0.0611)
- Original Ridge top-7: Spearman rho = -0.333 (p = 0.191)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.565-0.909)
- Patient bootstrap: adjusted HR 95% interval 0.329-1.01, HR > 1 in 3.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 50 pathways, adj HR 0.718 (p 0.082), signature: BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY, EARLY_PHASE_OF_HIV_LIFE_CYCLE, PYRIMIDINE_METABOLISM
  - z <= -1.2816: 35 pathways, adj HR 0.635 (p 0.0189), signature: PYRIMIDINE_CATABOLISM, BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY, EARLY_PHASE_OF_HIV_LIFE_CYCLE, PYRIMIDINE_METABOLISM
  - z <= -1.645: 22 pathways, adj HR 0.712 (p 0.096), signature: PYRIMIDINE_CATABOLISM, BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY, CYTOSOLIC_TRNA_AMINOACYLATION, SYNTHESIS_AND_INTERCONVERSION_OF_NUCLEOTIDE_DI_AND_TRIPHOSPHATES, EARLY_PHASE_OF_HIV_LIFE_CYCLE, PYRIMIDINE_METABOLISM

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| CYTOSOLIC_TRNA_AMINOACYLATION | robust | 5/5 | -2.36 | 0.581 | 0.983 | 0.706 | 0.195 | 1.01 |
| PYRIMIDINE_METABOLISM | supported | 4/5 | -4.53 | 0.689 | 1 | 0.882 | -0.29 | 1.07 |
| BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY | supported | 3/5 | -1.73 | 0.575 | 1 | 0.765 | 0.175 | 0.721 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | supported | 3/5 | -1.86 | 0.547 | 1 | 0.471 | 0.175 | 0.889 |
| SYNTHESIS_AND_INTERCONVERSION_OF_NUCLEOTIDE_DI_AND_TRIPHOSPHATES | fragile | 2/5 | -5.1 | 0.474 | 1 | 0.294 | 0.157 | 1.19 |
| PYRIMIDINE_CATABOLISM | fragile | 2/5 | -2.62 | 0.379 | 1 | 0.118 | 0.154 | 0.932 |
| GRB2_EVENTS_IN_ERBB2_SIGNALING | fragile | 2/5 | -1.21 | 0.151 | 1 | 0.0588 | 0.115 | 1.01 |
| SYNTHESIS_OF_SUBSTRATES_IN_N_GLYCAN_BIOSYTHESIS | fragile | 2/5 | -1.05 | 0.085 | 1 | 0 | 0.0999 | 1.57 |
| PURINE_METABOLISM | fragile | 2/5 | -2.64 | 0.0688 | 0.2 | 0 | -0.00333 | 0.919 |
| SHC1_EVENTS_IN_EGFR_SIGNALING | fragile | 2/5 | -1.34 | 0.0325 | 1 | 0 | 0.0569 | 1.27 |
| METABOLISM_OF_CARBOHYDRATES | fragile | 2/5 | -1.37 | 0.0325 | 1 | 0 | 0.0428 | 1.12 |
| E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION | fragile | 2/5 | -1.29 | 0.0213 | 1 | 0 | 0.0535 | 1 |
| SULFUR_AMINO_ACID_METABOLISM | fragile | 1/5 | -1.58 | 0.487 | 1 | 0.294 | 0.168 | 0.579 |
| RESOLUTION_OF_AP_SITES_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY | fragile | 1/5 | -1.52 | 0.438 | 1 | 0.0588 | 0.159 | 0.78 |
| PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS | fragile | 1/5 | -2.51 | 0.38 | 0.493 | 0.0588 | 0.0826 | 1.15 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
