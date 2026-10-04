# PAAD / GEMCITABINE: robust biomarker report

Pre-clinical model: **gdsc**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.778, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.475, events in treated patients = 47).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 15 | 181 | 100 | 47 | 0.0784 (Ridge, k=10) | 0.784 | 0.908 (0.678-1.22) | 0.516 | 0.516 | 0.475 | 0.886 |
| nit_nodeconf | 15 | 181 | 100 | 47 | 0.00898 (Ridge, k=10) | 0.35 | 0.83 (0.595-1.16) | 0.273 | 0.437 | 0.34 | 0.895 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.778 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0784, empirical p = 0.389

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.132 (p = 0.639)
- Original Ridge top-7: Spearman rho = -0.329 (p = 0.232)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.861-0.951)
- Patient bootstrap: adjusted HR 95% interval 0.699-1.42, HR > 1 in 40.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 50 pathways, adj HR 0.908 (p 0.517), signature: ABACAVIR_TRANSPORT_AND_METABOLISM, PURINE_METABOLISM
  - z <= -1.2816: 35 pathways, adj HR 0.913 (p 0.527), signature: ABACAVIR_TRANSPORT_AND_METABOLISM, PURINE_METABOLISM, GLUCOSE_METABOLISM
  - z <= -1.645: 22 pathways, adj HR 0.912 (p 0.526), signature: ABACAVIR_TRANSPORT_AND_METABOLISM, PURINE_METABOLISM, GLUCOSE_METABOLISM

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| ABACAVIR_TRANSPORT_AND_METABOLISM | supported | 4/5 | -3.3 | 0.978 | 1 | 1 | 0.4 | 0.872 |
| PURINE_METABOLISM | supported | 3/5 | -2.64 | 0.512 | 0.976 | 0.533 | 0.138 | 1.39 |
| VITAMIN_B5_PANTOTHENATE_METABOLISM | supported | 3/5 | -2.27 | 0.325 | 0.962 | 0.133 | 0.0706 | 1.08 |
| GLUCOSE_METABOLISM | fragile | 2/5 | -1.95 | 0.459 | 1 | 0.2 | 0.164 | 1.1 |
| CYTOSOLIC_TRNA_AMINOACYLATION | fragile | 2/5 | -2.36 | 0.449 | 0.989 | 0.2 | -0.15 | 1.09 |
| GLYCOLYSIS | fragile | 2/5 | -2.62 | 0.441 | 0.966 | 0.133 | 0.15 | 1.11 |
| PYRIMIDINE_CATABOLISM | fragile | 2/5 | -2.62 | 0.424 | 0.988 | 0.133 | 0.113 | 1.06 |
| PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS | fragile | 2/5 | -2.51 | 0.419 | 0.94 | 0.133 | 0.111 | 1.16 |
| FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS | fragile | 2/5 | -1.74 | 0.29 | 0.991 | 0 | 0.0843 | 0.868 |
| G1_S_SPECIFIC_TRANSCRIPTION | fragile | 2/5 | -2.64 | 0.229 | 0.913 | 0.0667 | 0.0566 | 1.31 |
| AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION | fragile | 2/5 | -1.28 | 0.198 | 1 | 0 | -0.115 | 0.972 |
| E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION | fragile | 2/5 | -2.13 | 0.164 | 1 | 0 | 0.0821 | 1.28 |
| GLUCONEOGENESIS | fragile | 2/5 | -1.42 | 0.0625 | 1 | 0 | 0.0822 | 1.13 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.58 | 0.06 | 1 | 0 | 0.0335 | 1.21 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | fragile | 2/5 | -1.86 | 0.035 | 1 | 0 | -0.0405 | 0.834 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
