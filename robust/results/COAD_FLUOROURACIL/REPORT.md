# COAD / FLUOROURACIL: robust biomarker report

**Signature status: promising, not validated** (permutation p = 0.0479, direction consistent across data versions: True, treatment-interaction p = 0.124, events in treated patients = 8 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| paper | 19 | 37 | 114 | 8 | 0.0762 (Ridge, k=7) | 0.129 | 2.08 (0.857-5.07) | 0.105 | 0.727 | 0.124 | 2.01 |
| current_ssgsea | 19 | 37 | 117 | 21 | 0.0122 (SVR, k=4) | 0.708 | 1.19 (0.726-1.96) | 0.484 | 0.596 | 0.159 | 1.09 |
| current_rank | 19 | 33 | 117 | 21 | 0.0669 (LinearRegression, k=5) | 0.341 | 1.22 (0.743-2.01) | 0.431 | 0.611 | 0.103 | 1.2 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.0479 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0762, empirical p = 0.479

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.296 (p = 0.218)
- Original Ridge top-7: Spearman rho = 0.0298 (p = 0.904)

## Fragility

- Leave-one-organoid-out: HR > 1 in 100% of refits (range 1.31-2.09)
- Patient bootstrap: adjusted HR 95% interval 0.447-1.79e+06, HR > 1 in 89.1%
- Proximity cut-off sensitivity:
  - z <= -1.0: 50 pathways, adj HR 2.09 (p 0.11), signature: REVERSIBLE_HYDRATION_OF_CARBON_DIOXIDE, HYALURONAN_UPTAKE_AND_DEGRADATION, ACTIVATION_OF_BH3_ONLY_PROTEINS
  - z <= -1.2816: 37 pathways, adj HR 2.08 (p 0.105), signature: REVERSIBLE_HYDRATION_OF_CARBON_DIOXIDE, HYALURONAN_UPTAKE_AND_DEGRADATION, ACTIVATION_OF_BH3_ONLY_PROTEINS
  - z <= -1.645: 23 pathways, adj HR 1.41 (p 0.444), signature: G1_S_SPECIFIC_TRANSCRIPTION, ABACAVIR_TRANSPORT_AND_METABOLISM, HYALURONAN_UPTAKE_AND_DEGRADATION

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| ACTIVATION_OF_BH3_ONLY_PROTEINS | robust | 5/5 | -1.29 | 0.851 | 1 | 0.947 | -0.347 | 0.798 |
| HYALURONAN_UPTAKE_AND_DEGRADATION | robust | 5/5 | -2.33 | 0.729 | 1 | 1 | 0.334 | 1.25 |
| REVERSIBLE_HYDRATION_OF_CARBON_DIOXIDE | robust | 5/5 | -1.63 | 0.701 | 1 | 0.895 | 0.263 | 1.35 |
| AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION | fragile | 2/5 | -2.06 | 0.455 | 1 | 0.158 | -0.149 | 1.22 |
| UNWINDING_OF_DNA | fragile | 2/5 | -1.34 | 0.378 | 1 | 0.0526 | -0.211 | 0.568 |
| PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS | fragile | 2/5 | -2 | 0.35 | 1 | 0.0526 | -0.178 | 0.652 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.61 | 0.278 | 1 | 0 | 0.106 | 1.35 |
| ENDOSOMAL_SORTING_COMPLEX_REQUIRED_FOR_TRANSPORT_ESCRT | fragile | 2/5 | -1.67 | 0.27 | 0.991 | 0 | -0.103 | 0.945 |
| REMOVAL_OF_THE_FLAP_INTERMEDIATE_FROM_THE_C_STRAND | fragile | 2/5 | -1.42 | 0.16 | 1 | 0 | -0.153 | 0.975 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | fragile | 2/5 | -1.65 | 0.124 | 1 | 0 | 0.0932 | 1.94 |
| PURINE_CATABOLISM | fragile | 2/5 | -2.13 | 0.121 | 1 | 0 | 0.131 | 1.98 |
| E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION | fragile | 2/5 | -2.56 | 0.113 | 1 | 0 | -0.136 | 0.561 |
| CDC6_ASSOCIATION_WITH_THE_ORC_ORIGIN_COMPLEX | fragile | 2/5 | -2.01 | 0.101 | 1 | 0 | -0.132 | 0.494 |
| G1_S_SPECIFIC_TRANSCRIPTION | fragile | 2/5 | -2.8 | 0.0938 | 1 | 0 | -0.135 | 0.573 |
| VITAMIN_B5_PANTOTHENATE_METABOLISM | fragile | 2/5 | -2.4 | 0.0725 | 1 | 0 | -0.052 | 0.677 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
