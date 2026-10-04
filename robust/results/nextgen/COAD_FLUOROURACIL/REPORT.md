# COAD / FLUOROURACIL: robust biomarker report

Pre-clinical model: **vdw**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.659, direction: inconsistent, treatment-interaction p = 0.463, events in treated patients = 20 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 19 | 174 | 115 | 20 | 0.0463 (Ridge, k=8) | 0.899 | 1.01 (0.567-1.8) | 0.97 | 0.464 | 0.463 | 1.06 |
| nit_nodeconf | 19 | 174 | 115 | 20 | 0.0248 (Ridge, k=10) | 0.881 | 0.967 (0.543-1.72) | 0.91 | 0.506 | 0.342 | 1.08 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.659 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0463, empirical p = 0.244

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.339 (p = 0.156)
- Original Ridge top-7: Spearman rho = 0.0684 (p = 0.781)

## Fragility

- Leave-one-organoid-out: HR > 1 in 84.2% of refits (range 0.753-1.36)
- Patient bootstrap: adjusted HR 95% interval 0.561-1.78, HR > 1 in 48.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 44 pathways, adj HR 0.993 (p 0.987), signature: REVERSIBLE_HYDRATION_OF_CARBON_DIOXIDE, G1_S_SPECIFIC_TRANSCRIPTION, E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION, AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION, HYALURONAN_UPTAKE_AND_DEGRADATION
  - z <= -1.2816: 33 pathways, adj HR 0.856 (p 0.679), signature: REVERSIBLE_HYDRATION_OF_CARBON_DIOXIDE, G1_S_SPECIFIC_TRANSCRIPTION, E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION, AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION, HYALURONAN_UPTAKE_AND_DEGRADATION, METABOLISM_OF_VITAMINS_AND_COFACTORS
  - z <= -1.645: 22 pathways, adj HR 0.986 (p 0.962), signature: G1_S_SPECIFIC_TRANSCRIPTION, CELL_CYCLE, E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION, AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION, HYALURONAN_UPTAKE_AND_DEGRADATION, METABOLISM_OF_VITAMINS_AND_COFACTORS

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION | robust | 5/5 | -2.56 | 0.686 | 1 | 0.842 | -0.278 | 0.969 |
| AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION | robust | 5/5 | -2.06 | 0.684 | 1 | 0.895 | -0.29 | 1 |
| SCFSKP2_MEDIATED_DEGRADATION_OF_P27_P21 | supported | 4/5 | -0.415 | 0.895 | 1 | 1 | 0.298 | 0.976 |
| HYALURONAN_UPTAKE_AND_DEGRADATION | supported | 4/5 | -2.33 | 0.762 | 1 | 0.895 | 0.318 | 0.973 |
| G1_S_SPECIFIC_TRANSCRIPTION | supported | 4/5 | -2.8 | 0.724 | 1 | 0.895 | -0.286 | 1 |
| PURINE_CATABOLISM | fragile | 2/5 | -2.13 | 0.426 | 0.994 | 0.263 | 0.19 | 1.52 |
| CELL_CYCLE | fragile | 2/5 | -3.41 | 0.146 | 1 | 0 | 0.137 | 1.11 |
| RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY | fragile | 2/5 | -1.53 | 0.0825 | 1 | 0 | 0.105 | 1.69 |
| BASE_EXCISION_REPAIR | fragile | 2/5 | -1.42 | 0.0725 | 1 | 0 | 0.129 | 1.49 |
| ACTIVATION_OF_BH3_ONLY_PROTEINS | fragile | 2/5 | -1.29 | 0.0625 | 1 | 0 | -0.126 | 0.881 |
| ENDOSOMAL_SORTING_COMPLEX_REQUIRED_FOR_TRANSPORT_ESCRT | fragile | 2/5 | -1.67 | 0.0475 | 0.947 | 0 | -0.0441 | 0.786 |
| CDC6_ASSOCIATION_WITH_THE_ORC_ORIGIN_COMPLEX | fragile | 2/5 | -2.01 | 0.0425 | 1 | 0 | -0.0755 | 0.959 |
| CELL_CYCLE_MITOTIC | fragile | 2/5 | -3.26 | 0.0312 | 1 | 0 | 0.0776 | 1.09 |
| MITOTIC_G1_G1_S_PHASES | fragile | 2/5 | -4.26 | 0.0275 | 1 | 0 | 0.00414 | 1.18 |
| SHC1_EVENTS_IN_ERBB4_SIGNALING | fragile | 2/5 | -1.06 | 0.0125 | 1 | 0 | -0.0931 | 0.559 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
