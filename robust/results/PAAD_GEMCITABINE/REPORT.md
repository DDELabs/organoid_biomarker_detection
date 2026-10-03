# PAAD / GEMCITABINE: robust biomarker report

Pre-clinical model: **cell-line proxy (GDSC 2012)**. Drug targets: CMPK1, RRM1, TYMS. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.573, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.122, events in treated patients = 50).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 15 | 35 | 103 | 50 | 0.198 (SVR, k=2) | 0.498 | 0.984 (0.733-1.32) | 0.917 | 0.494 | 0.122 | 1.07 |
| current_rank | 15 | 35 | 103 | 50 | 0.29 (SVR, k=7) | 0.471 | 0.861 (0.629-1.18) | 0.35 | 0.465 | 0.128 | 1.01 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.573 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.198, empirical p = 0.719

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.161 (p = 0.567)
- Original Ridge top-7: Spearman rho = 0.321 (p = 0.243)

## Fragility

- Leave-one-organoid-out: HR > 1 in 26.7% of refits (range 0.842-1.08)
- Patient bootstrap: adjusted HR 95% interval 0.722-1.3, HR > 1 in 39.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 50 pathways, adj HR 1.03 (p 0.83), signature: GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS, ABACAVIR_TRANSPORT_AND_METABOLISM, FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS, TRANSPORT_OF_VITAMINS_NUCLEOSIDES_AND_RELATED_MOLECULES, SYNTHESIS_OF_SUBSTRATES_IN_N_GLYCAN_BIOSYTHESIS
  - z <= -1.2816: 35 pathways, adj HR 0.984 (p 0.917), signature: GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS, ABACAVIR_TRANSPORT_AND_METABOLISM, CYTOSOLIC_TRNA_AMINOACYLATION, FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS, TRANSPORT_OF_VITAMINS_NUCLEOSIDES_AND_RELATED_MOLECULES
  - z <= -1.645: 22 pathways, adj HR 0.962 (p 0.798), signature: ABACAVIR_TRANSPORT_AND_METABOLISM, CYTOSOLIC_TRNA_AMINOACYLATION, FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS, TRANSPORT_OF_VITAMINS_NUCLEOSIDES_AND_RELATED_MOLECULES

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| ABACAVIR_TRANSPORT_AND_METABOLISM | supported | 4/5 | -3.3 | 0.814 | 1 | 0.933 | 0.323 | 0.91 |
| TRANSPORT_OF_VITAMINS_NUCLEOSIDES_AND_RELATED_MOLECULES | supported | 4/5 | -1.82 | 0.647 | 1 | 0.867 | 0.26 | 1.22 |
| FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS | supported | 4/5 | -1.74 | 0.641 | 1 | 0.733 | 0.228 | 0.973 |
| GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS | supported | 4/5 | -1.33 | 0.631 | 1 | 0.6 | -0.277 | 0.992 |
| CYTOSOLIC_TRNA_AMINOACYLATION | supported | 3/5 | -2.36 | 0.63 | 1 | 0.667 | -0.255 | 1.11 |
| VITAMIN_B5_PANTOTHENATE_METABOLISM | supported | 3/5 | -2.27 | 0.383 | 0.922 | 0.2 | 0.153 | 1.1 |
| E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION | supported | 3/5 | -1.29 | 0.316 | 1 | 0.0667 | 0.119 | 1.31 |
| PURINE_CATABOLISM | fragile | 2/5 | -1.58 | 0.191 | 1 | 0 | 0.116 | 1.44 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | fragile | 2/5 | -1.86 | 0.163 | 0.985 | 0 | -0.0724 | 0.828 |
| CDC6_ASSOCIATION_WITH_THE_ORC_ORIGIN_COMPLEX | fragile | 2/5 | -1.42 | 0.15 | 1 | 0 | 0.0687 | 1.22 |
| RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY | fragile | 2/5 | -2.08 | 0.14 | 1 | 0 | 0.0986 | 1.23 |
| BASE_EXCISION_REPAIR | fragile | 2/5 | -1.99 | 0.131 | 1 | 0 | 0.12 | 1.25 |
| E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION | fragile | 2/5 | -2.13 | 0.122 | 1 | 0 | 0.128 | 1.29 |
| ASSOCIATION_OF_LICENSING_FACTORS_WITH_THE_PRE_REPLICATIVE_COMPLEX | fragile | 2/5 | -2.1 | 0.05 | 1 | 0 | 0.0234 | 1.24 |
| G1_S_SPECIFIC_TRANSCRIPTION | fragile | 2/5 | -2.64 | 0.05 | 1 | 0 | 0.079 | 1.31 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
