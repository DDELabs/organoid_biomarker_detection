# BRCA / DOXORUBICIN: robust biomarker report

Pre-clinical model: **gdsc**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: promising, not validated** (permutation p = 0.285, direction: as expected (resistant score -> worse survival), treatment-interaction p = 0.388, events in treated patients = 22 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 37 | 126 | 311 | 22 | 0.248 (LinearRegression, k=9) | 0.288 | 1.19 (0.79-1.8) | 0.402 | 0.622 | 0.388 | 1.18 |
| nit_nodeconf | 37 | 126 | 311 | 22 | 0.0891 (Ridge, k=2) | 0.00446 | 1.35 (0.867-2.09) | 0.186 | 0.69 | 0.202 | 1.66 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.285 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.248, empirical p = 0.715

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.404 (p = 0.013)
- Original Ridge top-7: Spearman rho = -0.131 (p = 0.439)

## Fragility

- Leave-one-organoid-out: HR > 1 in 100% of refits (range 1.07-1.54)
- Patient bootstrap: adjusted HR 95% interval 0.819-1.88, HR > 1 in 83.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 16 pathways, adj HR 1.24 (p 0.41), signature: SIGNAL_TRANSDUCTION_BY_L1, INTEGRATION_OF_PROVIRUS, TELOMERE_MAINTENANCE
  - z <= -1.2816: 8 pathways, adj HR 1.18 (p 0.434), signature: SIGNAL_TRANSDUCTION_BY_L1
  - z <= -1.645: 0 pathways, adj HR nan (p nan), signature: 

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| SIGNAL_TRANSDUCTION_BY_L1 | robust | 5/5 | -1.59 | 0.821 | 0.982 | 0.973 | 0.184 | 1.2 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | supported | 4/5 | -1.36 | 0.541 | 0.986 | 0.811 | -0.121 | 0.933 |
| TELOMERE_MAINTENANCE | supported | 3/5 | -1.54 | 0.525 | 0.971 | 0.351 | 0.123 | 1.01 |
| INTEGRATION_OF_PROVIRUS | supported | 3/5 | -1.51 | 0.487 | 0.964 | 0.595 | 0.104 | 1.39 |
| INTERACTIONS_OF_VPR_WITH_HOST_CELLULAR_PROTEINS | fragile | 2/5 | -0.784 | 0.347 | 1 | 0.027 | -0.0848 | 0.865 |
| METABOLISM_OF_NON_CODING_RNA | fragile | 2/5 | -1.52 | 0.305 | 0.578 | 0.027 | -0.00954 | 0.922 |
| CHROMOSOME_MAINTENANCE | fragile | 2/5 | -1.01 | 0.158 | 0.984 | 0 | 0.0574 | 0.991 |
| SYNTHESIS_SECRETION_AND_INACTIVATION_OF_GLP1 | fragile | 2/5 | -1.04 | 0.146 | 0.744 | 0 | -0.0231 | 0.88 |
| ASSOCIATION_OF_TRIC_CCT_WITH_TARGET_PROTEINS_DURING_BIOSYNTHESIS | fragile | 2/5 | -0.695 | 0.135 | 0.981 | 0 | 0.0413 | 1.77 |
| REMOVAL_OF_THE_FLAP_INTERMEDIATE_FROM_THE_C_STRAND | fragile | 2/5 | -1.01 | 0.122 | 1 | 0 | 0.0458 | 1.1 |
| MYOGENESIS | fragile | 2/5 | -0.544 | 0.0963 | 1 | 0 | 0.0579 | 1.08 |
| RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY | fragile | 2/5 | -0.752 | 0.0838 | 1 | 0 | -0.0488 | 0.961 |
| CELL_CELL_JUNCTION_ORGANIZATION | fragile | 2/5 | -0.58 | 0.0225 | 1 | 0 | 0.0377 | 1.03 |
| REPAIR_SYNTHESIS_FOR_GAP_FILLING_BY_DNA_POL_IN_TC_NER | fragile | 2/5 | -0.821 | 0.0187 | 1 | 0 | 0.0306 | 1.03 |
| TRANSPORT_OF_RIBONUCLEOPROTEINS_INTO_THE_HOST_NUCLEUS | fragile | 2/5 | -0.813 | 0.0163 | 1 | 0 | -0.0276 | 0.843 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
