# PAAD / GEMCITABINE: robust biomarker report

Pre-clinical model: **tiriac**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.844, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.0164, events in treated patients = 47).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 44 | 183 | 100 | 47 | 0.0296 (SVR, k=9) | 0.152 | 0.824 (0.627-1.08) | 0.165 | 0.406 | 0.0164 | 1.21 |
| nit_nodeconf | 44 | 183 | 100 | 47 | 0.0468 (SVR, k=6) | 0.152 | 0.824 (0.627-1.08) | 0.165 | 0.406 | 0.0164 | 1.21 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.844 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0296, empirical p = 0.299

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.275 (p = 0.0707)
- Original Ridge top-7: Spearman rho = -0.292 (p = 0.0544)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.794-0.827)
- Patient bootstrap: adjusted HR 95% interval 0.635-1.08, HR > 1 in 11.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 51 pathways, adj HR 0.824 (p 0.165), signature: EARLY_PHASE_OF_HIV_LIFE_CYCLE
  - z <= -1.2816: 36 pathways, adj HR 0.824 (p 0.165), signature: EARLY_PHASE_OF_HIV_LIFE_CYCLE
  - z <= -1.645: 23 pathways, adj HR 0.767 (p 0.157), signature: BASE_EXCISION_REPAIR, ABACAVIR_TRANSPORT_AND_METABOLISM, METABOLISM_OF_NUCLEOTIDES, EARLY_PHASE_OF_HIV_LIFE_CYCLE

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | supported | 4/5 | -1.86 | 0.626 | 1 | 1 | 0.285 | 0.827 |
| ACTIVATION_OF_CHAPERONE_GENES_BY_ATF6_ALPHA | fragile | 2/5 | -0.301 | 0.289 | 1 | 0 | 0.195 | 1.16 |
| GLUCOSE_METABOLISM | fragile | 2/5 | -1.95 | 0.256 | 0.99 | 0 | 0.164 | 1.08 |
| APOBEC3G_MEDIATED_RESISTANCE_TO_HIV1_INFECTION | fragile | 2/5 | -1.78 | 0.249 | 1 | 0 | 0.055 | 1.03 |
| CDC6_ASSOCIATION_WITH_THE_ORC_ORIGIN_COMPLEX | fragile | 2/5 | -1.42 | 0.196 | 1 | 0 | 0.0743 | 1.23 |
| RAF_MAP_KINASE_CASCADE | fragile | 2/5 | -0.645 | 0.151 | 1 | 0 | 0.204 | 1.61 |
| FACILITATIVE_NA_INDEPENDENT_GLUCOSE_TRANSPORTERS | fragile | 2/5 | -1.74 | 0.13 | 1 | 0 | -0.037 | 0.829 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.58 | 0.122 | 1 | 0 | -0.182 | 1.19 |
| SHC1_EVENTS_IN_EGFR_SIGNALING | fragile | 2/5 | -1.34 | 0.0737 | 1 | 0 | 0.0991 | 1.62 |
| FORMATION_OF_TRANSCRIPTION_COUPLED_NER_TC_NER_REPAIR_COMPLEX | fragile | 2/5 | -0.836 | 0.0312 | 1 | 0 | 0.14 | 1.09 |
| P2Y_RECEPTORS | fragile | 2/5 | -1.54 | 0.0187 | 1 | 0 | -0.0576 | 0.832 |
| SIGNALLING_TO_RAS | fragile | 2/5 | -1.01 | 0.00375 | 1 | 0 | 0.142 | 1.08 |
| SYNTHESIS_OF_SUBSTRATES_IN_N_GLYCAN_BIOSYTHESIS | fragile | 2/5 | -1.05 | 0.0025 | 1 | 0 | 0.0129 | 1.04 |
| BASE_EXCISION_REPAIR | fragile | 1/5 | -1.99 | 0.334 | 1 | 0.0227 | -0.208 | 1.3 |
| APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B | fragile | 1/5 | -1.29 | 0.314 | 1 | 0 | -0.226 | 1.3 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
