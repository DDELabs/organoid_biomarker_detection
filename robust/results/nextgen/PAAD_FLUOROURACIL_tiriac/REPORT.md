# PAAD / FLUOROURACIL: robust biomarker report

Pre-clinical model: **tiriac**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: promising, not validated** (permutation p = 0.18, direction: as expected (resistant score -> worse survival), treatment-interaction p = 0.33, events in treated patients = 18 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 44 | 187 | 36 | 18 | 0.0936 (SVR, k=2) | 0.414 | 1.1 (0.624-1.93) | 0.749 | 0.501 | 0.33 | 1.02 |
| nit_nodeconf | 44 | 187 | 36 | 18 | 0.0232 (SVR, k=4) | 0.378 | 1.07 (0.649-1.76) | 0.791 | 0.545 | 0.395 | 0.967 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.18 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0936, empirical p = 0.517

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.109 (p = 0.482)
- Original Ridge top-7: Spearman rho = 0.312 (p = 0.0391)

## Fragility

- Leave-one-organoid-out: HR > 1 in 100% of refits (range 1.04-1.27)
- Patient bootstrap: adjusted HR 95% interval 0.499-2.31, HR > 1 in 49%
- Proximity cut-off sensitivity:
  - z <= -1.0: 44 pathways, adj HR 1.14 (p 0.712), signature: PYRIMIDINE_CATABOLISM, G1_S_TRANSITION, PURINE_SALVAGE, MITOTIC_G1_G1_S_PHASES, EARLY_PHASE_OF_HIV_LIFE_CYCLE
  - z <= -1.2816: 30 pathways, adj HR 1.26 (p 0.531), signature: G1_S_TRANSITION, PURINE_SALVAGE, MITOTIC_G1_G1_S_PHASES
  - z <= -1.645: 21 pathways, adj HR 1.27 (p 0.535), signature: G1_S_TRANSITION, PURINE_SALVAGE, MITOTIC_G1_G1_S_PHASES

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| PURINE_SALVAGE | robust | 5/5 | -2.02 | 0.981 | 1 | 1 | 0.431 | 1.17 |
| DNA_REPAIR | robust | 5/5 | 0.0527 | 0.52 | 1 | 0.864 | -0.207 | 0.938 |
| G1_S_TRANSITION | supported | 4/5 | -3.58 | 0.708 | 1 | 0.977 | -0.237 | 1.01 |
| FORMATION_OF_INCISION_COMPLEX_IN_GG_NER | supported | 4/5 | -0.982 | 0.583 | 1 | 0.818 | -0.215 | 1.03 |
| EARLY_PHASE_OF_HIV_LIFE_CYCLE | supported | 4/5 | -1.2 | 0.564 | 1 | 0.886 | 0.231 | 0.8 |
| NGF_SIGNALLING_VIA_TRKA_FROM_THE_PLASMA_MEMBRANE | supported | 4/5 | 0.0517 | 0.551 | 1 | 0.932 | 0.206 | 0.952 |
| MITOTIC_G1_G1_S_PHASES | fragile | 2/5 | -3.81 | 0.369 | 1 | 0.0227 | -0.167 | 0.975 |
| PYRIMIDINE_CATABOLISM | fragile | 2/5 | -1.93 | 0.3 | 1 | 0 | 0.148 | 1.24 |
| ABACAVIR_TRANSPORT_AND_METABOLISM | fragile | 2/5 | -2.35 | 0.285 | 1 | 0.0227 | -0.14 | 0.664 |
| ACTIVATION_OF_BH3_ONLY_PROTEINS | fragile | 2/5 | -1.13 | 0.224 | 1 | 0 | 0.142 | 1.15 |
| CELL_CYCLE_MITOTIC | fragile | 2/5 | -3.39 | 0.106 | 1 | 0 | -0.0292 | 0.955 |
| AUTODEGRADATION_OF_THE_E3_UBIQUITIN_LIGASE_COP1 | fragile | 2/5 | -0.545 | 0.0925 | 1 | 0 | -0.183 | 0.988 |
| SYNTHESIS_AND_INTERCONVERSION_OF_NUCLEOTIDE_DI_AND_TRIPHOSPHATES | fragile | 2/5 | -1.49 | 0.0675 | 1 | 0 | -0.0345 | 0.923 |
| RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | fragile | 2/5 | -1.84 | 0.0288 | 1 | 0 | 0.0291 | 1.2 |
| POL_SWITCHING | fragile | 2/5 | -1.47 | 0.0275 | 1 | 0 | -0.0613 | 0.939 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
