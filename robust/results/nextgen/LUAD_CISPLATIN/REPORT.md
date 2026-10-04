# LUAD / CISPLATIN: robust biomarker report

Pre-clinical model: **gdsc**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.968, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.921, events in treated patients = 25 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 65 | 279 | 80 | 25 | 0.0573 (Ridge, k=7) | 0.709 | 0.714 (0.445-1.15) | 0.163 | 0.492 | 0.921 | 0.767 |
| nit_nodeconf | 65 | 279 | 80 | 25 | 0.348 (SVR, k=10) | 0.882 | 0.902 (0.598-1.36) | 0.624 | 0.545 | 0.491 | 1.02 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.968 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0573, empirical p = 0.317

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.171 (p = 0.173)
- Original Ridge top-7: Spearman rho = 0.241 (p = 0.0535)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.674-0.877)
- Patient bootstrap: adjusted HR 95% interval 0.443-1.25, HR > 1 in 15%
- Proximity cut-off sensitivity:
  - z <= -1.0: 88 pathways, adj HR 0.716 (p 0.165), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS, EGFR_DOWNREGULATION
  - z <= -1.2816: 69 pathways, adj HR 0.839 (p 0.391), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS
  - z <= -1.645: 43 pathways, adj HR 0.785 (p 0.293), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS, HOST_INTERACTIONS_OF_HIV_FACTORS

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| EGFR_DOWNREGULATION | supported | 4/5 | -1.15 | 0.748 | 1 | 1 | -0.222 | 1.5 |
| RECYCLING_OF_BILE_ACIDS_AND_SALTS | supported | 4/5 | -2.73 | 0.642 | 1 | 0.923 | -0.221 | 1.11 |
| ABCA_TRANSPORTERS_IN_LIPID_HOMEOSTASIS | supported | 3/5 | -1.97 | 0.378 | 1 | 0.0923 | 0.144 | 1.25 |
| AMYLOIDS | supported | 3/5 | -2.34 | 0.306 | 1 | 0 | 0.127 | 1.03 |
| NFKB_IS_ACTIVATED_AND_SIGNALS_SURVIVAL | fragile | 2/5 | -2.51 | 0.3 | 1 | 0 | -0.115 | 0.802 |
| IL_6_SIGNALING | fragile | 2/5 | -1.65 | 0.294 | 1 | 0.0308 | 0.139 | 1.15 |
| HOST_INTERACTIONS_OF_HIV_FACTORS | fragile | 2/5 | -1.7 | 0.268 | 1 | 0.0154 | -0.129 | 1.45 |
| HEMOSTASIS | fragile | 2/5 | -2.13 | 0.254 | 1 | 0 | 0.116 | 1.08 |
| HDL_MEDIATED_LIPID_TRANSPORT | fragile | 2/5 | -2.73 | 0.195 | 1 | 0 | -0.106 | 0.957 |
| STEROID_HORMONES | fragile | 2/5 | -1.73 | 0.151 | 1 | 0 | 0.0735 | 1.01 |
| PEPTIDE_HORMONE_BIOSYNTHESIS | fragile | 2/5 | -1.34 | 0.105 | 1 | 0 | 0.0588 | 1.08 |
| METABOLISM_OF_STEROID_HORMONES_AND_VITAMINS_A_AND_D | fragile | 2/5 | -1.49 | 0.0813 | 1 | 0 | 0.0475 | 1.04 |
| ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC | fragile | 2/5 | -1.63 | 0.0688 | 1 | 0 | -0.0766 | 0.929 |
| NEF_MEDIATED_DOWNREGULATION_OF_MHC_CLASS_I_COMPLEX_CELL_SURFACE_EXPRESSION | fragile | 2/5 | -1.9 | 0.0612 | 1 | 0 | 0.0663 | 1.32 |
| P75NTR_SIGNALS_VIA_NFKB | fragile | 2/5 | -2.42 | 0.055 | 1 | 0 | -0.0761 | 0.963 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
