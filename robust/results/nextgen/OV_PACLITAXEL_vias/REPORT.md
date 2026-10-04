# OV / PACLITAXEL: robust biomarker report

Pre-clinical model: **vias**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 1, direction: inconsistent, treatment-interaction p = 0.244, events in treated patients = 191).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 14 | 322 | 307 | 191 | 0.00776 (Ridge, k=6) | 0.11 | 0.873 (0.757-1.01) | 0.0627 | 0.465 | 0.244 | 0.829 |
| nit_nodeconf | 14 | 322 | 307 | 191 | 0.158 (Ridge, k=3) | 0.134 | 1.16 (1.01-1.33) | 0.0382 | 0.537 | 0.63 | 0.928 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 1 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.00776, empirical p = 0.168

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.24 (p = 0.409)
- Original Ridge top-7: Spearman rho = 0.284 (p = 0.326)

## Fragility

- Leave-one-organoid-out: HR > 1 in 28.6% of refits (range 0.868-1.14)
- Patient bootstrap: adjusted HR 95% interval 0.746-0.979, HR > 1 in 1%
- Proximity cut-off sensitivity:
  - z <= -1.0: 126 pathways, adj HR 0.826 (p 0.00812), signature: HDL_MEDIATED_LIPID_TRANSPORT, CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_RAS, CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_CAMKII, RAP1_SIGNALLING, GLUCOSE_METABOLISM
  - z <= -1.2816: 86 pathways, adj HR 0.879 (p 0.0881), signature: ERKS_ARE_INACTIVATED, HDL_MEDIATED_LIPID_TRANSPORT, NRAGE_SIGNALS_DEATH_THROUGH_JNK, PREFOLDIN_MEDIATED_TRANSFER_OF_SUBSTRATE_TO_CCT_TRIC, PROTEIN_FOLDING, RAP1_SIGNALLING
  - z <= -1.645: 49 pathways, adj HR 0.821 (p 0.00877), signature: PREFOLDIN_MEDIATED_TRANSFER_OF_SUBSTRATE_TO_CCT_TRIC, PROTEIN_FOLDING, CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_RAS, CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_CAMKII, RAP1_SIGNALLING

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT | robust | 5/5 | -0.469 | 0.75 | 1 | 1 | 0.276 | 1.06 |
| MTORC1_MEDIATED_SIGNALLING | robust | 5/5 | -0.495 | 0.739 | 1 | 1 | 0.251 | 1.13 |
| ACTIVATED_TAK1_MEDIATES_P38_MAPK_ACTIVATION | supported | 4/5 | -0.304 | 0.495 | 1 | 0.929 | 0.187 | 1.01 |
| REGULATION_OF_INSULIN_SECRETION_BY_GLUCAGON_LIKE_PEPTIDE1 | supported | 4/5 | -0.609 | 0.46 | 1 | 0.786 | 0.172 | 1.19 |
| CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_RAS | supported | 3/5 | -2.07 | 0.606 | 1 | 0.714 | -0.366 | 1.25 |
| RAP1_SIGNALLING | fragile | 2/5 | -2.41 | 0.599 | 1 | 0.643 | -0.282 | 1.17 |
| CREB_PHOSPHORYLATION_THROUGH_THE_ACTIVATION_OF_CAMKII | fragile | 2/5 | -1.72 | 0.554 | 1 | 0.143 | -0.316 | 1.09 |
| SYNTHESIS_OF_PC | fragile | 2/5 | -2.65 | 0.181 | 1 | 0 | 0.182 | 1.02 |
| NUCLEAR_SIGNALING_BY_ERBB4 | fragile | 2/5 | -2.46 | 0.136 | 1 | 0 | -0.19 | 0.947 |
| PREFOLDIN_MEDIATED_TRANSFER_OF_SUBSTRATE_TO_CCT_TRIC | fragile | 2/5 | -2.87 | 0.0712 | 1 | 0 | 0.105 | 1.02 |
| ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING | fragile | 2/5 | -1.98 | 0.0712 | 1 | 0 | 0.147 | 1.09 |
| SHC1_EVENTS_IN_ERBB4_SIGNALING | fragile | 2/5 | -1.25 | 0.0563 | 1 | 0 | -0.192 | 0.992 |
| RECYCLING_OF_BILE_ACIDS_AND_SALTS | fragile | 2/5 | -2.43 | 0.0362 | 1 | 0 | 0.0721 | 0.916 |
| G2_M_DNA_DAMAGE_CHECKPOINT | fragile | 2/5 | -1.78 | 0.0312 | 1 | 0 | -0.0476 | 0.923 |
| MITOTIC_G2_G2_M_PHASES | fragile | 2/5 | -1.8 | 0.0288 | 1 | 0 | -0.118 | 0.987 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
