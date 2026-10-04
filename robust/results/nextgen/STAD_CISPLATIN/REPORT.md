# STAD / CISPLATIN: robust biomarker report

Pre-clinical model: **gdsc**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: not supported** (permutation p = 0.571, direction: inconsistent, treatment-interaction p = 0.758, events in treated patients = 16 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 16 | 279 | 45 | 16 | 0.279 (SVR, k=9) | 0.46 | 0.935 (0.546-1.6) | 0.806 | 0.526 | 0.758 | 1.11 |
| nit_nodeconf | 16 | 279 | 45 | 16 | 0.17 (SVR, k=3) | 0.524 | 1.1 (0.638-1.91) | 0.724 | 0.545 | 0.817 | 1.1 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.571 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.279, empirical p = 0.848

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = -0.285 (p = 0.284)
- Original Ridge top-7: Spearman rho = 0.0912 (p = 0.737)

## Fragility

- Leave-one-organoid-out: HR > 1 in 50% of refits (range 0.885-1.11)
- Patient bootstrap: adjusted HR 95% interval 0.199-1.99, HR > 1 in 20.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 88 pathways, adj HR 0.937 (p 0.811), signature: STEROID_HORMONES, TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS, LIPOPROTEIN_METABOLISM
  - z <= -1.2816: 69 pathways, adj HR 0.939 (p 0.82), signature: STEROID_HORMONES, TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS, LIPOPROTEIN_METABOLISM
  - z <= -1.645: 43 pathways, adj HR 1.01 (p 0.963), signature: NEF_MEDIATED_DOWNREGULATION_OF_MHC_CLASS_I_COMPLEX_CELL_SURFACE_EXPRESSION, STEROID_HORMONES, TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS, LIPOPROTEIN_METABOLISM

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| METAL_ION_SLC_TRANSPORTERS | supported | 4/5 | -3.77 | 0.996 | 1 | 1 | 0.372 | 0.859 |
| TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS | supported | 4/5 | -2.42 | 0.966 | 1 | 1 | 0.259 | 0.885 |
| LIPOPROTEIN_METABOLISM | supported | 4/5 | -2.22 | 0.63 | 1 | 0.562 | 0.137 | 1.19 |
| NEF_MEDIATED_DOWNREGULATION_OF_MHC_CLASS_I_COMPLEX_CELL_SURFACE_EXPRESSION | supported | 3/5 | -1.9 | 0.44 | 1 | 0.0625 | -0.0865 | 0.655 |
| STEROID_HORMONES | fragile | 2/5 | -1.73 | 0.635 | 1 | 0.5 | -0.168 | 1.22 |
| INTRINSIC_PATHWAY | fragile | 2/5 | -2.9 | 0.338 | 1 | 0.0625 | -0.119 | 0.965 |
| IRON_UPTAKE_AND_TRANSPORT | fragile | 2/5 | -3.12 | 0.14 | 0.839 | 0 | 0.0455 | 2.06 |
| REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS | fragile | 2/5 | -2.36 | 0.12 | 0.938 | 0 | -0.057 | 0.481 |
| HEMOSTASIS | fragile | 2/5 | -2.13 | 0.119 | 1 | 0 | 0.0699 | 1.05 |
| PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION | fragile | 2/5 | -1.99 | 0.116 | 1 | 0 | 0.0293 | 1.23 |
| LIPID_DIGESTION_MOBILIZATION_AND_TRANSPORT | fragile | 2/5 | -1.72 | 0.0925 | 1 | 0 | 0.0804 | 1.67 |
| SIGNALING_BY_NOTCH2 | fragile | 2/5 | -1.94 | 0.0775 | 1 | 0 | 0.0576 | 1.16 |
| SIGNALING_BY_RHO_GTPASES | fragile | 2/5 | -1.44 | 0.0762 | 1 | 0 | -0.0645 | 0.866 |
| PTM_GAMMA_CARBOXYLATION_HYPUSINE_FORMATION_AND_ARYLSULFATASE_ACTIVATION | fragile | 2/5 | -1.99 | 0.0413 | 1 | 0 | 0.0587 | 1.33 |
| CIRCADIAN_CLOCK | fragile | 2/5 | -1.52 | 0.035 | 1 | 0 | 0.0214 | 1.43 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
