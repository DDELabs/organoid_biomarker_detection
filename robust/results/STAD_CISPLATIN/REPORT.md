# STAD / CISPLATIN: robust biomarker report

Pre-clinical model: **cell-line proxy (GDSC 2012)**. Drug targets: A2M, ATOX1, MPG, TF. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.597, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.642, events in treated patients = 16 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 16 | 69 | 46 | 16 | 0.175 (LinearRegression, k=7) | 0.613 | 0.935 (0.546-1.6) | 0.806 | 0.51 | 0.642 | 0.919 |
| current_rank | 16 | 69 | 46 | 16 | 0.32 (Ridge, k=10) | 0.334 | 0.848 (0.48-1.5) | 0.569 | 0.508 | 0.894 | 0.852 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.597 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.175, empirical p = 0.693

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.424 (p = 0.102)
- Original Ridge top-7: Spearman rho = 0.182 (p = 0.499)

## Fragility

- Leave-one-organoid-out: HR > 1 in 6.25% of refits (range 0.827-1.12)
- Patient bootstrap: adjusted HR 95% interval 0.26-2.01, HR > 1 in 36.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 88 pathways, adj HR 0.837 (p 0.497), signature: TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS, TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX
  - z <= -1.2816: 69 pathways, adj HR 0.935 (p 0.806), signature: SIGNALING_BY_RHO_GTPASES, TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS
  - z <= -1.645: 43 pathways, adj HR 0.944 (p 0.839), signature: TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS, METAL_ION_SLC_TRANSPORTERS, LIPID_DIGESTION_MOBILIZATION_AND_TRANSPORT

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| METAL_ION_SLC_TRANSPORTERS | supported | 4/5 | -3.77 | 0.974 | 1 | 1 | 0.33 | 0.945 |
| TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS | supported | 4/5 | -2.42 | 0.906 | 1 | 1 | 0.28 | 0.849 |
| SIGNALING_BY_RHO_GTPASES | supported | 3/5 | -1.44 | 0.507 | 1 | 0.188 | -0.162 | 0.851 |
| TRAF6_MEDIATED_NFKB_ACTIVATION | fragile | 2/5 | -1.34 | 0.376 | 1 | 0.0625 | -0.163 | 0.9 |
| NRIF_SIGNALS_CELL_DEATH_FROM_THE_NUCLEUS | fragile | 2/5 | -1.51 | 0.279 | 1 | 0.0625 | 0.101 | 1.52 |
| DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY | fragile | 2/5 | -1.42 | 0.254 | 1 | 0 | -0.147 | 0.718 |
| LIPID_DIGESTION_MOBILIZATION_AND_TRANSPORT | fragile | 2/5 | -1.72 | 0.185 | 1 | 0 | 0.124 | 3.01 |
| LIPOPROTEIN_METABOLISM | fragile | 2/5 | -2.22 | 0.15 | 1 | 0 | 0.0878 | 1.54 |
| TRANSFERRIN_ENDOCYTOSIS_AND_RECYCLING | fragile | 2/5 | -2.47 | 0.142 | 1 | 0 | 0.123 | 1.61 |
| NEF_MEDIATED_DOWNREGULATION_OF_MHC_CLASS_I_COMPLEX_CELL_SURFACE_EXPRESSION | fragile | 2/5 | -1.9 | 0.102 | 1 | 0 | -0.0437 | 0.633 |
| HOST_INTERACTIONS_OF_HIV_FACTORS | fragile | 2/5 | -1.7 | 0.0638 | 1 | 0 | -0.0301 | 0.916 |
| HEMOSTASIS | fragile | 2/5 | -2.13 | 0.0475 | 1 | 0 | 0.0601 | 1.05 |
| REGULATED_PROTEOLYSIS_OF_P75NTR | fragile | 2/5 | -1.49 | 0.0462 | 1 | 0 | 0.0663 | 1.41 |
| ACTIVATION_OF_CHAPERONES_BY_ATF6_ALPHA | fragile | 2/5 | -1.31 | 0.0425 | 1 | 0 | -0.0761 | 0.954 |
| SIGNALING_BY_NOTCH2 | fragile | 2/5 | -1.94 | 0.04 | 1 | 0 | 0.0867 | 1.16 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
