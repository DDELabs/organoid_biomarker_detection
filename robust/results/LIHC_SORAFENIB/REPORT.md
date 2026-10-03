# LIHC / SORAFENIB: robust biomarker report

Pre-clinical model: **patient-derived organoids (LICOB)**. Drug targets: BRAF, FGFR1, FLT1, FLT3, FLT4, KDR, KIT, PDGFRB, RAF1, RET, ZHX2. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.95, direction: reversed (resistant score -> better survival), treatment-interaction p = 0.941, events in treated patients = 15 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 50 | 167 | 26 | 15 | 0.00328 (Ridge, k=6) | 0.025 | 0.389 (0.183-0.826) | 0.014 | 0.225 | 0.941 | 0.416 |
| current_rank | 50 | 167 | 26 | 15 | 0.00965 (Ridge, k=4) | 0.0118 | 0.389 (0.162-0.938) | 0.0354 | 0.257 | 0.212 | 0.38 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.95 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.00328, empirical p = 0.246

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.106 (p = 0.463)
- Original Ridge top-7: Spearman rho = 0.274 (p = 0.0541)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.37-0.837)
- Patient bootstrap: adjusted HR 95% interval 0.02-0.963, HR > 1 in 2.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 184 pathways, adj HR 0.837 (p 0.495), signature: CS_DS_DEGRADATION
  - z <= -1.2816: 167 pathways, adj HR 0.389 (p 0.014), signature: DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING, CS_DS_DEGRADATION
  - z <= -1.645: 147 pathways, adj HR 0.392 (p 0.0146), signature: DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING, CS_DS_DEGRADATION

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| CS_DS_DEGRADATION | supported | 4/5 | -2 | 0.656 | 1 | 1 | 0.202 | 0.82 |
| DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING | supported | 4/5 | -3.42 | 0.586 | 1 | 0.92 | 0.193 | 0.441 |
| TIE2_SIGNALING | supported | 3/5 | -3.92 | 0.236 | 1 | 0 | -0.175 | 0.411 |
| SIGNAL_REGULATORY_PROTEIN_SIRP_FAMILY_INTERACTIONS | supported | 3/5 | -4.48 | 0.165 | 1 | 0 | -0.0763 | 0.698 |
| REGULATION_OF_KIT_SIGNALING | fragile | 2/5 | -4.03 | 0.431 | 1 | 0.06 | -0.115 | 0.525 |
| AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE | fragile | 2/5 | -1.55 | 0.275 | 1 | 0 | -0.189 | 1.87 |
| HYALURONAN_UPTAKE_AND_DEGRADATION | fragile | 2/5 | -2.31 | 0.23 | 1 | 0.02 | 0.176 | 1.02 |
| DOWNSTREAM_SIGNAL_TRANSDUCTION | fragile | 2/5 | -2.92 | 0.171 | 1 | 0 | -0.127 | 0.262 |
| INSULIN_RECEPTOR_SIGNALLING_CASCADE | fragile | 2/5 | -3.34 | 0.152 | 1 | 0 | -0.126 | 0.324 |
| IL_7_SIGNALING | fragile | 2/5 | -3.18 | 0.149 | 1 | 0 | -0.0843 | 0.405 |
| PLATELET_SENSITIZATION_BY_LDL | fragile | 2/5 | -1.94 | 0.147 | 1 | 0 | -0.118 | 0.618 |
| NEURONAL_SYSTEM | fragile | 2/5 | -1.94 | 0.0875 | 1 | 0 | 0.128 | 1.42 |
| NETRIN1_SIGNALING | fragile | 2/5 | -2.68 | 0.0762 | 1 | 0 | -0.0538 | 0.512 |
| PI3K_CASCADE | fragile | 2/5 | -2.12 | 0.0762 | 1 | 0 | -0.114 | 0.358 |
| SIGNALING_BY_CONSTITUTIVELY_ACTIVE_EGFR | fragile | 2/5 | -3.31 | 0.07 | 1 | 0 | -0.152 | 0.752 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
