# LIHC / SORAFENIB: robust biomarker report

Pre-clinical model: **licob**. Drug targets: see drug_drugTarget.txt. Network: STRING > 700 (original precomputed proximity).

**Signature status: promising, not validated** (permutation p = 0.379, direction: as expected (resistant score -> worse survival), treatment-interaction p = 0.583, events in treated patients = 15 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nit | 50 | 278 | 25 | 15 | 0.0361 (SVR, k=10) | 0.837 | 1.42 (0.743-2.7) | 0.291 | 0.546 | 0.583 | 0.36 |
| nit_nodeconf | 50 | 278 | 25 | 15 | 0.0351 (SVR, k=2) | 0.94 | 1.62 (0.715-3.67) | 0.248 | 0.537 | 0.101 | 0.571 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.379 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.0361, empirical p = 0.383

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.0922 (p = 0.524)
- Original Ridge top-7: Spearman rho = 0.174 (p = 0.226)

## Fragility

- Leave-one-organoid-out: HR > 1 in 98% of refits (range 0.733-2.18)
- Patient bootstrap: adjusted HR 95% interval 0.463-4.81, HR > 1 in 88.5%
- Proximity cut-off sensitivity:
  - z <= -1.0: 185 pathways, adj HR 1.72 (p 0.105), signature: DOWNSTREAM_SIGNAL_TRANSDUCTION, SIGNAL_REGULATORY_PROTEIN_SIRP_FAMILY_INTERACTIONS, SOS_MEDIATED_SIGNALLING
  - z <= -1.2816: 167 pathways, adj HR 1.48 (p 0.235), signature: SIGNALING_BY_NOTCH4, DOWNSTREAM_SIGNAL_TRANSDUCTION, SIGNAL_REGULATORY_PROTEIN_SIRP_FAMILY_INTERACTIONS, SOS_MEDIATED_SIGNALLING
  - z <= -1.645: 147 pathways, adj HR 1.46 (p 0.248), signature: SIGNALING_BY_NOTCH4, DOWNSTREAM_SIGNAL_TRANSDUCTION, SIGNAL_REGULATORY_PROTEIN_SIRP_FAMILY_INTERACTIONS, SOS_MEDIATED_SIGNALLING

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| SIGNAL_REGULATORY_PROTEIN_SIRP_FAMILY_INTERACTIONS | robust | 5/5 | -4.48 | 0.843 | 1 | 1 | -0.338 | 0.645 |
| DOWNSTREAM_SIGNAL_TRANSDUCTION | robust | 5/5 | -2.92 | 0.664 | 1 | 0.98 | -0.245 | 0.281 |
| SIGNALING_BY_NOTCH4 | supported | 4/5 | -2.06 | 0.777 | 1 | 1 | 0.284 | 0.568 |
| SOS_MEDIATED_SIGNALLING | fragile | 2/5 | -4.13 | 0.555 | 1 | 0.68 | 0.242 | 0.796 |
| SEMA3A_PAK_DEPENDENT_AXON_REPULSION | fragile | 2/5 | -3.09 | 0.439 | 1 | 0 | 0.205 | 1.01 |
| SHC1_EVENTS_IN_EGFR_SIGNALING | fragile | 2/5 | -4.02 | 0.367 | 1 | 0.08 | -0.248 | 1.47 |
| DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING | fragile | 2/5 | -3.42 | 0.36 | 1 | 0.02 | 0.232 | 0.494 |
| REGULATION_OF_KIT_SIGNALING | fragile | 2/5 | -4.03 | 0.356 | 1 | 0 | -0.185 | 0.446 |
| SIGNALING_BY_PDGF | fragile | 2/5 | -2.82 | 0.34 | 1 | 0.04 | 0.173 | 0.319 |
| P130CAS_LINKAGE_TO_MAPK_SIGNALING_FOR_INTEGRINS | fragile | 2/5 | -2.99 | 0.247 | 1 | 0 | -0.197 | 0.547 |
| SIGNALING_BY_RHO_GTPASES | fragile | 2/5 | -2.78 | 0.24 | 1 | 0 | -0.169 | 0.321 |
| GRB2_SOS_PROVIDES_LINKAGE_TO_MAPK_SIGNALING_FOR_INTERGRINS_ | fragile | 2/5 | -3.35 | 0.24 | 1 | 0 | -0.186 | 0.404 |
| RAF_MAP_KINASE_CASCADE | fragile | 2/5 | -3.73 | 0.08 | 1 | 0 | 0.14 | 1.2 |
| PECAM1_INTERACTIONS | fragile | 2/5 | -3.96 | 0.0775 | 1 | 0 | -0.163 | 0.496 |
| TIE2_SIGNALING | fragile | 2/5 | -3.92 | 0.0712 | 1 | 0 | -0.145 | 0.342 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
