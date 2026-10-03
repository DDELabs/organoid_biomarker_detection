# LUAD / CISPLATIN: robust biomarker report

Pre-clinical model: **cell-line proxy (GDSC 2012)**. Drug targets: A2M, ATOX1, MPG, TF. Network: STRING v12 >= 700.

**Signature status: not supported** (permutation p = 0.88, direction: inconsistent, treatment-interaction p = 0.751, events in treated patients = 27 - underpowered, fewer than 30 events).

## Data versions

| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current_ssgsea | 65 | 69 | 83 | 27 | 0.147 (Ridge, k=6) | 0.797 | 0.81 (0.542-1.21) | 0.306 | 0.474 | 0.751 | 0.89 |
| current_rank | 65 | 69 | 83 | 27 | 0.0792 (LinearRegression, k=3) | 0.94 | 1.01 (0.629-1.61) | 0.98 | 0.506 | 0.497 | 0.914 |

## Honest significance (whole-pipeline permutation of organoid IC50)

- Robust signature, adjusted Cox HR: empirical p = 0.88 (500 permutations)
- Original approach (best of 3 models x k = 2..10): observed min p = 0.147, empirical p = 0.812

## Organoid cross-validation (leave one organoid out, selection inside folds)

- Robust learner: Spearman rho = 0.342 (p = 0.00535)
- Original Ridge top-7: Spearman rho = 0.104 (p = 0.411)

## Fragility

- Leave-one-organoid-out: HR > 1 in 0% of refits (range 0.714-0.861)
- Patient bootstrap: adjusted HR 95% interval 0.506-1.3, HR > 1 in 20%
- Proximity cut-off sensitivity:
  - z <= -1.0: 88 pathways, adj HR 0.854 (p 0.434), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS, EGFR_DOWNREGULATION, AMYLOIDS
  - z <= -1.2816: 69 pathways, adj HR 0.81 (p 0.306), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS, HOST_INTERACTIONS_OF_HIV_FACTORS, AMYLOIDS
  - z <= -1.645: 43 pathways, adj HR 0.795 (p 0.272), signature: RECYCLING_OF_BILE_ACIDS_AND_SALTS, HOST_INTERACTIONS_OF_HIV_FACTORS, AMYLOIDS

## Pathway scorecard (top 15)

| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |
|---|---|---|---|---|---|---|---|---|
| HOST_INTERACTIONS_OF_HIV_FACTORS | supported | 4/5 | -1.7 | 0.733 | 1 | 1 | -0.387 | 1.28 |
| RECYCLING_OF_BILE_ACIDS_AND_SALTS | supported | 4/5 | -2.73 | 0.713 | 1 | 1 | -0.459 | 1.14 |
| AMYLOIDS | supported | 3/5 | -2.34 | 0.55 | 1 | 0.677 | 0.299 | 1.05 |
| ABCA_TRANSPORTERS_IN_LIPID_HOMEOSTASIS | supported | 3/5 | -1.97 | 0.357 | 1 | 0.0154 | 0.239 | 1.38 |
| RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_ | fragile | 2/5 | -2.99 | 0.193 | 1 | 0 | -0.206 | 0.969 |
| IL_6_SIGNALING | fragile | 2/5 | -1.65 | 0.179 | 1 | 0 | 0.186 | 1.1 |
| GAMMA_CARBOXYLATION_TRANSPORT_AND_AMINO_TERMINAL_CLEAVAGE_OF_PROTEINS | fragile | 2/5 | -2.96 | 0.146 | 0.966 | 0 | 0.1 | 1.01 |
| PEPTIDE_HORMONE_BIOSYNTHESIS | fragile | 2/5 | -1.34 | 0.145 | 1 | 0 | 0.073 | 1.1 |
| SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX | fragile | 2/5 | -1.41 | 0.144 | 1 | 0 | 0.186 | 1.41 |
| SYNTHESIS_AND_INTERCONVERSION_OF_NUCLEOTIDE_DI_AND_TRIPHOSPHATES | fragile | 2/5 | -1.82 | 0.131 | 1 | 0 | 0.137 | 1.39 |
| PTM_GAMMA_CARBOXYLATION_HYPUSINE_FORMATION_AND_ARYLSULFATASE_ACTIVATION | fragile | 2/5 | -1.99 | 0.121 | 0.959 | 0 | 0.136 | 1.09 |
| SIGNALING_BY_RHO_GTPASES | fragile | 2/5 | -1.44 | 0.1 | 0.975 | 0 | -0.13 | 0.794 |
| HEMOSTASIS | fragile | 2/5 | -2.13 | 0.065 | 1 | 0 | 0.108 | 1.03 |
| CIRCADIAN_CLOCK | fragile | 2/5 | -1.52 | 0.0575 | 1 | 0 | -0.0383 | 0.894 |
| ABACAVIR_TRANSPORT_AND_METABOLISM | fragile | 2/5 | -1.37 | 0.0512 | 1 | 0 | 0.131 | 1.28 |

Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, consistent in every data version, patient HR direction matches organoid coefficient. Tier robust = 5/5, supported = 3-4/5.

Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`.
