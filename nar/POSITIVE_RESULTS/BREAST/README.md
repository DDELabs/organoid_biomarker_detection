# Breast cancer: chemotherapy response is predicted in cohorts the model has never seen

## What was tested
Can a model that reads 673 Reactome pathway scores from a breast tumour predict whether the patient responded to chemotherapy (responder = 1 as coded in each cohort; the response definitions differ between cohorts)? The model is trained on all other breast cohorts and tested on one cohort it never saw (leave-one-cohort-out). All choices (penalty, pathway selection) are made on the training cohorts only.

## Data
- 13 original cohorts, **3,511 treated patients, 1,202 responders** (BRIGHTNESS, GSE123845, GSE163882, GSE20194, GSE20271, GSE22093, GSE23988, GSE240671, GSE25066, GSE32646, GSE41998, ISPY2, TCGA-BRCA).
- 15 new cohorts added later (W2): 1,327 patients, 416 responders.
- Pathway scores are z-scored inside each cohort. Model: L2-penalised logistic regression, penalty fixed at 100 (the original model). 17 other models were tried in W3.

## Result 1. Original 13 cohorts (rule met in 13 of 13)
Mean AUC 0.752 (95% CI 0.728-0.775). Per cohort (fixed L2, DeLong CI, p vs 0.5):

| cohort | n | responders | AUC (95% CI) | p vs 0.5 |
|---|---|---|---|---|
| BRIGHTNESS | 482 | 236 | 0.674 (0.627-0.722) | 7.0e-13 |
| GSE123845 | 112 | 42 | 0.761 (0.670-0.851) | 1.6e-08 |
| GSE163882 | 222 | 80 | 0.715 (0.646-0.784) | 8.7e-10 |
| GSE20194 | 270 | 52 | 0.798 (0.736-0.860) | 5.2e-21 |
| GSE20271 | 178 | 26 | 0.801 (0.725-0.876) | 5.4e-15 |
| GSE22093 | 97 | 28 | 0.751 (0.646-0.855) | 2.4e-06 |
| GSE23988 | 61 | 20 | 0.846 (0.739-0.953) | 2.2e-10 |
| GSE240671 | 95 | 37 | 0.780 (0.679-0.881) | 5.9e-08 |
| GSE25066 | 438 | 82 | 0.817 (0.767-0.867) | 2.3e-35 |
| GSE32646 | 115 | 27 | 0.741 (0.634-0.848) | 1.0e-05 |
| GSE41998 | 253 | 69 | 0.715 (0.644-0.785) | 2.5e-09 |
| ISPY2 | 986 | 319 | 0.700 (0.667-0.733) | 4.1e-32 |
| TCGA-BRCA | 202 | 184 | 0.672 (0.554-0.790) | 0.004 |

The p column here is the DeLong test of the W3 re-run (identical AUCs). The original run in `robust/results/train_combined` uses a one-sided Mann-Whitney test; its largest p is 0.016 (TCGA-BRCA). BH q < 0.10 in 13 of 13 (`POSITIVE_INDEX.tsv` BR1). Sources: `tables/robust__part1_response_leave_one_cohort_out.tsv` (original run), `tables/19__breast_loco_all_models.tsv` (with CIs).

## Result 2. 15 new cohorts (rule met in 8 of 15 by p, 9 of 15 by BH q)
Trained on all breast cohorts of the v2 table except the held-out one:

| cohort | n | responders | AUC (95% CI) | p | q (BH over 15) |
|---|---|---|---|---|---|
| GSE34138 | 177 | 28 | 0.867 (0.780-0.955) | 7.3e-10 | 1.1e-08 |
| GSE42822 | 91 | 37 | 0.801 (0.705-0.898) | 1.2e-06 | 8.8e-06 |
| GSE4779 | 102 | 39 | 0.744 (0.641-0.846) | 3.8e-05 | 1.6e-04 |
| GSE6861 | 161 | 66 | 0.690 (0.605-0.775) | 4.3e-05 | 1.6e-04 |
| GSE32603 | 138 | 35 | 0.717 (0.612-0.822) | 1.3e-04 | 4.0e-04 |
| GSE173839 | 100 | 36 | 0.727 (0.620-0.835) | 1.7e-04 | 4.3e-04 |
| GSE37946 | 50 | 27 | 0.718 (0.577-0.859) | 0.009 | 0.018 |
| GSE50948 | 156 | 53 | 0.617 (0.522-0.712) | 0.017 | 0.032 |
| GSE21974 | 32 | 10 | 0.718 (0.514-0.922) | 0.053 | 0.089 |
| GSE66305 | 88 | 27 | 0.613 (0.482-0.744) | 0.092 | 0.139 |
| GSE18864 | 24 | 8 | 0.711 (0.478-0.944) | 0.106 | 0.144 |
| GSE18728 | 21 | 8 | 0.683 (0.437-0.929) | 0.185 | 0.231 |
| GSE16446 | 114 | 16 | 0.569 (0.413-0.725) | 0.381 | 0.439 |
| GSE76360 | 48 | 9 | 0.581 (0.366-0.796) | 0.460 | 0.492 |
| GSE59515 | 25 | 17 | 0.471 (0.222-0.719) | 0.842 | 0.842 |

What I think: the new cohorts that did not reach p < 0.05 are small (21-114 patients, 7-17 responders) or have a weak AUC (GSE59515 0.471, GSE66305 0.613). The mean over all 15 is 0.682, lower than the original 13 (0.752). So the effect is real but not uniform. Source: `tables/18__part1_v2.tsv` (table `v2_all_groups`).

## Result 3. All-model benchmark (W3): the signal does not depend on the model type
Same 13 held-out cohorts, 17 model types, settings tuned only inside the training cohorts.

| model | mean AUC (95% CI), 13 cohorts | difference vs fixed L2 (p, q) |
|---|---|---|
| SVM_linear | 0.764 (0.741-0.786) | +0.012 (p 0.023, q 0.100) |
| NetSmoothL2 | 0.759 (0.737-0.782) | +0.008 (p 0.074, q 0.141) |
| L2 | 0.759 (0.737-0.782) | +0.008 (p 0.060, q 0.141) |
| SVM_RBF | 0.756 (0.732-0.779) | +0.004 (p 0.577, q 0.613) |
| L1 | 0.755 (0.732-0.778) | +0.004 (p 0.348, q 0.394) |
| L2_fixed | 0.752 (0.728-0.775) | reference |
| ElasticNet | 0.751 (0.728-0.774) | -0.000 (p 0.943, q 0.943) |
| L2_shortlist | 0.744 (0.721-0.767) | -0.008 (p 0.297, q 0.389) |
| XGBoost | 0.743 (0.719-0.766) | -0.009 (p 0.327, q 0.394) |
| L2_fixed_shortlist | 0.741 (0.718-0.765) | -0.010 (p 0.152, q 0.236) |
| L1_shortlist | 0.741 (0.717-0.764) | -0.011 (p 0.150, q 0.236) |
| RandomForest | 0.740 (0.716-0.763) | -0.012 (p 0.215, q 0.305) |
| ElasticNet_shortlist | 0.738 (0.714-0.761) | -0.014 (p 0.070, q 0.141) |
| SVM_linear_shortlist | 0.735 (0.711-0.759) | -0.016 (p 0.031, q 0.105) |
| MLP | 0.734 (0.709-0.759) | -0.018 (p 0.039, q 0.109) |
| ER_HER2_only | 0.693 (0.666-0.720) | -0.065 (p 1.8e-05, q 1.0e-04) |
| GCN_prior | 0.677 (0.651-0.703) | -0.075 (p 2.3e-08, q 2.0e-07) |
| Proliferation_score | 0.646 (0.619-0.673) | -0.106 (p 1.6e-13, q 2.7e-12) |

- Baselines: ER+HER2 flags alone 0.693 (9 cohorts); proliferation score alone 0.646. The pathway model beats both (paired p < 0.001).
- 200 random sets of 183 pathways (the size of the network shortlist) average 0.742; the network shortlist gives 0.744. The shortlist adds nothing over random pathways.
- The best model (linear SVM, 0.764) is +0.012 over fixed L2 with p 0.023 but BH q 0.100, so it is a small and not clear gain. Source: `tables/19__model_comparison.tsv`, `tables/19__breast_loco_all_models.tsv`, `tables/19__random_pathway_sets_auc.tsv`; figures `19__fig_model_mean_auc`, `19__fig_model_heatmap`, `19__fig_roc_top3`.

## Result 4. Which pathways drive it (W4)
| k pathways | chosen on training cohorts by one-at-a-time strength | chosen by L2 weight | random sets |
|---|---|---|---|
| 1 | 0.680 | 0.540 | 0.561 |
| 5 | 0.689 | 0.666 | 0.638 |
| 20 | 0.710 | 0.722 | 0.696 |
| 100 | 0.734 | 0.744 | 0.730 |
| 200 | 0.740 | 0.750 | 0.743 |

- One pathway chosen on the training cohorts (EARLY_PHASE_OF_HIV_LIFE_CYCLE in 13 of 13 folds) already gives 0.680; five pathways chosen greedily plateau at about 0.72; about 100 or more are needed to be within 0.01 of the full model. Random sets of 100-200 pathways reach 0.730-0.743, so the information is spread over many correlated pathways.
- The top 30 weights are stable (28 of 30 selected in at least 80% of 100 bootstrap L1 fits) but are generic cell-machinery pathways (several are Reactome virus or drug pathways); weights of correlated pathways should not be read one by one.
- Not simply proliferation: proliferation alone gives 0.646; with proliferation regressed out of every pathway the model still gives 0.721 (0.694-0.747); dropping the 27 proliferation pathways gives 0.748. The model adds to proliferation + ER + HER2 in 11 of 13 cohorts (likelihood-ratio p < 0.05; Stouffer z 14.0).
- Subgroups (cohorts with enough responders and non-responders; receptor flags missing for GSE22093, GSE23988, TCGA-BRCA): ER+/HER2- AUC 0.670 (0.599-0.741), 1,111 patients, 8 cohorts; HER2+ 0.728 (0.669-0.788), 520 patients; triple-negative 0.720 (0.678-0.761), 1,470 patients. All three intervals exclude 0.5. In ER+/HER2- proliferation alone gives 0.627, so more of the signal there is proliferation.
Sources: `tables/20__auc_vs_k.tsv`, `20__top30_pathways.tsv`, `20__stability_selection.tsv`, `20__proliferation_tests.tsv`, `20__subgroup_auc_pooled.tsv`, `20__minimal_subset.tsv`, `20__family_weights.tsv`; figures `20__fig_*`.

## Caveats
- Cohorts mix drugs and response definitions (pCR, clinical response, RECIST-like). The test is "predict response in a new cohort", not "predict response to one named drug".
- TCGA-BRCA has 184 responders of 202, so its AUC (0.672) rests on only 18 non-responders.
- The network-proximity shortlist is from the paper's table, not a new STRING v12 run.
- Randomised HER2 trials (NOAH, CHER-LOB) do not show a drug-specific effect and are in `../BORDERLINE.md`, not here.

## What I think
This is the clearest positive in the project. A pathway model learned in other breast cohorts predicts chemotherapy response in 13 of 13 cohorts it never saw, and in 8 of 15 further cohorts, with AUC around 0.7-0.8 in the larger ones. It is a general chemotherapy-response signal that is not just proliferation or receptor status. It is not a drug-specific signal, and it does not need the network-proximity step: a random 183-pathway set does as well.

## Where things are
`figures/` (PNG + PDF + caption), `tables/` (copies of the source tables), `source_reports/` (copies of the agents' reports). Copy list with checksums: `../COPY_MANIFEST.tsv`.
