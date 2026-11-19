# W3 positives: all-model benchmark on every positive result

All numbers come from `tables/positives_all_models.tsv` (written by `code/pos_analyze.py`). Seeds fixed (20261119). No git commit made. Breast outputs untouched.

## What was done (plain words)
- Same model set as in the breast report: L2, L1, elastic net, linear SVM, RBF SVM, random forest, XGBoost, small MLP, network-smoothed L2, graph network (GCN with a network-proximity prior), and the four linear models again on the network-proximity shortlist only. Baselines: a proliferation score (mean of cell-cycle pathways, sign learned in training data), and 200 random pathway sets of the shortlist size (ridge).
- Every model is trained only on that result's original training data (lab models, or training cohorts) and tuned only by cross-validation inside it (leave-one-out for 6-19 lab models; cohort-grouped 4-fold for patient cohorts). The test patients are never used for training or tuning. The test sets, tests and pooling are the ones used originally.
- The models here are regressors. Targets: organoid drug response (organoid tasks); a martingale residual (how much worse than expected a patient did, from a Cox model) for patient survival cohorts; non-response for response cohorts.
- Main score per task: C-index (survival; 0.5 = chance) or AUC (mutation / response). The original test (log-rank, Cox HR, one-sided Mann-Whitney) is also in the table. 95% CI and the paired difference to the original method come from a bootstrap that resamples the same patients for every model (500-1000 resamples; lung 200; ovarian survival 300). q = Benjamini-Hochberg across models within a task.
- Original methods reproduced: bladder paper data log-rank p 0.044 (rank-1) and 0.031 (amino-acid pathway); lung pooled HR 1.363 (1.290-1.439) and shortlist 1.207 (17_lung: 1.36 and 1.21); rectal top-10 signature HR 1.87 (1.42-2.48; 05_A3_B1 reported 1.91, small difference from tie handling); olaparib correlation signature log-rank p 0.0045, HR 2.09 (22: 0.0044, 2.09). Ovarian joint original scores are read from the saved scores of 03_organoid_plus_patient.

## Result per positive (original method vs best added model)
| task | original method | n / events | original score (95% CI) | original test p | best added model | its score (95% CI) | paired difference to original |
|---|---|---|---|---|---|---|---|
| lung_untreated_prognosis | Cox_ridge_all | 2938 / 1320 | 0.601 (0.585-0.619) | 1.1e-28 | NetSmoothL2 | 0.622 (0.606-0.638) | +0.021 (+0.009 to +0.033), p=0.01, q=0.02 |
| bladder_ERCC2_currentdata | Paper_AA_pathway | 397 / 36 | 0.665 (0.571-0.752) | 0.00054 | MLP | 0.637 (0.551-0.724) | -0.028 (-0.136 to +0.080), p=0.63, q=0.63 |
| bladder_ERCC2_paperdata | Paper_AA_pathway | 402 / 36 | 0.639 (0.533-0.739) | 0.003 | RandomForest | 0.619 (0.532-0.708) | -0.019 (-0.139 to +0.109), p=0.73, q=0.78 |
| bladder_logrank_currentdata | Paper_AA_pathway | 79 / 28 | 0.572 (0.448-0.685) | 0.063 | L1_shortlist | 0.589 (0.466-0.712) | +0.017 (-0.159 to +0.180), p=0.79, q=0.9 |
| bladder_logrank_paperdata | Paper_AA_pathway | 79 / 20 | 0.678 (0.552-0.808) | 0.031 | L1_shortlist | 0.632 (0.472-0.795) | -0.046 (-0.190 to +0.092), p=0.51, q=0.52 |
| colon_cetuximab_BRAF | Paper_rank1 | 573 / 51 | 0.640 (0.563-0.719) | 0.00048 | L1_shortlist | 0.704 (0.632-0.770) | +0.064 (+0.011 to +0.113), p=0.018, q=0.036 |
| ovarian_paclitaxel_joint_survival | Joint_original | 564 / 426 | 0.559 (0.530-0.587) | 0.005 | SVM_RBF | 0.579 (0.551-0.607) | +0.020 (-0.006 to +0.046), p=0.13, q=0.27 |
| ovarian_paclitaxel_joint_response | Joint_original | 276 / 86 | 0.578 (0.483-0.660) | 0.034 | SVM_RBF | 0.659 (0.579-0.742) | +0.081 (-0.005 to +0.171), p=0.088, q=0.35 |
| ovarian_platinum_joint_survival | Joint_original | 683 / 523 | 0.583 (0.556-0.614) | 6e-08 | SVM_RBF | 0.585 (0.557-0.616) | +0.002 (-0.014 to +0.019), p=0.86, q=0.92 |
| ovarian_platinum_joint_response | Joint_original | 353 / 112 | 0.692 (0.625-0.745) | 1.5e-07 | SVM_RBF | 0.682 (0.617-0.748) | -0.009 (-0.054 to +0.042), p=0.78, q=0.81 |
| ovarian_olaparib_GSE63885 | Corr_signature | 75 / 66 | 0.586 (0.503-0.662) | 0.0045 | XGBoost | 0.580 (0.494-0.656) | -0.006 (-0.080 to +0.067), p=0.86, q=0.86 |
| ovarian_olaparib_GSE9891 | Corr_signature | 241 / 110 | 0.542 (0.478-0.603) | 0.54 | SVM_linear_shortlist | 0.538 (0.470-0.600) | -0.004 (-0.054 to +0.049), p=0.81, q=0.81 |
| ovarian_olaparib_TCGAOV | Corr_signature | 368 / 224 | 0.528 (0.485-0.573) | 0.65 | SVM_linear_shortlist | 0.527 (0.482-0.570) | -0.001 (-0.049 to +0.047), p=0.99, q=0.99 |
| rectal_GSE87211_5FU | Paper_top10_frozen | 189 / 42 | 0.684 (0.597-0.767) | 1.1e-05 | SVM_linear | 0.665 (0.578-0.752) | -0.019 (-0.081 to +0.041), p=0.51, q=0.54 |

What I think (per result):
- **Lung untreated prognosis.** Several models beat the original ridge Cox: NetSmoothL2 and plain ridge L2 reach C 0.622 and 0.620 vs 0.601, pooled HR per SD 1.44 vs 1.36; paired differences +0.021 and +0.019, bootstrap p 0.01, q 0.02. Random forest is also higher (+0.016, p 0.04, q 0.064). So here a model does better than the original, by a small but fairly clear margin. The shortlist version (original 1.21) and the graph network are worse. Random shortlist-size sets reach the original in 44% of draws, so the shortlist adds nothing.
- **Bladder cisplatin survival (79 patients, 20 or 28 deaths).** No model clearly beats the paper's amino-acid pathway; every paired p is above 0.5. Only the original-type pathway scores and L1/elastic-net on the shortlist come close; all-pathway models sit near chance (C 0.47-0.52). Random shortlist-size sets average 0.49 (paper data). With nine organoids none of these differences are reliable.
- **Bladder ERCC2.** The amino-acid pathway (AUC 0.64 paper data, 0.665 current) is not beaten; the best added models (proliferation score 0.634, MLP 0.637, random forest) are close and their differences have wide intervals (about -0.14 to +0.12). Note the proliferation baseline does almost as well, so the signal is not specific to the network method here. Random sets: 6.5% / 0% reach the original.
- **Colon cetuximab vs BRAF V600E.** L1 on the shortlist is clearly better than the original rank-1 pathway: AUC 0.704 vs 0.640, difference +0.064 (+0.011 to +0.113), p 0.018, q 0.036. Elastic net on the shortlist is similar (+0.058, p 0.038, q 0.061). Caution: BRAF V600E is tied to tumour subtype (MSI, right side), shown earlier in 13_gaps_bio, so a better score here may be a better subtype detector.
- **Rectal GSE87211 5-FU (189 patients, 42 events).** The frozen top-10 signature (C 0.684) is not beaten; linear models trained on all pathways are 0.665 (difference -0.019, p about 0.5). The rank-1 version alone is clearly worse (0.573). Random shortlist-size sets reach the original in 0.5% of draws, so it is not just any 10 pathways.
- **Ovarian olaparib in GSE63885 (75 patients, 66 deaths).** The correlation signature is the only model with a log-rank signal (p 0.0045; C 0.586). Ridge, SVM and L1 trained on the same six lines give log-rank p 0.27-0.38. XGBoost has C 0.580 and others are 0.52-0.57, none significantly different from the original. In the two context cohorts (TCGA-OV, GSE9891) the signature has p 0.65 and 0.54, and so does every model. With six lines this stays exploratory.
- **Ovarian joint, platinum and paclitaxel.** Platinum: joint original C 0.583 survival, AUC 0.692 response; no model differs clearly (best RBF SVM +0.002; NetSmoothL2 -0.018, p 0.09). Paclitaxel: joint original C 0.559; RBF SVM 0.579 (+0.020, p 0.13), so no clear gain. Paclitaxel response: the original joint (0.578) is below its own patient-only model (0.668, +0.089, p 0.012, q 0.19) and RBF SVM (0.659); after correction nothing is clear. Patient-only and joint are almost the same, so lab models add little (as in 03).

## Model summary over all tasks (paired difference to the original method, 14 tasks)
| model | mean score | mean difference | tasks better / worse | clearly (p<0.05) better / worse | sign-test p |
|---|---|---|---|---|---|
| RandomForest | 0.567 | -0.043 | 3/11 | 1/2 | 0.0574 |
| L1_shortlist | 0.565 | -0.045 | 3/11 | 1/5 | 0.0574 |
| SVM_RBF | 0.563 | -0.048 | 4/10 | 0/1 | 0.18 |
| ElasticNet_shortlist | 0.561 | -0.049 | 3/11 | 1/7 | 0.0574 |
| L2_shortlist | 0.561 | -0.049 | 1/13 | 0/5 | 0.00183 |
| NetSmoothL2 | 0.560 | -0.050 | 1/13 | 1/2 | 0.00183 |
| XGBoost | 0.560 | -0.051 | 1/13 | 0/1 | 0.00183 |
| L2 | 0.555 | -0.055 | 1/13 | 1/4 | 0.00183 |
| SVM_linear_shortlist | 0.555 | -0.056 | 0/14 | 0/4 | 0.000122 |
| SVM_linear | 0.555 | -0.056 | 2/12 | 0/2 | 0.0129 |
| RandomSets_median | 0.552 | -0.058 | 1/13 | 0/0 | 0.00183 |
| ElasticNet | 0.551 | -0.059 | 2/12 | 0/5 | 0.0129 |
| L1 | 0.538 | -0.072 | 0/14 | 0/6 | 0.000122 |
| MLP | 0.533 | -0.078 | 0/14 | 0/5 | 0.000122 |
| Proliferation_score | 0.530 | -0.081 | 0/14 | 0/4 | 0.000122 |
| GCN_prior | 0.523 | -0.088 | 0/14 | 0/6 | 0.000122 |

What I think: on average no added model beats the original methods; the originals are each tailored to their task (a hand-picked pathway or signature), and most added models are worse on average. Exceptions are lung (ridge-type and network-smoothed models) and colon cetuximab/BRAF (L1 and elastic net on the shortlist). Random forest has the best mean difference (-0.043) and is rarely clearly worse. Shortlist versions beat all-pathway versions in the small organoid tasks (bladder, BRAF) and lose in lung. The graph network and the MLP are the weakest, as in breast.

## Not done
- HER2 randomised trials (NOAH, CHER-LOB) from 18_other_cohorts: not added. Their saved result points the opposite way (NOAH score AUC 0.328 in the trastuzumab arm vs 0.589 in the chemo-alone arm, interaction p 0.041; CHER-LOB lapatinib arm p 0.379), so they are not positives, and a binary zoo trained on all breast cohorts per arm would need a new run.
- The network shortlists use the paper's proximity table (bundle values), not a new STRING v12 run. Grids are small (2-9 settings) to stay within 2 CPU processes. The GCN uses 150 steps.
- Lung, joint ovarian and olaparib models use regression on a surrogate target, so they are not the identical loss as the originals (Cox ridge, joint Cox/logistic). The original is always shown next to them.
- Small tasks (9 organoids, 6 lines) give very wide intervals; the olaparib and ERCC2/BLCA comparisons cannot separate models.

## Files
- Tables: `tables/positives_all_models.tsv` (task x model: n, events, metric, CI, p, HR, paired difference), `positives_per_cohort.tsv`, `positives_chosen_hyperparameters.tsv`, `positives_model_comparison.tsv`; held-out scores: `tables/preds_pos/*.pkl`.
- Figures (PNG 300 dpi + PDF + caption): `figures/fig_positives_heatmap`, `figures/fig_bar_<task>` (14 tasks).
- Code: `code/pzoo.py`, `pos_stats.py`, `pos_organoid.py`, `pos_lung.py`, `pos_ovarian.py`, `pos_analyze.py`, `pos_figures.py`. Logs: `logs/`.
- Input data: bundles in `nar/07_B3_models/_data/raw/bundles/`, `robust/results/training_table_v2/` (lung, olaparib patients), `training_table/` (rectal), TCGA mutation files used in 10_completeness and 02_organoid_only, `nar/03_organoid_plus_patient/cache/*scores.tsv.gz`.
