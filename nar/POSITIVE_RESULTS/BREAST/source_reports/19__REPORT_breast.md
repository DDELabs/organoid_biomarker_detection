# W3 breast: all-model benchmark (leave-one-cohort-out)

## What was done (plain words)
- Task: predict chemotherapy response (responder = 1) in a breast cohort the model has never seen. 13 cohorts, 3,511 treated patients, 673 pathway scores (z-scored inside each cohort).
- Step 1, check: I re-ran the existing L2 logistic model (penalty lambda = 100, no tuning). All 13 AUCs match `robust/results/train_combined/part1_response_leave_one_cohort_out.tsv` to 3 decimals (largest difference 0.0000). File: `tables/l2_reproduction_check.tsv`. This model is called **L2_fixed** and is the reference for every comparison.
- Step 2: every other model was trained on the 12 other cohorts. Settings (penalty, trees, depth, and so on) were picked by 4-fold cross-validation that keeps whole cohorts together, **inside the 12 training cohorts only** (score = mean AUC per cohort). The held-out cohort was never used for tuning. Chosen settings: `tables/preds/<model>_hp.tsv`. Seeds fixed (20261118).
- Models: L2 (tuned), L1, elastic net, linear SVM, RBF SVM, random forest, XGBoost, MLP (1 hidden layer of 32), network-smoothed L2 (each pathway score mixed with its 10 most gene-overlapping pathways), and a graph network (GCN from `nar/07_B3_models/code/nets.py`, nodes = pathways, edges = gene overlap, node prior = network proximity to the drug targets).
- Network-proximity shortlist: 183 pathways at z <= -1.2816 from the paper's own proximity table, for the breast drugs it contains (paclitaxel, docetaxel, doxorubicin, 5-FU; union). I could not use the STRING v12 recomputation because its cache is not in the repo. The shortlist uses no outcome data. Four linear models were run on the shortlist vs on all 673 pathways.
- Baselines: ER+HER2 only (logistic on the two receptor flags; 9 cohorts have them), proliferation score (mean of the 27 cell-cycle pathways; sign set from training cohorts), 200 random sets of 183 pathways (L2 each time).
- Statistics: per-cohort AUC with 95% CI and p versus 0.5 (DeLong). Paired comparison with L2_fixed on the same patients (DeLong test of correlated AUCs). Mean over cohorts = plain average; its CI and the paired mean difference use the sum of per-cohort variances (cohorts are independent). q = Benjamini-Hochberg across models.

## Main result
| model | mean AUC (95% CI) | mean difference vs L2_fixed (95% CI), p, q | cohorts better / worse | cohorts significantly (p<0.05) better / worse | mean AUC within ER/HER2 groups |
|---|---|---|---|---|---|
| SVM_linear | 0.764 (0.741-0.786) | +0.012 (+0.002 to +0.022), p=0.023, q=0.100 | 11/2 | 2/1 | 0.718 |
| NetSmoothL2 | 0.759 (0.737-0.782) | +0.008 (-0.001 to +0.016), p=0.074, q=0.141 | 10/3 | 2/1 | 0.713 |
| L2 | 0.759 (0.737-0.782) | +0.008 (-0.000 to +0.016), p=0.060, q=0.141 | 11/2 | 2/1 | 0.715 |
| SVM_RBF | 0.756 (0.732-0.779) | +0.004 (-0.010 to +0.019), p=0.577, q=0.613 | 6/7 | 1/0 | 0.708 |
| L1 | 0.755 (0.732-0.778) | +0.004 (-0.004 to +0.011), p=0.348, q=0.394 | 9/4 | 0/0 | 0.708 |
| L2_fixed | 0.752 (0.728-0.775) |  | 0/0 | 0/0 | 0.705 |
| ElasticNet | 0.751 (0.728-0.774) | -0.000 (-0.010 to +0.009), p=0.943, q=0.943 | 8/4 | 0/1 | 0.705 |
| L2_shortlist | 0.744 (0.721-0.767) | -0.008 (-0.022 to +0.007), p=0.297, q=0.389 | 7/6 | 0/2 | 0.699 |
| XGBoost | 0.743 (0.719-0.766) | -0.009 (-0.026 to +0.009), p=0.327, q=0.394 | 6/7 | 1/1 | 0.689 |
| L2_fixed_shortlist | 0.741 (0.718-0.765) | -0.010 (-0.024 to +0.004), p=0.152, q=0.236 | 5/8 | 0/1 | 0.694 |
| L1_shortlist | 0.741 (0.717-0.764) | -0.011 (-0.026 to +0.004), p=0.150, q=0.236 | 6/7 | 0/1 | 0.696 |
| RandomForest | 0.740 (0.716-0.763) | -0.012 (-0.031 to +0.007), p=0.215, q=0.305 | 5/8 | 0/1 | 0.683 |
| ElasticNet_shortlist | 0.738 (0.714-0.761) | -0.014 (-0.029 to +0.001), p=0.070, q=0.141 | 6/7 | 0/2 | 0.692 |
| SVM_linear_shortlist | 0.735 (0.711-0.759) | -0.016 (-0.031 to -0.002), p=0.031, q=0.105 | 6/7 | 0/2 | 0.693 |
| MLP | 0.734 (0.709-0.759) | -0.018 (-0.035 to -0.001), p=0.039, q=0.109 | 4/9 | 0/2 | 0.689 |
| ER_HER2_only | 0.693 (0.666-0.720) | -0.065 (-0.095 to -0.035), p=0.000, q=0.000 | 0/9 | 0/2 | 0.500 |
| GCN_prior | 0.677 (0.651-0.703) | -0.075 (-0.101 to -0.049), p=0.000, q=0.000 | 3/10 | 0/6 | 0.620 |
| Proliferation_score | 0.646 (0.619-0.673) | -0.106 (-0.134 to -0.078), p=0.000, q=0.000 | 0/13 | 0/7 | 0.613 |

What I think: the existing L2 is already near the top. Linear SVM has the best mean AUC (0.764 vs 0.752), a gain of 0.012. Its p is 0.023 but q is 0.100 after correcting for the 17 models compared, and the sign test across cohorts gives p=0.094. So it is a small gain that is borderline, not clear-cut. Tuned L2 and network-smoothed L2 are in the same range (+0.008, p about 0.06 to 0.07). Notice that the best linear SVM and L2 settings sat at the strongest-shrinkage edge of the grid (C=1e-4 for SVM; 0.001-0.003 for L2), so stronger shrinkage than lambda=100 helps a little. Random forest, XGBoost, RBF SVM and MLP did not beat the linear models (RF -0.012, XGBoost -0.009, MLP -0.018; MLP significant at p=0.039 but not after correction). The graph network was clearly worse (0.677); with only 150 training steps and a pathway-overlap graph it may be under-trained, so I would not read this as a verdict on graph models in general.

## Network-proximity shortlist vs all pathways
Using only the 183 shortlisted pathways lowers the mean AUC for every linear model (L2: 0.741 vs 0.752). The 200 random 183-pathway sets average 0.742 (range of set averages 0.723-0.758), and about 31% of random sets are as good as or better than the shortlist (0.744). What I think: the shortlist is no better than a random set of the same size here. Table `tables/random_pathway_sets_auc.tsv`.

## Baselines
ER+HER2 only: mean 0.693 (9 cohorts). Proliferation score alone: 0.646 (13 cohorts). Both are clearly below L2 (paired p<0.001), so the pathway model carries information beyond receptors and a simple proliferation score. AUC within ER/HER2 groups (receptor group differences removed): L2_fixed 0.705, SVM_linear 0.718. What I think: the proliferation score is a strong-but-incomplete part of the signal; W4 will dig into this.

## Best model, per cohort (linear SVM vs existing L2)
| cohort | n | responders | L2 AUC | SVM AUC (95% CI) | p vs 0.5 | difference | paired p |
|---|---|---|---|---|---|---|---|
| BRIGHTNESS | 482 | 236 | 0.674 | 0.683 (0.635-0.730) | 3.0e-14 | +0.008 | 0.502 |
| GSE123845 | 112 | 42 | 0.761 | 0.781 (0.694-0.869) | 3.1e-10 | +0.020 | 0.452 |
| GSE163882 | 222 | 80 | 0.715 | 0.762 (0.698-0.826) | 7.9e-16 | +0.047 | 0.008 |
| GSE20194 | 270 | 52 | 0.798 | 0.812 (0.752-0.872) | 2.2e-24 | +0.014 | 0.376 |
| GSE20271 | 178 | 26 | 0.801 | 0.808 (0.732-0.885) | 2.5e-15 | +0.008 | 0.699 |
| GSE22093 | 97 | 28 | 0.751 | 0.769 (0.668-0.870) | 1.8e-07 | +0.019 | 0.393 |
| GSE23988 | 61 | 20 | 0.846 | 0.857 (0.756-0.958) | 4.2e-12 | +0.011 | 0.636 |
| GSE240671 | 95 | 37 | 0.780 | 0.746 (0.641-0.851) | 4.2e-06 | -0.034 | 0.041 |
| GSE25066 | 438 | 82 | 0.817 | 0.794 (0.741-0.846) | 6.6e-28 | -0.023 | 0.051 |
| GSE32646 | 115 | 27 | 0.741 | 0.785 (0.687-0.882) | 9.9e-09 | +0.044 | 0.053 |
| GSE41998 | 253 | 69 | 0.715 | 0.729 (0.660-0.798) | 9.0e-11 | +0.015 | 0.226 |
| ISPY2 | 986 | 319 | 0.700 | 0.718 (0.686-0.750) | 1.9e-40 | +0.019 | 0.020 |
| TCGA-BRCA | 202 | 184 | 0.672 | 0.681 (0.568-0.794) | 1.7e-03 | +0.009 | 0.750 |

## Files
- Tables: `tables/breast_loco_all_models.tsv` (every model x cohort: n, responders, AUC, CI, p, AUC within ER/HER2, paired difference vs L2 with CI and p), `tables/model_comparison.tsv`, `tables/random_pathway_sets_auc.tsv`, `tables/pathway_flags.tsv`, `tables/proliferation_pathways_used.tsv`, test-set predictions in `tables/preds/`.
- Figures (PNG 300 dpi + PDF + `.caption.md`): `figures/fig_model_heatmap`, `figures/fig_model_mean_auc`, `figures/fig_roc_top3` (top 3 whole-pathway models: SVM_linear, NetSmoothL2, L2).
- Code: `code/` (common.py, prep.py, repro_l2.py, bench_breast.py, baselines.py, analyze.py, figures.py). Logs: `logs/`.
- Input data: `robust/results/training_table/` (patients.tsv.gz, pathway_scores_*.tsv.gz; 20,766 samples), duplicates dropped as in the original script.

## Limits
- Not done in this run: the other W3 settings (positives: ERCC2, cetuximab, rectal, lung); breast only, as asked. STRING v12 shortlist not used (see above).
- Hyper-parameter grids are small (2 to 9 settings) to keep CPU use to 2 processes. The GCN used 150 steps only.
- Breast cohorts mix drugs, so one shared shortlist is used for all cohorts rather than a drug-specific one.
