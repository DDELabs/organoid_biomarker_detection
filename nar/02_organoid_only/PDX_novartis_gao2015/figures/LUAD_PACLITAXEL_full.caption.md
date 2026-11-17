# Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD, paclitaxel

All numbers below come from `tables/LUAD_PACLITAXEL_*.tsv` and `LUAD_PACLITAXEL_summary.json`.
Lab models: 24 PDX models, response to paclitaxel measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (24 NSCLC models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 150 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `ionotropic activity of kainate receptors`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.18 (p = 0.393), SVR 0.27 (p = 0.197),
LinearRegression 0.19 (p = 0.377).
**e.** Kaplan-Meier curves for TCGA patients given paclitaxel, split at the median of the biomarker score.
n = 44, deaths = 21, log-rank p = 0.913, hazard ratio (HR, non-responders vs
responders) = 1.05 (95% CI 0.42-2.63). 5-year survival: predicted
responders 0.31, predicted non-responders 0.31.
**f.** The same split in patients with no drug record (comparison group): n = 327, deaths =
114, p = 0.00117, HR = 1.86 (1.27-2.73).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0815 (Ridge k=4); 0 of them have p < 0.05.
Ridge k = 2: p = 0.906. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 150 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
17; empirical p (random sets at least as good) = 0.524. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.681.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.952; 9.3% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.773; -1.1824 (user's LUAD choice): p = 0.913.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
