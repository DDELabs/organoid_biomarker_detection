# Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, paclitaxel

All numbers below come from `tables/BRCA_PACLITAXEL_*.tsv` and `BRCA_PACLITAXEL_summary.json`.
Lab models: 37 PDX models, response to paclitaxel measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (37 BRCA models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 150 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `the role of nef in hiv1 replication and disease pathogenesis`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.01 (p = 0.951), SVR 0.03 (p = 0.875),
LinearRegression -0.00 (p = 0.988).
**e.** Kaplan-Meier curves for TCGA patients given paclitaxel, split at the median of the biomarker score.
n = 231, deaths = 22, log-rank p = 0.641, hazard ratio (HR, non-responders vs
responders) = 1.23 (95% CI 0.52-2.92). 5-year survival: predicted
responders 0.91, predicted non-responders 0.90.
**f.** The same split in patients with no drug record (comparison group): n = 298, deaths =
91, p = 0.68, HR = 1.09 (0.72-1.66).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0371 (SVR k=9); 1 of them have p < 0.05.
Ridge k = 2: p = 0.194. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 150 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
19; empirical p (random sets at least as good) = 0.8. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.815.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.643; 2.1% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.316; -1.1824 (user's LUAD choice): p = 0.641.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
