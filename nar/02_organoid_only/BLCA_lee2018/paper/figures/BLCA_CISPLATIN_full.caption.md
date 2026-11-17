# Bladder organoids (Lee 2018), paper source-data NES -> TCGA-BLCA, cisplatin

All numbers below come from `tables/BLCA_CISPLATIN_*.tsv` and `BLCA_CISPLATIN_summary.json`.
Lab models: 9 organoids, response to cisplatin measured as AUC of fitted viability curve (CoderData)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). 

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 30 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `transport of glucose and other sugars bile salts and organic acids metal ions a…`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
-0.07 (p = 0.865), SVR -0.07 (p = 0.865),
LinearRegression -0.07 (p = 0.865).
**e.** Kaplan-Meier curves for TCGA patients given cisplatin, split at the median of the biomarker score.
n = 79, deaths = 20, log-rank p = 0.0435, hazard ratio (HR, non-responders vs
responders) = 2.63 (95% CI 0.99-6.96). 5-year survival: predicted
responders 0.75, predicted non-responders 0.44.
**f.** The same split in patients with no drug record (comparison group): n = 288, deaths =
79, p = 0.929, HR = 1.02 (0.66-1.59).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0171 (Ridge k=10); 2 of them have p < 0.05.
Ridge k = 2: p = 0.0604. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 10000 random sets of 30 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05 and isogenic t-test p<0.05 in the predicted direction"; the proximal set reaches it at rank
6; empirical p (random sets at least as good) = 0.0437. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.0663.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.297; 29.6% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.909; -1.1824 (user's LUAD choice): p = 0.0435.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
