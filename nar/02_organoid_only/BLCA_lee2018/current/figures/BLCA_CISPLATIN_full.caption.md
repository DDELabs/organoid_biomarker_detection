# Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, cisplatin

All numbers below come from `tables/BLCA_CISPLATIN_*.tsv` and `BLCA_CISPLATIN_summary.json`.
Lab models: 9 organoids, response to cisplatin measured as AUC of fitted viability curve (CoderData)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). 

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 47 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `transport of organic anions`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.23 (p = 0.546), SVR 0.23 (p = 0.546),
LinearRegression 0.23 (p = 0.546).
**e.** Kaplan-Meier curves for TCGA patients given cisplatin, split at the median of the biomarker score.
n = 79, deaths = 28, log-rank p = 0.0672, hazard ratio (HR, non-responders vs
responders) = 2.03 (95% CI 0.94-4.38). 5-year survival: predicted
responders 0.68, predicted non-responders 0.50.
**f.** The same split in patients with no drug record (comparison group): n = 283, deaths =
131, p = 0.107, HR = 0.75 (0.53-1.06).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.00821 (Ridge k=9); 10 of them have p < 0.05.
Ridge k = 2: p = 0.0707. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 10000 random sets of 47 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05 and isogenic t-test p<0.05 in the predicted direction"; the proximal set reaches it at rank
None; empirical p (random sets at least as good) = nan. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.145.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.036; 3.5% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.0672; -1.1824 (user's LUAD choice): p = 0.0672.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
