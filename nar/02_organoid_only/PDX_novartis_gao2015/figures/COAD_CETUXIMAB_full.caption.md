# Novartis PDX (CRC, Gao 2015) -> TCGA-COAD, cetuximab

All numbers below come from `tables/COAD_CETUXIMAB_*.tsv` and `COAD_CETUXIMAB_summary.json`.
Lab models: 41 PDX models, response to cetuximab measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (41 CRC models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 41 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `abc family proteins mediated transport`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.10 (p = 0.554), SVR 0.04 (p = 0.782),
LinearRegression -0.11 (p = 0.495).
**e.** Kaplan-Meier curves for TCGA patients given cetuximab, split at the median of the biomarker score.
n = 14, deaths = 6, log-rank p = 0.927, hazard ratio (HR, non-responders vs
responders) = 0.93 (95% CI 0.18-4.67). 5-year survival: predicted
responders 0.27, predicted non-responders 0.38.
**f.** The same split in patients with no drug record (comparison group): n = 282, deaths =
62, p = 0.161, HR = 0.70 (0.42-1.16).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.00238 (LinearRegression k=5); 11 of them have p < 0.05.
Ridge k = 2: p = 0.513. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 41 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
7; empirical p (random sets at least as good) = 0.711. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.914.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.887; 24.8% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.0276; -1.1824 (user's LUAD choice): p = 0.963.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
