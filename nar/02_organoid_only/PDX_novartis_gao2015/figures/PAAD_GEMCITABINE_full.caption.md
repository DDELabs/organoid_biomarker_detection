# Novartis PDX (PDAC, Gao 2015) -> TCGA-PAAD, gemcitabine

All numbers below come from `tables/PAAD_GEMCITABINE_*.tsv` and `PAAD_GEMCITABINE_summary.json`.
Lab models: 33 PDX models, response to gemcitabine measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (33 PDAC models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 45 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `hyaluronan uptake and degradation`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.09 (p = 0.632), SVR 0.06 (p = 0.745),
LinearRegression 0.11 (p = 0.547).
**e.** Kaplan-Meier curves for TCGA patients given gemcitabine, split at the median of the biomarker score.
n = 104, deaths = 51, log-rank p = 0.137, hazard ratio (HR, non-responders vs
responders) = 1.54 (95% CI 0.87-2.74). 5-year survival: predicted
responders 0.44, predicted non-responders 0.00.
**f.** The same split in patients with no drug record (comparison group): n = 60, deaths =
35, p = 0.0917, HR = 1.79 (0.90-3.54).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0241 (LinearRegression k=4); 2 of them have p < 0.05.
Ridge k = 2: p = 0.105. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 45 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
5; empirical p (random sets at least as good) = 0.306. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.262.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.591; 30.6% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.151; -1.1824 (user's LUAD choice): p = 0.0556.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
