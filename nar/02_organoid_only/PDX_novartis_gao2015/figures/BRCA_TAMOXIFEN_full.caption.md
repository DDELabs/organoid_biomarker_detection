# Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, tamoxifen

All numbers below come from `tables/BRCA_TAMOXIFEN_*.tsv` and `BRCA_TAMOXIFEN_summary.json`.
Lab models: 37 PDX models, response to tamoxifen measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (37 BRCA models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: computed on paper network, targets ESR1,ESR2. Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 150 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `spry regulation of fgf signaling`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.16 (p = 0.343), SVR 0.13 (p = 0.437),
LinearRegression 0.16 (p = 0.357).
**e.** Kaplan-Meier curves for TCGA patients given tamoxifen, split at the median of the biomarker score.
n = 253, deaths = 22, log-rank p = 0.786, hazard ratio (HR, non-responders vs
responders) = 1.13 (95% CI 0.48-2.65). 5-year survival: predicted
responders 0.92, predicted non-responders 0.86.
**f.** The same split in patients with no drug record (comparison group): n = 298, deaths =
91, p = 0.78, HR = 0.94 (0.62-1.43).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.107 (Ridge k=9); 0 of them have p < 0.05.
Ridge k = 2: p = 0.66. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 150 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
24; empirical p (random sets at least as good) = 0.603. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.782.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.866; 2.3% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.644; -1.1824 (user's LUAD choice): p = 0.557.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
