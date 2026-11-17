# Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, trastuzumab

All numbers below come from `tables/BRCA_TRASTUZUMAB_*.tsv` and `BRCA_TRASTUZUMAB_summary.json`.
Lab models: 36 PDX models, response to trastuzumab measured as BestAvgResponse (% tumour volume change)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). One mouse per PDX model (36 BRCA models).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: computed on paper network, targets ERBB2. Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 54 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `platelet adhesion to exposed collagen`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.11 (p = 0.51), SVR 0.13 (p = 0.437),
LinearRegression 0.09 (p = 0.615).
**e.** Kaplan-Meier curves for TCGA patients given trastuzumab, split at the median of the biomarker score.
n = 74, deaths = 4, log-rank p = 0.0916, hazard ratio (HR, non-responders vs
responders) = 30491286581.07 (95% CI 0.00-inf). 5-year survival: predicted
responders 1.00, predicted non-responders 0.88.
**f.** The same split in patients with no drug record (comparison group): n = 298, deaths =
91, p = 0.622, HR = 1.11 (0.73-1.68).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.088 (LinearRegression k=3); 0 of them have p < 0.05.
Ridge k = 2: p = 0.624. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 54 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
None; empirical p (random sets at least as good) = nan. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.129.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.131; 0.0% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.456; -1.1824 (user's LUAD choice): p = 0.51.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
