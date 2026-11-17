# Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, carboplatin

All numbers below come from `tables/LUAD_CARBOPLATIN_*.tsv` and `LUAD_CARBOPLATIN_summary.json`.
Lab models: 43 cell lines, response to carboplatin measured as AUC (CTRPv2)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). Carboplatin response from CTRPv2 on the same GDSC/CCLE LUAD lines; cisplatin proximity used (carboplatin has no protein targets in DrugBank).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 47 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `reversible hydration of carbon dioxide`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.00 (p = 0.985), SVR -0.07 (p = 0.649),
LinearRegression -0.17 (p = 0.274).
**e.** Kaplan-Meier curves for TCGA patients given carboplatin, split at the median of the biomarker score.
n = 77, deaths = 33, log-rank p = 0.579, hazard ratio (HR, non-responders vs
responders) = 0.82 (95% CI 0.41-1.65). 5-year survival: predicted
responders 0.20, predicted non-responders 0.29.
**f.** The same split in patients with no drug record (comparison group): n = 327, deaths =
114, p = 0.54, HR = 0.89 (0.62-1.29).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0666 (SVR k=3); 0 of them have p < 0.05.
Ridge k = 2: p = 0.379. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 47 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
None; empirical p (random sets at least as good) = nan. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.557.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.514; 0.0% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.761; -1.1824 (user's LUAD choice): p = 0.761.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
