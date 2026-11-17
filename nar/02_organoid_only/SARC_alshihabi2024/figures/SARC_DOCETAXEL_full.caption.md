# Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC, docetaxel

All numbers below come from `tables/SARC_DOCETAXEL_*.tsv` and `SARC_DOCETAXEL_summary.json`.
Lab models: 14 organoids, response to docetaxel measured as % viability at a fixed dose
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). Organoids span several sarcoma subtypes; TCGA-SARC is also mixed. Doxorubicin was not screened in these organoids.

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 150 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `signaling by fgfr1 mutants`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.34 (p = 0.239), SVR 0.27 (p = 0.345),
LinearRegression 0.34 (p = 0.239).
**e.** Kaplan-Meier curves for TCGA patients given docetaxel, split at the median of the biomarker score.
n = 37, deaths = 17, log-rank p = 0.751, hazard ratio (HR, non-responders vs
responders) = 1.17 (95% CI 0.45-3.04). 5-year survival: predicted
responders 0.56, predicted non-responders 0.53.
**f.** The same split in patients with no drug record (comparison group): n = 184, deaths =
66, p = 0.819, HR = 0.95 (0.58-1.53).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.094 (Ridge k=3); 0 of them have p < 0.05.
Ridge k = 2: p = 0.679. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 150 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
29; empirical p (random sets at least as good) = 0.916. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.98.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.83; 6.8% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.726; -1.1824 (user's LUAD choice): p = 0.751.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
