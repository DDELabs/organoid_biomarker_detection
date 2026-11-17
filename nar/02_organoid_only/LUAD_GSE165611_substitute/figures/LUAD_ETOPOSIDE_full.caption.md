# Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, etoposide

All numbers below come from `tables/LUAD_ETOPOSIDE_*.tsv` and `LUAD_ETOPOSIDE_summary.json`.
Lab models: 46 cell lines, response to etoposide measured as ln IC50 (GDSC1)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). Lab models are 48 GDSC lung adenocarcinoma cell lines (GDSC1), not patient-derived cells.

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: computed on paper network, targets TOP2A,TOP2B. Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 40 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `base excision repair`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.16 (p = 0.282), SVR 0.34 (p = 0.02),
LinearRegression -0.18 (p = 0.219).
**e.** Kaplan-Meier curves for TCGA patients given etoposide, split at the median of the biomarker score.
n = 15, deaths = 4, log-rank p = 0.8, hazard ratio (HR, non-responders vs
responders) = 1.30 (95% CI 0.17-9.74). 5-year survival: predicted
responders 0.73, predicted non-responders 0.00.
**f.** The same split in patients with no drug record (comparison group): n = 327, deaths =
114, p = 0.514, HR = 0.88 (0.61-1.28).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.048 (Ridge k=3); 1 of them have p < 0.05.
Ridge k = 2: p = 0.792. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 40 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
10; empirical p (random sets at least as good) = 0.554. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.724.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.619; 5.5% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.8; -1.1824 (user's LUAD choice): p = 0.815.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
