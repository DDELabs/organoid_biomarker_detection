# Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, vinorelbine

All numbers below come from `tables/LUAD_VINORELBINE_*.tsv` and `LUAD_VINORELBINE_summary.json`.
Lab models: 46 cell lines, response to vinorelbine measured as ln IC50 (GDSC2)
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). Lab models are 49 GDSC lung adenocarcinoma cell lines (GDSC2), not patient-derived cells.

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: computed on paper network, targets TUBB. Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 24 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `cell cycle mitotic`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.30 (p = 0.0408), SVR 0.35 (p = 0.0175),
LinearRegression 0.15 (p = 0.334).
**e.** Kaplan-Meier curves for TCGA patients given vinorelbine, split at the median of the biomarker score.
n = 24, deaths = 13, log-rank p = 0.337, hazard ratio (HR, non-responders vs
responders) = 1.74 (95% CI 0.56-5.42). 5-year survival: predicted
responders 0.53, predicted non-responders 0.15.
**f.** The same split in patients with no drug record (comparison group): n = 327, deaths =
114, p = 0.0129, HR = 0.62 (0.43-0.91).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.021 (Ridge k=10); 2 of them have p < 0.05.
Ridge k = 2: p = 0.553. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 24 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
7; empirical p (random sets at least as good) = 0.195. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.309.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.255; 5.1% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.337; -1.1824 (user's LUAD choice): p = 0.897.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
