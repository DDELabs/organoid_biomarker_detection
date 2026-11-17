# Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, doxorubicin

All numbers below come from `tables/OV_DOXORUBICIN_*.tsv` and `OV_DOXORUBICIN_summary.json`.
Lab models: 14 organoids, response to doxorubicin measured as AUC of 5-parameter viability curve
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). 14 organoids from 11 patients (3 patients give 2 organoids each).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 14 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `resolution of ap sites via the multiple nucleotide patch replacement pathway`. Same top pathway in all three models: False.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.02 (p = 0.935), SVR 0.38 (p = 0.185),
LinearRegression -0.26 (p = 0.375).
**e.** Kaplan-Meier curves for TCGA patients given doxorubicin, split at the median of the biomarker score.
n = 117, deaths = 98, log-rank p = 0.976, hazard ratio (HR, non-responders vs
responders) = 0.99 (95% CI 0.66-1.50). 5-year survival: predicted
responders 0.28, predicted non-responders 0.18.
**f.** The same split in patients with no drug record (comparison group): n = 37, deaths =
29, p = 0.497, HR = 0.77 (0.37-1.62).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.0473 (Ridge k=5); 1 of them have p < 0.05.
Ridge k = 2: p = 0.883. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 14 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
None; empirical p (random sets at least as good) = nan. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.983.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.895; 0.0% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.986; -1.1824 (user's LUAD choice): p = 0.604.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
