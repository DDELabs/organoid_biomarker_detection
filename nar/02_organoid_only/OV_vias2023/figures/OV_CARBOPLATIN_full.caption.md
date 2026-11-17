# Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, carboplatin

All numbers below come from `tables/OV_CARBOPLATIN_*.tsv` and `OV_CARBOPLATIN_summary.json`.
Lab models: 14 organoids, response to oxaliplatin measured as AUC of 5-parameter viability curve
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). No carboplatin/cisplatin screen in these organoids: oxaliplatin response (the only platinum screened) and the paper's cisplatin proximity are used as the platinum proxy. 14 organoids from 11 patients (3 patients give 2 organoids each).

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: paper table (MOESM3 Fig 2a,d). Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). 47 pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `metabolism of polyamines`. Same top pathway in all three models: True.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
0.55 (p = 0.0408), SVR 0.53 (p = 0.0537),
LinearRegression 0.54 (p = 0.047).
**e.** Kaplan-Meier curves for TCGA patients given carboplatin, split at the median of the biomarker score.
n = 347, deaths = 218, log-rank p = 0.543, hazard ratio (HR, non-responders vs
responders) = 1.09 (95% CI 0.83-1.42). 5-year survival: predicted
responders 0.35, predicted non-responders 0.31.
**f.** The same split in patients with no drug record (comparison group): n = 37, deaths =
29, p = 0.0383, HR = 2.19 (1.02-4.69).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the 27
treated tests with k = 2..10: 0.245 (LinearRegression k=2); 0 of them have p < 0.05.
Ridge k = 2: p = 0.53. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: 2000 random sets of 47 pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "survival p<0.05"; the proximal set reaches it at rank
16; empirical p (random sets at least as good) = 0.771. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = 0.532.
**j.** Permutation null: the lab-model responses were shuffled 1000 times and the whole rank-1 Ridge procedure rerun.
Empirical p = 0.628; 11.8% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = 0.252; -1.1824 (user's LUAD choice): p = 0.543.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
