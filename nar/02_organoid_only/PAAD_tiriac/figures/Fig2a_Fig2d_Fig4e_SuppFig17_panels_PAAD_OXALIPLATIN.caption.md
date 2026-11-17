# Supplementary figure: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | OXALIPLATIN

All panels follow the organoid-only strategy of Kong et al. 2020 (network proximity -> pathway scores ->
machine learning on lab models -> median split of drug-treated patients -> survival test).

- **a** Network proximity. Each bar is one Reactome pathway; the bar length is how close the pathway sits to the
  drug's protein targets in the STRING network (z-score; more negative = closer). Blue bars pass the paper's
  cut-off (z <= -1.2816); 47 of 668 pathways pass. Proximity values are the
  authors' published table for cisplatin (oxaliplatin has no protein target, so the platinum drug cisplatin is used as its proxy).
- **b** Same analysis with other cut-offs (the user's -1.182405 and -2.1925376, and -1.6449 = 5%). Bars show
  -log10 p of the rank-1 pathway in treated (blue) and untreated (grey) patients; the dotted line is p = 0.05.
- **c** Pathway activity (ssGSEA scores) of the proximal pathways in the 36 lab models, ordered by drug
  response (AUC; lower = more sensitive). Red = higher than average for that pathway, blue = lower.
- **d-f** Leave-one-out test inside the lab models: each model is predicted by a model trained on all the others
  (Ridge, SVR, LinearRegression, all proximal pathways). Spearman rho measures agreement (1 = perfect order).
- **g** p-value (as -log10 p; higher = smaller p) for 1 to 10 top pathways for each learner, in treated (solid)
  and untreated (dashed) patients.
- **h** Kaplan-Meier curves for drug-treated patients split at the median predicted response using the rank-1
  Ridge pathway (Synthesis of substrates in n glycan biosythesis). Log-rank p = 0.223; n = 21; deaths = 11
  (8 in predicted responders, 3 in predicted non-responders);
  HR 0.44 (95% CI 0.11-1.71) (non-responders vs responders; > 1 is the expected direction).
- **i** The same pathway and split rule in patients of the same cancer with no drug record ("untreated"), to see
  whether the split is drug-specific. Log-rank p = 0.072; n = 61; deaths = 35; HR 0.53 (95% CI 0.26-1.07).
- **j** The paper's bootstrap (Fig. 4e/f): 10,000 random pathway sets of the same size, each run through
  Ridge -> rank-1 -> log-rank. Share of random sets with p as small or smaller than the real one = 0.160.
- **k** The lab-model drug response shuffled 1,000 times and the whole rank-1 analysis repeated.
  Share of shuffles as good or better = 0.135.
- **l** The same shuffles, but keeping the best of all 30 tests (3 learners x 1-10 pathways), which mirrors picking
  the best result. Real best p = 0.201; share of shuffles as good or better = 0.645.
- **m** The paper's own bootstrap statistic (Fig. 4e): for each random set, the first rank (by |coefficient|) whose
  pathway alone gives p < 0.05 in treated patients.
  Real proximal pathways: first such rank = 41;
  share of random sets reaching that rank or better = 0.726 (paper, colon 5-FU: 0.0012).
- **n** Key numbers. The rank-1 pathway's own proximity z = -1.6744 (cut-off -1.2816).
- Panels d-f stand in for the paper's Supp. Fig. 17 (a 60/10/30 split of the organoids); leave-one-out uses every
  model once as test data, which is more stable at these sample sizes.

Data: current GDC STAR FPKM-UQ + TCGA-CDR OS + GDC biotab drug file. Lab-model response metric: AUC.
