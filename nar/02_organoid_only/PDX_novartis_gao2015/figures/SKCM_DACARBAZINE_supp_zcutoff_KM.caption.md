# Cut-off sensitivity, Novartis PDX (CM, Gao 2015) -> TCGA-SKCM dacarbazine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker metabolism of carbohydrates; n = 32, deaths = 21, log-rank p = 0.0011, HR = 0.21 (0.07-0.58).
- z <= -1.1824: biomarker metabolism of carbohydrates; n = 32, deaths = 21, log-rank p = 0.0011, HR = 0.21 (0.07-0.58).
- z <= -1.2816: biomarker metabolism of carbohydrates; n = 32, deaths = 21, log-rank p = 0.0011, HR = 0.21 (0.07-0.58).
- z <= -1.645: biomarker pol switching; n = 32, deaths = 21, log-rank p = 0.0917, HR = 0.43 (0.16-1.17).
- z <= -2.1925: biomarker lagging strand synthesis; n = 32, deaths = 21, log-rank p = 0.21, HR = 1.88 (0.69-5.09).
