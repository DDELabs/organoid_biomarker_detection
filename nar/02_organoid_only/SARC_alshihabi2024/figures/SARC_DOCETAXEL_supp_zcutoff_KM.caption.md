# Cut-off sensitivity, Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC docetaxel

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker signaling by fgfr1 mutants; n = 37, deaths = 17, log-rank p = 0.751, HR = 1.17 (0.45-3.04).
- z <= -1.1824: biomarker signaling by fgfr1 mutants; n = 37, deaths = 17, log-rank p = 0.751, HR = 1.17 (0.45-3.04).
- z <= -1.2816: biomarker signaling by fgfr1 mutants; n = 37, deaths = 17, log-rank p = 0.751, HR = 1.17 (0.45-3.04).
- z <= -1.645: biomarker p38mapk events; n = 37, deaths = 17, log-rank p = 0.975, HR = 0.99 (0.38-2.57).
- z <= -2.1925: biomarker recycling of bile acids and salts; n = 37, deaths = 17, log-rank p = 0.726, HR = 0.84 (0.31-2.27).
