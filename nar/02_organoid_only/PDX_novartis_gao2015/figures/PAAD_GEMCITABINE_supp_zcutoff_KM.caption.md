# Cut-off sensitivity, Novartis PDX (PDAC, Gao 2015) -> TCGA-PAAD gemcitabine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker hyaluronan uptake and degradation; n = 104, deaths = 51, log-rank p = 0.137, HR = 1.54 (0.87-2.74).
- z <= -1.1824: biomarker cdc6 association with the orc origin complex; n = 104, deaths = 51, log-rank p = 0.0556, HR = 1.73 (0.98-3.06).
- z <= -1.2816: biomarker hyaluronan uptake and degradation; n = 104, deaths = 51, log-rank p = 0.137, HR = 1.54 (0.87-2.74).
- z <= -1.645: biomarker cdc6 association with the orc origin complex; n = 104, deaths = 51, log-rank p = 0.0556, HR = 1.73 (0.98-3.06).
- z <= -2.1925: biomarker mitotic g1 g1 s phases; n = 104, deaths = 51, log-rank p = 0.151, HR = 0.66 (0.38-1.17).
