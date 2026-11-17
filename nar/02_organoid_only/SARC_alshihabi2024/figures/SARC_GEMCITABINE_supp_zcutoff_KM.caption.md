# Cut-off sensitivity, Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC gemcitabine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker early phase of hiv life cycle; n = 42, deaths = 20, log-rank p = 0.132, HR = 0.49 (0.19-1.26).
- z <= -1.1824: biomarker hyaluronan uptake and degradation; n = 42, deaths = 20, log-rank p = 0.0835, HR = 0.45 (0.18-1.14).
- z <= -1.2816: biomarker hyaluronan uptake and degradation; n = 42, deaths = 20, log-rank p = 0.0835, HR = 0.45 (0.18-1.14).
- z <= -1.645: biomarker early phase of hiv life cycle; n = 42, deaths = 20, log-rank p = 0.132, HR = 0.49 (0.19-1.26).
- z <= -2.1925: biomarker pyrimidine catabolism; n = 42, deaths = 20, log-rank p = 0.671, HR = 0.82 (0.33-2.03).
