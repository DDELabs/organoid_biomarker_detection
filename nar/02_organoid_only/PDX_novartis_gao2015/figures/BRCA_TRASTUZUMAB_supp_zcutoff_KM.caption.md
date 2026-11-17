# Cut-off sensitivity, Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA trastuzumab

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker activation of chaperone genes by xbp1s; n = 74, deaths = 4, log-rank p = 0.0916, HR = 30491293586.03 (0.00-inf).
- z <= -1.1824: biomarker insulin synthesis and processing; n = 74, deaths = 4, log-rank p = 0.51, HR = 0.45 (0.04-5.03).
- z <= -1.2816: biomarker platelet adhesion to exposed collagen; n = 74, deaths = 4, log-rank p = 0.0916, HR = 30491286581.07 (0.00-inf).
- z <= -1.645: biomarker basigin interactions; n = 74, deaths = 4, log-rank p = 0.578, HR = 1.95 (0.18-21.55).
- z <= -2.1925: biomarker sema4d induced cell migration and growth cone collapse; n = 74, deaths = 4, log-rank p = 0.456, HR = 2.43 (0.22-26.92).
