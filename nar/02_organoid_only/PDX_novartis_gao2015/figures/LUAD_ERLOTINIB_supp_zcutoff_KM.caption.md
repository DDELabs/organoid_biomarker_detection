# Cut-off sensitivity, Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD erlotinib

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker hs gag degradation; n = 21, deaths = 7, log-rank p = 0.0782, HR = 0.17 (0.02-1.56).
- z <= -1.1824: biomarker metabolism of porphyrins; n = 21, deaths = 7, log-rank p = 0.207, HR = 0.26 (0.03-2.41).
- z <= -1.2816: biomarker egfr downregulation; n = 21, deaths = 7, log-rank p = 0.0396, HR = 7.61 (0.81-71.03).
- z <= -1.645: biomarker egfr downregulation; n = 21, deaths = 7, log-rank p = 0.0396, HR = 7.61 (0.81-71.03).
- z <= -2.1925: biomarker nuclear receptor transcription pathway; n = 21, deaths = 7, log-rank p = 0.771, HR = 1.29 (0.23-7.28).
