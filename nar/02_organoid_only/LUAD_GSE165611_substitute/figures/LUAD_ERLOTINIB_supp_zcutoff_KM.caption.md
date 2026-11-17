# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD erlotinib

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker phase1 functionalization of compounds; n = 21, deaths = 7, log-rank p = 0.353, HR = 0.46 (0.08-2.50).
- z <= -1.1824: biomarker phase1 functionalization of compounds; n = 21, deaths = 7, log-rank p = 0.353, HR = 0.46 (0.08-2.50).
- z <= -1.2816: biomarker pi3k events in erbb2 signaling; n = 21, deaths = 7, log-rank p = 0.514, HR = 0.57 (0.10-3.17).
- z <= -1.645: biomarker pi3k events in erbb2 signaling; n = 21, deaths = 7, log-rank p = 0.514, HR = 0.57 (0.10-3.17).
- z <= -2.1925: biomarker recycling of bile acids and salts; n = 21, deaths = 7, log-rank p = 0.871, HR = 1.16 (0.19-7.24).
