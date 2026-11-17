# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD vinorelbine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker ras activation uopn ca2 infux through nmda receptor; n = 24, deaths = 13, log-rank p = 0.626, HR = 1.32 (0.43-4.13).
- z <= -1.1824: biomarker signaling by hippo; n = 24, deaths = 13, log-rank p = 0.897, HR = 1.08 (0.36-3.26).
- z <= -1.2816: biomarker cell cycle mitotic; n = 24, deaths = 13, log-rank p = 0.337, HR = 1.74 (0.56-5.42).
- z <= -1.645: biomarker cell cycle; n = 24, deaths = 13, log-rank p = 0.337, HR = 0.58 (0.18-1.80).
- z <= -2.1925: biomarker cell cycle mitotic; n = 24, deaths = 13, log-rank p = 0.337, HR = 1.74 (0.56-5.42).
