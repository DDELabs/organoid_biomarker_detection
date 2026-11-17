# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD docetaxel

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker endogenous sterols; n = 22, deaths = 13, log-rank p = 0.341, HR = 1.80 (0.53-6.09).
- z <= -1.1824: biomarker endogenous sterols; n = 22, deaths = 13, log-rank p = 0.341, HR = 1.80 (0.53-6.09).
- z <= -1.2816: biomarker endogenous sterols; n = 22, deaths = 13, log-rank p = 0.341, HR = 1.80 (0.53-6.09).
- z <= -1.645: biomarker signaling by hippo; n = 22, deaths = 13, log-rank p = 0.949, HR = 0.96 (0.31-3.01).
- z <= -2.1925: biomarker apoptotic cleavage of cell adhesion proteins; n = 22, deaths = 13, log-rank p = 0.398, HR = 1.64 (0.52-5.20).
