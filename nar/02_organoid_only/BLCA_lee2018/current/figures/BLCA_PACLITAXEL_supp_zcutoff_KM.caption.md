# Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA paclitaxel

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker endogenous sterols; n = 10, deaths = 6, log-rank p = 0.359, HR = 0.46 (0.08-2.53).
- z <= -1.1824: biomarker endogenous sterols; n = 10, deaths = 6, log-rank p = 0.359, HR = 0.46 (0.08-2.53).
- z <= -1.2816: biomarker endogenous sterols; n = 10, deaths = 6, log-rank p = 0.359, HR = 0.46 (0.08-2.53).
- z <= -1.645: biomarker activation of bh3 only proteins; n = 10, deaths = 6, log-rank p = 0.211, HR = 0.35 (0.06-1.95).
- z <= -2.1925: biomarker activation of bh3 only proteins; n = 10, deaths = 6, log-rank p = 0.211, HR = 0.35 (0.06-1.95).
