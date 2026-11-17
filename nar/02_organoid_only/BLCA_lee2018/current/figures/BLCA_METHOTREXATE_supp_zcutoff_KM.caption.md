# Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA methotrexate

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker signaling by tgf beta receptor complex; n = 12, deaths = 5, log-rank p = 0.614, HR = 0.63 (0.11-3.81).
- z <= -1.1824: biomarker signaling by tgf beta receptor complex; n = 12, deaths = 5, log-rank p = 0.614, HR = 0.63 (0.11-3.81).
- z <= -1.2816: biomarker signaling by tgf beta receptor complex; n = 12, deaths = 5, log-rank p = 0.614, HR = 0.63 (0.11-3.81).
- z <= -1.645: biomarker enos activation and regulation; n = 12, deaths = 5, log-rank p = 0.601, HR = 1.61 (0.27-9.67).
- z <= -2.1925: biomarker enos activation and regulation; n = 12, deaths = 5, log-rank p = 0.601, HR = 1.61 (0.27-9.67).
