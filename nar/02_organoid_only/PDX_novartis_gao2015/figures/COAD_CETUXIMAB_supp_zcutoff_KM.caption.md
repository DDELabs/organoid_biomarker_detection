# Cut-off sensitivity, Novartis PDX (CRC, Gao 2015) -> TCGA-COAD cetuximab

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker signaling by notch3; n = 14, deaths = 6, log-rank p = 0.963, HR = 0.96 (0.19-4.92).
- z <= -1.1824: biomarker signaling by constitutively active egfr; n = 14, deaths = 6, log-rank p = 0.963, HR = 0.96 (0.19-4.92).
- z <= -1.2816: biomarker abc family proteins mediated transport; n = 14, deaths = 6, log-rank p = 0.927, HR = 0.93 (0.18-4.67).
- z <= -1.645: biomarker signaling by constitutively active egfr; n = 14, deaths = 6, log-rank p = 0.963, HR = 0.96 (0.19-4.92).
- z <= -2.1925: biomarker l1cam interactions; n = 14, deaths = 6, log-rank p = 0.0276, HR = 0.00 (0.00-inf).
