# Cut-off sensitivity, Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA tamoxifen

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker formation of atp by chemiosmotic coupling; n = 253, deaths = 22, log-rank p = 0.557, HR = 1.30 (0.54-3.10).
- z <= -1.1824: biomarker formation of atp by chemiosmotic coupling; n = 253, deaths = 22, log-rank p = 0.557, HR = 1.30 (0.54-3.10).
- z <= -1.2816: biomarker spry regulation of fgf signaling; n = 253, deaths = 22, log-rank p = 0.786, HR = 1.13 (0.48-2.65).
- z <= -1.645: biomarker spry regulation of fgf signaling; n = 253, deaths = 22, log-rank p = 0.786, HR = 1.13 (0.48-2.65).
- z <= -2.1925: biomarker pecam1 interactions; n = 253, deaths = 22, log-rank p = 0.644, HR = 0.82 (0.34-1.94).
