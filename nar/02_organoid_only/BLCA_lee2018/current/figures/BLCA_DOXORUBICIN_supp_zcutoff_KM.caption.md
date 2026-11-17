# Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA doxorubicin

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker reversible hydration of carbon dioxide; n = 16, deaths = 8, log-rank p = 0.36, HR = 2.14 (0.40-11.31).
- z <= -1.1824: biomarker smooth muscle contraction; n = 16, deaths = 8, log-rank p = 0.612, HR = 0.68 (0.15-3.05).
- z <= -1.2816: biomarker cytosolic trna aminoacylation; n = 16, deaths = 8, log-rank p = 0.368, HR = 0.50 (0.11-2.32).
- z <= -1.645: biomarker cytosolic trna aminoacylation; n = 16, deaths = 8, log-rank p = 0.368, HR = 0.50 (0.11-2.32).
- z <= -2.1925: biomarker hyaluronan uptake and degradation; n = 16, deaths = 8, log-rank p = 0.597, HR = 0.67 (0.15-3.01).
