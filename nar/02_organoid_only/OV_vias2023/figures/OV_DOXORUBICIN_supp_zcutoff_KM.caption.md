# Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV doxorubicin

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker glycosaminoglycan metabolism; n = 117, deaths = 98, log-rank p = 0.347, HR = 0.82 (0.55-1.24).
- z <= -1.1824: biomarker smooth muscle contraction; n = 117, deaths = 98, log-rank p = 0.604, HR = 0.90 (0.60-1.34).
- z <= -1.2816: biomarker resolution of ap sites via the multiple nucleotide patch replacement …; n = 117, deaths = 98, log-rank p = 0.976, HR = 0.99 (0.66-1.50).
- z <= -1.645: biomarker synthesis and interconversion of nucleotide di and triphosphates; n = 117, deaths = 98, log-rank p = 0.615, HR = 1.11 (0.74-1.67).
- z <= -2.1925: biomarker mitochondrial trna aminoacylation; n = 117, deaths = 98, log-rank p = 0.986, HR = 1.00 (0.67-1.50).
