# Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV carboplatin

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker metabolism of polyamines; n = 347, deaths = 218, log-rank p = 0.543, HR = 1.09 (0.83-1.42).
- z <= -1.1824: biomarker metabolism of polyamines; n = 347, deaths = 218, log-rank p = 0.543, HR = 1.09 (0.83-1.42).
- z <= -1.2816: biomarker metabolism of polyamines; n = 347, deaths = 218, log-rank p = 0.543, HR = 1.09 (0.83-1.42).
- z <= -1.645: biomarker gamma carboxylation transport and amino terminal cleavage of proteins; n = 347, deaths = 218, log-rank p = 0.193, HR = 0.84 (0.64-1.09).
- z <= -2.1925: biomarker platelet activation signaling and aggregation; n = 347, deaths = 218, log-rank p = 0.252, HR = 1.17 (0.89-1.53).
