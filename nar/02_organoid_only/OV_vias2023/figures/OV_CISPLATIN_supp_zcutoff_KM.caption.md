# Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV cisplatin

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker metabolism of polyamines; n = 123, deaths = 77, log-rank p = 0.262, HR = 0.77 (0.49-1.21).
- z <= -1.1824: biomarker metabolism of polyamines; n = 123, deaths = 77, log-rank p = 0.262, HR = 0.77 (0.49-1.21).
- z <= -1.2816: biomarker metabolism of polyamines; n = 123, deaths = 77, log-rank p = 0.262, HR = 0.77 (0.49-1.21).
- z <= -1.645: biomarker gamma carboxylation transport and amino terminal cleavage of proteins; n = 123, deaths = 77, log-rank p = 0.545, HR = 0.87 (0.55-1.37).
- z <= -2.1925: biomarker platelet activation signaling and aggregation; n = 123, deaths = 77, log-rank p = 0.157, HR = 0.72 (0.46-1.14).
