# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD etoposide

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker repair synthesis for gap filling by dna pol in tc ner; n = 15, deaths = 4, log-rank p = 0.854, HR = 0.83 (0.11-6.25).
- z <= -1.1824: biomarker pyrimidine metabolism; n = 15, deaths = 4, log-rank p = 0.815, HR = 1.26 (0.18-9.05).
- z <= -1.2816: biomarker base excision repair; n = 15, deaths = 4, log-rank p = 0.8, HR = 1.30 (0.17-9.74).
- z <= -1.645: biomarker base excision repair; n = 15, deaths = 4, log-rank p = 0.8, HR = 1.30 (0.17-9.74).
- z <= -2.1925: biomarker resolution of ap sites via the single nucleotide replacement pathway; n = 15, deaths = 4, log-rank p = 0.8, HR = 1.30 (0.17-9.74).
