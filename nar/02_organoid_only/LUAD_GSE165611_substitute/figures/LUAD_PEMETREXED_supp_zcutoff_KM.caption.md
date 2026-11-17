# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD pemetrexed

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker metabolism of amino acids and derivatives; n = 59, deaths = 15, log-rank p = 0.00485, HR = 0.19 (0.05-0.69).
- z <= -1.1824: biomarker metabolism of amino acids and derivatives; n = 59, deaths = 15, log-rank p = 0.00485, HR = 0.19 (0.05-0.69).
- z <= -1.2816: biomarker metabolism of amino acids and derivatives; n = 59, deaths = 15, log-rank p = 0.00485, HR = 0.19 (0.05-0.69).
- z <= -1.645: biomarker g1 s transition; n = 59, deaths = 15, log-rank p = 0.0123, HR = 0.22 (0.06-0.81).
- z <= -2.1925: biomarker g1 s transition; n = 59, deaths = 15, log-rank p = 0.0123, HR = 0.22 (0.06-0.81).
