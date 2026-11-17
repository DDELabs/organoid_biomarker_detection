# Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD gemcitabine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker g1 s specific transcription; n = 25, deaths = 16, log-rank p = 0.0297, HR = 3.17 (1.06-9.44).
- z <= -1.1824: biomarker g1 s specific transcription; n = 25, deaths = 16, log-rank p = 0.0297, HR = 3.17 (1.06-9.44).
- z <= -1.2816: biomarker g1 s specific transcription; n = 25, deaths = 16, log-rank p = 0.0297, HR = 3.17 (1.06-9.44).
- z <= -1.645: biomarker g1 s transition; n = 25, deaths = 16, log-rank p = 0.221, HR = 0.53 (0.19-1.48).
- z <= -2.1925: biomarker g1 s transition; n = 25, deaths = 16, log-rank p = 0.221, HR = 0.53 (0.19-1.48).
