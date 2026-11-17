# Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA gemcitabine

Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.
-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.

- z <= -1.0: biomarker smad2 smad3 smad4 heterotrimer regulates transcription; n = 88, deaths = 34, log-rank p = 0.454, HR = 0.77 (0.39-1.52).
- z <= -1.1824: biomarker activation of bh3 only proteins; n = 88, deaths = 34, log-rank p = 0.509, HR = 0.80 (0.41-1.56).
- z <= -1.2816: biomarker activation of bh3 only proteins; n = 88, deaths = 34, log-rank p = 0.509, HR = 0.80 (0.41-1.56).
- z <= -1.645: biomarker activation of bh3 only proteins; n = 88, deaths = 34, log-rank p = 0.509, HR = 0.80 (0.41-1.56).
- z <= -2.1925: biomarker purine catabolism; n = 88, deaths = 34, log-rank p = 0.857, HR = 1.06 (0.54-2.09).
