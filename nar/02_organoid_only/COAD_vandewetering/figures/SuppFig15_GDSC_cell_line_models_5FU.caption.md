# Supp. Fig. 15 (colon): 5-FU models trained on cancer cell lines instead of organoids

Same pipeline, but the Ridge model is trained on GDSC1 cell-line 5-FU response (ln IC50; GDSC release 8.5) with
cell-line RNA-seq (Cell Model Passports, log2 TPM+1), either all screened lines with RNA-seq or colorectal lines
only. The top three pathways are each used alone to split 5-FU-treated TCGA-COAD patients (current GDC data and the
authors' 2019 files). Paper: all 893 lines, top pathway p = 0.17; colorectal lines, top pathway p = 0.86.

- all GDSC1 lines (903 lines) -> COAD, rank 1 Cdc6 association with the orc origin complex: n = 129, deaths = 25, p = 0.708, HR 0.86 (95% CI 0.39-1.90)
- all GDSC1 lines (903 lines) -> COAD, rank 2 E2f mediated regulation of dna replication: n = 129, deaths = 25, p = 0.219, HR 0.61 (95% CI 0.27-1.36)
- all GDSC1 lines (903 lines) -> COAD, rank 3 Association of licensing factors with the pre replicative complex: n = 129, deaths = 25, p = 0.420, HR 1.39 (95% CI 0.62-3.08)
- all GDSC1 lines (903 lines) -> COAD_authors, rank 1 Cdc6 association with the orc origin complex: n = 114, deaths = 8, p = 0.583, HR 1.49 (95% CI 0.36-6.26)
- all GDSC1 lines (903 lines) -> COAD_authors, rank 2 E2f mediated regulation of dna replication: n = 114, deaths = 8, p = 0.554, HR 0.65 (95% CI 0.16-2.73)
- all GDSC1 lines (903 lines) -> COAD_authors, rank 3 Association of licensing factors with the pre replicative complex: n = 114, deaths = 8, p = 0.593, HR 1.48 (95% CI 0.35-6.19)
- colorectal GDSC1 lines (43 lines) -> COAD, rank 1 E2f enabled inhibition of pre replication complex formation: n = 129, deaths = 25, p = 0.403, HR 0.71 (95% CI 0.32-1.58)
- colorectal GDSC1 lines (43 lines) -> COAD, rank 2 Pyrimidine metabolism: n = 129, deaths = 25, p = 0.257, HR 1.59 (95% CI 0.71-3.57)
- colorectal GDSC1 lines (43 lines) -> COAD, rank 3 Mitotic g1 g1 s phases: n = 129, deaths = 25, p = 0.774, HR 1.12 (95% CI 0.51-2.49)
- colorectal GDSC1 lines (43 lines) -> COAD_authors, rank 1 E2f enabled inhibition of pre replication complex formation: n = 114, deaths = 8, p = 0.506, HR 0.62 (95% CI 0.15-2.60)
- colorectal GDSC1 lines (43 lines) -> COAD_authors, rank 2 Pyrimidine metabolism: n = 114, deaths = 8, p = 0.251, HR 2.27 (95% CI 0.54-9.51)
- colorectal GDSC1 lines (43 lines) -> COAD_authors, rank 3 Mitotic g1 g1 s phases: n = 114, deaths = 8, p = 0.856, HR 1.14 (95% CI 0.28-4.56)
