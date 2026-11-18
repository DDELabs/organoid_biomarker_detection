# Is the paper's 5-FU marker really a 5-FU marker, or an oxaliplatin / FOLFOX marker?

Most TCGA "5-FU" patients also had oxaliplatin (FOLFOX). The paper's marker is the "Activation of BH3-only proteins"
pathway. Here we look at it against oxaliplatin in lab models and within regimen groups in patients.

- **a** 19 colorectal organoids: correlation of the BH3-only score with oxaliplatin (red) and 5-FU (orange) response.
  Filled = our recomputed score, open = authors' score. Negative = higher score goes with more sensitivity.
- **b** Cancer cell lines (GDSC1, GDSC2, CTRPv2, PRISM): colorectal lines only (filled) and all lines with cancer type
  removed (open).
- **c** TCGA colon 5-FU-treated patients: hazard ratio per standard deviation of the frozen BH3-only prediction
  (higher = predicted more resistant; HR > 1 = the paper's direction), in all patients and separately in patients
  with and without oxaliplatin. The title gives the score x oxaliplatin interaction p (does the marker work
  differently with oxaliplatin?).
- **d** GSE39582 (French colon cancer cohort, stage II-III): the same frozen score in each adjuvant regimen and in
  untreated patients, for relapse-free (RFS) and overall survival (OS), with interaction tests.

- TCGA colon, authors' 2019 files, OS, all 5-FU patients: n = 114, events = 8, HR per SD 1.52 (0.94-2.46), p = 0.091
- TCGA colon, authors' 2019 files, OS, 5-FU + oxaliplatin (FOLFOX-type): n = 94, events = 6, HR per SD 1.67 (0.94-2.98), p = 0.082
- TCGA colon, authors' 2019 files, OS, 5-FU without oxaliplatin: n = 20, events = 2, HR per SD 2.87 (0.25-33.19), p = 0.399
- TCGA colon, current data, OS, all 5-FU patients: n = 129, events = 25, HR per SD 1.20 (0.82-1.75), p = 0.357, interaction HR 0.88 (0.33-2.31), p = 0.789
- TCGA colon, current data, OS, 5-FU + oxaliplatin (FOLFOX-type): n = 106, events = 18, HR per SD 1.22 (0.77-1.92), p = 0.398
- TCGA colon, current data, OS, 5-FU without oxaliplatin: n = 23, events = 7, HR per SD 1.50 (0.65-3.48), p = 0.341
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, 5FU alone: n = 76, events = 28, HR per SD 0.87 (0.59-1.28), p = 0.481
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, FUFOL (5-FU + leucovorin): n = 50, events = 13, HR per SD 1.12 (0.66-1.88), p = 0.677
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, 5FU or FUFOL (no oxaliplatin): n = 126, events = 41, HR per SD 0.96 (0.70-1.31), p = 0.788
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, FOLFOX: n = 12, events = 7, HR per SD 0.69 (0.30-1.63), p = 0.400
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, no chemotherapy: n = 258, events = 62, HR per SD 1.20 (0.93-1.54), p = 0.153
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, FOLFOX vs 5FU/FUFOL (interaction): n = 138, events = 48, interaction HR 0.59 (0.15-2.29), p = 0.443
- GSE39582 stage II-III, frozen BH3-only (rank score), RFS, 5FU/FUFOL vs no chemotherapy (interaction): n = 384, events = 103, interaction HR 0.79 (0.53-1.19), p = 0.262
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, 5FU alone: n = 76, events = 19, HR per SD 1.00 (0.62-1.59), p = 0.987
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, FUFOL (5-FU + leucovorin): n = 51, events = 16, HR per SD 1.18 (0.73-1.89), p = 0.495
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, 5FU or FUFOL (no oxaliplatin): n = 127, events = 35, HR per SD 1.08 (0.77-1.51), p = 0.666
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, FOLFOX: n = 12, events = 5, HR per SD 1.15 (0.49-2.66), p = 0.749
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, no chemotherapy: n = 262, events = 89, HR per SD 0.98 (0.79-1.20), p = 0.814
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, FOLFOX vs 5FU/FUFOL (interaction): n = 139, events = 40, interaction HR 1.15 (0.31-4.27), p = 0.835
- GSE39582 stage II-III, frozen BH3-only (rank score), OS, 5FU/FUFOL vs no chemotherapy (interaction): n = 389, events = 124, interaction HR 1.11 (0.75-1.64), p = 0.618
