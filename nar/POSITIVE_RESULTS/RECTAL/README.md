# Rectal cancer (GSE87211): 5-FU organoid signatures predict disease-free survival

## What was tested
Frozen pathway signatures built from organoid or PDX 5-FU data (top-10 pathways by absolute ridge coefficient, nothing refitted) were scored in an independent cohort: GSE87211, locally advanced rectal cancer given 5-FU-based chemoradiotherapy (5-FU + radiotherapy, with or without oxaliplatin or cetuximab). Outcome: disease-free survival. Direction fixed beforehand: higher score = predicted more resistant, so HR above 1.

## Data
189 patients with rectal tumours and DFS, 42 recurrences. All patients received 5-FU; there is no untreated comparator arm, so this tests prognosis in treated patients.

## Results
| id | signature | HR per SD (95% CI) | p | BH q (196 tests) | random-signature null p |
|---|---|---|---|---|---|
| RE1 | van de Wetering organoids, current data, top-10 | **1.91 (1.43-2.56)** | 1.1e-05 | **0.002** | 0.001 |
| RE2 | Gao colorectal PDX, top-10 | **1.65 (1.22-2.22)** | 0.001 | **0.097** | 0.026 |
| RE3 | paper (authors' 2019 files), rank-1 BH3-only | 1.38 (1.00-1.89) | 0.048 | 0.650 | 0.326 |
| RE4 | paper (authors' 2019 files), top-10 | 1.44 (1.04-1.98) | 0.026 | 0.650 | 0.214 |
| RE5 | all-model benchmark: frozen top-10 C-index | C 0.684 (0.597-0.767) | 1.1e-05 | | |
| RE6 | GSE87211 oxaliplatin-response x joint-score interaction (90 vs 111 patients, 13 vs 22 responders) | OR 2.24 (1.03-4.86) | 0.041 | not corrected (1 of 32) | |
(`POSITIVE_INDEX.tsv` gives RE1 as 1.91; the 05 report text rounds the same fit to 1.92.) Sources: `tables/05__B1_external_results.tsv`, `tables/05__B1_frozen_signatures.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/03__interaction_all.tsv`.

## How robust is it
- RE1 is the only one of 196 external tests to survive BH at q < 0.05; RE2 passes q < 0.10. Both beat random 10-pathway signatures (null p 0.001 and 0.026).
- All-model benchmark: no added model beats the frozen top-10 (best linear SVM 0.665, difference -0.019, p 0.51). Random shortlist-size sets reach 0.684 in 0.5% of draws. The rank-1 pathway alone is much weaker (C 0.573).
- The signal is consistent across signatures built from different lab models (organoids, PDX).

## Caveats
- Prognostic in treated patients only. With no untreated arm it cannot be called drug-specific, and the same kind of score separates survival in untreated patients in other cohorts (nar/05, nar/14).
- The W2 response model for the same cohort (n 201, 35 responders, response = no tumour at surgery) gives AUC 0.427 (p 0.18), i.e. no response signal; the DFS signal is not a response signal.
- 42 events.

## What I think
This is the best external validation of the organoid-derived 5-FU pathways, and it is stable across signatures. It most likely reflects a general recurrence-risk signature that the organoid pathways pick up, not proof that organoid response transfers.

## Where things are
`figures/`, `tables/`. Checksums in `../COPY_MANIFEST.tsv`.
