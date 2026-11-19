# Bladder: cisplatin organoid model, ERCC2 mutation check, isogenic cells and the joint model

## What was tested
The paper's bladder cisplatin result (Kong et al. 2020, Fig 2e-h, 4c, 5b) reproduced with the paper's own inputs and with current data; then extra tests. Method: pathways close to cisplatin targets in the STRING network; ssGSEA scores; ridge on 9 organoid responses (7 patients; two lines are recurrences of the same patients; audit I07); patients scored with the top pathway and split at the median; log-rank on overall survival.

## Data
Lee 2018 bladder organoids (9 with cisplatin response; response = refit AUC, the paper's IC50 is not public); TCGA-BLCA cisplatin-treated 79 patients (20 deaths on paper data; 28 on current GDC data); 288 untreated; GSE31684, GSE48276 etc. for the joint model.

## Results
| id | test | n / events | effect (95% CI) | p | tier |
|---|---|---|---|---|---|
| BL1 | paper-data rank-1 pathway, treated TCGA-BLCA | 79 / 20 | HR 2.63 (0.99-6.96) | 0.044 | nominal (BH q 0.354 over 57; shuffle p 0.30) |
| BL1 | paper's amino-acid pathway, same patients | 79 / 20 | HR 2.66 (1.06-6.69) | 0.031 | nominal |
| BL2 | current GDC data, response-shuffle test | 79 / 28 | log-rank p 0.067, HR 2.03 (0.94-4.38) | shuffle p **0.036** | nominal (the patient log-rank itself is not < 0.05) |
| BL3 | ERCC2-mutant vs wild type, paper data, amino-acid pathway | 36 / 367 | one-sided Mann-Whitney | **0.0030** | strong |
| BL3 | same, current data | 36 / 366 | | **0.00055** | strong |
| BL4 | isogenic T24 sensitive vs resistant, amino-acid pathway | 3 vs 3 | t 12.7 | **0.00022** | strong (tiny n) |
| BL5 | joint organoid+patient model, 4 held-out cohorts | 157 / 99 | HR 1.24 (1.03-1.50) per SD risk | 0.024 | nominal |
| BL5 | patient-only model, same cohorts | 157 / 99 | HR 1.45 (1.18-1.79) | 5.1e-04 | strong |
| BL6 | joint model, gemcitabine, TCGA-BLCA 5-fold | 64 / 30 | HR 1.53 (1.09-2.14) | 0.014 | nominal |
| BL7 | joint model, GSE48276 MVAC x score interaction | 53 (8 treated events) | HR 5.30 (1.69-16.66) | 0.0043 | nominal (1 of 32) |
| BL8 | current-data rank-1 signature, GSE31684 interaction | 35 treated / 25 events | HR 1.82 (1.07-3.09) | 0.028 | nominal (q 0.650 over 196) |
| BL9 | ERCC2, all-model benchmark, paper amino-acid pathway AUC | 397 / 36 | AUC 0.665 (0.571-0.752) | 5.4e-04 | strong |
Sources: `tables/05__A3_master_results.tsv`, `tables/02__ERCC2_mutation_concordance.tsv`, `tables/02__isogenic_T24_tests.tsv`, `tables/02__paper_biomarker_amino_acid_KM.tsv`, `tables/03__summary_main.tsv`, `tables/03__interaction_all.tsv`, `tables/05__B1_external_results.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/13__G8_ERCC2_null_and_adjustment.tsv`.

## How robust is it
- ERCC2 is the sturdy one. Share of all pathways (both signs) as strong as the amino-acid pathway: 0.017 (paper data), 0.002 (current data); share of random proximal-size sets: 0.063 and 0.000. It holds after tumour-mutation-burden adjustment (p 0.024 paper, 0.008 current; OR per SD 0.63 and 0.54). The paper's rank-1 pathway does not carry it (p 0.88 / 0.33).
- Survival: with 20 deaths a p of 0.04 is easy to reach by chance (29.6% of response-shuffled runs also reach p < 0.05 on paper data). One line per patient (7 lines) leaves the treated statistics unchanged (p 0.044) but the shuffle p moves to 0.263. The amino-acid pathway is rank 6 in our refit (rank 30 of 47 on current data, treated p 0.063, untreated p 0.021).
- Isogenic T24: only the paper's NES for 3 vs 3 samples exist, so p 0.00022 rests on the paper's processed values. Our rank-1 pathway (p 0.017) points the wrong way there.
- Joint model: the organoid-only model is not significant (HR 1.07, p 0.45; gemcitabine 0.90, p 0.59). The gain comes from the patient part (patient-only 1.45); joint minus patient-only is not above zero in any of 20 unit x outcome comparisons (nar/14 G16). Fold-aligned joint HR 1.27, p 0.013.
- All-model benchmark: no added model beats the paper's amino-acid pathway for ERCC2 (best MLP 0.637); a simple proliferation score gets 0.634, so this signal is not specific to the network method.

## Caveats
Nine organoids from seven patients; AUC instead of IC50; TCGA treatment records are "ever given"; untreated groups differ in age and stage; all survival tests have 20-28 events. The five bladder cohorts' paper-signature tests lean the right way (OR/HR 1.16-1.57) but none is below 0.05.

## What I think
ERCC2 concordance is the strongest bladder result and it is specific. The cisplatin survival split is at the edge: p 0.03-0.07 on 20-28 deaths, helped by the choice of pathway. The joint model is a real held-out signal in bladder, but it is carried by patient data.

## Where things are
`figures/`, `tables/`, `source_reports/03__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
