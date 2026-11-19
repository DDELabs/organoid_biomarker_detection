# Ovarian cancer: response cohorts, joint models, olaparib signature, A2780 isogenic, BRCA

## What was tested
Ovarian chemotherapy response prediction from the breast/pan-ovarian pathway model, the joint organoid+patient model, the user's olaparib data, and sanity checks.

## Results
| id | test | n / events | effect (95% CI) | p | tier |
|---|---|---|---|---|---|
| OV1 | GSE30161 response (carboplatin/cisplatin +/- paclitaxel; CR vs PR/PD), trained on other ovarian cohorts | 55 / 32 CR | AUC 0.728 (0.596-0.860) | **0.0043**, q 0.014 | strong |
| OV2 | GSE9891 response (RFS >= 12 months), same | 229 / 146 | AUC 0.624 (0.551-0.697) | **0.0018**, q 0.006 | strong |
| OV3 | olaparib correlation signature (6 ovarian lines, user's GSE229119) in platinum-treated GSE63885 | 75 / 66 deaths | HR 2.09 (1.24-3.50) | 0.0045 (exact permutation 0.014) | nominal |
| OV4 | BRCA1/2 somatic mutated vs wild type, platinum-treated TCGA-OV, PFS | 253 / 184 (16 mutated) | HR 0.42 (0.21-0.82) | 0.0087 | nominal |
| OV4 | same, OS | 253 / 149 | HR 0.47 (0.22-1.01) | 0.049 | nominal |
| OV5 | A2780 vs cisplatin-resistant (GSE15709), ovarian-organoid rank-1 pathway (polyamines) | 5 vs 5 | t 3.03 | 0.017 | nominal |
| OV5 | same pair, lung cell-line rank-1 pathway | 5 vs 5 | t -2.65 | 0.029 | nominal |
| OV6 | joint model, paclitaxel with platinum, 3 held-out cohorts | 564 / 426 | HR 1.14 (1.03-1.26) risk scale | 0.0099 | strong |
| OV6 | patient-only | 564 / 426 | HR 1.22 (1.10-1.34) | 1.1e-04 | strong |
| OV7 | joint model, platinum (cell-line proxy), 3 cohorts | 683 / 523 | HR 1.17 (1.07-1.29) | 6.6e-04 | strong |
| OV7 | patient-only | 683 / 523 | HR 1.19 (1.08-1.30) | 2.5e-04 | strong |
| OV8 | Vias organoid rank-1 signature, GSE9891 RFS | 240 / 179 | HR 1.23 (1.06-1.43) | 0.0059 | nominal (q 0.379 over 196) |
Sources: `tables/18__part1_v2.tsv`, `tables/22__ovary_signature_survival.tsv`, `tables/22__ovary_BRCA_somatic_survival_platinum.tsv`, `tables/10__Fig4bc_isogenic_other_cancers.tsv`, `tables/03__summary_main.tsv`, `tables/05__B1_external_results.tsv`, `tables/19__positives_all_models_rest.tsv`.

## Best ovarian result and whether it meets the rule
The best ovarian result by size and evidence is OV2 (GSE9891, 229 patients, AUC 0.624, p 0.0018) and the joint platinum survival model (OV7, C 0.583, p 6e-08 in the benchmark file). Both meet the rule. The best single lab-model-driven result is OV3 (olaparib signature, GSE63885) and it meets the p < 0.05 arm only nominally.

## How robust is it
- OV1/OV2: from the W2 run (pool of ovarian cohorts). GSE9891's AUC was 0.565 with the smaller v1 pool; the other three ovarian response cohorts (GSE15622 0.483, GSE51373 0.615, GSE63885 0.611) are not significant.
- OV3 is exploratory. The same signature gives p 0.65 (TCGA-OV) and 0.54 (GSE9891); with 6 lines a perfect correlation occurs by chance in 69% of orderings; ridge, SVM and L1 trained on the same lines give p 0.27-0.38; XGBoost has C 0.580 (not different from the original). Separately, the paper-method olaparib model points the other way in GSE9891 (HR 0.58, p 0.005) and TCGA-OV: see `../BORDERLINE.md`.
- OV6/OV7: carried by patient data (organoid- or cell-line-only HR 0.97 and 1.06, not significant). Platinum benchmark: joint C 0.583 and AUC 0.692 are not beaten by any added model.
- OV4 is a clinical sanity check on the cohort (germline BRCA is not public; 16 somatic mutants), not a result of the pathway model.

## What I think
Ovarian gives a real but modest signal in the larger response cohorts and in survival, mostly from patient-trained models. The olaparib signature in GSE63885 is intriguing, but with six lines and one cohort I would keep it exploratory.

## Where things are
`figures/`, `tables/`, `source_reports/22__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
