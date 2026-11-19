# Lung: a pathway score that predicts survival in people who got no chemotherapy

## What was tested
Lung cohorts were added in W1 (20 builders; 3,878 new samples; table v2). The question was whether pathway scores learned in lung chemo/no-chemo cohorts predict survival (a) in patients who did not get chemotherapy (general prognosis) and (b) as a chemo-benefit score. The prognosis score is trained only on the other cohorts, never on the cohort it is tested in.

## Data
14 tested cohorts of untreated or surgery-only patients (GSE14814 observation arm, GSE37745, GSE42127, GSE68465 no-chemo arm, GSE157011, GSE162520, GSE30219, GSE31210, GSE3141, GSE41271, GSE50081, GSE72094, GSE81089, GSE8894); **2,938 patients, 1,320 events**. GSE29013 (6 events) was not tested. Endpoint overall survival (GSE8894: relapse-free). Score: ridge Cox on 673 pathways, adjusted for age (and stage where known).

## Result 1 (strong). General prognostic score in untreated patients
| score | pathways | pooled HR per SD (95% CI) | p |
|---|---|---|---|
| general (main-effect) prognostic score | all 673 | **1.36 (1.29-1.44)** | 1.1e-28 |
| general prognostic score | 47 cisplatin shortlist | 1.21 (1.14-1.28) | 2.0e-11 |
| chemo-benefit score | all 673 | 0.98 (0.93-1.04) | 0.47 |
| chemo-benefit score | shortlist | 0.99 (0.94-1.05) | 0.75 |
- HR above 1 in 14 of 14 cohorts; p < 0.05 in 9 of 14 (GSE157011 1.20, GSE30219 1.49, GSE31210 1.80, GSE41271 1.72, GSE42127 1.37, GSE68465 1.38 (n=340, 162 deaths, 1.18-1.61), GSE72094 1.41, GSE81089 1.46, GSE8894 1.41). The observation arm of the randomised GSE14814 gives 1.06 (0.69-1.63), 30 deaths.
- Source: `tables/17__lung_untreated_prognosis.tsv`; figure `17__lung_untreated_prognosis`.

## Result 2 (strong). All-model benchmark on the same data (W3)
The original ridge Cox gives pooled C-index 0.601 (0.585-0.619). Network-smoothed L2 gives 0.622 (0.606-0.638), +0.021 (+0.009 to +0.033), paired p 0.010, q 0.020; plain ridge 0.620. The shortlist version is worse, and random shortlist-size sets reach the original in 44% of draws. Source: `tables/19__positives_all_models_lung.tsv`; figure `19__fig_bar_lung_untreated_prognosis`.

## Result 3 (nominal only). Two lab-model-trained tests in the 57-test survey
| model -> TCGA-LUAD patients | n / events | HR (95% CI) | log-rank p | BH q over 57 | shuffle p |
|---|---|---|---|---|---|
| lung adenocarcinoma cell lines, gemcitabine | 25 / 16 | 3.17 (1.06-9.44) | 0.030 | 0.339 | 0.138 |
| NSCLC PDX, erlotinib | 21 / 7 | 7.61 (0.81-71.0) | 0.040 | 0.354 | 0.045 |
Both meet the p < 0.05 arm of the rule but not the q arm, and the intervals are wide. Source: `tables/05__A3_master_results.tsv`; figure `05__A3_Fig1_master_forest`.

## Caveats
- The score is prognostic: it separates survival in untreated patients, but the chemo-specific part does not (HR 0.98), and in the randomised JBR.10-type trial GSE14814 the interaction is HR 1.02 (p 0.95). It is not evidence of chemo benefit prediction.
- Stage was approximated for GSE68465; GSE161537 and GSE162520 are targeted panels (185 of 673 pathways cannot be scored).
- Immunotherapy response (GSE218989 AUC 0.566 p 0.031, GSE161537 and GSE126044 with pan-cancer models) and the GSE29013 interaction are in `../BORDERLINE.md`.

## What I think
Lung gives one solid result: pathway scores carry prognostic information across 14 independent cohorts, and several models do slightly better than the original. It is a prognosis result, not a drug-response result. The two lab-model tests are weak leads.

## Where things are
`figures/`, `tables/`, `source_reports/17__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
