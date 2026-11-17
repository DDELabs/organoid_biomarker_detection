# Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.223, n = 21, deaths = 11, HR 0.44 (95% CI 0.11-1.71)
- treated, 2 pathway(s): p = 0.469, n = 21, deaths = 11, HR 1.60 (95% CI 0.44-5.82)
- treated, 3 pathway(s): p = 0.945, n = 21, deaths = 11, HR 0.96 (95% CI 0.27-3.42)
- treated, 4 pathway(s): p = 0.945, n = 21, deaths = 11, HR 0.96 (95% CI 0.27-3.42)
- treated, 5 pathway(s): p = 0.833, n = 21, deaths = 11, HR 0.87 (95% CI 0.24-3.12)
- treated, 6 pathway(s): p = 0.701, n = 21, deaths = 11, HR 0.78 (95% CI 0.22-2.79)
- treated, 7 pathway(s): p = 0.882, n = 21, deaths = 11, HR 1.10 (95% CI 0.31-3.95)
- treated, 8 pathway(s): p = 0.975, n = 21, deaths = 11, HR 1.02 (95% CI 0.28-3.67)
- treated, 9 pathway(s): p = 0.631, n = 21, deaths = 11, HR 0.73 (95% CI 0.21-2.62)
- treated, 10 pathway(s): p = 0.223, n = 21, deaths = 11, HR 0.44 (95% CI 0.11-1.71)
- untreated, 1 pathway(s): p = 0.072, n = 61, deaths = 35, HR 0.53 (95% CI 0.26-1.07)
- untreated, 2 pathway(s): p = 0.001, n = 61, deaths = 35, HR 0.30 (95% CI 0.14-0.64)
- untreated, 3 pathway(s): p = 0.139, n = 61, deaths = 35, HR 0.60 (95% CI 0.30-1.19)
- untreated, 4 pathway(s): p = 0.004, n = 61, deaths = 35, HR 0.37 (95% CI 0.18-0.75)
- untreated, 5 pathway(s): p = 0.111, n = 61, deaths = 35, HR 0.58 (95% CI 0.29-1.14)
- untreated, 6 pathway(s): p = 0.081, n = 61, deaths = 35, HR 0.55 (95% CI 0.28-1.09)
- untreated, 7 pathway(s): p = 0.026, n = 61, deaths = 35, HR 0.46 (95% CI 0.23-0.92)
- untreated, 8 pathway(s): p = 0.030, n = 61, deaths = 35, HR 0.47 (95% CI 0.24-0.94)
- untreated, 9 pathway(s): p = 0.036, n = 61, deaths = 35, HR 0.48 (95% CI 0.24-0.97)
- untreated, 10 pathway(s): p = 0.247, n = 61, deaths = 35, HR 0.67 (95% CI 0.34-1.32)
