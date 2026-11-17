# Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.497, n = 21, deaths = 11, HR 0.65 (95% CI 0.18-2.30)
- treated, 2 pathway(s): p = 0.203, n = 21, deaths = 11, HR 0.42 (95% CI 0.11-1.66)
- treated, 3 pathway(s): p = 0.434, n = 21, deaths = 11, HR 1.69 (95% CI 0.45-6.39)
- treated, 4 pathway(s): p = 0.434, n = 21, deaths = 11, HR 1.69 (95% CI 0.45-6.39)
- treated, 5 pathway(s): p = 0.434, n = 21, deaths = 11, HR 1.69 (95% CI 0.45-6.39)
- treated, 6 pathway(s): p = 0.459, n = 21, deaths = 11, HR 0.62 (95% CI 0.17-2.22)
- treated, 7 pathway(s): p = 0.818, n = 21, deaths = 11, HR 0.86 (95% CI 0.24-3.11)
- treated, 8 pathway(s): p = 0.452, n = 21, deaths = 11, HR 0.60 (95% CI 0.15-2.33)
- treated, 9 pathway(s): p = 0.948, n = 21, deaths = 11, HR 1.05 (95% CI 0.28-3.95)
- treated, 10 pathway(s): p = 0.434, n = 21, deaths = 11, HR 1.69 (95% CI 0.45-6.39)
- untreated, 1 pathway(s): p = 0.242, n = 61, deaths = 35, HR 1.49 (95% CI 0.76-2.91)
- untreated, 2 pathway(s): p = 0.328, n = 61, deaths = 35, HR 1.40 (95% CI 0.71-2.73)
- untreated, 3 pathway(s): p = 0.194, n = 61, deaths = 35, HR 1.56 (95% CI 0.79-3.05)
- untreated, 4 pathway(s): p = 0.702, n = 61, deaths = 35, HR 0.88 (95% CI 0.45-1.71)
- untreated, 5 pathway(s): p = 0.607, n = 61, deaths = 35, HR 1.19 (95% CI 0.61-2.31)
- untreated, 6 pathway(s): p = 0.798, n = 61, deaths = 35, HR 1.09 (95% CI 0.56-2.12)
- untreated, 7 pathway(s): p = 0.128, n = 61, deaths = 35, HR 0.60 (95% CI 0.31-1.17)
- untreated, 8 pathway(s): p = 0.013, n = 61, deaths = 35, HR 0.43 (95% CI 0.22-0.85)
- untreated, 9 pathway(s): p = 0.005, n = 61, deaths = 35, HR 0.38 (95% CI 0.19-0.77)
- untreated, 10 pathway(s): p = 0.399, n = 61, deaths = 35, HR 0.75 (95% CI 0.38-1.47)
