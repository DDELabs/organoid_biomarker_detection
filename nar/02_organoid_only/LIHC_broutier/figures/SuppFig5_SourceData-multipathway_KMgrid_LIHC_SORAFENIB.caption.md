# Kaplan-Meier grid: Liver cancer organoids (Broutier 2017) | LIHC (current GDC) | SORAFENIB

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_LIHC_SORAFENIB.tsv`.

- treated, 1 pathway(s): p = 0.139, n = 30, deaths = 15, HR 0.43 (95% CI 0.14-1.35)
- treated, 2 pathway(s): p = 0.259, n = 30, deaths = 15, HR 0.54 (95% CI 0.18-1.60)
- treated, 3 pathway(s): p = 0.481, n = 30, deaths = 15, HR 1.52 (95% CI 0.47-4.88)
- treated, 4 pathway(s): p = 0.556, n = 30, deaths = 15, HR 1.46 (95% CI 0.41-5.13)
- treated, 5 pathway(s): p = 0.470, n = 30, deaths = 15, HR 1.56 (95% CI 0.46-5.30)
- treated, 6 pathway(s): p = 0.411, n = 30, deaths = 15, HR 1.66 (95% CI 0.49-5.60)
- treated, 7 pathway(s): p = 0.305, n = 30, deaths = 15, HR 1.83 (95% CI 0.57-5.90)
- treated, 8 pathway(s): p = 0.143, n = 30, deaths = 15, HR 2.50 (95% CI 0.71-8.81)
- treated, 9 pathway(s): p = 0.143, n = 30, deaths = 15, HR 2.50 (95% CI 0.71-8.81)
- treated, 10 pathway(s): p = 0.143, n = 30, deaths = 15, HR 2.50 (95% CI 0.71-8.81)
- untreated, 1 pathway(s): p = 0.374, n = 328, deaths = 108, HR 1.19 (95% CI 0.81-1.73)
- untreated, 2 pathway(s): p = 0.757, n = 328, deaths = 108, HR 0.94 (95% CI 0.64-1.38)
- untreated, 3 pathway(s): p = 0.825, n = 328, deaths = 108, HR 1.04 (95% CI 0.71-1.53)
- untreated, 4 pathway(s): p = 0.214, n = 328, deaths = 108, HR 0.78 (95% CI 0.53-1.15)
- untreated, 5 pathway(s): p = 0.434, n = 328, deaths = 108, HR 0.86 (95% CI 0.59-1.26)
- untreated, 6 pathway(s): p = 0.814, n = 328, deaths = 108, HR 0.96 (95% CI 0.65-1.40)
- untreated, 7 pathway(s): p = 0.876, n = 328, deaths = 108, HR 1.03 (95% CI 0.71-1.50)
- untreated, 8 pathway(s): p = 0.727, n = 328, deaths = 108, HR 1.07 (95% CI 0.73-1.56)
- untreated, 9 pathway(s): p = 0.911, n = 328, deaths = 108, HR 1.02 (95% CI 0.70-1.49)
- untreated, 10 pathway(s): p = 0.865, n = 328, deaths = 108, HR 1.03 (95% CI 0.71-1.51)
