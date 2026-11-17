# Kaplan-Meier grid: Colorectal PDX (Isella 2017, cetuximab) | COAD (current GDC) | CETUXIMAB

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_CETUXIMAB.tsv`.

- treated, 1 pathway(s): p = 0.648, n = 15, deaths = 7, HR 0.70 (95% CI 0.15-3.22)
- treated, 2 pathway(s): p = 0.755, n = 15, deaths = 7, HR 1.27 (95% CI 0.28-5.84)
- treated, 3 pathway(s): p = 0.755, n = 15, deaths = 7, HR 1.27 (95% CI 0.28-5.84)
- treated, 4 pathway(s): p = 0.755, n = 15, deaths = 7, HR 1.27 (95% CI 0.28-5.84)
- treated, 5 pathway(s): p = 0.755, n = 15, deaths = 7, HR 1.27 (95% CI 0.28-5.84)
- treated, 6 pathway(s): p = 0.755, n = 15, deaths = 7, HR 1.27 (95% CI 0.28-5.84)
- treated, 7 pathway(s): p = 0.916, n = 15, deaths = 7, HR 0.92 (95% CI 0.20-4.22)
- treated, 8 pathway(s): p = 0.484, n = 15, deaths = 7, HR 0.58 (95% CI 0.13-2.70)
- treated, 9 pathway(s): p = 0.484, n = 15, deaths = 7, HR 0.58 (95% CI 0.13-2.70)
- treated, 10 pathway(s): p = 0.393, n = 15, deaths = 7, HR 0.49 (95% CI 0.09-2.62)
- untreated, 1 pathway(s): p = 0.046, n = 302, deaths = 68, HR 1.63 (95% CI 1.01-2.65)
- untreated, 2 pathway(s): p = 0.019, n = 302, deaths = 68, HR 1.77 (95% CI 1.09-2.88)
- untreated, 3 pathway(s): p = 0.008, n = 302, deaths = 68, HR 1.92 (95% CI 1.18-3.12)
- untreated, 4 pathway(s): p = 0.005, n = 302, deaths = 68, HR 1.97 (95% CI 1.21-3.21)
- untreated, 5 pathway(s): p = 0.001, n = 302, deaths = 68, HR 2.21 (95% CI 1.35-3.61)
- untreated, 6 pathway(s): p = 0.002, n = 302, deaths = 68, HR 2.17 (95% CI 1.33-3.55)
- untreated, 7 pathway(s): p = 0.004, n = 302, deaths = 68, HR 2.03 (95% CI 1.24-3.31)
- untreated, 8 pathway(s): p = 0.026, n = 302, deaths = 68, HR 1.72 (95% CI 1.06-2.79)
- untreated, 9 pathway(s): p = 0.021, n = 302, deaths = 68, HR 1.75 (95% CI 1.08-2.84)
- untreated, 10 pathway(s): p = 0.043, n = 302, deaths = 68, HR 1.63 (95% CI 1.01-2.64)
