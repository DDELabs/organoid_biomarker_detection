# Kaplan-Meier grid: Liver cancer organoids (Ji 2023) | LIHC (current GDC) | SORAFENIB

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_LIHC_SORAFENIB.tsv`.

- treated, 1 pathway(s): p = 0.188, n = 30, deaths = 15, HR 0.48 (95% CI 0.15-1.47)
- treated, 2 pathway(s): p = 0.702, n = 30, deaths = 15, HR 1.23 (95% CI 0.42-3.60)
- treated, 3 pathway(s): p = 0.160, n = 30, deaths = 15, HR 0.43 (95% CI 0.13-1.44)
- treated, 4 pathway(s): p = 0.150, n = 30, deaths = 15, HR 0.42 (95% CI 0.12-1.41)
- treated, 5 pathway(s): p = 0.028, n = 30, deaths = 15, HR 0.28 (95% CI 0.08-0.94)
- treated, 6 pathway(s): p = 0.146, n = 30, deaths = 15, HR 0.43 (95% CI 0.13-1.38)
- treated, 7 pathway(s): p = 0.618, n = 30, deaths = 15, HR 0.76 (95% CI 0.26-2.23)
- treated, 8 pathway(s): p = 0.123, n = 30, deaths = 15, HR 0.39 (95% CI 0.12-1.34)
- treated, 9 pathway(s): p = 0.320, n = 30, deaths = 15, HR 0.59 (95% CI 0.20-1.70)
- treated, 10 pathway(s): p = 0.574, n = 30, deaths = 15, HR 0.74 (95% CI 0.26-2.12)
- untreated, 1 pathway(s): p = 0.129, n = 328, deaths = 108, HR 1.34 (95% CI 0.92-1.96)
- untreated, 2 pathway(s): p = 0.315, n = 328, deaths = 108, HR 1.21 (95% CI 0.83-1.77)
- untreated, 3 pathway(s): p = 0.598, n = 328, deaths = 108, HR 0.90 (95% CI 0.62-1.32)
- untreated, 4 pathway(s): p = 0.084, n = 328, deaths = 108, HR 0.72 (95% CI 0.49-1.05)
- untreated, 5 pathway(s): p = 0.163, n = 328, deaths = 108, HR 0.76 (95% CI 0.52-1.12)
- untreated, 6 pathway(s): p = 0.557, n = 328, deaths = 108, HR 0.89 (95% CI 0.61-1.30)
- untreated, 7 pathway(s): p = 0.528, n = 328, deaths = 108, HR 1.13 (95% CI 0.77-1.65)
- untreated, 8 pathway(s): p = 0.943, n = 328, deaths = 108, HR 0.99 (95% CI 0.68-1.44)
- untreated, 9 pathway(s): p = 0.606, n = 328, deaths = 108, HR 1.10 (95% CI 0.76-1.61)
- untreated, 10 pathway(s): p = 0.968, n = 328, deaths = 108, HR 1.01 (95% CI 0.69-1.47)
