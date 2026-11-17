# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | IRINOTECAN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_IRINOTECAN.tsv`.

- treated, 1 pathway(s): p = 0.257, n = 32, deaths = 15, HR 0.55 (95% CI 0.19-1.56)
- treated, 2 pathway(s): p = 0.937, n = 32, deaths = 15, HR 1.04 (95% CI 0.36-2.99)
- treated, 3 pathway(s): p = 0.549, n = 32, deaths = 15, HR 0.72 (95% CI 0.25-2.10)
- treated, 4 pathway(s): p = 0.509, n = 32, deaths = 15, HR 1.43 (95% CI 0.49-4.16)
- treated, 5 pathway(s): p = 0.607, n = 32, deaths = 15, HR 1.32 (95% CI 0.45-3.87)
- treated, 6 pathway(s): p = 0.607, n = 32, deaths = 15, HR 1.32 (95% CI 0.45-3.87)
- treated, 7 pathway(s): p = 0.978, n = 32, deaths = 15, HR 0.99 (95% CI 0.36-2.73)
- treated, 8 pathway(s): p = 0.607, n = 32, deaths = 15, HR 1.32 (95% CI 0.45-3.87)
- treated, 9 pathway(s): p = 0.607, n = 32, deaths = 15, HR 1.32 (95% CI 0.45-3.87)
- treated, 10 pathway(s): p = 0.842, n = 32, deaths = 15, HR 0.90 (95% CI 0.33-2.50)
- untreated, 1 pathway(s): p = 0.043, n = 302, deaths = 68, HR 0.60 (95% CI 0.36-0.99)
- untreated, 2 pathway(s): p = 0.987, n = 302, deaths = 68, HR 1.00 (95% CI 0.62-1.61)
- untreated, 3 pathway(s): p = 0.857, n = 302, deaths = 68, HR 0.96 (95% CI 0.59-1.55)
- untreated, 4 pathway(s): p = 0.519, n = 302, deaths = 68, HR 0.85 (95% CI 0.53-1.39)
- untreated, 5 pathway(s): p = 0.813, n = 302, deaths = 68, HR 0.94 (95% CI 0.59-1.52)
- untreated, 6 pathway(s): p = 0.572, n = 302, deaths = 68, HR 1.15 (95% CI 0.71-1.85)
- untreated, 7 pathway(s): p = 0.655, n = 302, deaths = 68, HR 1.12 (95% CI 0.69-1.80)
- untreated, 8 pathway(s): p = 0.724, n = 302, deaths = 68, HR 1.09 (95% CI 0.68-1.76)
- untreated, 9 pathway(s): p = 0.862, n = 302, deaths = 68, HR 0.96 (95% CI 0.59-1.55)
- untreated, 10 pathway(s): p = 0.884, n = 302, deaths = 68, HR 1.04 (95% CI 0.64-1.67)
