# Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | IRINOTECAN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_IRINOTECAN.tsv`.

- treated, 1 pathway(s): p = 0.302, n = 13, deaths = 5, HR 0.33 (95% CI 0.04-3.00)
- treated, 2 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 3 pathway(s): p = 0.302, n = 13, deaths = 5, HR 0.33 (95% CI 0.04-3.00)
- treated, 4 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 5 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 6 pathway(s): p = 0.832, n = 13, deaths = 5, HR 1.22 (95% CI 0.19-7.65)
- treated, 7 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 8 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 9 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- treated, 10 pathway(s): p = 0.552, n = 13, deaths = 5, HR 1.84 (95% CI 0.24-14.10)
- untreated, 1 pathway(s): p = 0.021, n = 61, deaths = 35, HR 2.22 (95% CI 1.11-4.43)
- untreated, 2 pathway(s): p = 0.004, n = 61, deaths = 35, HR 2.69 (95% CI 1.35-5.38)
- untreated, 3 pathway(s): p = 0.009, n = 61, deaths = 35, HR 2.46 (95% CI 1.23-4.92)
- untreated, 4 pathway(s): p = 0.193, n = 61, deaths = 35, HR 1.56 (95% CI 0.80-3.03)
- untreated, 5 pathway(s): p = 0.039, n = 61, deaths = 35, HR 2.01 (95% CI 1.02-3.96)
- untreated, 6 pathway(s): p = 0.059, n = 61, deaths = 35, HR 1.89 (95% CI 0.97-3.71)
- untreated, 7 pathway(s): p = 0.013, n = 61, deaths = 35, HR 2.34 (95% CI 1.18-4.64)
- untreated, 8 pathway(s): p = 0.251, n = 61, deaths = 35, HR 1.48 (95% CI 0.76-2.89)
- untreated, 9 pathway(s): p = 0.594, n = 61, deaths = 35, HR 1.20 (95% CI 0.62-2.33)
- untreated, 10 pathway(s): p = 0.594, n = 61, deaths = 35, HR 1.20 (95% CI 0.62-2.33)
