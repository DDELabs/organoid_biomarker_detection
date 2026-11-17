# Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | SN-38 (models) / IRINOTECAN (patients)

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_SN-38.tsv`.

- treated, 1 pathway(s): p = 0.570, n = 13, deaths = 5, HR 0.59 (95% CI 0.10-3.63)
- treated, 2 pathway(s): p = 0.570, n = 13, deaths = 5, HR 0.59 (95% CI 0.10-3.63)
- treated, 3 pathway(s): p = 0.570, n = 13, deaths = 5, HR 0.59 (95% CI 0.10-3.63)
- treated, 4 pathway(s): p = 0.411, n = 13, deaths = 5, HR 0.47 (95% CI 0.08-2.91)
- treated, 5 pathway(s): p = 0.570, n = 13, deaths = 5, HR 0.59 (95% CI 0.10-3.63)
- treated, 6 pathway(s): p = 0.570, n = 13, deaths = 5, HR 0.59 (95% CI 0.10-3.63)
- treated, 7 pathway(s): p = 0.210, n = 13, deaths = 5, HR 0.26 (95% CI 0.03-2.52)
- treated, 8 pathway(s): p = 0.252, n = 13, deaths = 5, HR 0.28 (95% CI 0.03-2.80)
- treated, 9 pathway(s): p = 0.023, n = 13, deaths = 5, HR not estimable (0 deaths in one group)
- treated, 10 pathway(s): p = 0.023, n = 13, deaths = 5, HR not estimable (0 deaths in one group)
- untreated, 1 pathway(s): p = 0.186, n = 61, deaths = 35, HR 1.58 (95% CI 0.80-3.11)
- untreated, 2 pathway(s): p = 0.201, n = 61, deaths = 35, HR 1.55 (95% CI 0.79-3.05)
- untreated, 3 pathway(s): p = 0.562, n = 61, deaths = 35, HR 1.22 (95% CI 0.63-2.37)
- untreated, 4 pathway(s): p = 0.797, n = 61, deaths = 35, HR 1.09 (95% CI 0.56-2.13)
- untreated, 5 pathway(s): p = 0.463, n = 61, deaths = 35, HR 1.28 (95% CI 0.66-2.50)
- untreated, 6 pathway(s): p = 0.952, n = 61, deaths = 35, HR 1.02 (95% CI 0.52-1.99)
- untreated, 7 pathway(s): p = 0.996, n = 61, deaths = 35, HR 1.00 (95% CI 0.51-1.94)
- untreated, 8 pathway(s): p = 0.583, n = 61, deaths = 35, HR 1.21 (95% CI 0.62-2.35)
- untreated, 9 pathway(s): p = 0.996, n = 61, deaths = 35, HR 1.00 (95% CI 0.51-1.94)
- untreated, 10 pathway(s): p = 0.606, n = 61, deaths = 35, HR 0.84 (95% CI 0.43-1.64)
