# Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.958, n = 37, deaths = 18, HR 0.98 (95% CI 0.38-2.48)
- treated, 2 pathway(s): p = 0.205, n = 37, deaths = 18, HR 1.82 (95% CI 0.71-4.66)
- treated, 3 pathway(s): p = 0.097, n = 37, deaths = 18, HR 2.24 (95% CI 0.84-5.92)
- treated, 4 pathway(s): p = 0.470, n = 37, deaths = 18, HR 1.42 (95% CI 0.55-3.69)
- treated, 5 pathway(s): p = 0.142, n = 37, deaths = 18, HR 2.04 (95% CI 0.77-5.40)
- treated, 6 pathway(s): p = 0.130, n = 37, deaths = 18, HR 2.08 (95% CI 0.79-5.50)
- treated, 7 pathway(s): p = 0.450, n = 37, deaths = 18, HR 1.43 (95% CI 0.56-3.65)
- treated, 8 pathway(s): p = 0.130, n = 37, deaths = 18, HR 2.08 (95% CI 0.79-5.50)
- treated, 9 pathway(s): p = 0.141, n = 37, deaths = 18, HR 2.04 (95% CI 0.77-5.40)
- treated, 10 pathway(s): p = 0.909, n = 37, deaths = 18, HR 1.06 (95% CI 0.41-2.75)
- untreated, 1 pathway(s): p = 0.412, n = 61, deaths = 35, HR 0.75 (95% CI 0.39-1.48)
- untreated, 2 pathway(s): p = 0.216, n = 61, deaths = 35, HR 0.65 (95% CI 0.32-1.30)
- untreated, 3 pathway(s): p = 0.259, n = 61, deaths = 35, HR 1.47 (95% CI 0.75-2.89)
- untreated, 4 pathway(s): p = 0.222, n = 61, deaths = 35, HR 1.52 (95% CI 0.77-2.99)
- untreated, 5 pathway(s): p = 0.039, n = 61, deaths = 35, HR 2.05 (95% CI 1.02-4.10)
- untreated, 6 pathway(s): p = 0.363, n = 61, deaths = 35, HR 1.37 (95% CI 0.69-2.70)
- untreated, 7 pathway(s): p = 0.617, n = 61, deaths = 35, HR 1.19 (95% CI 0.60-2.33)
- untreated, 8 pathway(s): p = 0.762, n = 61, deaths = 35, HR 1.11 (95% CI 0.56-2.18)
- untreated, 9 pathway(s): p = 0.141, n = 61, deaths = 35, HR 1.68 (95% CI 0.84-3.39)
- untreated, 10 pathway(s): p = 0.319, n = 61, deaths = 35, HR 1.41 (95% CI 0.72-2.79)
