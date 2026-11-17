# Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | GEMCITABINE

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_GEMCITABINE.tsv`.

- treated, 1 pathway(s): p = 0.023, n = 104, deaths = 51, HR 0.52 (95% CI 0.30-0.92)
- treated, 2 pathway(s): p = 0.137, n = 104, deaths = 51, HR 0.66 (95% CI 0.38-1.15)
- treated, 3 pathway(s): p = 0.113, n = 104, deaths = 51, HR 0.64 (95% CI 0.36-1.12)
- treated, 4 pathway(s): p = 0.431, n = 104, deaths = 51, HR 0.80 (95% CI 0.46-1.39)
- treated, 5 pathway(s): p = 0.229, n = 104, deaths = 51, HR 0.71 (95% CI 0.41-1.24)
- treated, 6 pathway(s): p = 0.257, n = 104, deaths = 51, HR 0.72 (95% CI 0.41-1.27)
- treated, 7 pathway(s): p = 0.040, n = 104, deaths = 51, HR 0.55 (95% CI 0.31-0.98)
- treated, 8 pathway(s): p = 0.005, n = 104, deaths = 51, HR 0.44 (95% CI 0.25-0.79)
- treated, 9 pathway(s): p = 0.012, n = 104, deaths = 51, HR 0.49 (95% CI 0.28-0.86)
- treated, 10 pathway(s): p = 0.022, n = 104, deaths = 51, HR 0.52 (95% CI 0.30-0.92)
- untreated, 1 pathway(s): p = 0.017, n = 61, deaths = 35, HR 0.45 (95% CI 0.23-0.88)
- untreated, 2 pathway(s): p = 0.040, n = 61, deaths = 35, HR 0.49 (95% CI 0.25-0.98)
- untreated, 3 pathway(s): p = 0.026, n = 61, deaths = 35, HR 0.45 (95% CI 0.22-0.92)
- untreated, 4 pathway(s): p = 0.002, n = 61, deaths = 35, HR 0.32 (95% CI 0.15-0.68)
- untreated, 5 pathway(s): p = 0.019, n = 61, deaths = 35, HR 0.43 (95% CI 0.21-0.89)
- untreated, 6 pathway(s): p = 0.073, n = 61, deaths = 35, HR 0.53 (95% CI 0.26-1.07)
- untreated, 7 pathway(s): p = 0.032, n = 61, deaths = 35, HR 0.48 (95% CI 0.24-0.95)
- untreated, 8 pathway(s): p = 0.024, n = 61, deaths = 35, HR 0.45 (95% CI 0.22-0.91)
- untreated, 9 pathway(s): p = 0.013, n = 61, deaths = 35, HR 0.42 (95% CI 0.21-0.85)
- untreated, 10 pathway(s): p = 0.007, n = 61, deaths = 35, HR 0.40 (95% CI 0.20-0.79)
