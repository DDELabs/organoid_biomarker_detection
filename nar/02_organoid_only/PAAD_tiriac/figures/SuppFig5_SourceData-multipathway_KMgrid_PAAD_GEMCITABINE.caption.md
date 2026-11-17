# Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | GEMCITABINE

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_GEMCITABINE.tsv`.

- treated, 1 pathway(s): p = 0.223, n = 104, deaths = 51, HR 1.41 (95% CI 0.81-2.46)
- treated, 2 pathway(s): p = 0.382, n = 104, deaths = 51, HR 1.28 (95% CI 0.74-2.23)
- treated, 3 pathway(s): p = 0.241, n = 104, deaths = 51, HR 1.39 (95% CI 0.80-2.42)
- treated, 4 pathway(s): p = 0.423, n = 104, deaths = 51, HR 1.26 (95% CI 0.72-2.21)
- treated, 5 pathway(s): p = 0.506, n = 104, deaths = 51, HR 1.21 (95% CI 0.69-2.13)
- treated, 6 pathway(s): p = 0.762, n = 104, deaths = 51, HR 0.92 (95% CI 0.52-1.62)
- treated, 7 pathway(s): p = 0.446, n = 104, deaths = 51, HR 0.80 (95% CI 0.45-1.42)
- treated, 8 pathway(s): p = 0.023, n = 104, deaths = 51, HR 0.52 (95% CI 0.29-0.92)
- treated, 9 pathway(s): p = 0.031, n = 104, deaths = 51, HR 0.54 (95% CI 0.31-0.95)
- treated, 10 pathway(s): p = 0.263, n = 104, deaths = 51, HR 0.72 (95% CI 0.41-1.28)
- untreated, 1 pathway(s): p = 0.897, n = 61, deaths = 35, HR 1.05 (95% CI 0.54-2.05)
- untreated, 2 pathway(s): p = 0.191, n = 61, deaths = 35, HR 1.56 (95% CI 0.80-3.05)
- untreated, 3 pathway(s): p = 0.711, n = 61, deaths = 35, HR 1.13 (95% CI 0.58-2.21)
- untreated, 4 pathway(s): p = 0.354, n = 61, deaths = 35, HR 1.37 (95% CI 0.70-2.69)
- untreated, 5 pathway(s): p = 0.672, n = 61, deaths = 35, HR 1.15 (95% CI 0.59-2.25)
- untreated, 6 pathway(s): p = 0.459, n = 61, deaths = 35, HR 1.29 (95% CI 0.66-2.54)
- untreated, 7 pathway(s): p = 0.427, n = 61, deaths = 35, HR 0.76 (95% CI 0.39-1.49)
- untreated, 8 pathway(s): p = 0.197, n = 61, deaths = 35, HR 0.65 (95% CI 0.33-1.26)
- untreated, 9 pathway(s): p = 0.121, n = 61, deaths = 35, HR 0.59 (95% CI 0.30-1.16)
- untreated, 10 pathway(s): p = 0.364, n = 61, deaths = 35, HR 0.74 (95% CI 0.38-1.43)
