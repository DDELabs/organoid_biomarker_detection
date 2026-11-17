# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_READ_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.093, n = 42, deaths = 5, HR 5.37 (95% CI 0.60-48.34)
- treated, 2 pathway(s): p = 0.494, n = 42, deaths = 5, HR 1.85 (95% CI 0.31-11.12)
- treated, 3 pathway(s): p = 0.571, n = 42, deaths = 5, HR 0.60 (95% CI 0.10-3.61)
- treated, 4 pathway(s): p = 0.124, n = 42, deaths = 5, HR 0.21 (95% CI 0.02-1.88)
- treated, 5 pathway(s): p = 0.550, n = 42, deaths = 5, HR 0.58 (95% CI 0.10-3.51)
- treated, 6 pathway(s): p = 0.329, n = 42, deaths = 5, HR 0.35 (95% CI 0.04-3.17)
- treated, 7 pathway(s): p = 0.161, n = 42, deaths = 5, HR 0.23 (95% CI 0.03-2.12)
- treated, 8 pathway(s): p = 0.267, n = 42, deaths = 5, HR 0.30 (95% CI 0.03-2.79)
- treated, 9 pathway(s): p = 0.334, n = 42, deaths = 5, HR 0.35 (95% CI 0.04-3.20)
- treated, 10 pathway(s): p = 0.097, n = 42, deaths = 5, HR 0.17 (95% CI 0.02-1.69)
- untreated, 1 pathway(s): p = 0.326, n = 89, deaths = 18, HR 0.62 (95% CI 0.23-1.65)
- untreated, 2 pathway(s): p = 0.262, n = 89, deaths = 18, HR 0.59 (95% CI 0.23-1.51)
- untreated, 3 pathway(s): p = 0.106, n = 89, deaths = 18, HR 0.46 (95% CI 0.18-1.21)
- untreated, 4 pathway(s): p = 0.171, n = 89, deaths = 18, HR 0.52 (95% CI 0.20-1.35)
- untreated, 5 pathway(s): p = 0.375, n = 89, deaths = 18, HR 0.64 (95% CI 0.23-1.76)
- untreated, 6 pathway(s): p = 0.193, n = 89, deaths = 18, HR 0.50 (95% CI 0.17-1.45)
- untreated, 7 pathway(s): p = 0.167, n = 89, deaths = 18, HR 0.48 (95% CI 0.17-1.39)
- untreated, 8 pathway(s): p = 0.118, n = 89, deaths = 18, HR 0.45 (95% CI 0.16-1.26)
- untreated, 9 pathway(s): p = 0.192, n = 89, deaths = 18, HR 0.53 (95% CI 0.20-1.41)
- untreated, 10 pathway(s): p = 0.263, n = 89, deaths = 18, HR 0.57 (95% CI 0.21-1.54)
