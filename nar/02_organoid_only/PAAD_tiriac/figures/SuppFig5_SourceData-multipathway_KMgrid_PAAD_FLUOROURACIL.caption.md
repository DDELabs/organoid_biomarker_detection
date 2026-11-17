# Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_PAAD_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.932, n = 37, deaths = 18, HR 0.96 (95% CI 0.38-2.44)
- treated, 2 pathway(s): p = 0.415, n = 37, deaths = 18, HR 0.67 (95% CI 0.26-1.75)
- treated, 3 pathway(s): p = 0.296, n = 37, deaths = 18, HR 0.59 (95% CI 0.22-1.59)
- treated, 4 pathway(s): p = 0.277, n = 37, deaths = 18, HR 0.58 (95% CI 0.22-1.56)
- treated, 5 pathway(s): p = 0.532, n = 37, deaths = 18, HR 0.74 (95% CI 0.28-1.92)
- treated, 6 pathway(s): p = 0.689, n = 37, deaths = 18, HR 0.82 (95% CI 0.32-2.14)
- treated, 7 pathway(s): p = 0.637, n = 37, deaths = 18, HR 0.79 (95% CI 0.30-2.07)
- treated, 8 pathway(s): p = 0.637, n = 37, deaths = 18, HR 0.79 (95% CI 0.30-2.07)
- treated, 9 pathway(s): p = 0.649, n = 37, deaths = 18, HR 0.80 (95% CI 0.31-2.09)
- treated, 10 pathway(s): p = 0.964, n = 37, deaths = 18, HR 0.98 (95% CI 0.38-2.50)
- untreated, 1 pathway(s): p = 0.281, n = 61, deaths = 35, HR 0.69 (95% CI 0.34-1.37)
- untreated, 2 pathway(s): p = 0.770, n = 61, deaths = 35, HR 0.90 (95% CI 0.46-1.76)
- untreated, 3 pathway(s): p = 0.026, n = 61, deaths = 35, HR 0.46 (95% CI 0.23-0.93)
- untreated, 4 pathway(s): p = 0.643, n = 61, deaths = 35, HR 0.85 (95% CI 0.44-1.67)
- untreated, 5 pathway(s): p = 0.715, n = 61, deaths = 35, HR 1.13 (95% CI 0.58-2.20)
- untreated, 6 pathway(s): p = 0.015, n = 61, deaths = 35, HR 2.27 (95% CI 1.15-4.46)
- untreated, 7 pathway(s): p = 0.008, n = 61, deaths = 35, HR 2.44 (95% CI 1.24-4.79)
- untreated, 8 pathway(s): p = 0.008, n = 61, deaths = 35, HR 2.44 (95% CI 1.24-4.79)
- untreated, 9 pathway(s): p = 0.023, n = 61, deaths = 35, HR 2.15 (95% CI 1.10-4.20)
- untreated, 10 pathway(s): p = 0.865, n = 61, deaths = 35, HR 1.06 (95% CI 0.54-2.06)
