# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.866, n = 129, deaths = 25, HR 1.07 (95% CI 0.49-2.35)
- treated, 2 pathway(s): p = 0.426, n = 129, deaths = 25, HR 0.72 (95% CI 0.33-1.61)
- treated, 3 pathway(s): p = 0.429, n = 129, deaths = 25, HR 0.73 (95% CI 0.33-1.61)
- treated, 4 pathway(s): p = 0.346, n = 129, deaths = 25, HR 0.68 (95% CI 0.31-1.52)
- treated, 5 pathway(s): p = 0.728, n = 129, deaths = 25, HR 0.87 (95% CI 0.40-1.91)
- treated, 6 pathway(s): p = 0.753, n = 129, deaths = 25, HR 1.13 (95% CI 0.52-2.49)
- treated, 7 pathway(s): p = 0.688, n = 129, deaths = 25, HR 1.18 (95% CI 0.53-2.61)
- treated, 8 pathway(s): p = 0.386, n = 129, deaths = 25, HR 1.42 (95% CI 0.64-3.14)
- treated, 9 pathway(s): p = 0.494, n = 129, deaths = 25, HR 1.32 (95% CI 0.60-2.90)
- treated, 10 pathway(s): p = 0.407, n = 129, deaths = 25, HR 1.40 (95% CI 0.63-3.08)
- untreated, 1 pathway(s): p = 0.419, n = 302, deaths = 68, HR 1.22 (95% CI 0.75-1.96)
- untreated, 2 pathway(s): p = 0.061, n = 302, deaths = 68, HR 1.58 (95% CI 0.97-2.57)
- untreated, 3 pathway(s): p = 0.442, n = 302, deaths = 68, HR 1.21 (95% CI 0.75-1.95)
- untreated, 4 pathway(s): p = 0.855, n = 302, deaths = 68, HR 1.05 (95% CI 0.65-1.69)
- untreated, 5 pathway(s): p = 0.459, n = 302, deaths = 68, HR 1.20 (95% CI 0.74-1.93)
- untreated, 6 pathway(s): p = 0.399, n = 302, deaths = 68, HR 1.23 (95% CI 0.76-1.98)
- untreated, 7 pathway(s): p = 0.810, n = 302, deaths = 68, HR 1.06 (95% CI 0.66-1.71)
- untreated, 8 pathway(s): p = 0.767, n = 302, deaths = 68, HR 1.07 (95% CI 0.67-1.73)
- untreated, 9 pathway(s): p = 0.967, n = 302, deaths = 68, HR 1.01 (95% CI 0.63-1.62)
- untreated, 10 pathway(s): p = 0.567, n = 302, deaths = 68, HR 1.15 (95% CI 0.71-1.85)
