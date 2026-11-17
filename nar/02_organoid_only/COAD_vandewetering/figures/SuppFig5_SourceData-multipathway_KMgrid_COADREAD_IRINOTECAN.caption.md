# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | IRINOTECAN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COADREAD_IRINOTECAN.tsv`.

- treated, 1 pathway(s): p = 0.392, n = 39, deaths = 17, HR 0.65 (95% CI 0.25-1.74)
- treated, 2 pathway(s): p = 0.530, n = 39, deaths = 17, HR 1.37 (95% CI 0.51-3.70)
- treated, 3 pathway(s): p = 0.906, n = 39, deaths = 17, HR 1.06 (95% CI 0.40-2.85)
- treated, 4 pathway(s): p = 0.401, n = 39, deaths = 17, HR 1.53 (95% CI 0.56-4.18)
- treated, 5 pathway(s): p = 0.466, n = 39, deaths = 17, HR 1.45 (95% CI 0.53-3.98)
- treated, 6 pathway(s): p = 0.466, n = 39, deaths = 17, HR 1.45 (95% CI 0.53-3.98)
- treated, 7 pathway(s): p = 0.877, n = 39, deaths = 17, HR 1.08 (95% CI 0.41-2.82)
- treated, 8 pathway(s): p = 0.466, n = 39, deaths = 17, HR 1.45 (95% CI 0.53-3.98)
- treated, 9 pathway(s): p = 0.466, n = 39, deaths = 17, HR 1.45 (95% CI 0.53-3.98)
- treated, 10 pathway(s): p = 0.978, n = 39, deaths = 17, HR 1.01 (95% CI 0.39-2.65)
- untreated, 1 pathway(s): p = 0.028, n = 391, deaths = 86, HR 0.61 (95% CI 0.39-0.95)
- untreated, 2 pathway(s): p = 0.962, n = 391, deaths = 86, HR 0.99 (95% CI 0.65-1.51)
- untreated, 3 pathway(s): p = 0.818, n = 391, deaths = 86, HR 1.05 (95% CI 0.69-1.61)
- untreated, 4 pathway(s): p = 0.807, n = 391, deaths = 86, HR 0.95 (95% CI 0.62-1.46)
- untreated, 5 pathway(s): p = 0.841, n = 391, deaths = 86, HR 0.96 (95% CI 0.63-1.47)
- untreated, 6 pathway(s): p = 0.649, n = 391, deaths = 86, HR 1.10 (95% CI 0.72-1.69)
- untreated, 7 pathway(s): p = 0.714, n = 391, deaths = 86, HR 1.08 (95% CI 0.71-1.66)
- untreated, 8 pathway(s): p = 0.955, n = 391, deaths = 86, HR 1.01 (95% CI 0.66-1.55)
- untreated, 9 pathway(s): p = 0.848, n = 391, deaths = 86, HR 0.96 (95% CI 0.63-1.47)
- untreated, 10 pathway(s): p = 0.889, n = 391, deaths = 86, HR 0.97 (95% CI 0.64-1.49)
