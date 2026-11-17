# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.440, n = 112, deaths = 21, HR 0.71 (95% CI 0.30-1.70)
- treated, 2 pathway(s): p = 0.134, n = 112, deaths = 21, HR 0.51 (95% CI 0.20-1.26)
- treated, 3 pathway(s): p = 0.016, n = 112, deaths = 21, HR 0.31 (95% CI 0.11-0.85)
- treated, 4 pathway(s): p = 0.001, n = 112, deaths = 21, HR 0.17 (95% CI 0.05-0.57)
- treated, 5 pathway(s): p = 0.003, n = 112, deaths = 21, HR 0.23 (95% CI 0.08-0.67)
- treated, 6 pathway(s): p = 0.006, n = 112, deaths = 21, HR 0.24 (95% CI 0.08-0.72)
- treated, 7 pathway(s): p = 0.012, n = 112, deaths = 21, HR 0.27 (95% CI 0.09-0.81)
- treated, 8 pathway(s): p = 0.014, n = 112, deaths = 21, HR 0.31 (95% CI 0.11-0.84)
- treated, 9 pathway(s): p = 0.072, n = 112, deaths = 21, HR 0.43 (95% CI 0.17-1.11)
- treated, 10 pathway(s): p = 0.672, n = 112, deaths = 21, HR 0.83 (95% CI 0.35-1.98)
- untreated, 1 pathway(s): p = 0.470, n = 302, deaths = 68, HR 1.19 (95% CI 0.74-1.93)
- untreated, 2 pathway(s): p = 0.920, n = 302, deaths = 68, HR 0.98 (95% CI 0.61-1.57)
- untreated, 3 pathway(s): p = 0.332, n = 302, deaths = 68, HR 0.79 (95% CI 0.49-1.28)
- untreated, 4 pathway(s): p = 0.737, n = 302, deaths = 68, HR 0.92 (95% CI 0.57-1.48)
- untreated, 5 pathway(s): p = 0.605, n = 302, deaths = 68, HR 0.88 (95% CI 0.55-1.42)
- untreated, 6 pathway(s): p = 0.912, n = 302, deaths = 68, HR 0.97 (95% CI 0.61-1.57)
- untreated, 7 pathway(s): p = 0.143, n = 302, deaths = 68, HR 0.70 (95% CI 0.43-1.13)
- untreated, 8 pathway(s): p = 0.108, n = 302, deaths = 68, HR 0.68 (95% CI 0.42-1.09)
- untreated, 9 pathway(s): p = 0.096, n = 302, deaths = 68, HR 0.67 (95% CI 0.41-1.08)
- untreated, 10 pathway(s): p = 0.091, n = 302, deaths = 68, HR 0.66 (95% CI 0.41-1.07)
