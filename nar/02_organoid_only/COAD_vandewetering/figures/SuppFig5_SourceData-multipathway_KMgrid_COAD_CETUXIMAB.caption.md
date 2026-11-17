# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | CETUXIMAB

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_CETUXIMAB.tsv`.

- treated, 1 pathway(s): p = 0.329, n = 15, deaths = 7, HR 0.44 (95% CI 0.08-2.42)
- treated, 2 pathway(s): p = 0.226, n = 15, deaths = 7, HR 0.36 (95% CI 0.06-2.02)
- treated, 3 pathway(s): p = 0.226, n = 15, deaths = 7, HR 0.36 (95% CI 0.06-2.02)
- treated, 4 pathway(s): p = 0.139, n = 15, deaths = 7, HR 0.29 (95% CI 0.05-1.64)
- treated, 5 pathway(s): p = 0.226, n = 15, deaths = 7, HR 0.36 (95% CI 0.06-2.02)
- treated, 6 pathway(s): p = 0.789, n = 15, deaths = 7, HR 0.81 (95% CI 0.18-3.68)
- treated, 7 pathway(s): p = 0.789, n = 15, deaths = 7, HR 0.81 (95% CI 0.18-3.68)
- treated, 8 pathway(s): p = 0.139, n = 15, deaths = 7, HR 0.29 (95% CI 0.05-1.64)
- treated, 9 pathway(s): p = 0.139, n = 15, deaths = 7, HR 0.29 (95% CI 0.05-1.64)
- treated, 10 pathway(s): p = 0.018, n = 15, deaths = 7, HR not estimable (0 deaths in one group)
- untreated, 1 pathway(s): p = 0.193, n = 302, deaths = 68, HR 1.37 (95% CI 0.85-2.21)
- untreated, 2 pathway(s): p = 0.329, n = 302, deaths = 68, HR 1.27 (95% CI 0.79-2.05)
- untreated, 3 pathway(s): p = 0.744, n = 302, deaths = 68, HR 1.08 (95% CI 0.67-1.75)
- untreated, 4 pathway(s): p = 0.694, n = 302, deaths = 68, HR 1.10 (95% CI 0.68-1.77)
- untreated, 5 pathway(s): p = 0.738, n = 302, deaths = 68, HR 1.08 (95% CI 0.67-1.75)
- untreated, 6 pathway(s): p = 0.775, n = 302, deaths = 68, HR 0.93 (95% CI 0.58-1.51)
- untreated, 7 pathway(s): p = 0.718, n = 302, deaths = 68, HR 0.92 (95% CI 0.57-1.48)
- untreated, 8 pathway(s): p = 0.754, n = 302, deaths = 68, HR 1.08 (95% CI 0.67-1.74)
- untreated, 9 pathway(s): p = 0.881, n = 302, deaths = 68, HR 0.96 (95% CI 0.60-1.56)
- untreated, 10 pathway(s): p = 0.932, n = 302, deaths = 68, HR 1.02 (95% CI 0.63-1.65)
