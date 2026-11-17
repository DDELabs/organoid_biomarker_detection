# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COADREAD_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.695, n = 154, deaths = 26, HR 0.86 (95% CI 0.40-1.86)
- treated, 2 pathway(s): p = 0.049, n = 154, deaths = 26, HR 0.44 (95% CI 0.19-1.02)
- treated, 3 pathway(s): p = 0.033, n = 154, deaths = 26, HR 0.40 (95% CI 0.17-0.96)
- treated, 4 pathway(s): p = 0.001, n = 154, deaths = 26, HR 0.23 (95% CI 0.09-0.61)
- treated, 5 pathway(s): p = 0.014, n = 154, deaths = 26, HR 0.33 (95% CI 0.13-0.83)
- treated, 6 pathway(s): p = 0.004, n = 154, deaths = 26, HR 0.26 (95% CI 0.10-0.70)
- treated, 7 pathway(s): p = 0.005, n = 154, deaths = 26, HR 0.27 (95% CI 0.10-0.71)
- treated, 8 pathway(s): p = 0.0007, n = 154, deaths = 26, HR 0.19 (95% CI 0.07-0.56)
- treated, 9 pathway(s): p = 0.025, n = 154, deaths = 26, HR 0.38 (95% CI 0.16-0.92)
- treated, 10 pathway(s): p = 0.258, n = 154, deaths = 26, HR 0.64 (95% CI 0.29-1.40)
- untreated, 1 pathway(s): p = 0.975, n = 391, deaths = 86, HR 0.99 (95% CI 0.65-1.52)
- untreated, 2 pathway(s): p = 0.718, n = 391, deaths = 86, HR 0.93 (95% CI 0.61-1.41)
- untreated, 3 pathway(s): p = 0.171, n = 391, deaths = 86, HR 0.74 (95% CI 0.48-1.14)
- untreated, 4 pathway(s): p = 0.530, n = 391, deaths = 86, HR 0.87 (95% CI 0.57-1.33)
- untreated, 5 pathway(s): p = 0.174, n = 391, deaths = 86, HR 0.75 (95% CI 0.49-1.14)
- untreated, 6 pathway(s): p = 0.306, n = 391, deaths = 86, HR 0.80 (95% CI 0.52-1.23)
- untreated, 7 pathway(s): p = 0.035, n = 391, deaths = 86, HR 0.63 (95% CI 0.41-0.97)
- untreated, 8 pathway(s): p = 0.016, n = 391, deaths = 86, HR 0.59 (95% CI 0.38-0.91)
- untreated, 9 pathway(s): p = 0.024, n = 391, deaths = 86, HR 0.61 (95% CI 0.39-0.94)
- untreated, 10 pathway(s): p = 0.031, n = 391, deaths = 86, HR 0.62 (95% CI 0.40-0.96)
