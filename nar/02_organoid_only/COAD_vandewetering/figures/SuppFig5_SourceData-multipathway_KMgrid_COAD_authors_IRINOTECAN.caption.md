# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | IRINOTECAN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_authors_IRINOTECAN.tsv`.

- treated, 1 pathway(s): p = 0.182, n = 22, deaths = 3, HR not estimable (0 deaths in one group)
- treated, 2 pathway(s): p = 0.182, n = 22, deaths = 3, HR not estimable (0 deaths in one group)
- treated, 3 pathway(s): p = 0.083, n = 22, deaths = 3, HR not estimable (0 deaths in one group)
- treated, 4 pathway(s): p = 0.259, n = 22, deaths = 3, HR 0.27 (95% CI 0.02-3.08)
- treated, 5 pathway(s): p = 0.083, n = 22, deaths = 3, HR not estimable (0 deaths in one group)
- treated, 6 pathway(s): p = 0.433, n = 22, deaths = 3, HR 0.39 (95% CI 0.03-4.44)
- treated, 7 pathway(s): p = 0.083, n = 22, deaths = 3, HR not estimable (0 deaths in one group)
- treated, 8 pathway(s): p = 0.433, n = 22, deaths = 3, HR 0.39 (95% CI 0.03-4.44)
- treated, 9 pathway(s): p = 0.433, n = 22, deaths = 3, HR 0.39 (95% CI 0.03-4.44)
- treated, 10 pathway(s): p = 0.433, n = 22, deaths = 3, HR 0.39 (95% CI 0.03-4.44)
- untreated, 1 pathway(s): p = 0.718, n = 298, deaths = 41, HR 0.89 (95% CI 0.48-1.67)
- untreated, 2 pathway(s): p = 0.926, n = 298, deaths = 41, HR 1.03 (95% CI 0.56-1.90)
- untreated, 3 pathway(s): p = 0.987, n = 298, deaths = 41, HR 1.00 (95% CI 0.54-1.86)
- untreated, 4 pathway(s): p = 0.606, n = 298, deaths = 41, HR 1.17 (95% CI 0.64-2.17)
- untreated, 5 pathway(s): p = 0.843, n = 298, deaths = 41, HR 1.07 (95% CI 0.57-2.00)
- untreated, 6 pathway(s): p = 0.523, n = 298, deaths = 41, HR 0.81 (95% CI 0.43-1.54)
- untreated, 7 pathway(s): p = 0.878, n = 298, deaths = 41, HR 1.05 (95% CI 0.57-1.94)
- untreated, 8 pathway(s): p = 0.882, n = 298, deaths = 41, HR 1.05 (95% CI 0.56-1.97)
- untreated, 9 pathway(s): p = 0.528, n = 298, deaths = 41, HR 1.22 (95% CI 0.65-2.29)
- untreated, 10 pathway(s): p = 0.773, n = 298, deaths = 41, HR 0.91 (95% CI 0.48-1.71)
