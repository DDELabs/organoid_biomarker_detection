# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_authors_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.014, n = 114, deaths = 8, HR not estimable (0 deaths in one group)
- treated, 2 pathway(s): p = 0.882, n = 114, deaths = 8, HR 0.90 (95% CI 0.22-3.61)
- treated, 3 pathway(s): p = 0.373, n = 114, deaths = 8, HR 1.90 (95% CI 0.45-8.04)
- treated, 4 pathway(s): p = 0.314, n = 114, deaths = 8, HR 2.24 (95% CI 0.45-11.20)
- treated, 5 pathway(s): p = 0.129, n = 114, deaths = 8, HR 3.24 (95% CI 0.65-16.18)
- treated, 6 pathway(s): p = 0.253, n = 114, deaths = 8, HR 2.47 (95% CI 0.50-12.33)
- treated, 7 pathway(s): p = 0.076, n = 114, deaths = 8, HR 3.96 (95% CI 0.78-20.17)
- treated, 8 pathway(s): p = 0.216, n = 114, deaths = 8, HR 2.67 (95% CI 0.53-13.43)
- treated, 9 pathway(s): p = 0.314, n = 114, deaths = 8, HR 2.24 (95% CI 0.45-11.20)
- treated, 10 pathway(s): p = 0.314, n = 114, deaths = 8, HR 2.24 (95% CI 0.45-11.20)
- untreated, 1 pathway(s): p = 0.161, n = 298, deaths = 41, HR 0.64 (95% CI 0.34-1.20)
- untreated, 2 pathway(s): p = 0.656, n = 298, deaths = 41, HR 0.87 (95% CI 0.47-1.61)
- untreated, 3 pathway(s): p = 0.968, n = 298, deaths = 41, HR 0.99 (95% CI 0.53-1.83)
- untreated, 4 pathway(s): p = 0.393, n = 298, deaths = 41, HR 0.76 (95% CI 0.41-1.42)
- untreated, 5 pathway(s): p = 0.422, n = 298, deaths = 41, HR 0.78 (95% CI 0.42-1.45)
- untreated, 6 pathway(s): p = 0.844, n = 298, deaths = 41, HR 0.94 (95% CI 0.51-1.74)
- untreated, 7 pathway(s): p = 0.632, n = 298, deaths = 41, HR 0.86 (95% CI 0.46-1.59)
- untreated, 8 pathway(s): p = 0.619, n = 298, deaths = 41, HR 0.85 (95% CI 0.46-1.59)
- untreated, 9 pathway(s): p = 0.907, n = 298, deaths = 41, HR 0.96 (95% CI 0.52-1.78)
- untreated, 10 pathway(s): p = 0.960, n = 298, deaths = 41, HR 0.98 (95% CI 0.53-1.82)
