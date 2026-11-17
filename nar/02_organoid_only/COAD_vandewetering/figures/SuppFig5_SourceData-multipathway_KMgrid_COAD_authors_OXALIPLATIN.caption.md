# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | OXALIPLATIN

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COAD_authors_OXALIPLATIN.tsv`.

- treated, 1 pathway(s): p = 0.844, n = 100, deaths = 7, HR 1.18 (95% CI 0.23-6.11)
- treated, 2 pathway(s): p = 0.790, n = 100, deaths = 7, HR 0.82 (95% CI 0.18-3.65)
- treated, 3 pathway(s): p = 0.469, n = 100, deaths = 7, HR 1.82 (95% CI 0.35-9.43)
- treated, 4 pathway(s): p = 0.491, n = 100, deaths = 7, HR 1.77 (95% CI 0.34-9.11)
- treated, 5 pathway(s): p = 0.281, n = 100, deaths = 7, HR 2.40 (95% CI 0.46-12.39)
- treated, 6 pathway(s): p = 0.521, n = 100, deaths = 7, HR 1.70 (95% CI 0.33-8.77)
- treated, 7 pathway(s): p = 0.531, n = 100, deaths = 7, HR 1.62 (95% CI 0.35-7.36)
- treated, 8 pathway(s): p = 0.221, n = 100, deaths = 7, HR 2.71 (95% CI 0.52-14.16)
- treated, 9 pathway(s): p = 0.714, n = 100, deaths = 7, HR 1.33 (95% CI 0.29-5.99)
- treated, 10 pathway(s): p = 0.322, n = 100, deaths = 7, HR 2.25 (95% CI 0.43-11.69)
- untreated, 1 pathway(s): p = 0.660, n = 298, deaths = 41, HR 0.87 (95% CI 0.47-1.61)
- untreated, 2 pathway(s): p = 0.504, n = 298, deaths = 41, HR 0.81 (95% CI 0.43-1.51)
- untreated, 3 pathway(s): p = 0.847, n = 298, deaths = 41, HR 1.06 (95% CI 0.57-1.97)
- untreated, 4 pathway(s): p = 0.720, n = 298, deaths = 41, HR 1.12 (95% CI 0.61-2.07)
- untreated, 5 pathway(s): p = 0.132, n = 298, deaths = 41, HR 1.61 (95% CI 0.86-3.01)
- untreated, 6 pathway(s): p = 0.723, n = 298, deaths = 41, HR 1.12 (95% CI 0.60-2.08)
- untreated, 7 pathway(s): p = 0.058, n = 298, deaths = 41, HR 1.82 (95% CI 0.97-3.40)
- untreated, 8 pathway(s): p = 0.102, n = 298, deaths = 41, HR 1.67 (95% CI 0.90-3.09)
- untreated, 9 pathway(s): p = 0.253, n = 298, deaths = 41, HR 1.43 (95% CI 0.77-2.64)
- untreated, 10 pathway(s): p = 0.336, n = 298, deaths = 41, HR 1.35 (95% CI 0.73-2.50)
