# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_READ_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.409, n = 59, deaths = 5, HR 2.12 (95% CI 0.34-13.17)
- treated, 2 pathway(s): p = 0.146, n = 59, deaths = 5, HR 4.62 (95% CI 0.49-43.22)
- treated, 3 pathway(s): p = 0.419, n = 59, deaths = 5, HR 2.09 (95% CI 0.34-13.03)
- treated, 4 pathway(s): p = 0.495, n = 59, deaths = 5, HR 1.88 (95% CI 0.30-11.89)
- treated, 5 pathway(s): p = 0.675, n = 59, deaths = 5, HR 0.68 (95% CI 0.11-4.12)
- treated, 6 pathway(s): p = 0.437, n = 59, deaths = 5, HR 0.50 (95% CI 0.08-2.99)
- treated, 7 pathway(s): p = 0.610, n = 59, deaths = 5, HR 0.63 (95% CI 0.10-3.81)
- treated, 8 pathway(s): p = 0.663, n = 59, deaths = 5, HR 0.67 (95% CI 0.11-4.06)
- treated, 9 pathway(s): p = 0.630, n = 59, deaths = 5, HR 0.64 (95% CI 0.11-3.90)
- treated, 10 pathway(s): p = 0.655, n = 59, deaths = 5, HR 0.67 (95% CI 0.11-4.02)
- untreated, 1 pathway(s): p = 0.024, n = 89, deaths = 18, HR 3.34 (95% CI 1.09-10.20)
- untreated, 2 pathway(s): p = 0.032, n = 89, deaths = 18, HR 3.17 (95% CI 1.04-9.68)
- untreated, 3 pathway(s): p = 0.185, n = 89, deaths = 18, HR 1.92 (95% CI 0.72-5.14)
- untreated, 4 pathway(s): p = 0.114, n = 89, deaths = 18, HR 2.17 (95% CI 0.81-5.81)
- untreated, 5 pathway(s): p = 0.154, n = 89, deaths = 18, HR 2.01 (95% CI 0.75-5.40)
- untreated, 6 pathway(s): p = 0.077, n = 89, deaths = 18, HR 2.46 (95% CI 0.87-6.93)
- untreated, 7 pathway(s): p = 0.245, n = 89, deaths = 18, HR 1.74 (95% CI 0.67-4.51)
- untreated, 8 pathway(s): p = 0.340, n = 89, deaths = 18, HR 1.58 (95% CI 0.61-4.08)
- untreated, 9 pathway(s): p = 0.146, n = 89, deaths = 18, HR 2.03 (95% CI 0.76-5.44)
- untreated, 10 pathway(s): p = 0.142, n = 89, deaths = 18, HR 2.05 (95% CI 0.77-5.48)
