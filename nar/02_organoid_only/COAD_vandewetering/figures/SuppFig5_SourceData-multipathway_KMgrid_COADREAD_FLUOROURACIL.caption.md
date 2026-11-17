# Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | FLUOROURACIL

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).
SVR and LinearRegression results for every k are in the table `all_tests_COADREAD_FLUOROURACIL.tsv`.

- treated, 1 pathway(s): p = 0.603, n = 188, deaths = 30, HR 1.21 (95% CI 0.59-2.48)
- treated, 2 pathway(s): p = 0.794, n = 188, deaths = 30, HR 0.91 (95% CI 0.44-1.86)
- treated, 3 pathway(s): p = 0.505, n = 188, deaths = 30, HR 0.78 (95% CI 0.38-1.61)
- treated, 4 pathway(s): p = 0.696, n = 188, deaths = 30, HR 0.87 (95% CI 0.42-1.78)
- treated, 5 pathway(s): p = 0.654, n = 188, deaths = 30, HR 0.85 (95% CI 0.41-1.74)
- treated, 6 pathway(s): p = 0.838, n = 188, deaths = 30, HR 0.93 (95% CI 0.45-1.90)
- treated, 7 pathway(s): p = 0.821, n = 188, deaths = 30, HR 1.09 (95% CI 0.53-2.24)
- treated, 8 pathway(s): p = 0.701, n = 188, deaths = 30, HR 1.15 (95% CI 0.56-2.36)
- treated, 9 pathway(s): p = 0.702, n = 188, deaths = 30, HR 1.15 (95% CI 0.56-2.37)
- treated, 10 pathway(s): p = 0.666, n = 188, deaths = 30, HR 1.17 (95% CI 0.57-2.40)
- untreated, 1 pathway(s): p = 0.104, n = 391, deaths = 86, HR 1.42 (95% CI 0.93-2.19)
- untreated, 2 pathway(s): p = 0.022, n = 391, deaths = 86, HR 1.65 (95% CI 1.07-2.54)
- untreated, 3 pathway(s): p = 0.315, n = 391, deaths = 86, HR 1.24 (95% CI 0.81-1.90)
- untreated, 4 pathway(s): p = 0.323, n = 391, deaths = 86, HR 1.24 (95% CI 0.81-1.90)
- untreated, 5 pathway(s): p = 0.063, n = 391, deaths = 86, HR 1.49 (95% CI 0.97-2.29)
- untreated, 6 pathway(s): p = 0.138, n = 391, deaths = 86, HR 1.38 (95% CI 0.90-2.11)
- untreated, 7 pathway(s): p = 0.585, n = 391, deaths = 86, HR 1.12 (95% CI 0.74-1.72)
- untreated, 8 pathway(s): p = 0.385, n = 391, deaths = 86, HR 1.21 (95% CI 0.79-1.84)
- untreated, 9 pathway(s): p = 0.368, n = 391, deaths = 86, HR 1.21 (95% CI 0.79-1.85)
- untreated, 10 pathway(s): p = 0.279, n = 391, deaths = 86, HR 1.26 (95% CI 0.83-1.93)
