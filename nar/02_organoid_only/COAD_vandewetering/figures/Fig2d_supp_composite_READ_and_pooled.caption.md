# Main figure: Colorectal organoids (van de Wetering 2015)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, FLUOROURACIL in READ (current GDC): 19 models, 30 proximal pathways; rank-1 pathway 'Purine ribonucleoside monophosphate biosynthesis'; treated n = 59, deaths = 5, p = 0.409, HR 2.12 (95% CI 0.34-13.17); untreated n = 89, deaths = 18, p = 0.024, HR 3.34 (95% CI 1.09-10.20); Ridge leave-one-out rho = 0.58; shuffle share = 0.219.
- Row 2, OXALIPLATIN in READ (current GDC): 19 models, 36 proximal pathways; rank-1 pathway 'Metabolism of polyamines'; treated n = 42, deaths = 5, p = 0.093, HR 5.37 (95% CI 0.60-48.34); untreated n = 89, deaths = 18, p = 0.326, HR 0.62 (95% CI 0.23-1.65); Ridge leave-one-out rho = 0.09; shuffle share = 0.126.
- Row 3, FLUOROURACIL in COAD+READ pooled (current GDC): 19 models, 30 proximal pathways; rank-1 pathway 'Purine ribonucleoside monophosphate biosynthesis'; treated n = 188, deaths = 30, p = 0.603, HR 1.21 (95% CI 0.59-2.48); untreated n = 391, deaths = 86, p = 0.104, HR 1.42 (95% CI 0.93-2.19); Ridge leave-one-out rho = 0.58; shuffle share = 0.775.
- Row 4, IRINOTECAN in COAD+READ pooled (current GDC): 19 models, 66 proximal pathways; rank-1 pathway 'Hyaluronan uptake and degradation'; treated n = 39, deaths = 17, p = 0.392, HR 0.65 (95% CI 0.25-1.74); untreated n = 391, deaths = 86, p = 0.028, HR 0.61 (95% CI 0.39-0.95); Ridge leave-one-out rho = -0.22; shuffle share = 0.703.
- Row 5, OXALIPLATIN in COAD+READ pooled (current GDC): 19 models, 36 proximal pathways; rank-1 pathway 'Metabolism of polyamines'; treated n = 154, deaths = 26, p = 0.695, HR 0.86 (95% CI 0.40-1.86); untreated n = 391, deaths = 86, p = 0.975, HR 0.99 (95% CI 0.65-1.52); Ridge leave-one-out rho = 0.09; shuffle share = 0.822.
