# Main figure: Colorectal PDX (Isella 2017, cetuximab)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, CETUXIMAB in COAD (current GDC): 192 models, 41 proximal pathways; rank-1 pathway 'Signaling by notch3'; treated n = 15, deaths = 7, p = 0.648, HR 0.70 (95% CI 0.15-3.22); untreated n = 302, deaths = 68, p = 0.046, HR 1.63 (95% CI 1.01-2.65); Ridge leave-one-out rho = 0.50; shuffle share = 0.760.
