# Main figure: Liver cancer organoids (Broutier 2017)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, SORAFENIB in LIHC (current GDC): 5 models, 15 proximal pathways; rank-1 pathway 'Reversible hydration of carbon dioxide'; treated n = 30, deaths = 15, p = 0.139, HR 0.43 (95% CI 0.14-1.35); untreated n = 328, deaths = 108, p = 0.374, HR 1.19 (95% CI 0.81-1.73); Ridge leave-one-out rho = -0.70; shuffle share = 0.521.
