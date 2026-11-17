# Main figure: Liver cancer organoids (Ji 2023)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, SORAFENIB in LIHC (current GDC): 61 models, 22 proximal pathways; rank-1 pathway 'Signaling by activated point mutants of fgfr1'; treated n = 30, deaths = 15, p = 0.188, HR 0.48 (95% CI 0.15-1.47); untreated n = 328, deaths = 108, p = 0.129, HR 1.34 (95% CI 0.92-1.96); Ridge leave-one-out rho = -0.16; shuffle share = 0.583.
