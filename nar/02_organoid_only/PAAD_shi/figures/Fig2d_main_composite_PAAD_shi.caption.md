# Main figure: Pancreatic organoids (Shi 2022)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, FLUOROURACIL in PAAD (current GDC): 39 models, 38 proximal pathways; rank-1 pathway 'Inhibition of replication initiation of damaged dna by rb1 e2f1'; treated n = 37, deaths = 18, p = 0.958, HR 0.98 (95% CI 0.38-2.48); untreated n = 61, deaths = 35, p = 0.412, HR 0.75 (95% CI 0.39-1.48); Ridge leave-one-out rho = -0.01; shuffle share = 0.945.
- Row 2, GEMCITABINE in PAAD (current GDC): 39 models, 45 proximal pathways; rank-1 pathway 'G1 s specific transcription'; treated n = 104, deaths = 51, p = 0.023, HR 0.52 (95% CI 0.30-0.92); untreated n = 61, deaths = 35, p = 0.017, HR 0.45 (95% CI 0.23-0.88); Ridge leave-one-out rho = -0.08; shuffle share = 0.100.
- Row 3, IRINOTECAN in PAAD (current GDC): 39 models, 78 proximal pathways; rank-1 pathway 'Fanconi anemia pathway'; treated n = 13, deaths = 5, p = 0.302, HR 0.33 (95% CI 0.04-3.00); untreated n = 61, deaths = 35, p = 0.021, HR 2.22 (95% CI 1.11-4.43); Ridge leave-one-out rho = -0.27; shuffle share = 0.351.
- Row 4, OXALIPLATIN in PAAD (current GDC): 39 models, 46 proximal pathways; rank-1 pathway 'Intrinsic pathway'; treated n = 21, deaths = 11, p = 0.497, HR 0.65 (95% CI 0.18-2.30); untreated n = 61, deaths = 35, p = 0.242, HR 1.49 (95% CI 0.76-2.91); Ridge leave-one-out rho = -0.35; shuffle share = 0.277.
