# Main figure: Pancreatic organoids (Tiriac 2018)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, FLUOROURACIL in PAAD (current GDC): 36 models, 38 proximal pathways; rank-1 pathway 'Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation'; treated n = 37, deaths = 18, p = 0.932, HR 0.96 (95% CI 0.38-2.44); untreated n = 61, deaths = 35, p = 0.281, HR 0.69 (95% CI 0.34-1.37); Ridge leave-one-out rho = 0.01; shuffle share = 0.873.
- Row 2, GEMCITABINE in PAAD (current GDC): 36 models, 45 proximal pathways; rank-1 pathway 'Homologous recombination repair of replication independent double strand breaks'; treated n = 104, deaths = 51, p = 0.223, HR 1.41 (95% CI 0.81-2.46); untreated n = 61, deaths = 35, p = 0.897, HR 1.05 (95% CI 0.54-2.05); Ridge leave-one-out rho = 0.37; shuffle share = 0.805.
- Row 3, OXALIPLATIN in PAAD (current GDC): 36 models, 47 proximal pathways; rank-1 pathway 'Synthesis of substrates in n glycan biosythesis'; treated n = 21, deaths = 11, p = 0.223, HR 0.44 (95% CI 0.11-1.71); untreated n = 61, deaths = 35, p = 0.072, HR 0.53 (95% CI 0.26-1.07); Ridge leave-one-out rho = 0.14; shuffle share = 0.135.
- Row 4, SN-38 in PAAD (current GDC): 36 models, 78 proximal pathways; rank-1 pathway 'Regulation of the fanconi anemia pathway'; treated n = 13, deaths = 5, p = 0.570, HR 0.59 (95% CI 0.10-3.63); untreated n = 61, deaths = 35, p = 0.186, HR 1.58 (95% CI 0.80-3.11); Ridge leave-one-out rho = 0.02; shuffle share = 0.626.
