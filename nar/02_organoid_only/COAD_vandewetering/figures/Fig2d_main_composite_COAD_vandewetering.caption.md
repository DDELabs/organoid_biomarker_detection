# Main figure: Colorectal organoids (van de Wetering 2015)

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= -1.2816); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

- Row 1, FLUOROURACIL in COAD (authors' 2019 files): 19 models, 37 proximal pathways; rank-1 pathway 'Activation of bh3 only proteins'; treated n = 114, deaths = 8, p = 0.014, HR not estimable (0 deaths in one group); untreated n = 298, deaths = 41, p = 0.161, HR 0.64 (95% CI 0.34-1.20); Ridge leave-one-out rho = 0.44; shuffle share = 0.034.
- Row 2, IRINOTECAN in COAD (authors' 2019 files): 19 models, 76 proximal pathways; rank-1 pathway 'Hyaluronan uptake and degradation'; treated n = 22, deaths = 3, p = 0.182, HR not estimable (0 deaths in one group); untreated n = 298, deaths = 41, p = 0.718, HR 0.89 (95% CI 0.48-1.67); Ridge leave-one-out rho = -0.21; shuffle share = 0.317.
- Row 3, OXALIPLATIN in COAD (authors' 2019 files): 19 models, 43 proximal pathways; rank-1 pathway 'Bile acid and bile salt metabolism'; treated n = 100, deaths = 7, p = 0.844, HR 1.18 (95% CI 0.23-6.11); untreated n = 298, deaths = 41, p = 0.660, HR 0.87 (95% CI 0.47-1.61); Ridge leave-one-out rho = 0.15; shuffle share = 0.946.
- Row 4, CETUXIMAB in COAD (current GDC): 19 models, 32 proximal pathways; rank-1 pathway 'Gastrin creb signalling pathway via pkc and mapk'; treated n = 15, deaths = 7, p = 0.329, HR 0.44 (95% CI 0.08-2.42); untreated n = 302, deaths = 68, p = 0.193, HR 1.37 (95% CI 0.85-2.21); Ridge leave-one-out rho = -0.40; shuffle share = 0.539.
- Row 5, FLUOROURACIL in COAD (current GDC): 19 models, 30 proximal pathways; rank-1 pathway 'Purine ribonucleoside monophosphate biosynthesis'; treated n = 129, deaths = 25, p = 0.866, HR 1.07 (95% CI 0.49-2.35); untreated n = 302, deaths = 68, p = 0.419, HR 1.22 (95% CI 0.75-1.96); Ridge leave-one-out rho = 0.58; shuffle share = 0.941.
- Row 6, IRINOTECAN in COAD (current GDC): 19 models, 66 proximal pathways; rank-1 pathway 'Hyaluronan uptake and degradation'; treated n = 32, deaths = 15, p = 0.257, HR 0.55 (95% CI 0.19-1.56); untreated n = 302, deaths = 68, p = 0.043, HR 0.60 (95% CI 0.36-0.99); Ridge leave-one-out rho = -0.22; shuffle share = 0.238.
- Row 7, OXALIPLATIN in COAD (current GDC): 19 models, 36 proximal pathways; rank-1 pathway 'Metabolism of polyamines'; treated n = 112, deaths = 21, p = 0.440, HR 0.71 (95% CI 0.30-1.70); untreated n = 302, deaths = 68, p = 0.470, HR 1.19 (95% CI 0.74-1.93); Ridge leave-one-out rho = 0.09; shuffle share = 0.448.
