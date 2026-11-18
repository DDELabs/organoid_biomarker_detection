# Organoids + patients together: one model per cancer

Folder: `nar/03_organoid_plus_patient/`. Code in `code/`, tables in `tables/`, figures in `figures/` (main) and
`figures/supp/` (one per drug). Every number below was written by the code into a file in `tables/`.

## The question in one sentence
Does a model that learns drug response from lab models (organoids, PDX or cell lines) **and** from treated patients
at the same time predict outcome in patients it has never seen better than a model that learns from only one of them?

## What was done

1. **Data reused.** The combined training table (`robust/results/training_table`, 20,766 samples, 673 Reactome pathway
   rank scores) for patients; lab-model sets from `data/curated/preclinical/`. Two additions:
   - STORM liver trial (GSE109211, sorafenib vs placebo, 139 samples with expression) was not in the table; I scored it
     the same way (same rank score, same gene sets).
   - Ovarian organoids (Vias et al. 2023, 14 organoids) in `data/vias2023_raw/`. Both files have exactly the sha256
     values in `robust/data_builders/provenance/organoids_ovarian_vias2023.md`, i.e. the same data the organoid-only
     agent uses in `nar/02_organoid_only/`.
2. **11 cancer x drug units** (`code/units.py`), each with lab models screened with the drug and treated patients:

| Unit | Lab models (n with response) | Patient cohorts |
|---|---|---|
| CRC_5FU | van de Wetering organoids, 5-FU (19); Isella PDX in the alignment figure | GSE39582, GSE103479, GSE14333, TCGA-COAD+READ, GSE72970, GSE87211 (survival); GSE28702, GSE104645, GSE19860, GSE62080 (response) |
| CRC_OXALIPLATIN | van de Wetering, oxaliplatin | TCGA-COADREAD, GSE72970, GSE39582, GSE28702, GSE19860, GSE104645, GSE87211 |
| CRC_IRINOTECAN | van de Wetering, irinotecan | GSE72970, GSE39582, TCGA-COADREAD, GSE62080, GSE104645 |
| BLCA_CISPLATIN | Lee 2018 organoids, cisplatin (9) | TCGA-BLCA, GSE31684, GSE48276, GSE70691 (survival); GSE169455, GSE52219 (response) |
| BLCA_GEMCITABINE | Lee 2018, gemcitabine | TCGA-BLCA |
| PAAD_GEMCITABINE | Tiriac 2018 (36) + Shi 2022 (39) organoids | TCGA-PAAD |
| LIHC_SORAFENIB | Ji 2023 (61) + Broutier 2017 (5) organoids | GSE109211 STORM (response; sorafenib vs placebo) |
| OV_PACLITAXEL | Vias 2023 organoids, paclitaxel (14) | TCGA-OV, GSE9891, GSE63885 (survival); GSE51373, GSE15622 (response) |
| OV_PLATINUM | ovarian cell lines, carboplatin (CTRPv2 + PRISM); no platinum organoid screen exists | same ovarian cohorts |
| LUNG_PLATINUM | lung cell lines: carboplatin (CTRPv2, 138), cisplatin (PRISM, 93) | TCGA-LUAD, TCGA-LUSC, GSE42127 |
| SARC_GEMCITABINE | Al Shihabi 2024 sarcoma organoids (14) | TCGA-SARC (16 treated) |

   TCGA-LIHC has only 3 sorafenib-treated patients, so it was not used.

3. **Common scale (alignment).** Step 1: every pathway standardised inside each patient cohort and inside each lab set
   (mean 0, SD 1). Step 2 (Celligner-like): directions with much more spread in lab models than in tumours (missing
   immune/stroma, culture effects) are removed (contrastive PCA reused from `robust/obd/bridge.py`), then the lab-model
   covariance is matched to the tumour covariance (CORAL, shrunk towards identity).
4. **Three models on the same features.** Features = the drug's network-proximal pathway shortlist (STRING v12 lists from
   `robust/results/survlearn`; sorafenib from the Kong et al. 2020 table, 22 pathways). Everything was also refitted on
   all 673 pathways.
   - *Organoid-only*: ridge regression on lab-model drug sensitivity (rank-normalised within each set).
   - *Patient-only*: penalised Cox regression in treated patients, stratified by cohort (response cohorts enter as a
     logistic task).
   - *Joint*: one shared set of pathway weights fitted to both tasks at once
     (loss = w x lab-model loss + (1 - w) x patient loss + ridge penalty). w (0.25/0.5/0.75) and the penalty were chosen
     by inner cross-validation on the training patients only.
5. **Patients enter training only through cross-fitting.** With 2 or more usable cohorts: leave one cohort out, train
   on the rest, test on the cohort never seen. With one cohort (bladder gemcitabine, pancreas, liver, sarcoma): 5-fold
   cross-fitting by patient, 5 repeats. Untreated / placebo / other-drug patients are never used for training.
6. **Held-out tests:** C-index (survival) or AUROC (response); HR per SD with 95% CI; log-rank of median split;
   Kaplan-Meier curves; treatment x score interaction where a comparison group exists.
7. **Learning curves:** 25-100% of lab models (20 random draws) and 10-100% of training patients (10 draws),
   hyper-parameters fixed at the median of the main run.

Re-run (about 2 hours, 2 processes): `cd code; python3 run_joint.py --jobs 2; python3 align_figs.py;
python3 learning_curves.py --jobs 2; python3 make_figures.py; python3 report_numbers.py`. Seed 20240611.

## Results

### Alignment (`tables/alignment_metrics.tsv`)
"Model-vs-tumour AUROC" = how well a classifier tells lab models from tumours (1 = always; 0.5 or below = cannot).
"Neighbour share" = fraction of a lab model's 15 nearest neighbours that are tumours (perfect mixing in brackets).

| Cancer | AUROC raw -> standardised -> aligned | Neighbour share raw -> aligned (perfect) |
|---|---|---|
| colorectal | 1.00 -> 0.36 -> 0.35 | 0.00 -> 0.87 (0.89) |
| bladder | 1.00 -> 0.15 -> 0.12 | 0.02 -> 0.98 (0.97) |
| pancreas | 1.00 -> 0.26 -> 0.27 | 0.00 -> 0.54 (0.57) |
| liver | 1.00 -> 0.25 -> 0.22 | 0.01 -> 0.56 (0.67) |
| ovarian | 1.00 -> 0.25 -> 0.24 | 0.03 -> 0.86 (0.91) |
| lung | 1.00 -> 0.34 -> 0.32 | 0.00 -> 0.84 (0.82) |
| sarcoma | 1.00 -> 0.14 -> 0.11 | 0.05 -> 0.92 (0.94) |

What I think: before alignment, lab models and tumours are perfectly separable in every cancer. Standardising inside
each cohort does most of the work; the Celligner-like step brings neighbour mixing close to perfect. Lab models and
tumours are on a common scale before any model is trained.

### Pooled held-out results, network shortlist (`tables/summary_main.tsv`)
HR per SD of the benefit score in held-out treated patients (Cox stratified by cohort; HR < 1 = higher score, longer
event-free survival). Concordance pairs are compared only within a cohort.

| Unit | Model | Survival n / events | C-index | HR per SD (95% CI), p | Response: n / responders, AUROC |
|---|---|---|---|---|---|
| CRC_5FU | organoid-only | 792 / 301 | 0.511 | 0.93 (0.83-1.05), 0.23 | 723 / 340, 0.490 |
| | patient-only | | 0.546 | 0.89 (0.80-0.99), 0.026 | 0.508 |
| | joint | | 0.532 | 0.89 (0.80-1.00), 0.042 | 0.502 |
| CRC_OXALIPLATIN | organoid-only | 264 / 92 | 0.453 | 1.04 (0.85-1.28), 0.71 | 427 / 205, 0.551 |
| | patient-only | | 0.457 | 1.24 (0.99-1.54), 0.060 | 0.515 |
| | joint | | 0.468 | 1.09 (0.89-1.33), 0.43 | 0.540 |
| CRC_IRINOTECAN | organoid-only | 118 / 100 | 0.467 | 0.98 (0.79-1.23), 0.86 | 178 / 96, 0.491 |
| | patient-only | | 0.425 | 1.07 (0.87-1.32), 0.50 | 0.450 |
| | joint | | 0.473 | 0.95 (0.76-1.19), 0.64 | 0.493 |
| BLCA_CISPLATIN | organoid-only | 157 / 99 | 0.507 | 0.93 (0.77-1.12), 0.45 | 229 / 92, 0.499 |
| | patient-only | | 0.560 | 0.69 (0.56-0.85), 0.001 | 0.580 |
| | joint | | 0.528 | 0.81 (0.67-0.97), 0.024 | 0.536 |
| BLCA_GEMCITABINE | organoid-only | 64 / 30 | 0.423 | 1.11 (0.77-1.59), 0.59 | 54 / 36, 0.356 |
| | patient-only | | 0.632 | 0.65 (0.46-0.90), 0.010 | 0.721 |
| | joint | | 0.594 | 0.65 (0.47-0.92), 0.014 | 0.642 |
| PAAD_GEMCITABINE | organoid-only | 92 / 57 | 0.371 | 1.39 (1.04-1.87), 0.027 | 59 / 29, 0.432 |
| | patient-only | | 0.577 | 0.92 (0.70-1.22), 0.56 | 0.520 |
| | joint | | 0.521 | 1.14 (0.85-1.52), 0.38 | 0.491 |
| LIHC_SORAFENIB | organoid-only | no public survival times | | | 67 / 21, 0.965 |
| | patient-only | | | | 0.947 |
| | joint | | | | 0.947 |
| OV_PACLITAXEL | organoid-only | 564 / 426 | 0.502 | 1.04 (0.94-1.14), 0.49 | 272 / 187, 0.492 |
| | patient-only | | 0.558 | 0.82 (0.74-0.91), <0.001 | 0.631 |
| | joint | | 0.553 | 0.88 (0.80-0.97), 0.010 | 0.584 |
| OV_PLATINUM | lab-only (cell lines) | 683 / 523 | 0.520 | 0.95 (0.87-1.03), 0.20 | 347 / 236, 0.583 |
| | patient-only | | 0.550 | 0.84 (0.77-0.92), <0.001 | 0.610 |
| | joint | | 0.550 | 0.85 (0.78-0.93), 0.001 | 0.618 |
| LUNG_PLATINUM | lab-only (cell lines) | 277 / 97 | 0.482 | 1.10 (0.90-1.34), 0.35 | 126 / 99, 0.490 |
| | patient-only | | 0.552 | 0.89 (0.73-1.07), 0.22 | 0.594 |
| | joint | | 0.524 | 0.98 (0.81-1.18), 0.84 | 0.571 |
| SARC_GEMCITABINE | organoid-only | 16 / 9 | 0.312 | 1.61 (0.71-3.62), 0.25 | — |
| | patient-only | | 0.638 | 1.13 (0.59-2.19), 0.71 | — |
| | joint | | 0.450 | 1.64 (0.74-3.67), 0.23 | — |

What I think:
- The **organoid-only** score reaches held-out survival C-index 0.31-0.52. Its best results are colorectal oxaliplatin
  response (AUROC 0.551) and ovarian platinum response from cell lines (0.583); liver is discussed below. Lab-model
  drug response on its own does not carry over to patient survival in these cohorts.
- The **patient-only** model is the strongest single source in most cancers (bladder, ovarian, colorectal 5-FU).
- The **joint** model sits between the two. Its pooled HR is below 1 with a CI excluding 1 in five units (CRC_5FU 0.89,
  BLCA_CISPLATIN 0.81, BLCA_GEMCITABINE 0.65, OV_PACLITAXEL 0.88, OV_PLATINUM 0.85). It is above patient-only on
  concordance in colorectal oxaliplatin (response 0.540 vs 0.515; survival 0.468 vs 0.457), irinotecan (survival 0.473
  vs 0.425; response 0.493 vs 0.450) and ovarian platinum response (0.618 vs 0.610).
- Inner cross-validation chose the smallest lab-model weight on the grid (w = 0.25) in 55 of 61 joint fits
  (`tables/joint_weight_chosen.tsv`): the patients "ask" for the lab-model task to act as a light prior.
- All 673 pathways vs the shortlist (`tables/summary_all.tsv`): the pooled joint survival C-index was higher with all
  pathways in 6 of 10 units (e.g. CRC_5FU 0.560 vs 0.532; OV_PLATINUM 0.583 vs 0.550; BLCA_GEMCITABINE 0.689 vs
  0.594). The shortlist is easier to interpret; all pathways carry a little more prognostic information.

### Drug vs comparison (interaction) tests, joint model, shortlist (`tables/interaction_all.tsv`)
Drug-specific direction: interaction HR below 1 (survival), interaction OR above 1 (response).

| Cohort (comparison) | Unit | Treated / comparison (events or responders) | Interaction (95% CI), p |
|---|---|---|---|
| GSE39582 (5-FU vs no chemo) | CRC_5FU | 200 / 288 (75 / 62) | HR 0.90 (0.66-1.24), 0.52 |
| GSE103479 (5-FU vs surgery) | CRC_5FU | 66 / 84 (23 / 18) | HR 0.81 (0.43-1.52), 0.51; organoid-only 0.70 (0.38-1.29), 0.25 |
| GSE14333 (chemo vs none) | CRC_5FU | 87 / 139 (27 / 23) | HR 1.03 (0.65-1.61), 0.92 |
| TCGA-COADREAD (chemo vs none) | CRC_5FU | 140 / 369 (28 / 71) | HR 1.26 (0.80-1.97), 0.32 |
| GSE72970 (FOLFOX vs FOLFIRI) | CRC_OXALIPLATIN | 40 / 84 (34 / 80) | PFS HR 1.24 (0.84-1.84), 0.29; response OR 0.62 (0.26-1.45), 0.27 |
| GSE72970 (FOLFIRI vs FOLFOX) | CRC_IRINOTECAN | 88 / 36 (83 / 31) | PFS HR 0.70 (0.44-1.11), 0.13; organoid-only 0.72 (0.46-1.13), 0.16 |
| GSE87211 (5-FU-RT + oxaliplatin vs without) | CRC_OXALIPLATIN | 90 / 111 (13 / 22) | response OR 2.24 (1.04-4.82), 0.040; organoid-only 2.31 (1.06-5.03), 0.035 |
| TCGA-BLCA (cisplatin vs none) | BLCA_CISPLATIN | 68 / 322 (33 / 132) | HR 0.77 (0.54-1.10), 0.15 |
| GSE48276 (MVAC vs none) | BLCA_CISPLATIN | 16 / 37 (8 / 20) | HR 0.19 (0.06-0.59), 0.004; organoid-only 0.33 (0.13-0.83), 0.018 |
| GSE31684 (chemo vs none) | BLCA_CISPLATIN | 35 / 58 (26 / 13) | HR 1.11 (0.56-2.20), 0.76 |
| TCGA-BLCA (gemcitabine vs none) | BLCA_GEMCITABINE | 64 / 322 (30 / 132) | HR 0.73 (0.52-1.02), 0.068 |
| TCGA-PAAD (gemcitabine vs none) | PAAD_GEMCITABINE | 92 / 76 (57 / 41) | HR 2.13 (1.39-3.27), 0.0006 |
| GSE109211 (sorafenib vs placebo, randomised) | LIHC_SORAFENIB | 67 / 72 (21 / 21) | OR 1.34 (0.15-11.9), 0.80 |
| GSE9891 (platinum + paclitaxel vs platinum) | OV_PACLITAXEL | 191 / 47 (138 / 41) | RFS HR 0.92 (0.66-1.27), 0.60 |
| GSE63885 (+paclitaxel vs +cyclophosphamide) | OV_PACLITAXEL | 41 / 34 (36 / 30) | OS HR 1.27 (0.72-2.27), 0.41 |
| GSE15622 (carboplatin vs paclitaxel, randomised) | OV_PLATINUM | 15 / 20 (9 / 13) | OR 3.17 (0.69-14.7), 0.14; cell-line-only OR 5.92 (1.04-33.6), 0.045 |
| GSE42127 (platinum vs none) | LUNG_PLATINUM | 49 / 124 (9 / 52) | HR 1.25 (0.67-2.31), 0.49 |
| TCGA-LUAD / TCGA-LUSC (platinum vs none) | LUNG_PLATINUM | 124 / 369; 104 / 386 | HR 1.19 (0.86-1.63), 0.29; 0.99 (0.67-1.46), 0.94 |
| TCGA-SARC (gemcitabine vs none) | SARC_GEMCITABINE | 16 / 230 (9 / 123) | HR 1.50 (0.81-2.78), 0.20 |

What I think: 33 joint-model interaction tests were run; 3 have p < 0.05 (about 1.7 expected by chance). Two point the
drug-specific way: GSE48276 MVAC bladder (HR 0.19, but only 8 treated events) and GSE87211 oxaliplatin response
(OR 2.24). The lab-only score also gives GSE15622 carboplatin vs paclitaxel OR 5.92 (35 patients). In all three the
organoid/cell-line-only score points the same way as the joint score, which is the most encouraging sign that lab
models add drug-specific information. TCGA-PAAD points the other way (HR 2.13): the score separates untreated
patients (HR 0.53 per SD) more than gemcitabine-treated ones, i.e. it behaves as a prognostic score there. None of these
would survive correction for 33 tests on its own; they are leads for a larger, randomised cohort.

### Liver (STORM) note (`tables/storm_auroc_by_arm.tsv`)
The sorafenib-arm AUROC is 0.947-0.965, but the same scores give AUROC 0.894-0.957 in the **placebo** arm. What I
think: this is not sorafenib-specific. STORM's "response" is recurrence-based and the cohort has one dominant
expression axis (the earlier audit found its first principal component alone gives AUROC 0.96), so any score that loads
on that axis looks excellent in both arms. The interaction test (OR 1.34, p = 0.80) is the honest reading.

### Learning curves (`tables/learning_curves_mean.tsv`; panel c of each main figure)
Joint model, 25% -> 100% of lab models: CRC_5FU 0.519 -> 0.526; BLCA_CISPLATIN 0.557 -> 0.548; PAAD 0.514 -> 0.513;
OV_PACLITAXEL 0.551 -> 0.558; OV_PLATINUM 0.542 -> 0.550; LUNG 0.549 -> 0.550; LIHC 0.900 -> 0.945.
Patient-only, 10% -> 100% of training patients: CRC_5FU 0.519 -> 0.543; BLCA_GEMCITABINE 0.557 -> 0.668;
OV_PACLITAXEL 0.527 -> 0.560; OV_PLATINUM 0.530 -> 0.550; PAAD 0.497 -> 0.556.
What I think: adding patients raises held-out performance steadily; adding lab models changes it little, because the
joint fit keeps the lab task as a light prior. With only 10% of the patients the joint model is at or slightly above
patient-only in CRC_5FU (0.522 vs 0.519) and CRC_OXALIPLATIN (0.512 vs 0.501) and about equal in lung (0.514 vs
0.517): small patient sets are where organoids help most.
Sarcoma: no learning curve could be drawn; each held-out fold has only 3 treated patients, below the minimum of 10.

## Figures (each with a plain-language `.caption.md`; PNG 300 dpi + PDF)
- Main, one per cancer: `figures/main_{colorectal,bladder,pancreas,liver,ovarian,lung,sarcoma}` — alignment
  before/after, learning curves, concordance per held-out cohort, HR forest, KM curves (joint and organoid-only),
  interaction forest.
- Alignment per cancer (PCA + UMAP, raw / standardised / aligned): `figures/alignment_<cancer>`.
- Cross-cancer forest plot: `figures/forest_cross_cancer`.
- Supplementary per drug (shortlist vs all pathways, KM in every held-out cohort for every model): `figures/supp/unit_<UNIT>`.

## Tables
`summary_main.tsv` (pooled), `summary_all.tsv` (every cohort), `interaction_all.tsv`, `alignment_metrics.tsv`,
`learning_curves.tsv`, `learning_curves_mean.tsv`, `joint_weight_chosen.tsv`, `storm_auroc_by_arm.tsv`,
`<unit>__<features>__{eval,interaction,folds,meta}.tsv`, `data_versions.tsv` (size and sha256 of every input file).

## Limits (one line each)
- Treatment was randomised only in STORM and GSE15622; "untreated" groups elsewhere differ in stage and age.
- TCGA "treated" means the drug appears in the record at any time; some patients started after progression.
- Lung and ovarian platinum use cell lines: no public lung or platinum-screened ovarian organoid panel with expression was available.
- Small lab sets: bladder cisplatin 9 organoids, ovarian 14, sarcoma 14.
- Single-cohort units (pancreas, bladder gemcitabine, liver, sarcoma) are tested on held-out patients of the same cohort, not a new cohort.
- The joint weight grid started at 0.25 and inner CV nearly always chose it; lighter weights were not tried.
