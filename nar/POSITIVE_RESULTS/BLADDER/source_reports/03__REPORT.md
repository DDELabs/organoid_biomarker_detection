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
| BLCA_CISPLATIN | Lee 2018 organoids, cisplatin (9 lines from 7 patients; SCBO-3_2 and SCBO-11_2 are recurrence lines) | TCGA-BLCA, GSE31684, GSE48276, GSE70691 (survival); GSE169455, GSE52219 (response) |
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

Re-run (about 3 hours, 2 processes): `cd code; python3 run_joint.py --jobs 2; python3 align_figs.py;
python3 learning_curves.py --jobs 2; python3 run_joint.py --features shortlist --fold-align --jobs 2; python3 make_figures.py; python3 report_numbers.py`. Seed 20240611.

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
Convention (audit fix I05, same as every other folder): the model's benefit score is turned into a **risk score**
(risk = minus benefit). HR per SD of the risk score in held-out treated patients (Cox stratified by cohort);
**HR above 1 = expected direction** (higher risk score, shorter event-free survival). Concordance pairs are compared
only within a cohort. Old tables in this report used the benefit score (HR below 1 = expected); the new HR is 1/old
HR for the same fit. GSE39582 now uses the shared definition (audit fix I01, see "Patient definitions" below), so the
CRC_5FU rows use 729 patients instead of 792.

| Unit | Model | Survival n / events | C-index | HR per SD of risk score (95% CI), p | Response: n / responders, AUROC |
|---|---|---|---|---|---|
| CRC_5FU | organoid-only | 729 / 273 | 0.503 | 1.07 (0.94-1.20), 0.31 | 723 / 340, 0.490 |
|  | patient-only |  | 0.536 | 1.11 (0.99-1.24), 0.083 | 0.504 |
|  | joint |  | 0.525 | 1.11 (0.99-1.25), 0.078 | 0.501 |
| CRC_OXALIPLATIN | organoid-only | 255 / 88 | 0.453 | 0.97 (0.79-1.20), 0.78 | 419 / 199, 0.552 |
|  | patient-only |  | 0.449 | 0.80 (0.63-1.00), 0.048 | 0.521 |
|  | joint |  | 0.464 | 0.92 (0.75-1.13), 0.41 | 0.540 |
| CRC_IRINOTECAN | organoid-only | 109 / 96 | 0.465 | 1.02 (0.82-1.28), 0.84 | 170 / 90, 0.493 |
|  | patient-only |  | 0.459 | 0.84 (0.68-1.05), 0.13 | 0.489 |
|  | joint |  | 0.436 | 0.93 (0.75-1.15), 0.50 | 0.480 |
| BLCA_CISPLATIN | organoid-only | 157 / 99 | 0.507 | 1.07 (0.89-1.29), 0.45 | 229 / 92, 0.499 |
|  | patient-only |  | 0.560 | 1.45 (1.18-1.79), <0.001 | 0.580 |
|  | joint |  | 0.528 | 1.24 (1.03-1.50), 0.024 | 0.536 |
| BLCA_GEMCITABINE | organoid-only | 64 / 30 | 0.423 | 0.90 (0.63-1.30), 0.59 | 54 / 36, 0.356 |
|  | patient-only |  | 0.632 | 1.55 (1.11-2.15), 0.01 | 0.721 |
|  | joint |  | 0.594 | 1.53 (1.09-2.14), 0.014 | 0.642 |
| PAAD_GEMCITABINE | organoid-only | 92 / 57 | 0.371 | 0.72 (0.53-0.96), 0.027 | 59 / 29, 0.432 |
|  | patient-only |  | 0.577 | 1.09 (0.82-1.44), 0.56 | 0.520 |
|  | joint |  | 0.521 | 0.88 (0.66-1.17), 0.38 | 0.491 |
| LIHC_SORAFENIB | organoid-only | no public survival times |  |  | 67 / 21, 0.965 |
|  | patient-only |  |  |  | 0.947 |
|  | joint |  |  |  | 0.947 |
| OV_PACLITAXEL | organoid-only | 564 / 426 | 0.502 | 0.97 (0.88-1.06), 0.49 | 272 / 187, 0.492 |
|  | patient-only |  | 0.558 | 1.22 (1.10-1.34), <0.001 | 0.631 |
|  | joint |  | 0.553 | 1.14 (1.03-1.26), 0.0099 | 0.584 |
| OV_PLATINUM | lab-only (cell lines) | 683 / 523 | 0.520 | 1.06 (0.97-1.15), 0.20 | 347 / 236, 0.583 |
|  | patient-only |  | 0.550 | 1.19 (1.08-1.30), <0.001 | 0.610 |
|  | joint |  | 0.550 | 1.17 (1.07-1.29), <0.001 | 0.618 |
| LUNG_PLATINUM | lab-only (cell lines) | 277 / 97 | 0.482 | 0.91 (0.74-1.11), 0.35 | 126 / 99, 0.490 |
|  | patient-only |  | 0.552 | 1.13 (0.93-1.36), 0.22 | 0.594 |
|  | joint |  | 0.524 | 1.02 (0.84-1.23), 0.83 | 0.571 |
| SARC_GEMCITABINE | organoid-only | 16 / 9 | 0.312 | 0.62 (0.28-1.40), 0.25 | — |
|  | patient-only |  | 0.637 | 0.88 (0.46-1.70), 0.71 |  |
|  | joint |  | 0.450 | 0.61 (0.27-1.36), 0.23 |  |

What I think:
- The **organoid-only** score reaches held-out survival C-index 0.31-0.52. Its best results are colorectal oxaliplatin
  response (AUROC 0.552) and ovarian platinum response from cell lines (0.583); liver is discussed below. Lab-model
  drug response on its own does not carry over to patient survival in these cohorts.
- The **patient-only** model is the strongest single source in most cancers (bladder, ovarian). In colorectal 5-FU its
  HR is now 1.11 (0.99-1.24), p = 0.083 (before the fix: HR 0.89 in benefit scale, p = 0.026).
- The **joint** model sits between the two. Its pooled HR is above 1 with a CI excluding 1 in four units
  (BLCA_CISPLATIN 1.24, BLCA_GEMCITABINE 1.53, OV_PACLITAXEL 1.14, OV_PLATINUM 1.17). CRC_5FU is 1.11 (0.99-1.25),
  p = 0.078: before the fix it was counted as a fifth unit (p = 0.042). The joint model is above patient-only on
  concordance in colorectal oxaliplatin (response 0.540 vs 0.521; survival 0.464 vs 0.449) and ovarian platinum response
  (0.618 vs 0.610). In irinotecan it is now below patient-only (survival 0.436 vs 0.459; response 0.480 vs 0.489);
  before the FOLFIRINOX fix (I02) it was above (0.473 vs 0.425; 0.493 vs 0.450).
- Inner cross-validation chose the smallest lab-model weight on the grid (w = 0.25) in 56 of 61 joint fits
  (`tables/joint_weight_chosen.tsv`): the patients "ask" for the lab-model task to act as a light prior.
- All 673 pathways vs the shortlist (`tables/summary_all.tsv`): the pooled joint survival C-index was higher with all
  pathways in 6 of 10 units (e.g. CRC_5FU 0.546 vs 0.525; OV_PLATINUM 0.583 vs 0.550; BLCA_GEMCITABINE 0.689 vs
  0.594). The shortlist is easier to interpret; all pathways carry a little more prognostic information.

### Patient definitions used here (audit fixes I01, I02; `nar/15_fixes/definitions/GSE39582_treatment.tsv`)
- GSE39582 "5-FU-treated" = arm 5FU, FUFOL or FOLFOX, stage II-III only: 139 patients (138 with RFS). Patients whose
  chemotherapy type is not given (66) are no longer counted as 5-FU treated; FOLFIRI is excluded. Comparison group =
  "no chemotherapy", stage II-III: 264 (258 with RFS). Before: 200 treated / 288 comparison.
- FOLFIRINOX records (oxaliplatin and irinotecan together) are removed from the oxaliplatin and irinotecan units and
  from both arms of the FOLFOX-vs-FOLFIRI comparison in GSE72970 (4 patients; 36 vs 84 and 84 vs 36 instead of 40 vs 84
  and 88 vs 36). They stay in 5-FU-backbone analyses.

### Alignment inside each fold (audit fix I06; `tables/*__shortlist_foldalign__*.tsv`, `nar/15_fixes/tables/03_foldalign_comparison.tsv`)
The main runs align lab models to the tumours using the covariance of all patients of the unit, including the held-out
cohort (no outcome labels are used, but it is a form of transductive step). Seven units were re-run with the alignment
redone inside each leave-one-cohort-out fold from training cohorts only. Held-out C-index / AUROC change by less than
0.02 for almost every row (largest: OV_PACLITAXEL joint response 0.584 -> 0.643, organoid-only response
0.492 -> 0.531). The joint pooled HR for BLCA_CISPLATIN is 1.24 -> 1.27 (p 0.024 -> 0.013), OV_PACLITAXEL 1.14 -> 1.22
(p 0.010 -> <0.001), CRC_5FU 1.11 -> 1.11 (p 0.078 -> 0.078). What I think: the alignment is not what produces the
held-out signal. The learning-curve hyper-parameters are still fixed at the median of the main run, which were tuned
with those cohorts included (stated limit).

| Unit | Model | Outcome | C / AUROC transductive | C / AUROC fold-aligned | HR (risk) transductive, p | HR (risk) fold-aligned, p |
|---|---|---|---|---|---|---|
| BLCA_CISPLATIN | organoid-only | survival | 0.507 | 0.516 | 1.07, 0.45 | 1.10, 0.32 |
| BLCA_CISPLATIN | organoid-only | response | 0.499 | 0.509 |  |  |
| BLCA_CISPLATIN | patient-only | survival | 0.560 | 0.560 | 1.45, <0.001 | 1.45, <0.001 |
| BLCA_CISPLATIN | patient-only | response | 0.580 | 0.580 |  |  |
| BLCA_CISPLATIN | joint | survival | 0.528 | 0.540 | 1.24, 0.024 | 1.27, 0.013 |
| BLCA_CISPLATIN | joint | response | 0.536 | 0.542 |  |  |
| CRC_5FU | organoid-only | survival | 0.503 | 0.502 | 1.07, 0.31 | 1.06, 0.32 |
| CRC_5FU | organoid-only | response | 0.490 | 0.492 |  |  |
| CRC_5FU | patient-only | survival | 0.536 | 0.536 | 1.11, 0.083 | 1.11, 0.083 |
| CRC_5FU | patient-only | response | 0.504 | 0.504 |  |  |
| CRC_5FU | joint | survival | 0.525 | 0.525 | 1.11, 0.078 | 1.11, 0.078 |
| CRC_5FU | joint | response | 0.501 | 0.502 |  |  |
| LIHC_SORAFENIB | organoid-only | response | 0.965 | 0.965 |  |  |
| LIHC_SORAFENIB | patient-only | response | 0.947 | 0.947 |  |  |
| LIHC_SORAFENIB | joint | response | 0.947 | 0.948 |  |  |
| LUNG_PLATINUM | organoid-only | survival | 0.482 | 0.479 | 0.91, 0.35 | 0.96, 0.65 |
| LUNG_PLATINUM | organoid-only | response | 0.490 | 0.479 |  |  |
| LUNG_PLATINUM | patient-only | survival | 0.552 | 0.552 | 1.13, 0.22 | 1.13, 0.22 |
| LUNG_PLATINUM | patient-only | response | 0.594 | 0.594 |  |  |
| LUNG_PLATINUM | joint | survival | 0.524 | 0.521 | 1.02, 0.83 | 1.03, 0.74 |
| LUNG_PLATINUM | joint | response | 0.571 | 0.573 |  |  |
| OV_PACLITAXEL | organoid-only | survival | 0.502 | 0.520 | 0.97, 0.49 | 1.04, 0.46 |
| OV_PACLITAXEL | organoid-only | response | 0.492 | 0.531 |  |  |
| OV_PACLITAXEL | patient-only | survival | 0.558 | 0.558 | 1.22, <0.001 | 1.22, <0.001 |
| OV_PACLITAXEL | patient-only | response | 0.631 | 0.631 |  |  |
| OV_PACLITAXEL | joint | survival | 0.553 | 0.566 | 1.14, 0.0099 | 1.22, <0.001 |
| OV_PACLITAXEL | joint | response | 0.584 | 0.643 |  |  |
| PAAD_GEMCITABINE | organoid-only | survival | 0.371 | 0.371 | 0.72, 0.027 | 0.72, 0.027 |
| PAAD_GEMCITABINE | organoid-only | response | 0.432 | 0.432 |  |  |
| PAAD_GEMCITABINE | patient-only | survival | 0.577 | 0.577 | 1.09, 0.56 | 1.09, 0.56 |
| PAAD_GEMCITABINE | patient-only | response | 0.520 | 0.520 |  |  |
| PAAD_GEMCITABINE | joint | survival | 0.521 | 0.523 | 0.88, 0.38 | 0.89, 0.42 |
| PAAD_GEMCITABINE | joint | response | 0.491 | 0.494 |  |  |
| SARC_GEMCITABINE | organoid-only | survival | 0.312 | 0.312 | 0.62, 0.25 | 0.62, 0.25 |
| SARC_GEMCITABINE | patient-only | survival | 0.637 | 0.637 | 0.88, 0.71 | 0.88, 0.71 |
| SARC_GEMCITABINE | joint | survival | 0.450 | 0.438 | 0.61, 0.23 | 0.62, 0.24 |

### Drug vs comparison (interaction) tests, joint model, shortlist (`tables/interaction_all.tsv`)
Drug-specific direction: interaction HR above 1 (survival) or interaction OR above 1 (response), in the risk-score
convention. Old tables listed survival interactions in the benefit scale (HR below 1 = drug-specific); the new HR is
1/old HR (e.g. GSE48276: 0.19 -> 5.30; TCGA-PAAD: 2.13 -> 0.47). Response ORs are unchanged.

| Cohort | Unit | Outcome | Treated / comparison (events or responders) | Joint interaction (95% CI), p | Organoid/lab-only interaction, p |
|---|---|---|---|---|---|
| GSE39582 | CRC_5FU | survival (RFS) | 137 / 252 (47 / 61) | HR 1.19 (0.83-1.69), 0.34 | 1.16, 0.43 |
| GSE103479 | CRC_5FU | survival (PFS) | 66 / 84 (23 / 18) | HR 1.22 (0.65-2.29), 0.55 | 1.43, 0.25 |
| GSE14333 | CRC_5FU | survival (DFS) | 87 / 139 (27 / 23) | HR 0.98 (0.62-1.54), 0.92 | 1.06, 0.82 |
| TCGA-COADREAD | CRC_5FU | survival (PFS) | 140 / 369 (28 / 71) | HR 0.80 (0.51-1.25), 0.33 | 0.78, 0.29 |
| TCGA-COADREAD | CRC_OXALIPLATIN | survival (PFS) | 128 / 400 (38 / 93) | HR 0.93 (0.63-1.36), 0.70 | 0.98, 0.92 |
| GSE72970 | CRC_OXALIPLATIN | survival (PFS) | 36 / 84 (31 / 80) | HR 0.80 (0.53-1.19), 0.26 | 0.90, 0.61 |
| GSE72970 | CRC_OXALIPLATIN | response | 36 / 84 (21 / 38) | OR 0.62 (0.26-1.50), 0.29 | 0.70, 0.43 |
| GSE39582 | CRC_OXALIPLATIN | survival (RFS) | 11 / 252 (6 / 61) | HR 2.69 (0.58-12.53), 0.21 | 1.72, 0.52 |
| GSE104645 | CRC_OXALIPLATIN | response | 113 / 54 (66 / 38) | OR 0.97 (0.49-1.91), 0.92 | 1.02, 0.94 |
| GSE87211 | CRC_OXALIPLATIN | survival (DFS) | 80 / 95 (13 / 21) | HR 0.87 (0.44-1.71), 0.69 | 0.92, 0.80 |
| GSE87211 | CRC_OXALIPLATIN | response | 90 / 111 (13 / 22) | OR 2.24 (1.03-4.86), 0.041 | 2.31, 0.035 |
| GSE72970 | CRC_IRINOTECAN | survival (PFS) | 84 / 36 (80 / 31) | HR 1.25 (0.81-1.92), 0.32 | 1.43, 0.13 |
| GSE72970 | CRC_IRINOTECAN | response | 84 / 36 (38 / 21) | OR 1.30 (0.55-3.10), 0.55 | 1.52, 0.38 |
| GSE104645 | CRC_IRINOTECAN | response | 54 / 113 (38 / 66) | OR 0.89 (0.43-1.86), 0.76 | 1.18, 0.64 |
| TCGA-BLCA | BLCA_CISPLATIN | survival (PFS) | 68 / 322 (33 / 132) | HR 1.30 (0.91-1.87), 0.15 | 1.12, 0.57 |
| GSE31684 | BLCA_CISPLATIN | survival (RFS) | 35 / 58 (26 / 13) | HR 0.90 (0.45-1.77), 0.76 | 1.10, 0.79 |
| GSE48276 | BLCA_CISPLATIN | survival (OS) | 16 / 37 (8 / 20) | HR 5.30 (1.69-16.66), 0.0043 | 3.05, 0.018 |
| TCGA-BLCA | BLCA_GEMCITABINE | survival (PFS) | 64 / 322 (30 / 132) | HR 1.38 (0.98-1.94), 0.068 | 1.21, 0.32 |
| TCGA-PAAD | PAAD_GEMCITABINE | survival (PFS) | 92 / 76 (57 / 41) | HR 0.47 (0.31-0.72), <0.001 | 0.66, 0.065 |
| GSE109211 | LIHC_SORAFENIB | response | 67 / 72 (21 / 21) | OR 1.34 (0.15-11.94), 0.80 | 2.69, 0.33 |
| TCGA-OV | OV_PACLITAXEL | survival (PFS) | 332 / 57 (252 / 32) | HR 0.92 (0.59-1.43), 0.70 | 0.97, 0.90 |
| GSE9891 | OV_PACLITAXEL | survival (RFS) | 191 / 47 (138 / 41) | HR 1.09 (0.79-1.52), 0.60 | 1.07, 0.66 |
| GSE9891 | OV_PACLITAXEL | response | 183 / 46 (120 / 26) | OR 0.80 (0.41-1.55), 0.51 | 1.00, 0.99 |
| GSE63885 | OV_PACLITAXEL | survival (OS) | 41 / 34 (36 / 30) | HR 0.79 (0.44-1.40), 0.41 | 0.76, 0.30 |
| GSE63885 | OV_PACLITAXEL | response | 41 / 34 (38 / 27) | OR 0.33 (0.07-1.52), 0.15 | 0.36, 0.18 |
| GSE15622 | OV_PACLITAXEL | response | 20 / 15 (13 / 9) | OR 1.39 (0.32-6.07), 0.66 | 0.15, 0.062 |
| TCGA-OV | OV_PLATINUM | survival (PFS) | 370 / 46 (278 / 23) | HR 1.05 (0.68-1.64), 0.81 | 1.42, 0.17 |
| GSE15622 | OV_PLATINUM | response | 15 / 20 (9 / 13) | OR 3.17 (0.69-14.68), 0.14 | 5.91, 0.045 |
| TCGA-LUAD | LUNG_PLATINUM | survival (PFS) | 124 / 369 (56 / 143) | HR 0.84 (0.61-1.16), 0.29 | 1.01, 0.94 |
| TCGA-LUSC | LUNG_PLATINUM | survival (PFS) | 104 / 386 (32 / 110) | HR 1.01 (0.68-1.50), 0.94 | 0.85, 0.45 |
| GSE42127 | LUNG_PLATINUM | survival (OS) | 49 / 124 (9 / 52) | HR 0.80 (0.43-1.49), 0.49 | 0.69, 0.35 |
| TCGA-SARC | SARC_GEMCITABINE | survival (PFS) | 16 / 230 (9 / 123) | HR 0.67 (0.36-1.24), 0.20 | 1.05, 0.89 |

What I think: 32 joint-model interaction tests were run; 3 have p < 0.05 (about 1.6 expected by chance). Two point the
drug-specific way: GSE48276 MVAC bladder (HR 5.30, but only 8 treated events) and GSE87211 oxaliplatin response
(OR 2.24). The lab-only score also gives GSE15622 carboplatin vs paclitaxel OR 5.91 (35 patients). In all three the
organoid/cell-line-only score points the same way as the joint score, which is the most encouraging sign that lab
models add drug-specific information. TCGA-PAAD points the other way (HR 0.47): the score separates untreated
patients more than gemcitabine-treated ones, i.e. it behaves as a prognostic score there. After the I01/I02 fixes the GSE72970 oxaliplatin-vs-irinotecan survival test is HR 0.80 (0.53-1.19), p = 0.26
(36 vs 84 patients; before: 40 vs 84, HR 1.24 in the benefit scale, p = 0.29), so it points away from the drug-specific
direction in both versions. None of these would survive correction for 32 tests on its own; they are leads for a larger,
randomised cohort. KM curves in the main figures are now drawn per cohort (one endpoint each); the pooled number is the
stratified Cox HR (audit fix I13).

### Liver (STORM) note (`tables/storm_auroc_by_arm.tsv`)
The sorafenib-arm AUROC is 0.947-0.965, but the same scores give AUROC 0.894-0.957 in the **placebo** arm. What I
think: this is not sorafenib-specific. STORM's "response" is recurrence-based and the cohort has one dominant
expression axis (the earlier audit found its first principal component alone gives AUROC 0.96), so any score that loads
on that axis looks excellent in both arms. The interaction test (OR 1.34, p = 0.80) is the honest reading.

### Learning curves (`tables/learning_curves_mean.tsv`; panel c of each main figure)
Joint model, 25% -> 100% of lab models: CRC_5FU 0.514 -> 0.520; BLCA_CISPLATIN 0.557 -> 0.548; PAAD 0.514 -> 0.513;
OV_PACLITAXEL 0.551 -> 0.558; OV_PLATINUM 0.542 -> 0.550; LUNG 0.549 -> 0.550; LIHC 0.900 -> 0.945.
Patient-only, 10% -> 100% of training patients: CRC_5FU 0.508 -> 0.533; BLCA_GEMCITABINE 0.557 -> 0.668;
OV_PACLITAXEL 0.527 -> 0.560; OV_PLATINUM 0.530 -> 0.550; PAAD 0.497 -> 0.556.
(CRC_5FU values changed with the new GSE39582 definition: lab models 0.519 -> 0.526 before, patients 0.519 -> 0.543 before.)
What I think: adding patients raises held-out performance steadily; adding lab models changes it little, because the
joint fit keeps the lab task as a light prior. With only 10% of the patients the joint model is at or slightly above
patient-only in CRC_5FU (0.516 vs 0.508) and CRC_OXALIPLATIN (0.517 vs 0.509) and about equal in lung (0.514 vs
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
- Small lab sets: bladder cisplatin 9 organoid lines from 7 patients (audit I07; a one-line-per-patient sensitivity for the nar/02 bladder reproduction is in `nar/02_organoid_only/BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv`), ovarian 14, sarcoma 14.
- Single-cohort units (pancreas, bladder gemcitabine, liver, sarcoma) are tested on held-out patients of the same cohort, not a new cohort.
- The joint weight grid started at 0.25 and inner CV nearly always chose it; lighter weights were not tried.
- Fold-aligned runs (I06) exist for seven units only (CRC_5FU, BLCA_CISPLATIN, PAAD, LIHC, OV_PACLITAXEL, LUNG, SARC); the other units use the label-free transductive alignment.
