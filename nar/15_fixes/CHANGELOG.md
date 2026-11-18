# CHANGELOG: audit fixes applied (nar/15_fixes)

Date: 2026-11-18. Source of issues: `nar/12_audit/issues.tsv` (24 issues). Backups of every file before editing:
`nar/15_fixes/backup/`. Nothing was deleted. Nothing was committed. `nar/07_B3_models/` and `nar/13_gaps_bio/` were not
touched (their items are in `PENDING_FOR_RUNNING_AGENTS.md`).

## Status by issue
| ID | Status | One line |
|---|---|---|
| I01 | fixed (03, 05); 09 sensitivity added | one GSE39582 definition: 139 treated / 264 untreated (stage II-III) |
| I02 | fixed (03, 05) | FOLFIRINOX (4 patients) out of both GSE72970 arms: 36 vs 84 / 84 vs 36 |
| I03 | fixed (02 text + new table) | best-of-30 sentence now says direction |
| I04 | partly fixed | 03 main figures now 7.5 x 9.6 in, 05 Fig 5 split in two, 10 SuppFig4 redrawn; 02 large figures NOT re-rendered |
| I05 | fixed (03) | all HR now "HR > 1 = expected" |
| I06 | fixed (7 of 11 units re-run) | alignment inside each fold; changes small |
| I07 | fixed as sensitivity + text | 9 lines = 7 patients; one-line-per-patient run; SMBO-3 identity still to be checked in curation |
| I08 | fixed (text) | binomial p flagged as assuming independence; dependence-aware version is G9 (14, pending) |
| I09 | fixed | EPV flag + reduced model: 20/22 and 48/54 keep direction |
| I10 | fixed | navitoclax sentence corrected |
| I11 | fixed | docs/PROJECT_SUMMARY.md and MASTER_NUMBERS |
| I12 | pending | needs a Broutier bundle in 07 (not touched) and re-run of 10 c01 |
| I13 | fixed (03) | KM drawn per cohort only |
| I14 | partly fixed | difference stated in 02 text and MASTER_NUMBERS; 02 not re-parsed |
| I15 | fixed | "top-rank bootstrap (paper statistic)" vs "rank-1 p bootstrap" labels |
| I16 | fixed | wording in 09 |
| I17 | fixed (text) | 79/288 vs 77/294 stated |
| I18 | fixed | template wording |
| I19 | fixed | 09 sensitivity table |
| I20 | fixed | sorafenib IPW p = NA |
| I21 | pending | 07 |
| I22 | partly | 02 GI sha256 table written by the audit; 07 data_versions pending |
| I23 | partly | 08 Fig3 duplicate letter fixed, 10 SuppFig4 redrawn, A3 labels edited; font in 05 B1_Fig7, Fig1c row label, at-risk tables not all done |
| I24 | fixed (text) | samples and time points stated in 02 and 08 |

## Files changed and numbers that changed

### nar/03_organoid_plus_patient (I01, I02, I05, I06, I13, I04)
Files: `code/units.py`, `code/run_joint.py`, `code/make_figures.py` (edited earlier); re-run outputs `cache/CRC_*`,
`tables/CRC_*`, `tables/summary_main.tsv`, `summary_all.tsv`, `interaction_all.tsv`, `learning_curves*.tsv`,
`joint_weight_chosen.tsv`, new `tables/*__shortlist_foldalign__*` (CRC_5FU, BLCA_CISPLATIN, PAAD, LIHC, OV_PACLITAXEL,
LUNG, SARC), all `figures/main_*`, `figures/forest_cross_cancer*`, `figures/supp/*`, `REPORT.md`.
New helper tables: `nar/15_fixes/tables/03_pooled_table.md`, `03_interaction_table.md`, `03_foldalign_comparison.tsv/.md`.

Pooled held-out results, shortlist (HR old = benefit-score scale, HR new = risk-score scale; for the same fit HR_new = 1/HR_old).
Rows without a change in n and C differ only by the scale flip.

| Unit | Model | Outcome | n old -> new | events old -> new | C / AUROC old -> new | HR per SD: old (benefit scale) -> new (risk scale), p old -> new |
|---|---|---|---|---|---|---|
| CRC_5FU | organoid-only | survival | 792 -> 729 | 301 -> 273 | 0.511 -> 0.503 | 0.93 -> 1.07 (p 0.229 -> 0.314) |
| CRC_5FU | organoid-only | response | 723 -> 723 | 340 -> 340 | 0.490 -> 0.490 |  |
| CRC_5FU | patient-only | survival | 792 -> 729 | 301 -> 273 | 0.546 -> 0.536 | 0.89 -> 1.11 (p 0.0263 -> 0.0829) |
| CRC_5FU | patient-only | response | 723 -> 723 | 340 -> 340 | 0.508 -> 0.504 |  |
| CRC_5FU | joint | survival | 792 -> 729 | 301 -> 273 | 0.532 -> 0.525 | 0.89 -> 1.11 (p 0.0418 -> 0.0779) |
| CRC_5FU | joint | response | 723 -> 723 | 340 -> 340 | 0.502 -> 0.501 |  |
| CRC_OXALIPLATIN | organoid-only | survival | 264 -> 255 | 92 -> 88 | 0.453 -> 0.453 | 1.04 -> 0.97 (p 0.714 -> 0.777) |
| CRC_OXALIPLATIN | organoid-only | response | 427 -> 419 | 205 -> 199 | 0.551 -> 0.552 |  |
| CRC_OXALIPLATIN | patient-only | survival | 264 -> 255 | 92 -> 88 | 0.457 -> 0.449 | 1.24 -> 0.80 (p 0.0598 -> 0.0476) |
| CRC_OXALIPLATIN | patient-only | response | 427 -> 419 | 205 -> 199 | 0.515 -> 0.521 |  |
| CRC_OXALIPLATIN | joint | survival | 264 -> 255 | 92 -> 88 | 0.468 -> 0.464 | 1.08 -> 0.92 (p 0.432 -> 0.409) |
| CRC_OXALIPLATIN | joint | response | 427 -> 419 | 205 -> 199 | 0.540 -> 0.540 |  |
| CRC_IRINOTECAN | organoid-only | survival | 118 -> 109 | 100 -> 96 | 0.467 -> 0.465 | 0.98 -> 1.02 (p 0.863 -> 0.838) |
| CRC_IRINOTECAN | organoid-only | response | 178 -> 170 | 96 -> 90 | 0.491 -> 0.493 |  |
| CRC_IRINOTECAN | patient-only | survival | 118 -> 109 | 100 -> 96 | 0.425 -> 0.459 | 1.07 -> 0.84 (p 0.497 -> 0.129) |
| CRC_IRINOTECAN | patient-only | response | 178 -> 170 | 96 -> 90 | 0.450 -> 0.489 |  |
| CRC_IRINOTECAN | joint | survival | 118 -> 109 | 100 -> 96 | 0.473 -> 0.436 | 0.95 -> 0.93 (p 0.637 -> 0.504) |
| CRC_IRINOTECAN | joint | response | 178 -> 170 | 96 -> 90 | 0.493 -> 0.480 |  |
| BLCA_CISPLATIN | organoid-only | survival | 157 -> 157 | 99 -> 99 | 0.507 -> 0.507 | 0.93 -> 1.07 (p 0.451 -> 0.451) |
| BLCA_CISPLATIN | organoid-only | response | 229 -> 229 | 92 -> 92 | 0.499 -> 0.499 |  |
| BLCA_CISPLATIN | patient-only | survival | 157 -> 157 | 99 -> 99 | 0.560 -> 0.560 | 0.69 -> 1.45 (p 0.000513 -> 0.000513) |
| BLCA_CISPLATIN | patient-only | response | 229 -> 229 | 92 -> 92 | 0.580 -> 0.580 |  |
| BLCA_CISPLATIN | joint | survival | 157 -> 157 | 99 -> 99 | 0.528 -> 0.528 | 0.81 -> 1.24 (p 0.0238 -> 0.0238) |
| BLCA_CISPLATIN | joint | response | 229 -> 229 | 92 -> 92 | 0.536 -> 0.536 |  |
| BLCA_GEMCITABINE | organoid-only | survival | 64 -> 64 | 30 -> 30 | 0.423 -> 0.423 | 1.11 -> 0.90 (p 0.587 -> 0.587) |
| BLCA_GEMCITABINE | organoid-only | response | 54 -> 54 | 36 -> 36 | 0.356 -> 0.356 |  |
| BLCA_GEMCITABINE | patient-only | survival | 64 -> 64 | 30 -> 30 | 0.632 -> 0.632 | 0.65 -> 1.55 (p 0.00996 -> 0.00996) |
| BLCA_GEMCITABINE | patient-only | response | 54 -> 54 | 36 -> 36 | 0.721 -> 0.721 |  |
| BLCA_GEMCITABINE | joint | survival | 64 -> 64 | 30 -> 30 | 0.594 -> 0.594 | 0.65 -> 1.53 (p 0.0136 -> 0.0136) |
| BLCA_GEMCITABINE | joint | response | 54 -> 54 | 36 -> 36 | 0.642 -> 0.642 |  |
| PAAD_GEMCITABINE | organoid-only | survival | 92 -> 92 | 57 -> 57 | 0.371 -> 0.371 | 1.39 -> 0.72 (p 0.0268 -> 0.0268) |
| PAAD_GEMCITABINE | organoid-only | response | 59 -> 59 | 29 -> 29 | 0.432 -> 0.432 |  |
| PAAD_GEMCITABINE | patient-only | survival | 92 -> 92 | 57 -> 57 | 0.577 -> 0.577 | 0.92 -> 1.09 (p 0.557 -> 0.557) |
| PAAD_GEMCITABINE | patient-only | response | 59 -> 59 | 29 -> 29 | 0.520 -> 0.520 |  |
| PAAD_GEMCITABINE | joint | survival | 92 -> 92 | 57 -> 57 | 0.521 -> 0.521 | 1.14 -> 0.88 (p 0.381 -> 0.381) |
| PAAD_GEMCITABINE | joint | response | 59 -> 59 | 29 -> 29 | 0.491 -> 0.491 |  |
| LIHC_SORAFENIB | organoid-only | response | 67 -> 67 | 21 -> 21 | 0.965 -> 0.965 |  |
| LIHC_SORAFENIB | patient-only | response | 67 -> 67 | 21 -> 21 | 0.947 -> 0.947 |  |
| LIHC_SORAFENIB | joint | response | 67 -> 67 | 21 -> 21 | 0.947 -> 0.947 |  |
| OV_PACLITAXEL | organoid-only | survival | 564 -> 564 | 426 -> 426 | 0.502 -> 0.502 | 1.04 -> 0.97 (p 0.487 -> 0.487) |
| OV_PACLITAXEL | organoid-only | response | 272 -> 272 | 187 -> 187 | 0.492 -> 0.492 |  |
| OV_PACLITAXEL | patient-only | survival | 564 -> 564 | 426 -> 426 | 0.558 -> 0.558 | 0.82 -> 1.22 (p 0.000109 -> 0.000109) |
| OV_PACLITAXEL | patient-only | response | 272 -> 272 | 187 -> 187 | 0.631 -> 0.631 |  |
| OV_PACLITAXEL | joint | survival | 564 -> 564 | 426 -> 426 | 0.553 -> 0.553 | 0.88 -> 1.14 (p 0.00989 -> 0.00989) |
| OV_PACLITAXEL | joint | response | 272 -> 272 | 187 -> 187 | 0.584 -> 0.584 |  |
| OV_PLATINUM | organoid-only | survival | 683 -> 683 | 523 -> 523 | 0.520 -> 0.520 | 0.95 -> 1.06 (p 0.201 -> 0.201) |
| OV_PLATINUM | organoid-only | response | 347 -> 347 | 236 -> 236 | 0.583 -> 0.583 |  |
| OV_PLATINUM | patient-only | survival | 683 -> 683 | 523 -> 523 | 0.550 -> 0.550 | 0.84 -> 1.19 (p 0.000254 -> 0.000254) |
| OV_PLATINUM | patient-only | response | 347 -> 347 | 236 -> 236 | 0.610 -> 0.610 |  |
| OV_PLATINUM | joint | survival | 683 -> 683 | 523 -> 523 | 0.550 -> 0.550 | 0.85 -> 1.17 (p 0.000661 -> 0.000661) |
| OV_PLATINUM | joint | response | 347 -> 347 | 236 -> 236 | 0.618 -> 0.618 |  |
| LUNG_PLATINUM | organoid-only | survival | 277 -> 277 | 97 -> 97 | 0.482 -> 0.482 | 1.10 -> 0.91 (p 0.349 -> 0.349) |
| LUNG_PLATINUM | organoid-only | response | 126 -> 126 | 99 -> 99 | 0.490 -> 0.490 |  |
| LUNG_PLATINUM | patient-only | survival | 277 -> 277 | 97 -> 97 | 0.552 -> 0.552 | 0.89 -> 1.13 (p 0.217 -> 0.217) |
| LUNG_PLATINUM | patient-only | response | 126 -> 126 | 99 -> 99 | 0.594 -> 0.594 |  |
| LUNG_PLATINUM | joint | survival | 277 -> 277 | 97 -> 97 | 0.524 -> 0.524 | 0.98 -> 1.02 (p 0.835 -> 0.835) |
| LUNG_PLATINUM | joint | response | 126 -> 126 | 99 -> 99 | 0.571 -> 0.571 |  |
| SARC_GEMCITABINE | organoid-only | survival | 16 -> 16 | 9 -> 9 | 0.312 -> 0.312 | 1.61 -> 0.62 (p 0.251 -> 0.251) |
| SARC_GEMCITABINE | patient-only | survival | 16 -> 16 | 9 -> 9 | 0.637 -> 0.637 | 1.13 -> 0.88 (p 0.708 -> 0.708) |
| SARC_GEMCITABINE | joint | survival | 16 -> 16 | 9 -> 9 | 0.450 -> 0.450 | 1.64 -> 0.61 (p 0.226 -> 0.226) |

Interaction tests, joint model (survival HR flipped to risk scale; response OR unchanged):

| Cohort | Unit | Outcome | treated / comparison n old -> new | interaction old (benefit scale) -> new (risk scale) | p old -> new |
|---|---|---|---|---|---|
| GSE39582 | CRC_5FU | survival (RFS) | 200/288 -> 137/252 | 0.90 -> 1.19 | 0.523 -> 0.34 |
| GSE103479 | CRC_5FU | survival (PFS) | 66/84 -> 66/84 | 0.81 -> 1.22 | 0.514 -> 0.546 |
| GSE14333 | CRC_5FU | survival (DFS) | 87/139 -> 87/139 | 1.02 -> 0.98 | 0.917 -> 0.917 |
| TCGA-COADREAD | CRC_5FU | survival (PFS) | 140/369 -> 140/369 | 1.26 -> 0.80 | 0.322 -> 0.327 |
| TCGA-COADREAD | CRC_OXALIPLATIN | survival (PFS) | 133/400 -> 128/400 | 1.07 -> 0.93 | 0.715 -> 0.699 |
| GSE72970 | CRC_OXALIPLATIN | survival (PFS) | 40/84 -> 36/84 | 1.24 -> 0.80 | 0.286 -> 0.264 |
| GSE72970 | CRC_OXALIPLATIN | response | 40/84 -> 36/84 | 0.62 -> 0.62 | 0.269 -> 0.289 |
| GSE39582 | CRC_OXALIPLATIN | survival (RFS) | 11/288 -> 11/252 | 0.38 -> 2.69 | 0.253 -> 0.208 |
| GSE104645 | CRC_OXALIPLATIN | response | 113/54 -> 113/54 | 1.00 -> 0.97 | 0.993 -> 0.924 |
| GSE87211 | CRC_OXALIPLATIN | survival (DFS) | 80/95 -> 80/95 | 1.17 -> 0.87 | 0.646 -> 0.686 |
| GSE87211 | CRC_OXALIPLATIN | response | 90/111 -> 90/111 | 2.24 -> 2.24 | 0.0402 -> 0.0408 |
| GSE72970 | CRC_IRINOTECAN | survival (PFS) | 88/36 -> 84/36 | 0.70 -> 1.25 | 0.129 -> 0.32 |
| GSE72970 | CRC_IRINOTECAN | response | 88/36 -> 84/36 | 1.46 -> 1.30 | 0.407 -> 0.547 |
| GSE104645 | CRC_IRINOTECAN | response | 54/113 -> 54/113 | 0.85 -> 0.89 | 0.671 -> 0.756 |
| TCGA-BLCA | BLCA_CISPLATIN | survival (PFS) | 68/322 -> 68/322 | 0.77 -> 1.30 | 0.152 -> 0.152 |
| GSE31684 | BLCA_CISPLATIN | survival (RFS) | 35/58 -> 35/58 | 1.11 -> 0.90 | 0.758 -> 0.758 |
| GSE48276 | BLCA_CISPLATIN | survival (OS) | 16/37 -> 16/37 | 0.19 -> 5.30 | 0.00428 -> 0.00428 |
| TCGA-BLCA | BLCA_GEMCITABINE | survival (PFS) | 64/322 -> 64/322 | 0.73 -> 1.38 | 0.0678 -> 0.0678 |
| TCGA-PAAD | PAAD_GEMCITABINE | survival (PFS) | 92/76 -> 92/76 | 2.13 -> 0.47 | 0.000557 -> 0.000557 |
| GSE109211 | LIHC_SORAFENIB | response | 67/72 -> 67/72 | 1.34 -> 1.34 | 0.795 -> 0.795 |
| TCGA-OV | OV_PACLITAXEL | survival (PFS) | 332/57 -> 332/57 | 1.09 -> 0.92 | 0.702 -> 0.702 |
| GSE9891 | OV_PACLITAXEL | survival (RFS) | 191/47 -> 191/47 | 0.92 -> 1.09 | 0.599 -> 0.599 |
| GSE9891 | OV_PACLITAXEL | response | 183/46 -> 183/46 | 0.80 -> 0.80 | 0.51 -> 0.51 |
| GSE63885 | OV_PACLITAXEL | survival (OS) | 41/34 -> 41/34 | 1.27 -> 0.79 | 0.411 -> 0.411 |
| GSE63885 | OV_PACLITAXEL | response | 41/34 -> 41/34 | 0.33 -> 0.33 | 0.153 -> 0.153 |
| GSE15622 | OV_PACLITAXEL | response | 20/15 -> 20/15 | 1.39 -> 1.39 | 0.664 -> 0.664 |
| TCGA-OV | OV_PLATINUM | survival (PFS) | 370/46 -> 370/46 | 0.95 -> 1.05 | 0.813 -> 0.813 |
| GSE15622 | OV_PLATINUM | response | 15/20 -> 15/20 | 3.17 -> 3.17 | 0.139 -> 0.139 |
| TCGA-LUAD | LUNG_PLATINUM | survival (PFS) | 124/369 -> 124/369 | 1.19 -> 0.84 | 0.291 -> 0.291 |
| TCGA-LUSC | LUNG_PLATINUM | survival (PFS) | 104/386 -> 104/386 | 0.99 -> 1.01 | 0.942 -> 0.942 |
| GSE42127 | LUNG_PLATINUM | survival (OS) | 49/124 -> 49/124 | 1.25 -> 0.80 | 0.485 -> 0.485 |
| TCGA-SARC | SARC_GEMCITABINE | survival (PFS) | 16/230 -> 16/230 | 1.50 -> 0.67 | 0.2 -> 0.2 |

Other 03 numbers: joint inner-CV weight w = 0.25 chosen in 56 of 61 fits (was 55 of 61); joint model HR above 1 with CI
excluding 1 in 4 units (was "below 1 in 5 units" in the old scale; CRC_5FU dropped out, p 0.042 -> 0.078);
all-pathway joint C-index higher than shortlist in 6 of 10 units (unchanged count; CRC_5FU 0.560 vs 0.532 -> 0.546 vs 0.525);
learning curves CRC_5FU joint 25% -> 100% lab models 0.519 -> 0.526 became 0.514 -> 0.520, patient-only 10% -> 100%
0.519 -> 0.543 became 0.508 -> 0.533; joint vs patient-only at 10% patients CRC_5FU 0.522 vs 0.519 -> 0.516 vs 0.508,
CRC_OXALIPLATIN 0.512 vs 0.501 -> 0.517 vs 0.509. Interaction tests with p < 0.05: 3 of 32 (text said 33 tests; the table has 32 joint rows; same 3
tests). GSE39582 treated/comparison 200/288 -> 137/252 (survival, with RFS and expression); GSE72970 oxaliplatin 40/84 -> 36/84,
irinotecan 88/36 -> 84/36 (PFS interaction oxaliplatin: HR 1.24 old scale -> 0.80, p 0.29 -> 0.26; irinotecan 0.70 -> 1.25, p 0.13 -> 0.32).
Fold-aligned (I06) vs transductive, pooled C / AUROC: largest change OV_PACLITAXEL joint response 0.584 -> 0.643, organoid-only
response 0.492 -> 0.531; HR (risk) joint BLCA_CISPLATIN 1.24 -> 1.27 (p 0.024 -> 0.013), OV_PACLITAXEL 1.14 -> 1.22 (p 0.010 -> <0.001),
CRC_5FU 1.11 -> 1.11; full list in `15_fixes/tables/03_foldalign_comparison.tsv`.

### nar/05_A3_B1 (I01, I02, I04, I08, I09, I15)
Files: `code/b1_external.py`, `b1_figures.py`, `a3_figures.py`, `REPORT.md`, `tables/B1_external_results.tsv`,
`B1_summary_counts.tsv`, `B1_frozen_signatures.tsv`, new `tables/A3_adjusted_cox_epv*.tsv`, `A3_master_results.tsv` (flag column),
`figures/B1_Fig5*` (now two figures, Fig5 and Fig5b), `A3_Fig*`, `B1_Fig1`, `B1_Fig7`.
- B1 GSE39582 treated n: RFS 202 -> 138, OS 203 -> 139; comparator RFS 294 -> 258. GSE72970 FOLFOX 40 -> 36.
- B1 counts: treated p < 0.05 17 -> 14; in expected direction 13 -> 10; CRC 5-FU 10 -> 7; null p < 0.05 CRC 5-FU 6 -> 5, oxaliplatin 0 -> 1. Interactions (5, 3 expected), null total (10), BH q total (1) unchanged.
- Paper signature GSE39582 RFS HR 1.01 (0.81-1.27) p 0.93 -> 0.95 (0.72-1.25) p 0.70; interaction 0.80 p 0.19 -> 0.79 p 0.22; null p 0.18 -> 0.24; OS HR 1.05 p 0.73 -> 1.08 p 0.60.
- vdW 5-FU current GSE39582 OS rank-1 HR 1.43 p 0.006 -> 1.33 p 0.054; interaction 1.28 p 0.13 -> 1.19 p 0.34; null p 0.27 -> 0.45.
- Novartis PDX top-10 GSE39582 OS HR 1.64 p 0.00014 -> 1.48 p 0.008; interaction 1.60 p 0.006 -> 1.46 p 0.042; null p 0.021 -> 0.090; BH q 0.28 -> 0.65; RFS interaction 1.11 p 0.54 -> 1.13 p 0.52.
- GSE72970 oxaliplatin vdW current rank-1 PFS HR 1.22 p 0.23 -> 1.10 p 0.58; interaction 1.18 p 0.41 -> 1.09 p 0.67; FOLFOX median-split log-rank p 0.015 -> 0.17.
- A3 (I09): adjusted models with < 3 deaths per term: 27 of 49 (flagged); 20 of 22 unflagged keep direction; reduced model 48 of 54 (was 42 of 49 for all full models); rectum oxaliplatin adjusted HR 23.6 -> reduced model 5.37 (p 0.13).
- A3 (I08, I15): text only; binomial p 0.023 (treated, 7/57) and 0.00002 (untreated, 12/57) unchanged but labelled "assumes independence".

### nar/02_organoid_only (I03, I07, I14, I15, I17, I18, I24)
Files: `REPORT_GI.md`, `REPORT_BLCA_OV_LUAD_OTHER.md`, `code/report_sections.py`, `report_sections.md`, new
`tables/cross_set_best_direction.tsv`, `BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv`,
`BLCA_lee2018/sensitivity_one_line_per_patient/`, `tables/data_versions_02GI_sha256.tsv` (in 12_audit).
- No headline number changed. Text changes: COAD oxaliplatin best setting (Ridge k=4, p 0.001) "beats 99.2% of shuffles but HR 0.17 (opposite direction)" (was "not enough on its own"); bootstrap rows relabelled (0.0074 top-rank; 0.0072 / 0.007 rank-1 p).
- Bladder, one line per patient (7 lines instead of 9, new section 3.1b): treated p 0.0435, HR 2.63 (0.99-6.96) unchanged; shuffled-response empirical p 0.295 -> 0.263; bootstrap p 0.078 -> 0.084; rank-1 p bootstrap 0.064 -> 0.083.
- HCT-8 isogenic: 02 p = 0.039 is parental 0 h (3) vs resistant 0 h + 24 h + 48 h (9); stated.

### nar/08_depmap (I10, I23, I24)
`REPORT.md`: "every CI crosses 0" -> "all but one (PRISM navitoclax rho +0.405, CI 0.030-0.681, n 27, p 0.036)". HCT-8 note
(3 vs 3 per time point: 0 h p 0.18, 24 h p 0.007, 48 h p 0.019 vs 02 p 0.039). `code/02_partA.py`, `figures/Fig3*`: duplicate panel letter "e" fixed (re-rendered).

### nar/09_registered_test (I01, I16, I19, I20)
`REPORT.md`, new `code/../15_fixes/code/f09_sensitivity.py`, `tables/stage1_sensitivity_no_missing_start.tsv`, `tables/stage2_sensitivity_GSE39582_harmonised.tsv`.
- Sorafenib IPW p 1.000 -> NA (CI 0.82 (0.19-2.89) unchanged).
- Drop patients without a drug start day (PFI interaction HR old -> new): 5-FU 0.93 -> 0.93, oxaliplatin 0.93 -> 0.91, irinotecan unchanged, gemcitabine 1.12 (p 0.199) -> 1.09 (p 0.326), paclitaxel 1.00 -> 1.01, pemetrexed 0.98 -> 0.95, sorafenib 0.82 -> 0.71.
- GSE39582 harmonised sensitivity: RFS interaction 1.06 (0.74-1.52), p 0.37 (156 vs 298) -> 1.17 (0.80-1.71), one-sided p 0.21 (137 vs 252); OS 1.20 (0.81-1.80), p 0.18. Primary result unchanged (frozen).
- "pre-registered" -> "internally time-stamped analysis plan".

### nar/10_completeness (I23)
`figures/SuppFig4_cross_test_matrix_all_sets.*` redrawn by `15_fixes/code/f10_suppfig4.py` (fonts larger). I12 not fixed.

### docs/PROJECT_SUMMARY.md (I11)
Colon 5-FU current-data p 0.62 -> 0.866 (129 patients, 25 deaths); bladder cisplatin p 0.063 -> 0.067 (79 patients, 28 deaths).

### nar/12_audit (reconciliation)
`code/a02_master_numbers.py`, `a03_reconcile.py` edited and re-run; `tables/MASTER_NUMBERS.tsv` regenerated (copy in
`15_fixes/tables/`; previous in `15_fixes/backup/MASTER_NUMBERS_before.tsv`). New/updated rows: Joint_CRC_5FU C = 0.525, HR 1.11 (0.99-1.25), p = 0.078, n = 729/273 (was C 0.532, HR 0.89, p 0.042, n 792/301);
B1_counts 17/13 -> 14/10; Reg_S2_ox_GSE72970 notes; G11_GSE39582, G11_GSE72970 marked FIXED.
Remaining small n differences: 03 137 vs 252 (needs expression and RFS), 05 138 vs 258.

### Running logs
`15_fixes/logs/` (lc_crc, foldalign, b1, bladder, fig3_08, mf03, rn03). The machine restarted during the chain: everything had
finished except the 05 B1 step (b1_external.py, b1_summary.py, b1_figures.py, a3_figures.py), which was relaunched and finished.
