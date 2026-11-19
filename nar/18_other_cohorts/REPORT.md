# W2: other-cancer cohort harvest

All numbers below come from files written by `code/run_v2.py`, `code/arm_check.py` and `code/make_figures.py` (tables in `tables/`, figures in `figures/`). Seeds are fixed inside the scripts. Data: `robust/results/training_table_v2/` (patients.tsv.gz, pathway scores) and the unchanged v1 table `robust/results/training_table/`.

## 1. What was built

- The v1 table is byte-identical to before (sha256 check: every file OK, no differences).
- The v2 table has 149 cohorts (v1 had 73). 76 are new; 20 are the lung cohorts from W1 and 56 are the non-lung cohorts from this work (W2).
- Samples in v2: 32226 (v1: 20766).

What I think: the harvest worked as planned; nothing in v1 was changed, so every old result can still be reproduced.

### New non-lung cohorts built, by cancer

| cancer_group | cohorts | samples | with_response | responders | with_survival | events | randomised |
|---|---|---|---|---|---|---|---|
| bladder | 3 | 479 | 89 | 16 | 479 | 173 | 0 |
| breast | 26 | 3215 | 1348 | 423 | 1884 | 502 | 3 |
| colorectal | 10 | 683 | 285 | 123 | 421 | 146 | 0 |
| head and neck | 2 | 367 | 0 | 0 | 367 | 184 | 0 |
| melanoma | 3 | 41 | 41 | 14 | 25 | 21 | 0 |
| myeloma | 1 | 264 | 239 | 113 | 0 | 0 | 0 |
| oesophagus | 1 | 41 | 41 | 10 | 0 | 0 | 0 |
| ovarian | 6 | 1187 | 55 | 32 | 1187 | 807 | 1 |
| prostate | 1 | 248 | 0 | 0 | 248 | 56 | 0 |
| stomach | 3 | 1057 | 0 | 0 | 1057 | 481 | 0 |

### Every new non-lung cohort

| cohort | cancer_group | samples | treated | untreated_comparison | with_response | responders | with_survival | events | randomised |
|---|---|---|---|---|---|---|---|---|---|
| GSE16446 | breast | 120 | 120 | 0 | 114 | 16 | 114 | 25 | 0 |
| GSE50948 | breast | 156 | 156 | 0 | 156 | 53 | 0 | 0 | 1 |
| GSE34138 | breast | 178 | 178 | 0 | 177 | 28 | 0 | 0 | 0 |
| GSE17705 | breast | 298 | 298 | 0 | 0 | 0 | 298 | 71 | 0 |
| GSE66305 | breast | 88 | 88 | 0 | 88 | 27 | 0 | 0 | 1 |
| GSE42822 | breast | 91 | 91 | 0 | 91 | 37 | 0 | 0 | 0 |
| GSE6861 | breast | 161 | 161 | 0 | 161 | 66 | 0 | 0 | 0 |
| GSE4779 | breast | 102 | 102 | 0 | 102 | 39 | 0 | 0 | 0 |
| GSE32603 | breast | 141 | 141 | 0 | 138 | 35 | 141 | 40 | 0 |
| GSE6532_GPL96 | breast | 327 | 190 | 137 | 0 | 0 | 119 | 44 | 0 |
| GSE6532_GPL570 | breast | 87 | 87 | 0 | 0 | 0 | 87 | 28 | 0 |
| GSE58644 | breast | 321 | 191 | 10 | 0 | 0 | 314 | 69 | 0 |
| GSE7390 | breast | 198 | 0 | 198 | 0 | 0 | 198 | 91 | 0 |
| GSE11121 | breast | 200 | 0 | 200 | 0 | 0 | 200 | 46 | 0 |
| GSE37946 | breast | 50 | 50 | 0 | 50 | 27 | 0 | 0 | 0 |
| GSE76360 | breast | 50 | 50 | 0 | 48 | 9 | 0 | 0 | 0 |
| GSE21974 | breast | 32 | 32 | 0 | 32 | 10 | 0 | 0 | 0 |
| GSE18728 | breast | 21 | 21 | 0 | 21 | 8 | 0 | 0 | 0 |
| GSE173839 | breast | 105 | 105 | 0 | 100 | 36 | 0 | 0 | 1 |
| GSE18864 | breast | 24 | 24 | 0 | 24 | 8 | 0 | 0 | 0 |
| GSE103668 | breast | 21 | 21 | 0 | 21 | 7 | 0 | 0 | 0 |
| GSE9195 | breast | 77 | 77 | 0 | 0 | 0 | 77 | 13 | 0 |
| GSE12093 | breast | 136 | 136 | 0 | 0 | 0 | 136 | 20 | 0 |
| GSE45255 | breast | 139 | 121 | 10 | 0 | 0 | 136 | 32 | 0 |
| GSE31519 | breast | 67 | 58 | 9 | 0 | 0 | 64 | 23 | 0 |
| GSE59515 | breast | 25 | 25 | 0 | 25 | 17 | 0 | 0 | 0 |
| GSE140082 | ovarian | 380 | 380 | 0 | 0 | 0 | 380 | 235 | 1 |
| GSE32062 | ovarian | 260 | 260 | 0 | 0 | 0 | 260 | 193 | 0 |
| GSE26712 | ovarian | 185 | 185 | 0 | 0 | 0 | 185 | 129 | 0 |
| GSE49997 | ovarian | 194 | 194 | 0 | 0 | 0 | 194 | 124 | 0 |
| GSE17260 | ovarian | 110 | 110 | 0 | 0 | 0 | 110 | 76 | 0 |
| GSE30161 | ovarian | 58 | 58 | 0 | 55 | 32 | 58 | 50 | 0 |
| GSE35452 | colorectal | 46 | 46 | 0 | 46 | 24 | 0 | 0 | 0 |
| GSE119409 | colorectal | 66 | 66 | 0 | 56 | 15 | 0 | 0 | 0 |
| GSE150082 | colorectal | 39 | 39 | 0 | 39 | 16 | 0 | 0 | 0 |
| GSE45404_GPL4133 | colorectal | 38 | 38 | 0 | 38 | 16 | 0 | 0 | 0 |
| GSE45404_GPL570 | colorectal | 42 | 42 | 0 | 42 | 19 | 0 | 0 | 0 |
| GSE145037 | colorectal | 31 | 31 | 0 | 31 | 20 | 0 | 0 | 0 |
| GSE133057 | colorectal | 33 | 33 | 0 | 33 | 13 | 33 | 11 | 0 |
| GSE115821 | melanoma | 14 | 14 | 0 | 14 | 2 | 0 | 0 | 0 |
| GSE168204 | melanoma | 2 | 2 | 0 | 2 | 2 | 0 | 0 | 0 |
| GSE100797 | melanoma | 25 | 25 | 0 | 25 | 10 | 25 | 21 | 0 |
| GSE13507 | bladder | 165 | 27 | 138 | 0 | 0 | 165 | 69 | 0 |
| GSE176307 | bladder | 90 | 90 | 0 | 89 | 16 | 90 | 79 | 0 |
| GSE32894 | bladder | 224 | 0 | 0 | 0 | 0 | 224 | 25 | 0 |
| GSE65858 | head and neck | 270 | 0 | 0 | 0 | 0 | 270 | 133 | 0 |
| GSE41613 | head and neck | 97 | 0 | 0 | 0 | 0 | 97 | 51 | 0 |
| GSE17536 | colorectal | 177 | 0 | 0 | 0 | 0 | 177 | 83 | 0 |
| GSE33113 | colorectal | 89 | 0 | 0 | 0 | 0 | 89 | 18 | 0 |
| GSE38832 | colorectal | 122 | 0 | 0 | 0 | 0 | 122 | 34 | 0 |
| GSE26253 | stomach | 432 | 432 | 0 | 0 | 0 | 432 | 177 | 0 |
| GSE84437 | stomach | 433 | 0 | 0 | 0 | 0 | 433 | 209 | 0 |
| GSE15459 | stomach | 192 | 0 | 0 | 0 | 0 | 192 | 95 | 0 |
| GSE104958 | oesophagus | 41 | 41 | 0 | 41 | 10 | 0 | 0 | 0 |
| GSE9782 | myeloma | 264 | 264 | 0 | 239 | 113 | 0 | 0 | 0 |
| GSE116918 | prostate | 248 | 0 | 0 | 0 | 0 | 248 | 56 | 0 |

### Samples before and after, per cancer group

| cancer_group | samples_v1 | samples_v2 | with_response_v1 | with_response_v2 | with_survival_v1 | with_survival_v2 |
|---|---|---|---|---|---|---|
| ACC | 79 | 79 | 1 | 1 | 79 | 79 |
| CHOL | 35 | 35 | 6 | 6 | 35 | 35 |
| DLBC | 47 | 47 | 0 | 0 | 47 | 47 |
| MESO | 83 | 83 | 28 | 28 | 83 | 83 |
| PCPG | 178 | 178 | 1 | 1 | 178 | 178 |
| THCA | 504 | 504 | 1 | 1 | 504 | 504 |
| THYM | 119 | 119 | 2 | 2 | 119 | 119 |
| UVM | 79 | 79 | 0 | 0 | 79 | 79 |
| bladder | 772 | 1251 | 243 | 332 | 600 | 1079 |
| brain | 1429 | 1429 | 114 | 114 | 1399 | 1399 |
| breast | 7772 | 10987 | 3569 | 4917 | 4861 | 6745 |
| cervix | 291 | 291 | 74 | 74 | 291 | 291 |
| colorectal | 2272 | 2955 | 753 | 1038 | 1861 | 2282 |
| head and neck | 519 | 886 | 77 | 77 | 519 | 886 |
| kidney | 1254 | 1254 | 348 | 348 | 1243 | 1243 |
| liver | 513 | 513 | 149 | 149 | 366 | 366 |
| lung | 1174 | 5052 | 134 | 645 | 1174 | 4570 |
| melanoma | 528 | 569 | 114 | 155 | 528 | 553 |
| oesophagus | 212 | 253 | 56 | 97 | 184 | 184 |
| ovarian | 850 | 2037 | 373 | 428 | 758 | 1945 |
| pancreas | 177 | 177 | 65 | 65 | 177 | 177 |
| prostate | 497 | 745 | 29 | 29 | 497 | 745 |
| sarcoma | 259 | 259 | 24 | 24 | 259 | 259 |
| stomach | 390 | 1447 | 117 | 117 | 390 | 1447 |
| testis | 134 | 134 | 43 | 43 | 134 | 134 |
| uterus | 599 | 599 | 86 | 86 | 599 | 599 |
| myeloma | 0 | 264 | 0 | 239 | 0 | 0 |

### Cohorts not built, with the reason

These come from `skipped.tsv` of the v2 table. Plain reason: the public record lacks the arm, response or survival data that the test needs.

| accession | reason |
|---|---|
| GSE279390 | VESPER bladder trial: gene data only; arm, response and survival are not in the public record |
| GSE156915 | FOCUS colorectal trial: gene data only; arm and survival are not in the public record |
| GSE319641 | NeoTRIP breast trial: gene data only; arm and response are not in the public record |
| GSE234519 | NSABP B-41 breast trial: gene data only; arm and response are not in the public record |
| GSE116335 | CALGB 40601 breast trial: gene data only; arm and response are not in the public record |
| GSE87304 | bladder chemo cohort: subtype only; response and survival are not in the public record |
| GSE235067 | bladder chemo cohort: 66-gene panel and no response in the public record |
| GSE186102 | eribulin breast cohort: targeted gene panel, only 10 complete responses |
| GSE232569 | Piedmont pemetrexed-platinum lung study: gene data only (per-sample files); response and survival are not in the public record |
| CPTAC-3 | lung part of CPTAC-3 (GDC): not a GEO series; drug names are missing in the GDC treatment facet and outcome was not checked, so no treated group can be defined |
| GSE28796 | dose-dense docetaxel breast cohort: only 2 complete responders (plan: exclude) |
| GSE20181 | letrozole breast cohort: response is not in the public record (only in the paper) |
| GSE87455 | breast chemotherapy + bevacizumab cohort: response is not in the public record (only in the paper) |
| GSE2034 | breast relapse-free survival cohort: only a bone-relapse flag is in the public record, not survival time |
| GSE16391 | BIG 1-98 breast cohort: the 0/1 code for letrozole versus tamoxifen is not explained in the record, and it is a matched case-control sample |
| GSE3149 | ovarian cohort: no clinical data in the public record |
| GSE13876 | ovarian cohort: two-colour custom array with several samples per patient (415 samples, 157 patients); survival status is in a free-text field per channel; not built |
| GSE94104 | rectal chemoradiotherapy cohort: the 3-level regression grade has no stated direction (which grade is best) in the record |
| GSE111177 | prostate androgen-deprivation cohort: no time-to-event or response definition in the record (metastasis yes/no only) |
| GSE5851 | cetuximab colorectal cohort: no response or survival in the record or the supplementary files (raw CEL files only) |
| GSE14210 | stomach cisplatin/5-FU cohort: no per-sample response or time-to-progression in the record or the supplementary files |
| GSE62254 | ACRG stomach cohort: clinical table is not in the GEO record or supplementary files |

What I think: every skipped cohort has a one-line reason. Most are gene-data-only deposits where the clinical table sits in a paper, not in the public record.

### Duplicate check

`tables/duplicate_check.tsv` has 385 comparisons (new cohort versus existing cohorts of the same cancer, per-sample correlation of the 673 pathway scores). Highest correlation seen: 0.865; samples above 0.95: 0. In the v2 table, `possible_duplicates` flags: GSE25066 (50 samples), GSE70691 (11 samples). Flagged samples are removed before analysis (the loader drops them).

What I think: GSE25066 (50 samples) and GSE70691 (11 samples) have samples that look like copies of samples in other cohorts; removing them avoids counting a patient twice.

## 2. Arm sizes and events versus the verify files

`tables/arm_check.tsv`: 69 of 90 comparisons match exactly. The arm and event checks for the key trials all match (NOAH 63/51/42 with 53 pCR; CHER-LOB 23/31/34 with 27 pCR; ICON7 199/181 with 235 PFS and 96 OS events; GSE14814 71/62 with 60 deaths; I-SPY2 durvalumab 71/34 with 36 pCR; APEX bortezomib arm 188). The rows that show False in the sample-count block are cohorts where the builder kept only part of the series (for example normal tissue, other time points or samples with no label dropped); the cohort-level reason is in the cohort build notes.

One real issue was found and fixed: the first build of GSE9782 (APEX) had empty progression and death flags because of a column-name mismatch in `robust/build_training_cohorts.py`. I corrected the two lines in that builder (progression 189 and death 156 events now read correctly; the GEO file leaves 0 blank). The v2 table on disk was built before the fix, so `run_v2.py` reads the corrected GSE9782 values from the builder. A full rebuild of the v2 table will give the same numbers.

What I think: arm sizes and events agree with the source records; the APEX fix is what makes the bortezomib survival test possible.

## 3. Old results reproduced

Leave-one-cohort-out AUC from the v1 pooled run was recomputed with the lung-only v2 table: 31 of 33 cohorts match to 9 decimals. All 13 breast cohorts match. The two that differ are TCGA-LUAD (v1 0.552 vs 0.545) and TCGA-LUSC (v1 0.651 vs 0.667); this is expected because the W1 lung cohorts were added to the lung training pool.

What I think: the builder did not disturb the old breast numbers, so differences later are due to added data only.

## 4. Response prediction, one cohort held out (Figures 1 and 2)

Table: AUC of the pathway model trained on the other cohorts of the same cancer. The 95% interval is the Hanley-McNeil interval. `v1` is the original run where the cohort existed; `v2` is with the new training pool. Rectal chemoradiotherapy cohorts are a separate group here.

| cancer | cohort | n | responders | AUC v1 | AUC v2 | 95% CI v2 | p | new cohort |
|---|---|---|---|---|---|---|---|---|
| bladder | GSE169455 | 149 | 48 | 0.501 | 0.536 | 0.436 to 0.636 | 0.4784 |  |
| bladder | GSE176307 | 89 | 16 |  | 0.600 | 0.441 to 0.760 | 0.2132 | yes |
| bladder | TCGA-BLCA | 71 | 45 | 0.579 | 0.566 | 0.429 to 0.703 | 0.3612 |  |
| breast | BRIGHTNESS | 482 | 236 | 0.674 | 0.689 | 0.642 to 0.736 | 0.0000 |  |
| breast | GSE123845 | 112 | 42 | 0.761 | 0.760 | 0.664 to 0.856 | 0.0000 |  |
| breast | GSE163882 | 222 | 80 | 0.715 | 0.706 | 0.632 to 0.779 | 0.0000 |  |
| breast | GSE16446 | 114 | 16 |  | 0.569 | 0.413 to 0.725 | 0.3805 | yes |
| breast | GSE173839 | 100 | 36 |  | 0.727 | 0.620 to 0.835 | 0.0002 | yes |
| breast | GSE18728 | 21 | 8 |  | 0.683 | 0.437 to 0.929 | 0.1847 | yes |
| breast | GSE18864 | 24 | 8 |  | 0.711 | 0.478 to 0.944 | 0.1056 | yes |
| breast | GSE20194 | 270 | 52 | 0.798 | 0.822 | 0.749 to 0.895 | 0.0000 |  |
| breast | GSE20271 | 178 | 26 | 0.801 | 0.757 | 0.645 to 0.870 | 0.0000 |  |
| breast | GSE21974 | 32 | 10 |  | 0.718 | 0.514 to 0.922 | 0.0535 | yes |
| breast | GSE22093 | 97 | 28 | 0.751 | 0.743 | 0.627 to 0.859 | 0.0002 |  |
| breast | GSE23988 | 61 | 20 | 0.846 | 0.866 | 0.757 to 0.975 | 0.0000 |  |
| breast | GSE240671 | 95 | 37 | 0.780 | 0.748 | 0.643 to 0.853 | 0.0000 |  |
| breast | GSE25066 | 438 | 82 | 0.817 | 0.837 | 0.781 to 0.893 | 0.0000 |  |
| breast | GSE32603 | 138 | 35 |  | 0.717 | 0.612 to 0.822 | 0.0001 | yes |
| breast | GSE32646 | 115 | 27 | 0.741 | 0.760 | 0.646 to 0.873 | 0.0000 |  |
| breast | GSE34138 | 177 | 28 |  | 0.867 | 0.780 to 0.955 | 0.0000 | yes |
| breast | GSE37946 | 50 | 27 |  | 0.718 | 0.577 to 0.859 | 0.0086 | yes |
| breast | GSE41998 | 253 | 69 | 0.715 | 0.719 | 0.644 to 0.794 | 0.0000 |  |
| breast | GSE42822 | 91 | 37 |  | 0.801 | 0.705 to 0.898 | 0.0000 | yes |
| breast | GSE4779 | 102 | 39 |  | 0.744 | 0.641 to 0.846 | 0.0000 | yes |
| breast | GSE50948 | 156 | 53 |  | 0.617 | 0.522 to 0.712 | 0.0168 | yes |
| breast | GSE59515 | 25 | 17 |  | 0.471 | 0.222 to 0.719 | 0.8419 | yes |
| breast | GSE66305 | 88 | 27 |  | 0.613 | 0.482 to 0.744 | 0.0924 | yes |
| breast | GSE6861 | 161 | 66 |  | 0.690 | 0.605 to 0.775 | 0.0000 | yes |
| breast | GSE76360 | 48 | 9 |  | 0.581 | 0.366 to 0.796 | 0.4595 | yes |
| breast | ISPY2 | 986 | 319 | 0.700 | 0.716 | 0.680 to 0.751 | 0.0000 |  |
| breast | TCGA-BRCA | 202 | 184 | 0.672 | 0.678 | 0.563 to 0.793 | 0.0128 |  |
| colorectal | GSE104645 | 181 | 112 | 0.522 | 0.522 | 0.436 to 0.608 | 0.6205 |  |
| colorectal | GSE19860 | 40 | 15 | 0.565 | 0.565 | 0.378 to 0.752 | 0.5025 |  |
| colorectal | GSE28702 | 83 | 42 | 0.577 | 0.577 | 0.454 to 0.700 | 0.2310 |  |
| colorectal | GSE62080 | 21 | 9 | 0.806 | 0.806 | 0.605 to 1.000 | 0.0209 |  |
| colorectal | GSE72970 | 124 | 63 | 0.493 | 0.493 | 0.391 to 0.595 | 0.8887 |  |
| colorectal | GSE87211 | 201 | 35 | 0.427 | 0.427 | 0.326 to 0.528 | 0.1767 |  |
| colorectal | TCGA-COAD | 74 | 53 | 0.418 | 0.418 | 0.270 to 0.566 | 0.2752 |  |
| kidney | BRAUN2020_CHECKMATE | 281 | 44 | 0.472 | 0.472 | 0.381 to 0.564 | 0.5587 |  |
| kidney | EMTAB3267 | 43 | 19 | 0.482 | 0.482 | 0.307 to 0.658 | 0.8545 |  |
| lung | GSE161537 | 76 | 20 |  | 0.604 | 0.455 to 0.754 | 0.1694 | yes |
| lung | GSE207422 | 24 | 9 |  | 0.407 | 0.173 to 0.642 | 0.4743 | yes |
| lung | GSE218989 | 355 | 168 |  | 0.533 | 0.473 to 0.593 | 0.2877 | yes |
| lung | GSE33072 | 40 | 24 |  | 0.492 | 0.307 to 0.677 | 0.9450 | yes |
| lung | TCGA-LUAD | 79 | 57 | 0.552 | 0.545 | 0.405 to 0.686 | 0.5366 |  |
| lung | TCGA-LUSC | 55 | 45 | 0.651 | 0.667 | 0.495 to 0.839 | 0.1040 |  |
| melanoma | GSE100797 | 25 | 10 |  | 0.640 | 0.412 to 0.868 | 0.2555 | yes |
| melanoma | GSE78220 | 28 | 15 | 0.472 | 0.477 | 0.259 to 0.695 | 0.8538 |  |
| melanoma | GSE91061 | 49 | 10 | 0.582 | 0.590 | 0.384 to 0.795 | 0.3921 |  |
| melanoma | TCGA-SKCM | 37 | 23 | 0.525 | 0.562 | 0.371 to 0.753 | 0.5414 |  |
| oesophagus | GSE104958 | 41 | 10 |  | 0.655 | 0.449 to 0.861 | 0.1493 | yes |
| oesophagus | GSE45670 | 28 | 11 |  | 0.695 | 0.488 to 0.903 | 0.0904 | yes |
| ovarian | GSE15622 | 35 | 22 | 0.465 | 0.483 | 0.281 to 0.684 | 0.8779 |  |
| ovarian | GSE30161 | 55 | 32 |  | 0.728 | 0.596 to 0.860 | 0.0043 | yes |
| ovarian | GSE51373 | 28 | 16 | 0.635 | 0.615 | 0.404 to 0.825 | 0.3182 |  |
| ovarian | GSE63885 | 75 | 65 | 0.623 | 0.611 | 0.434 to 0.788 | 0.2651 |  |
| ovarian | GSE9891 | 229 | 146 | 0.565 | 0.624 | 0.551 to 0.697 | 0.0018 |  |
| rectal chemoradiotherapy | GSE119409 | 56 | 15 |  | 0.444 | 0.276 to 0.612 | 0.5293 | yes |
| rectal chemoradiotherapy | GSE133057 | 33 | 13 |  | 0.685 | 0.492 to 0.877 | 0.0801 | yes |
| rectal chemoradiotherapy | GSE145037 | 31 | 20 |  | 0.473 | 0.256 to 0.689 | 0.8204 | yes |
| rectal chemoradiotherapy | GSE150082 | 39 | 16 |  | 0.492 | 0.306 to 0.678 | 0.9431 | yes |
| rectal chemoradiotherapy | GSE35452 | 46 | 24 |  | 0.557 | 0.390 to 0.724 | 0.5165 | yes |
| rectal chemoradiotherapy | GSE45404_GPL4133 | 38 | 16 |  | 0.602 | 0.417 to 0.788 | 0.2939 | yes |
| rectal chemoradiotherapy | GSE45404_GPL570 | 42 | 19 |  | 0.611 | 0.438 to 0.784 | 0.2251 | yes |

Summary per cancer (cohorts that existed in v1):

| cancer | cohorts | mean_auc_v1 | mean_auc_v2 | new_cohorts | mean_auc_new |
|---|---|---|---|---|---|
| bladder | 2.0 | 0.54 | 0.551 | 1.0 | 0.6 |
| breast | 13.0 | 0.752 | 0.754 | 15.0 | 0.682 |
| colorectal | 7.0 | 0.544 | 0.544 |  |  |
| kidney | 2.0 | 0.477 | 0.477 |  |  |
| lung | 2.0 | 0.601 | 0.606 | 4.0 | 0.509 |
| melanoma | 3.0 | 0.526 | 0.543 | 1.0 | 0.64 |
| oesophagus |  |  |  | 2.0 | 0.675 |
| ovarian | 4.0 | 0.572 | 0.583 | 1.0 | 0.728 |
| rectal chemoradiotherapy |  |  |  | 7.0 | 0.552 |

What I think: for breast, the old 13 cohorts move little when the pool grows, and the new breast cohorts mostly fall in the same 0.6 to 0.87 range, with the small cohorts (n of 21 to 50) having wide intervals. Ovarian, rectal chemoradiotherapy and the other new groups mostly have intervals that include 0.5; the cohorts with a clear AUC are listed above with p. I leave the reading to you.

## 5. Randomised, drug-specific response tests (Figure 3)

| test | drug arm n (responders) | control arm n (responders) | AUC drug arm | AUC control arm | p (AUC difference) | interaction coefficient | p (interaction) |
|---|---|---|---|---|---|---|---|
| NOAH (randomised): trastuzumab + chemo vs chemo alone | 63 (31) | 51 (13) | 0.328 | 0.589 | 0.025 | -0.888 | 0.041 |
| CHER-LOB (randomised): chemotherapy + lapatinib vs chemotherapy + trastuzumab | 31 (8) | 23 (6) | 0.370 | 0.373 | 0.986 | -0.604 | 0.379 |
| CHER-LOB (randomised): chemotherapy + trastuzumab + lapatinib vs chemotherapy + trastuzumab | 34 (13) | 23 (6) | 0.242 | 0.578 | 0.039 | -1.125 | 0.128 |
| APEX (randomised): bortezomib vs dexamethasone, response | 71 (38) | 70 (28) | 0.511 | 0.527 | 0.871 | -0.038 | 0.911 |
| TRANSFER trastuzumab: trained in other trastuzumab-treated breast cohorts, tested in NOAH | 63 (31) | 51 (13) | 0.592 | 0.518 | 0.535 | 0.369 | 0.404 |

What I think: in NOAH the score separates responders better in the chemotherapy-alone arm (AUC 0.589) than in the trastuzumab arm (0.328), so there is no sign of a trastuzumab-specific effect, and the arm difference has p=0.025 in the opposite direction to the hope. The same direction shows in the CHER-LOB dual arm. APEX shows no arm difference. Arm sizes are 23 to 71 patients, so intervals are wide.

### Survival interaction in randomised trials (Figure 4, orange rows)

| test | endpoint | n | treated | events | interaction HR | 95% CI | p |
|---|---|---|---|---|---|---|---|
| ICON7 (randomised): bevacizumab x score | PFS | 380 | 199 | 235 | 0.870 | 0.672 to 1.126 | 0.291 |
| ICON7 (randomised): bevacizumab x score | OS | 380 | 199 | 96 | 0.793 | 0.519 to 1.210 | 0.282 |
| APEX (randomised): bortezomib x score | PFS | 140 | 71 | 111 | 0.857 | 0.572 to 1.283 | 0.454 |
| APEX (randomised): bortezomib x score | OS | 156 | 80 | 77 | 0.673 | 0.427 to 1.062 | 0.089 |

What I think: ICON7 and APEX both give an interaction HR below 1 (higher score, more benefit from the drug) but every interval includes 1. APEX overall survival is the closest (HR 0.673, p=0.089). The scores were cross-fitted, so no patient was scored by a model that saw them.

## 6. Survival interaction in observational cohorts (Figure 4, blue rows)

| setting | test cohort | n | treated | events | train cohorts | interaction HR | 95% CI | p |
|---|---|---|---|---|---|---|---|---|
| colorectal 5-FU after surgery | GSE103479 | 152 | 66 | 42 | 4 | 0.866 | 0.474 to 1.582 | 0.640 |
| colorectal 5-FU after surgery | GSE14333 | 226 | 87 | 50 | 4 | 0.681 | 0.390 to 1.189 | 0.176 |
| colorectal 5-FU after surgery | GSE39582 | 489 | 201 | 142 | 4 | 1.165 | 0.837 to 1.621 | 0.365 |
| colorectal 5-FU after surgery | TCGA-COAD | 366 | 97 | 74 | 4 | 1.195 | 0.721 to 1.980 | 0.489 |
| colorectal 5-FU after surgery | TCGA-READ | 123 | 37 | 21 | 4 | 0.773 | 0.217 to 2.758 | 0.692 |
| bladder chemotherapy | GSE13507 | 165 | 27 | 69 | 3 | 1.841 | 0.965 to 3.514 | 0.064 |
| bladder chemotherapy | GSE31684 | 93 | 35 | 39 | 3 | 0.907 | 0.434 to 1.898 | 0.796 |
| bladder chemotherapy | GSE48276 | 53 | 16 | 28 | 3 | 1.381 | 0.627 to 3.041 | 0.423 |
| bladder chemotherapy | TCGA-BLCA | 401 | 79 | 170 | 3 | 1.375 | 0.979 to 1.930 | 0.066 |
| breast tamoxifen vs none | TCGA-BRCA | 570 | 242 | 86 | 9 | 1.366 | 0.829 to 2.249 | 0.221 |
| breast chemotherapy vs none | TCGA-BRCA | 869 | 541 | 125 | 5 | 1.149 | 0.808 to 1.634 | 0.440 |
| ovarian platinum vs none | TCGA-OV | 414 | 368 | 300 | 7 | 1.256 | 0.811 to 1.945 | 0.307 |

What I think: the new tamoxifen and ovarian platinum cohorts are in the training pool (9 and 7 training cohorts) but only TCGA cohorts have at least 15 treated, 15 untreated and 15 events inside one cohort, so they are the only test rows for breast and ovarian. The smallest p is 0.064 (bladder GSE13507) and 0.066 (TCGA-BLCA), both with HR above 1.

## 7. Files

- Code: `code/arm_check.py`, `code/run_v2.py`, `code/make_figures.py`, `code/make_report.py`, `code/build_all.sh`.
- Tables: `tables/arm_check.tsv`, `duplicate_check.tsv`, `reproduce_v1_part1.tsv`, `part1_v2.tsv`, `part2_v2.tsv`, `part3_v2.tsv`, `new_randomised_tests.tsv`, `new_randomised_survival_tests.tsv`.
- Figures (PNG, PDF, caption): `figures/fig1_breast_loco_forest`, `fig2_other_cancers_loco_forest`, `fig3_randomised_drug_specific_forest`, `fig4_survival_interaction_forest`.
- Not overwritten: `robust/results/train_combined/` and `robust/results/training_table/`.
