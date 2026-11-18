# W1 Lung cohorts: harvest, table rebuild and analyses

All numbers come from the files in `tables/` (written by `code/`). Seeds are fixed (20260101). p < 0.001 is shown as "<0.001".

## 1. What was done
- Added a `"lung"` list to `robust/fetch_geo_matrices.py` (and an optional cancer filter on the command line). Added 20 builders plus `LUNG_BUILDERS`/`SKIPPED_LUNG` to `robust/build_training_cohorts.py`, and a `main_v2` (run with `--out training_table_v2`), `("lung","lung")` in `GROUP`, and a `histology` column to `robust/build_training_table.py`. Only additions; old cohorts untouched.
- `robust/results/training_table/` sha256 checked before and after (`tables/training_table_v1.sha256`): all 7 files identical. In `training_table_v2` all 20,766 old rows and all old pathway scores are identical to v1.
- v2 has 24,644 samples (v1 + 3,878 lung samples). Old cohorts were not rebuilt from raw data (their raw files are not on this machine); v2 extends the finished v1 table, which gives identical old rows.
- Duplicate check (same test as v1: pathway-profile correlation > 0.95) of every new lung cohort against GSE42127, TCGA-LUAD, TCGA-LUSC: none flagged. Highest correlation 0.91 (`lung_duplicate_check.tsv`). What I think: no sign of shared patients, but this test is crude across platforms.
- Data file names, sizes and sha256: `tables/lung_data_files.tsv`.

## 2. Built and skipped, per cohort (`tables/lung_cohorts_built.tsv`, `lung_cohorts_check.tsv`, `lung_arm_counts.tsv`)
| Cohort | Status | n | Arms / outcome as built | Matches `_verify`? |
|---|---|---|---|---|
| GSE14814 | built (randomised) | 133 | chemo 71 / observation 62; 60 deaths (OS years x12) | yes |
| GSE68465 | built | 443 | chemo 89 / none 341 / unknown 13; 236 deaths; 205 relapses | yes |
| GSE37745 | built | 196 | chemo 29 / none 71 / unknown 96; 145 deaths | yes |
| GSE29013 | built | 55 | chemo 34 / none 21; 18 deaths | yes |
| GSE31210 | built | 204 | 22 patients flagged by the authors (incomplete resection/adjuvant) removed from 226; 30 deaths, 54 relapses | n differs on purpose |
| GSE50081 | built | 181 | treatment not recorded; 75 deaths | yes |
| GSE8894 | built | 138 | relapse-free survival only; 69 events | yes |
| GSE3141 | built | 111 | survival parsed from text; 58 deaths | n/a |
| GSE30219 | built | 293 tumours | 200 deaths, 114 relapses | 293 vs 296 in plan (14 normals excluded; I count 293 tumour samples) |
| GSE72094 | built | 442 | 122 deaths (420 with survival) | yes |
| GSE41271 | built | 275 | 124 deaths, 122 relapses from dates | yes |
| GSE157011 | built | 484 | 264 deaths | yes |
| GSE33072 (BATTLE) | built | 131 | sorafenib 39, erlotinib 25, erlotinib+bexarotene 1; 65 with no treatment recorded. 8-week disease control exists for 40 (sorafenib 39, 23 controlled; erlotinib arm has none in GEO); PFS for 66 | yes |
| GSE161537 | built | 82 | anti-PD-1/PD-L1; response known 76 (CR+PR 20); 61 deaths | yes |
| GSE162520 | built | 92 | design says surgery, treatment not recorded; 21 deaths, 35 PFS events | yes |
| GSE135222 | built | 27 | anti-PD-1; 21 PFS events | yes |
| GSE126044 | built | 16 | anti-PD-1; 5 responders / 11 | yes |
| GSE218989 | built | 355 | **bulk** RNA-seq (TPM matrix, 355 samples x 19,916 genes), not single-cell; 168 responders / 187 | yes |
| GSE207422 | built (bulk part only) | 24 | pre-treatment bulk biopsies; anti-PD-1 + platinum chemo; major pathological response 9 / 24 | n=24 of 39 (the other samples are single-cell or post-treatment) |
| GSE81089 | built | 196 | 93 deaths; time from surgery to last contact | yes (199 tumours, 3 without status) |

Skipped:
- GSE232569 (priority 2): only per-sample normalised expression files; no response or survival in the record.
- CPTAC-3 lung (priority 2): not a GEO series; drug names missing in the GDC treatment facet; not downloaded.
- POPLAR/OAK: EGA, not reachable.
- Priority 4 and excluded rows not in the W1 steps (GSE11969, GSE13213, GSE26939, GSE19188, GSE136961, GSE93157, GSE190265, GSE166449, GSE61676, GSE43580): not attempted.

Step 3 notes: GSE162520 and GSE81089 supplements hold only expression tables; their clinical data are in the series matrix and were used. Treatment was not recorded for either, so the drug field is empty.
Caveats: GSE161537 and GSE162520 are targeted panels (~2,560 genes), so 185 of 673 pathways cannot be scored (set to 0 after within-cohort standardising). Where chemo is "PLATINUM", the regimen was not given. GSE68465 stage is approximated from pT/pN text. GSE81089 uses -1 for a few missing FPKM values, set to 0.

## 3. Analysis a: response, leave-one-cohort-out (`tables/lung_loco_response.tsv`, Figure `lung_auc_by_cohort`)
AUC (95% CI, Hanley-McNeil; p Mann-Whitney) for the model trained on other cohorts of the same group. Other models in the table.
| Test cohort | n | responders | same-group LOCO | breast model | pan-cancer (all other) | pan-cancer without lung |
|---|---|---|---|---|---|---|
| GSE218989 (IO) | 355 | 168 | 0.566 (0.507-0.626) p=0.031 | 0.545 (0.485-0.606) | 0.549 (0.489-0.609) | 0.537 (0.476-0.597) |
| GSE161537 (IO) | 76 | 20 | 0.539 (0.390-0.689) p=0.61 | 0.588 (0.439-0.738) | 0.681 (0.537-0.825) p=0.017 | 0.663 (0.518-0.809) p=0.031 |
| GSE207422 (IO+chemo, MPR) | 24 | 9 | 0.541 (0.297-0.785) | 0.593 (0.351-0.835) | 0.719 (0.496-0.941) p=0.084 | 0.681 (0.451-0.912) |
| GSE126044 (IO) | 16 | 5 | 0.800 (0.539-1.061) p=0.069 | 0.800 | 0.836 (0.595-1.078) p=0.038 | 0.745 |
| TCGA-LUAD (chemo) | 79 | 57 | 0.552 (0.412-0.692) p=0.48 | 0.563 | 0.486 | 0.489 |
| TCGA-LUSC (chemo) | 55 | 45 | 0.651 (0.475-0.827) p=0.14 | 0.422 | 0.518 | 0.469 |
| GSE33072 (sorafenib labels) | 40 | 24 | no partner | 0.482 (0.297-0.667) | 0.492 | 0.518 |
What I think: the largest cohort (GSE218989) gives AUC about 0.55, a small signal. The pan-cancer model does best in GSE161537 (0.68), but with 20 responders the interval is wide. Upper limits above 1 for n=16 are the approximation; the figure cuts them at 1. Several tests, so single p-values should be read with care.

## 4. Analysis b: adjuvant chemo x score interaction (`tables/lung_adjuvant_interaction.tsv`, Figures `lung_forest_interaction`, `lung_gse14814_km`)
Overall survival, Cox stratified by cohort, adjusted for age and stage (stage IV removed), trained on the other 6 lung cohorts (ridge 300), score = predicted benefit (high = more benefit), 1 SD. Interaction HR below 1 means chemo works better when the score is high. "Treated" in TCGA means any recorded chemo.
| Held-out cohort | n | chemo/none | events | HR (95% CI) | p |
|---|---|---|---|---|---|
| GSE14814 (randomised) | 133 | 71/62 | 60 | 1.02 (0.59-1.75) | 0.95 |
| GSE29013 | 55 | 34/21 | 18 | 0.34 (0.13-0.95) | 0.038 |
| GSE37745 | 99 | 28/71 | 73 | 1.28 (0.80-2.06) | 0.30 |
| GSE42127 | 172 | 49/123 | 60 | 1.67 (0.82-3.40) | 0.16 |
| GSE68465 | 429 | 89/340 | 226 | 0.90 (0.69-1.19) | 0.47 |
| TCGA-LUAD | 478 | 122/356 | 166 | 0.79 (0.55-1.14) | 0.21 |
| TCGA-LUSC | 487 | 106/381 | 207 | 0.77 (0.53-1.12) | 0.17 |
| Pooled (inverse variance) | 1,853 | | 810 | 0.91 (0.77-1.06) | 0.23 |
In GSE14814, chemo vs observation overall HR 0.72 (p=0.20); in the predicted high-benefit half 0.55 (p=0.11, 30 deaths), low half 0.91 (p=0.80). What I think: in the randomised trial the score does not separate who benefits (interaction HR 1.02). The one nominal hit (GSE29013, 18 events) is one of 14 tests, and the pooled estimate is close to 1. GSE29013 aside, the TCGA cohorts lean below 1 but their intervals include 1; TCGA chemo is not randomised. The pooled p is approximate because training sets overlap.

## 5. Analysis d: same with the paper's cisplatin network shortlist (47 of 673 pathways; `lung_cisplatin_shortlist.tsv`)
The shortlist uses the authors' own proximity table (z <= -1.2816), not a new STRING run.
| Held-out cohort | HR (95% CI) | p | share of 100 random 47-pathway sets with a more favourable estimate |
|---|---|---|---|
| GSE14814 | 1.42 (0.86-2.34) | 0.17 | 0.75 |
| GSE29013 | 0.25 (0.06-1.12) | 0.071 | 0.02 |
| GSE37745 | 1.21 (0.76-1.94) | 0.42 | 0.74 |
| GSE42127 | 2.07 (0.97-4.41) | 0.059 | 0.84 |
| GSE68465 | 0.92 (0.71-1.19) | 0.53 | 0.05 |
| TCGA-LUAD | 0.87 (0.61-1.23) | 0.43 | 0.36 |
| TCGA-LUSC | 0.83 (0.57-1.21) | 0.34 | 0.15 |
| Pooled | 0.98 (0.84-1.15) | 0.81 | |
What I think: the shortlist does no better than all pathways (pooled 0.98 vs 0.91) and in the randomised trial points the other way (1.42), and sits near the middle of random sets. Selecting by network closeness to cisplatin gave no gain here.

## 6. Analysis c: untreated patients only (`tables/lung_untreated_prognosis.tsv`, Figure `lung_untreated_prognosis`)
Patients with no chemo recorded, or surgery-only cohorts with treatment not recorded; 14 cohorts tested, GSE29013 (6 events) not tested, 2,938 patients, 1,320 events. Endpoint OS (GSE8894 relapse-free). Scores learned on the chemo/no-chemo cohorts, never on the cohort tested.
| Score | pathways | pooled HR per SD (95% CI) | p |
|---|---|---|---|
| general (main-effect) prognostic score | all 673 | 1.36 (1.29-1.44) | <0.001 |
| general prognostic score | cisplatin shortlist | 1.21 (1.14-1.28) | <0.001 |
| chemo-benefit score | all 673 | 0.98 (0.93-1.04) | 0.47 |
| chemo-benefit score | shortlist | 0.99 (0.94-1.05) | 0.75 |
The general score was above 1 in 14 of 14 cohorts (all pathways); individual cohort rows are in the table (e.g. GSE68465 untreated n=340, 162 deaths, HR 1.38, 1.18-1.61; GSE14814 observation arm n=62, 30 deaths, HR 1.06, 0.69-1.63).
What I think: the general part of the pathway model predicts survival in people who got no chemo, so it carries prognostic information. The chemo-specific part does not predict survival in untreated people (HR near 1), which is what a drug-specific score should do, but it also did not show benefit in the randomised trial.

## 7. Files
Code: `code/` (`fetch_suppl.py`, `check_built.py`, `provenance.py`, `arm_counts.py`, `run_lung.py`, `make_figures.py`). Tables: `tables/`. Figures (PNG+PDF+caption): `figures/lung_forest_interaction`, `lung_auc_by_cohort`, `lung_untreated_prognosis`, `lung_gse14814_km`. Logs: `logs/`. Built cohort files: `data/external/train_cohorts/` (git-ignored). No file over 50 MB was added to the repo.
