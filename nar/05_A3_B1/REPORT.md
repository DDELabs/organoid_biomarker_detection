# A3 statistical controls + B1 external validation

Folder `nar/05_A3_B1/`. Code in `code/`, tables in `tables/`, figures in `figures/` (PNG 300 dpi + PDF + a
plain-language caption file `*.caption.md` for each). Seed 20261117. At most 2 CPU processes. Nothing was retrained:
every patient score uses the coefficients saved by the organoid-only agents (`nar/02_organoid_only/`).

Re-run: `cd code; python a3_build_gi.py; python a3_build_other.py; python a3_analyze.py; python a3_figures.py;
python b1_external.py; python b1_summary.py; python b1_figures.py` (b1_external.py takes about 30 minutes because of
the random-signature nulls).

## In short

- **A3, 57 organoid-only analyses.** 7 of 57 have a treated-patient log-rank p < 0.05; chance alone gives 2.85
  (binomial p = 0.023 for 7 or more). Only 4 of the 7 point the expected way (predicted non-responders die sooner).
  After correction for 57 tests, none passes: smallest Benjamini-Hochberg q = 0.063 (dacarbazine, Novartis melanoma
  PDX, HR 0.21, the opposite direction), smallest Holm p = 0.063. The paper's colon 5-FU result: p = 0.014,
  q = 0.27, Holm p = 0.79.
- Untreated patients show **more** small p-values than treated patients (12 of 57 vs 2.85 expected, binomial
  p = 0.00002). Note (audit I08): the binomial p assumes the 57 tests are independent. They are not: untreated
  tests share patients (11 LUAD tests use the same 327 untreated patients; ovarian carboplatin and cisplatin give
  identical p). The binomial p values here are therefore too small; a permutation version that keeps the
  dependence is gap G9 (`nar/14_gaps_stats/`). The same pathway splits often mark general prognosis, not drug response.
- Score x treatment interaction: 5 of 57 p < 0.05 (2.85 expected); 3 of the 5 are in the opposite direction; smallest
  BH q = 0.36.
- Adjusting for age, sex, stage, purity, MSI and regimen keeps the direction in 42 of 49 estimable analyses
  (correlation of log HR 0.77). But 27 of those 49 full models have fewer than 3 deaths per estimated term
  (flagged in `tables/A3_adjusted_cox_epv.tsv`, audit I09); among the 22 with >= 3 deaths per term, 20 keep the
  direction. A reduced model (at most deaths/3 - 1 covariates: stage, then age, purity, MSI, sex, regimen) keeps
  the direction in 48 of 54 estimable analyses.
- Power: median 17 deaths per analysis; median minimal detectable HR (80% power) = 3.9. Only 2 of 57 analyses had
  80% power for HR 1.5 (ovarian carboplatin and paclitaxel), 5 of 57 for HR 2.0. A significant HR from these
  analyses is expected to be 1.7-8 times too large (type-M exaggeration) if the true HR is 1.25-2.
- **B1, external cohorts (196 tests, 15 frozen signatures, 21 cohort/endpoint/comparison settings).** 14 treated-arm
  tests had p < 0.05 (9.8 expected), 10 in the expected direction (before the GSE39582 / FOLFIRINOX definition fix, audit
  I01/I02: 17 and 13), mostly the 5-FU signatures in GSE87211 DFS and the Novartis PDX signature in GSE39582. Of 114 tests with a comparator arm, 5 interactions had p < 0.05 (3 expected direction). 10 of 196
  beat the random-signature null at 0.05 (9.8 expected). After BH across the 196 primary tests, 1 remains
  (vdW 5-FU current top-10 in GSE87211 rectal DFS, HR 1.92 per SD, q = 0.002; no comparator arm there).
- The paper's own signature (BH3-only pathway) in GSE39582: RFS HR 0.95 (0.72-1.25), p = 0.70 (stage II-III, n = 138; before: 1.01, p = 0.93, n = 202); interaction with
  untreated 0.79, p = 0.22 (before 0.80, p = 0.19).

What I think: the organoid-only signals look like the number of small p-values one gets from 57 small tests plus a
real prognostic (not drug-specific) component. The 5-FU multi-pathway scores carry prognostic information in
external colorectal cohorts, but they do not separate treated from untreated patients more than random pathway sets
do. The study is mainly limited by the small number of deaths per test.

---

## Part A3: statistical controls across all 57 organoid-only analyses

### What was done

1. **Exact rebuild.** For all 57 lab-model set x drug x TCGA cohort analyses (23 GI from
   `nar/02_organoid_only/_data/results/*.pkl`; 34 bladder / ovarian / lung / sarcoma / PDX from the agent's cached
   ssGSEA matrices and saved coefficient tables), I rebuilt the per-patient score, the median split and survival.
   All 114 log-rank p-values (57 treated + 57 untreated) are identical to the stored ones
   (`_data/a3_check_gi.tsv`, `_data/a3_check_other.tsv`).
2. **Master table** `tables/A3_master_results.tsv`: n, deaths (also per predicted group), log-rank p, HR (95% CI),
   untreated n / p / HR, response-shuffle (permutation) p, random-pathway bootstrap p, direction, BH q, Holm p,
   adjusted HRs, interaction HRs, power, minimal detectable HR, and the multi-pathway (top-10) log-rank p.
3. **Multiple testing** across the whole family of 57: Benjamini-Hochberg q (false-discovery rate) and Holm
   (chance of any false positive). Counts: `tables/A3_expected_by_chance.tsv`.
4. **Adjusted Cox** (treated patients): age (per 10 years), sex (if both present), stage (I-IV as a number), tumour
   purity (ABSOLUTE), MSI (colon/rectum: Xena `MSI_updated_Oct62011`, then BCR field, then the 2019 COAD biotab), and
   the regimen = indicators for the two most common partner drugs given to 10-90% of the treated patients (leucovorin
   excluded). For 5-FU this is oxaliplatin (FOLFOX vs 5-FU without oxaliplatin) and bevacizumab. Missing values =
   median plus a "missing" flag. Binary covariates with < 3 patients in a level or no deaths in a level were dropped
   (cannot be estimated). All model terms: `tables/A3_adjusted_cox_all_terms.tsv`.
5. **Interaction** score x treatment: Cox on treated + untreated (no drug record) patients: continuous score per SD,
   group, and their product; also adjusted (age, sex, stage, purity, MSI) and a median-split version.
6. **Power** (Schoenfeld formula, 1:1 split, two-sided 0.05): power for true HR 1.5 and 2.0, minimal detectable HR
   at 80% power.
7. **Winner's curse** (Gelman & Carlin 2014): for the 7 analyses with p < 0.05, power, type S (wrong sign) and type M
   (exaggeration) for true HRs 1.25, 1.5, 2.0 (`tables/A3_winners_curse.tsv`). SE of log HR = 2/sqrt(deaths), so it
   also works when the Cox HR is not estimable.

### Results

**Expected by chance (Fig. 2).**

| Family (57 tests) | p < 0.05 observed | expected | binomial p | BH q < 0.05 | Holm < 0.05 |
|---|---|---|---|---|---|
| treated log-rank (main test) | 7 | 2.85 | 0.023 (assumes independence) | 0 | 0 |
| response-shuffle p | 3 | 2.85 | 0.55 | 0 | 0 |
| random-pathway bootstrap p (rank-1 p statistic, not the paper's top-rank statistic) | 5 | 2.85 | 0.16 | 0 | 0 |
| untreated log-rank | 12 | 2.85 | 0.00002 (assumes independence; too small, see G9) | 0 | 0 |
| score x treatment interaction | 5 | 2.85 | 0.16 | 0 | 0 |
| multi-pathway (top-10) log-rank | 7 | 2.85 | 0.023 | 0 | 0 |

What I think: there are somewhat more small treated p-values than chance (7 vs 2.85), but the untreated test, which
cannot reflect drug response, has even more (12). So the excess is better explained by general prognosis than by
drug-specific prediction. The 57 tests share patients (e.g. colon and colon + rectum; 11 LUAD untreated tests use the
same 327 patients), so the binomial p assumes independence that does not hold and is too small. Use it as a guide
only; the dependence-aware permutation version is gap G9 (`nar/14_gaps_stats/`). [Audit I08.]

**The 7 nominal hits (master table; Fig. 1).**

| Analysis | n / deaths | HR (95% CI) | p | BH q | Holm | perm. p | rank-1 p bootstrap | direction |
|---|---|---|---|---|---|---|---|---|
| Novartis PDX melanoma, dacarbazine (SKCM) | 32 / 21 | 0.21 (0.07-0.58) | 0.0011 | 0.063 | 0.063 | 0.065 | 0.013 | opposite |
| Lung cell lines, pemetrexed (LUAD) | 59 / 15 | 0.19 (0.05-0.69) | 0.0048 | 0.14 | 0.27 | 0.12 | 0.11 | opposite |
| vdW colon, 5-FU (authors' 2019 files) | 114 / 8 | not estimable (0 vs 8 deaths) | 0.014 | 0.27 | 0.79 | 0.034 | 0.007 | expected |
| Shi pancreas, gemcitabine | 104 / 51 | 0.52 | 0.023 | 0.32 | 1.0 | 0.10 | 0.022 | opposite |
| Lung cell lines, gemcitabine | 25 / 16 | 3.17 | 0.030 | 0.34 | 1.0 | 0.14 | 0.048 | expected |
| Novartis PDX NSCLC, erlotinib | 21 / 7 | 7.61 | 0.040 | 0.35 | 1.0 | 0.045 | 0.027 | expected |
| Lee bladder, cisplatin (paper NES) | 79 / 20 | 2.63 (0.99-6.96) | 0.044 | 0.35 | 1.0 | 0.30 | 0.066 | expected |

Overall, 25 analyses have HR > 1 (expected), 29 HR < 1, 3 not estimable (2 lean expected, 1 opposite).
What I think: no single analysis survives correction for 57 tests. Three of the 7 hits are in the wrong direction.

**Adjusted Cox (Fig. 3).** 49 analyses gave an estimable adjusted HR for the median split; 42 of 49 kept the same
direction; correlation of log HR adjusted vs unadjusted = 0.77. Analyses with < 10 deaths become unstable (very wide
intervals, e.g. rectum oxaliplatin adjusted HR 23.6). Audit I09: 27 of the 49 full models have < 3 deaths per
estimated term (split + covariates) and are flagged (`low_epv_full` in `tables/A3_master_results.tsv`, details in
`tables/A3_adjusted_cox_epv.tsv`, summary `tables/A3_adjusted_cox_epv_summary.tsv`). Among the 22 unflagged models,
20 keep the direction. A reduced model with at most deaths/3 - 1 covariates (priority: stage, age, purity, MSI, sex,
regimen; >= 3 deaths per term) keeps the direction in 48 of 54 estimable analyses; in analyses with < 6 deaths it has
no covariates and equals the unadjusted Cox. Reduced model, rectum oxaliplatin: HR 5.37, p = 0.13 (was 23.6). Notable: paper bladder cisplatin 2.63 -> 2.27 (0.80-6.46),
p = 0.13; dacarbazine 0.21 -> 0.09 (0.02-0.46), p = 0.004; pemetrexed 0.19 -> 0.09 (0.02-0.44), p = 0.003; sarcoma
docetaxel 1.17 -> 6.86 (1.37-34.4), p = 0.019. For the paper's colon 5-FU result the median-split HR is not
estimable (no deaths in predicted responders); the continuous score gives HR 1.52 per SD (p = 0.091) unadjusted and
1.97 (p = 0.15) adjusted.
What I think: adjustment mostly keeps the picture; it does not create a robust new expected-direction signal.

**Interaction (Fig. 4).** 25 of 57 interaction HRs are > 1. p < 0.05 in 5: bladder cisplatin (current data) 1.64
(1.04-2.57), p = 0.033 (expected); and 4 below 1: Novartis PDX pancreas gemcitabine 0.55 (p = 0.006), lung pemetrexed
0.46 (p = 0.015), sarcoma gemcitabine 0.56 (p = 0.031), melanoma dacarbazine 0.58 (p = 0.037). Smallest BH q = 0.36.
The paper's colon 5-FU: interaction 1.66 (0.92-2.98), p = 0.093.
What I think: only bladder cisplatin (current data) behaves like a drug-specific marker, at nominal level only.

**Power and winner's curse (Fig. 5).** Median 17 deaths per analysis; minimal detectable HR at 80% power from 1.46
(ovarian, 216-218 deaths) to 25 (3 deaths); median 3.9. Only 2 of 57 analyses had 80% power for HR 1.5, 5 for HR 2.0.
For the 7 hits, power for a true HR 1.5 is 0.08-0.31 and a significant estimate would exaggerate the true HR by a
factor 1.8-4.5 (type M); for true HR 2.0, 1.2-2.7. For the paper's colon 5-FU (8 deaths): power 0.09 at HR 1.5,
exaggeration 4.2, chance of wrong sign 6%.
What I think: most analyses could only detect very large effects; the significant ones, if real, are probably much
smaller than they look.

### Figures (A3)

| Figure | Content |
|---|---|
| `figures/A3_Fig1_master_forest` | HR (95% CI) of all 57 analyses, coloured by direction, with p and BH q |
| `figures/A3_Fig2_QQ_expected_by_chance` | (a) QQ plot vs uniform; (b) observed vs expected count of p < 0.05 |
| `figures/A3_Fig3_adjusted_vs_unadjusted` | adjusted vs unadjusted HR: (a) median split, (b) continuous |
| `figures/A3_Fig4_interaction_forest` | score x treatment interaction HR per analysis |
| `figures/A3_Fig5_power_vs_events` | power vs deaths and minimal detectable HR |

Tables: `A3_master_results.tsv`, `A3_expected_by_chance.tsv`, `A3_adjusted_cox_all_terms.tsv`, `A3_winners_curse.tsv`.

---

## Part B1: frozen organoid signatures in external patient cohorts

### What was done

- **Signatures** (`tables/B1_frozen_signatures.tsv`): for each organoid-only analysis whose drug matches a cohort,
  (i) the rank-1 pathway alone and (ii) the multi-pathway model = top-10 pathways by |Ridge coefficient|, with the
  saved coefficients. Identical signatures were used once (ovarian cisplatin = carboplatin; ovarian docetaxel =
  paclitaxel: same organoid screen, same top pathways). 15 distinct signatures x 2 forms: colorectal 5-FU (vdW 2019
  files, vdW current, Novartis CRC PDX), oxaliplatin and irinotecan (vdW 2019 and current), bladder cisplatin (paper
  NES, current), ovarian platinum and paclitaxel (Vias), liver sorafenib (Ji, Broutier), lung platinum (cell-line
  cisplatin, carboplatin).
- **Scoring**: within-sample rank scores (singscore-like, platform-independent; training table, and STORM scored by the
  nar/03 loader); each pathway standardised within cohort; score = sum(coef x z), then standardised within cohort.
  Pre-specified direction: higher = predicted resistant, so HR or OR(non-response) > 1 is expected.
- **Tests** (`tables/B1_external_results.tsv`, 196 rows): survival -> Cox HR per SD (also median-split log-rank);
  response -> logistic OR of non-response per SD. With a comparator arm, the effect in the comparator and the
  score x arm interaction. Primary statistic = interaction when a comparator exists, else the treated-arm effect.
- **Random-signature null**: rank-1 -> all 673 pathways used alone (exact null); top-10 -> 1,000 random 10-pathway sets
  with random normal weights. Null p = share of random signatures with |z| at least as large.
- Cohorts: GSE39582 (5-FU-based = 5FU + FUFOL + FOLFOX, 139 patients (138 with RFS), vs untreated 264 (258), stage II-III, shared definition of audit I01; GSE72970 FOLFIRINOX patients removed from both arms, audit I02; subsets 5FU 76, FUFOL 50, FOLFOX 12 treated; RFS and OS),
  GSE103479, GSE14333, GSE72970 (FOLFOX vs FOLFIRI), GSE87211, GSE28702, GSE19860, GSE62080, GSE104645 (added:
  colorectal, in the training table); GSE169455, GSE52219, GSE48276, GSE70691 and GSE31684 (added) for bladder;
  GSE9891, GSE15622, GSE63885, GSE51373 for ovarian; GSE109211 (STORM) for liver; GSE42127 for lung.

### Results

Counts (`tables/B1_summary_counts.tsv`):

| Group | tests | treated p < 0.05 (expected dir.) | interaction p < 0.05 (expected dir.) of tests with comparator | null p < 0.05 | BH q < 0.05 |
|---|---|---|---|---|---|
| all | 196 | 14 (10) | 5 (3) of 114 | 10 | 1 |
| colorectal 5-FU | 96 | 7 (7) | 1 (1) of 48 | 5 | 1 |
| colorectal oxaliplatin | 32 | 1 (0) | 0 of 24 | 1 | 0 |
| colorectal irinotecan | 16 | 0 | 1 (1) of 12 | 0 | 0 |
| bladder cisplatin | 20 | 1 (0) | 1 (1) of 8 | 2 | 0 |
| ovarian platinum | 12 | 1 (1) | 0 of 4 | 0 | 0 |
| ovarian taxane | 12 | 0 | 2 (0) of 10 | 2 | 0 |
| liver sorafenib | 4 | 4 (2) | 0 of 4 | 0 | 0 |
| lung platinum | 4 | 0 | 0 of 4 | 0 | 0 |

Expected by chance at 0.05: 9.8 of 196.

**Colorectal (Fig. 1, Fig. 5a-f).**
- Paper signature (vdW 5-FU 2019 files, rank-1 BH3-only) in GSE39582: RFS HR 0.95 (0.72-1.25), p = 0.70,
  interaction 0.79, p = 0.22, null p = 0.24 (GSE39582 now = 5FU, FUFOL, FOLFOX, stage II-III, n = 138 with RFS; before
  n = 202: HR 1.01, p = 0.93, interaction 0.80, p = 0.19, null p = 0.18); OS HR 1.08, p = 0.60 (before 1.05, p = 0.73). Subsets (unchanged): 5FU alone RFS 0.88 (p = 0.48); FUFOL 1.12
  (p = 0.68); FOLFOX 0.57 (p = 0.40, 12 patients, 7 events). GSE103479 OS 0.98 (p = 0.93); GSE14333 DFS 1.18
  (p = 0.40); GSE72970 PFS 0.94 (p = 0.53); GSE87211 DFS 1.38 (1.00-1.89), p = 0.048 (null p = 0.33).
  What I think: the paper's 5-FU pathway does not carry over to independent 5-FU cohorts.
- vdW 5-FU current, GSE39582 OS: rank-1 HR 1.33 (0.99-1.78), p = 0.054, and HR 1.13 in untreated; interaction
  1.19, p = 0.34; null p = 0.45 (before the fix: HR 1.43, p = 0.006, interaction 1.28, p = 0.13, null p = 0.27). Top-10 in GSE87211 DFS: HR 1.92 (1.43-2.56), p = 0.00001, null p = 0.001, the one
  BH-surviving test (no comparator arm).
- Novartis CRC PDX 5-FU top-10, GSE39582 OS: HR 1.48 (1.11-1.98), p = 0.008; untreated HR 1.00; interaction 1.46
  (1.01-2.10), p = 0.042; random-signature null p = 0.090; BH q = 0.65 (before the fix: HR 1.64, p = 0.00014,
  interaction 1.60, p = 0.006, null p = 0.021, q = 0.28).
  What I think: this is the most "predictive-looking" result in B1 (treated-only effect, beats random sets), but it is
  OS only (RFS interaction 1.13, p = 0.52), no longer beats the random sets at 0.05 (null p 0.090) and does not survive correction for 196 tests. It is a candidate for a
  pre-registered test, not a finding.
- Oxaliplatin, GSE72970 FOLFOX vs FOLFIRI: vdW current rank-1 PFS HR 1.10 (0.78-1.55), p = 0.58, interaction 1.09,
  p = 0.67 (36 FOLFOX patients after removing 4 FOLFIRINOX; before: 40, HR 1.22, p = 0.23, interaction 1.18, p = 0.41;
  median-split log-rank in FOLFOX p = 0.17, before 0.015, Fig. 5e). Response cohorts GSE28702 / GSE19860
  ORs 0.60-0.78 (opposite).
- Irinotecan: GSE72970 vdW current top-10 PFS interaction 1.65 (1.02-2.67), p = 0.042, null p = 0.070 (unchanged: the 84 FOLFIRI patients did not change).

**Bladder (Fig. 2, Fig. 5g-h).** Paper cisplatin rank-1: GSE169455 OR non-pCR 1.18 (0.84-1.67), p = 0.34; GSE52219
1.57, p = 0.37; GSE48276 OS 1.28, p = 0.40 (16 treated); GSE70691 1.22, p = 0.25; GSE31684 1.16, p = 0.52. All in the
expected direction, none below 0.05. Current-data cisplatin rank-1: GSE31684 interaction 1.82 (1.07-3.09), p = 0.028,
null p = 0.033; but GSE169455 OR 0.73 (p = 0.078) and GSE48276 HR 0.43 (p = 0.021) go the other way.
What I think: the paper bladder signature leans the right way in all 5 cohorts, each with a small
effect; that consistency is worth a pooled pre-specified test.

**Ovarian (Fig. 3, Fig. 5i-j).** Platinum (Vias, oxaliplatin screen) rank-1: GSE9891 RFS HR 1.23 (1.06-1.43),
p = 0.006, null p = 0.092; OS 1.21, p = 0.056. Paclitaxel rank-1: GSE9891 OS platinum+taxane 1.14 (p = 0.31) vs platinum
alone, interaction 0.88, p = 0.52. Paclitaxel top-10, GSE63885: interactions 0.51 (OS, p = 0.024) and 0.12 (response,
p = 0.011): opposite direction.

**Liver (Fig. 4, Fig. 6).** STORM (sorafenib vs placebo): all four scores are strongly linked to response in both arms
(e.g. Broutier rank-1 OR non-response 3.94, p = 0.001 on sorafenib and 5.42 on placebo); Ji rank-1 OR 0.025 on
sorafenib and 0.091 on placebo (opposite). Interactions 0.27-1.28, all p > 0.07.
What I think: these pathways are prognostic for recurrence in STORM, not sorafenib-specific.

**Lung (Fig. 4, Fig. 5k-l).** GSE42127 adjuvant platinum vs observation (49 treated, 9 deaths): HRs 0.60-1.14, all
p >= 0.095; interactions 0.69-1.12, p >= 0.26.

**Pancreas:** no external pancreatic cohort with drug and outcome exists in `data/curated/trials/` or the training
table, so pancreas could only be tested in TCGA (part A3).

### Figures (B1)

| Figure | Content |
|---|---|
| `figures/B1_Fig1_forest_colorectal` | forest: treated effect and interaction, every colorectal signature x cohort |
| `figures/B1_Fig2_forest_bladder` | same, bladder |
| `figures/B1_Fig3_forest_ovarian` | same, ovarian |
| `figures/B1_Fig4_forest_liver_and_lung` | same, liver (STORM) and lung |
| `figures/B1_Fig5_KM_main_cohorts` | KM high vs low score, treated vs comparator: GSE39582, GSE103479, GSE72970, GSE48276, GSE9891, GSE42127 |
| `figures/B1_Fig6_STORM_sorafenib` | responder share by score and arm, STORM |
| `figures/B1_Fig7_heatmap_signature_by_cohort` | signed z of the primary test, signature x cohort |

Tables: `B1_frozen_signatures.tsv`, `B1_external_results.tsv` (with BH q of the primary test), `B1_summary_counts.tsv`;
per-patient scores of the rank-1 signatures in `_data/b1_patient_scores.tsv.gz`.

---

## Data versions (`tables/data_versions.tsv`)

| File | Bytes | sha256 (first 16) |
|---|---|---|
| nar/02_organoid_only/tables/cross_set_summary.tsv | 15,663 | 4e39a711d34389a6 |
| nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv | 22,746 | eaff50b122a7ed88 |
| nar/02_organoid_only/_data/raw/Survival_SupplementalTable_S1_20171025_xena_sp.tsv | 2,419,504 | a5e704158bb5c51c |
| robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt | 901,812 | f430a975433d82e0 |
| nar/05_A3_B1/_data/raw/COADREAD_clinicalMatrix | 607,684 | 4742eb0d07a5c40e |
| robust/results/training_table/patients.tsv.gz | 345,036 | a68da5ad8a0aa0ae |
| robust/results/training_table/pathway_scores_geo.tsv.gz | 11,879,734 | ce3e32a577cd2999 |
| robust/results/training_table/pathway_scores_tcga.tsv.gz | 11,150,748 | 0949e9fbe202dfd3 |
| data/curated/trials/GSE109211/expression.tsv.gz | 7,313,144 | 0f577406d247f2f3 |
| nar/03_organoid_plus_patient/cache/scores_GSE109211.tsv.gz | 823,935 | 2b31c76dfdf5ce62 |

Other read-only inputs: `nar/02_organoid_only/_data/results/*.pkl`, `nar/02_organoid_only/*/tables/*_coefficients.tsv`,
`nar/02_organoid_only/BLCA_OV_LUAD_OTHER/_data/cache/ssgsea_*_tcga.tsv.gz` and `treat_*.tsv`, MOESM3 sheet
`Figure_2h` (bladder paper NES), BCR biotab drug files. `COADREAD_clinicalMatrix` (MSI) was downloaded from the UCSC
Xena TCGA hub into `_data/raw/`.

## Limits (one line each)

- MSI is known for about half of colon/rectum patients; the rest enter with a "missing" flag.
- Adjusted HR not estimable in 8 analyses (no deaths in one predicted group or < 8 deaths): left as NA.
- The "untreated" TCGA group = no drug record, not a randomised control; the interaction is observational.
- GSE39582 FOLFIRI vs untreated (11 patients) was skipped (< 8 treated with RFS after excluding stage IV).
- No pancreatic external cohort exists; sorafenib only in STORM, whose "response" is recurrence-based.
- GSE14333 and GSE31684 record "chemotherapy" without the drug name; used as 5-FU / cisplatin proxies, as in nar/03.
