# Statistics gap-filling: G9, G10, G12, G15, G16, G17, G18

Folder `nar/14_gaps_stats/`. Code in `code/`, tables in `tables/`, figures in `figures/` (PNG 300 dpi + PDF + a
plain-language `*.caption.md` for each). Seed 20261118. At most 2 CPU processes. Nothing was retrained for G9, G10, G12,
G17 and G18: every patient score is the frozen score rebuilt by `nar/05_A3_B1` (`_data/a3_patients_*.tsv.gz`, which
reproduce all 114 stored log-rank p-values exactly). Input files with size and sha256: `tables/data_versions.tsv`.

Re-run (about 15 minutes in all):
`cd code; python g9_permutation_null.py; python g9_figure.py; python g10_meta_equivalence.py; python g12_baseline_table.py;
python g17_pfi_survey.py; python g18_score_agreement.py; python g16_joint_extended_grid.py --jobs 2;
python g16_bootstrap_cindex.py; python g15_citation_table.py; python data_versions.py`

## In short

| Gap | Main number | What I think |
|---|---|---|
| G9 shuffle null | 7 of 57 treated p < 0.05; shuffle null mean 2.87, 95th percentile 6, shuffle p = 0.030. Untreated 12 of 57, shuffle p = 0.0006 | The shared-patient structure barely changes the binomial picture (0.023 -> 0.030). The untreated excess is still larger. |
| G10 meta-analysis | Pooled HR per SD 0.96 (95% CI 0.90-1.03), tau^2 0.006, I^2 20%, 95% prediction interval 0.81-1.15 | Across 57 analyses the average effect is close to 1; an average HR above 1.02 per SD is ruled out (one-sided 95%). |
| G12 baseline table | Colon: treated median age 60 vs 73 untreated; stage III-IV 84% vs 23%; MSI-H 6% vs 19% | The "untreated" group is a very different population, so untreated splits measure prognosis in a different case mix. |
| G15 citations | OpenAlex, Crossref and Semantic Scholar refused by the network proxy; PubMed has no "cited by" link here | In PubMed I found no paper that re-tests the BH3-only or amino-acid biomarker. This is not a full citation search. |
| G16 joint model | Joint minus patient-only: 0 of 20 unit x outcome intervals above 0 (original grid: 2 below 0); joint minus organoid-only: 6 of 20 above 0 | Adding organoid data gives no measurable gain over patient-only models; the gain over organoid-only comes from the patient data. |
| G17 PFI | PFI: 8 of 57 treated p < 0.05 (binomial p = 0.007), only 1 of 8 in the expected direction; untreated 4 of 57 | The PFI hits are mostly in the opposite direction and are mostly other analyses than the OS hits (2 overlap). |
| G18 score agreement | Median Spearman ssGSEA vs rank score per pathway 0.98 (5,688 pathway x cohort pairs); with the rank score 3 of the 7 hits stay below 0.05 | The scores agree closely. Borderline hits move with the method, including the paper's colon 5-FU result (p 0.014 -> 0.32). |

---

## G9. Dependence-aware null for the 57-analysis survey

**What was done.** Survival (time and event together) was shuffled between patients inside each TCGA survival pool,
5,000 times. Colon, rectum and colon+rectum analyses share patients, so they form one pool with one shared shuffle.
The paper-era 2019 colon files and the paper's bladder biotab survival are their own pools. Shuffling was done
separately for "treated with any survey drug" and "no drug record" patients. A patient in several analyses always gets
the same shuffled survival. Frozen scores and median splits were not touched. All 57 treated and 57 untreated log-rank
tests were recomputed each time (own numpy log-rank, identical to statsmodels `survdiff` on all 57 analyses).
Code `code/g9_permutation_null.py`, `code/g9_figure.py`. Tables `tables/G9_permutation_null_summary.tsv`,
`tables/G9_null_count_distribution.tsv`, `tables/G9_per_analysis_permutation.tsv`. Figure `figures/G9_permutation_null`.

| Statistic | Observed | Shuffle null mean | Null 95th pct | Shuffle p |
|---|---|---|---|---|
| treated: count p < 0.05 (of 57) | 7 | 2.87 | 6 | 0.030 |
| treated: count p < 0.05 in expected direction | 4 | 1.47 | 4 | 0.066 |
| treated: count p < 0.01 | 2 | 0.53 | 2 | 0.100 |
| treated: smallest p | 0.0011 | 0.018 (mean) | 0.054 (5th pct) | 0.050 |
| untreated: count p < 0.05 | 12 | 2.96 | 7 | 0.0006 |
| untreated: smallest p | 0.0012 | 0.019 (mean) | 0.059 (5th pct) | 0.063 |
| groups of overlapping patients (of 11) with >= 1 treated p < 0.05 | 6 | 2.34 | 5 | 0.012 |

The 11 patient groups are clusters of analyses whose treated patients overlap by at least half (for example all
colorectal 5-FU analyses form one cluster). Family-wise adjusted p (min-p over the 57 shuffled tests): melanoma PDX
dacarbazine 0.050 (opposite direction), lung cell-line pemetrexed 0.22 (opposite), paper's colon 5-FU 0.54, all others
>= 0.72.

What I think: the shuffle null is almost the same as the binomial (null mean 2.87 vs 2.85; both have 95th percentile
6; p 0.030 vs 0.023), because most of the 57 analyses use different patients or different drugs. So the earlier
binomial summary was fair. The best single analysis is at the edge of family-wise significance (0.050) and points the
opposite way. The untreated excess (12, p = 0.0006) remains the stronger signal, which supports reading the excess as
prognosis rather than drug response.

## G10. Evidence-of-absence: random-effects meta-analysis and equivalence bounds

**What was done.** Effect = log HR per SD of the frozen organoid score in treated patients (continuous Cox, from
`A3_master_results.tsv`), signed so that > 0 is the expected direction. REML tau^2 (between-analysis variance),
Hartung-Knapp-Sidik-Jonkman 95% CI (a CI that also accounts for the uncertainty in tau^2), 95% prediction interval
(where the true effect of a new, similar analysis would fall). Also DerSimonian-Laird, subsets by drug class and TCGA
cohort, one analysis per patient cluster (the one with most deaths), median-split HR, adjusted HR and the score x
treatment interaction. Equivalence: a bound b is "ruled out" when the one-sided 95% upper limit (upper end of the 90% CI)
is below b. Caveat: analyses sharing patients are not independent; the one-per-cluster row is the check for that.
Code `code/g10_meta_equivalence.py`. Tables `tables/G10_random_effects_meta.tsv`, `tables/G10_per_analysis_equivalence.tsv`,
`tables/G10_equivalence_counts.tsv`. Figure `figures/G10_meta_forest`.

| Effect, subset | k | Pooled HR (95% CI) | one-sided 95% upper | tau^2 | I^2 | 95% prediction interval |
|---|---|---|---|---|---|---|
| continuous score, all | 57 | 0.96 (0.90-1.03), p = 0.27 | 1.02 | 0.006 | 20% | 0.81-1.15 |
| one analysis per patient cluster | 11 | 0.88 (0.72-1.07) | 1.03 | 0.031 | 57% | 0.56-1.38 |
| excluding the two paper-era file sets | 53 | 0.95 (0.89-1.02) | 1.01 | 0.007 | 17% | 0.80-1.14 |
| antimetabolites | 16 | 0.95 (0.81-1.10) | 1.07 | 0.020 | 45% | 0.67-1.33 |
| platinum | 13 | 1.02 (0.89-1.18) | 1.15 | 0.007 | 23% | 0.81-1.29 |
| microtubule (taxanes, vinca) | 10 | 1.02 (0.90-1.15) | 1.13 | 0 | 0% | 0.90-1.15 |
| topoisomerase / anthracycline | 8 | 0.90 (0.72-1.11) | 1.07 | 0 | 0% | 0.71-1.12 |
| targeted / other | 9 | 0.89 (0.66-1.20) | 1.14 | 0 | 0% | 0.66-1.21 |
| median split, all | 54 | 0.98 (0.87-1.11) | 1.09 | 0.011 | 19% | 0.77-1.25 |
| adjusted continuous score, all | 49 | 0.99 (0.92-1.06) | 1.05 | 0 | 20% | 0.92-1.06 |
| score x treatment interaction, all | 57 | 0.95 (0.87-1.04) | 1.02 | 0.018 | 11% | 0.71-1.26 |

Per TCGA cohort the pooled HR runs from 0.79 (colon+rectum, k = 3) to 1.11 (bladder, k = 8, 95% CI 0.88-1.40); all
include 1 (table). Per analysis: the one-sided 95% upper limit is below 1.25 in 18 of 57 analyses, below 1.5 in 33,
below 2.0 in 44 and below 3.0 in 52. Equivalence in both directions (HR between 1/b and b) holds for 4, 13, 27 and 44
analyses at b = 1.25, 1.5, 2.0 and 3.0.

What I think: the pooled estimate is the honest summary. On average, one SD higher organoid-predicted resistance goes
with an HR of 0.96 in treated patients, and an average effect of HR 1.1 per SD or more is ruled out. Single analyses are
weak: 24 of 57 cannot exclude an HR of 1.5 per SD. Within drug classes the upper limits are 1.07-1.15. A small true
effect in a few drug-cohort pairs remains possible (prediction interval up to 1.15; 1.38 with one analysis per patient
cluster). Part (b) of the original G10 item (replication probability of the two headline results) is covered by the
power / type-M numbers in `nar/05_A3_B1/tables/A3_winners_curse.tsv` and was not repeated.

## G12. Treated vs "no drug record" patients per TCGA cohort

**What was done.** Same clinical variables as the A3 adjusted Cox model (`a3_analyze.clinical()`): age, sex, stage
(TCGA-CDR), ABSOLUTE purity, MSI (Xena COADREAD + 2019 COAD biotab). One row per cohort (treated with any survey drug)
and one per cohort x drug. Mann-Whitney test for age and purity; Fisher's exact test for sex, stage III-IV and MSI-H.
Code `code/g12_baseline_table.py`. Tables `tables/G12_baseline_treated_vs_untreated.tsv` (all numbers, incl. % male),
`tables/G12_baseline_readable.tsv` (formatted). Figure `figures/G12_baseline_treated_vs_untreated`.

| Cohort | n treated / untreated | median age (p) | % stage III-IV (p) | median purity (p) | % MSI-H (p) |
|---|---|---|---|---|---|
| BLCA | 107 / 289 | 65 vs 70 (8e-05) | 71 vs 66 (0.33) | 0.61 vs 0.60 (0.91) | |
| BRCA | 435 / 298 | 53 vs 60 (1e-08) | 27 vs 25 (0.49) | 0.57 vs 0.61 (0.004) | |
| COAD | 136 / 302 | 60 vs 73 (4e-14) | 84 vs 23 (2e-32) | 0.65 vs 0.66 (0.41) | 6 vs 19 (0.015) |
| COAD (2019 paper files) | 121 / 298 | 60 vs 73 (8e-13) | 83 vs 23 (8e-30) | 0.66 vs 0.65 (0.80) | 7 vs 19 (0.036) |
| COADREAD | 199 / 391 | 60 vs 71 (2e-17) | 81 vs 25 (1e-39) | 0.65 vs 0.66 (0.29) | 4 vs 15 (0.006) |
| READ | 61 / 89 | 62 vs 68 (7e-04) | 78 vs 30 (2e-08) | 0.64 vs 0.70 (0.26) | 0 vs 0 |
| LIHC | 30 / 328 | 64 vs 61 (0.67) | 45 vs 24 (0.024) | 0.69 vs 0.72 (0.60) | |
| LUAD | 164 / 327 | 64 vs 67 (0.006) | 38 vs 12 (8e-11) | 0.45 vs 0.45 (0.92) | |
| OV | 379 / 37 | 58 vs 64 (0.003) | 94 vs 97 (0.71) | 0.81 vs 0.74 (0.14) | |
| PAAD | 109 / 61 | 65 vs 65 (0.50) | 6 vs 2 (0.42) | 0.40 vs 0.40 (0.24) | |
| SARC | 43 / 184 | 57 vs 63 (0.010) | no stage in TCGA-CDR | 0.86 vs 0.75 (0.005) | |
| SKCM | 32 / 313 | 60 vs 61 (0.72) | 46 vs 47 (1.0) | 0.66 vs 0.72 (0.66) | |

What I think: in colon, rectum and lung the untreated patients are older and mostly early-stage (chemotherapy is given
for stage III and advanced disease). In colon they are also about three times more often MSI-high. So "untreated" is not
a control arm; it is a different, mostly lower-risk population. Interaction tests that use it should be read with this
in mind, and the 12 of 57 untreated splits below 0.05 speak to prognosis in a different case mix.

## G15. Papers citing Kong et al. 2020

**What was done.** I tested the citation APIs first: `api.openalex.org`, `api.crossref.org` and
`api.semanticscholar.org` were all refused by this machine's network proxy (CONNECT 403) on 2026-11-18, so no cited-by
list could be downloaded. Substitute: PubMed "similar articles" for PMID 33127883 (top 20), targeted PubMed keyword
searches (organoid + network-based + drug response; organoid + transfer learning / drug-response prediction, 2021-2026)
and a same-group search (Kim S, POSTECH, network-based). Raw PubMed records `_data/pubmed/*.json`; code
`code/g15_citation_table.py`; table `tables/G15_literature_related_to_Kong2020.tsv` (41 records, each with DOI and my
one-line reading). Records relevant to the novelty claim (according to PubMed):

| PMID | DOI | Relation |
|---|---|---|
| 35764641 | [10.1038/s41467-022-31535-6](https://doi.org/10.1038/s41467-022-31535-6) | Same group (Kong et al. 2022, Nat Commun): network-based ML for immunotherapy response; no re-test of the organoid biomarkers |
| 41999762 | [10.1016/j.crmeth.2026.101411](https://doi.org/10.1016/j.crmeth.2026.101411) | Same group (Kim S lab, 2026, Cell Rep Methods): SubNetDL, mutation network propagation over 10 TCGA cancer-drug pairs; different method, no re-test |
| 36675023 | [10.3390/ijms24021508](https://doi.org/10.3390/ijms24021508) | New transcriptomic 5-FU resistance biomarkers; does not test the BH3-only pathway |
| 35805961 | [10.3390/ijms23136956](https://doi.org/10.3390/ijms23136956) | Biomarkers in drug-treated bladder organoids; does not test the amino-acid pathway |
| 33882314 | [10.1016/j.celrep.2021.109026](https://doi.org/10.1016/j.celrep.2021.109026) | Targeted RNA-seq screening in colorectal organoids; no patient re-test |
| 34154646 | [10.1186/s13073-021-00920-z](https://doi.org/10.1186/s13073-021-00920-z) | Network systems pharmacogenomics (medulloblastoma); related approach |
| 33709991 | [10.1097/DCR.0000000000001971](https://doi.org/10.1097/DCR.0000000000001971) | Blinded PDO chemotherapy-response study, stage IV colorectal (functional avatar) |
| 35696669 | [10.1021/acsbiomaterials.2c00354](https://doi.org/10.1021/acsbiomaterials.2c00354) | PDO chemotherapy-response prediction, colorectal (functional avatar) |

What I think: in PubMed I found no paper that re-tests the BH3-only (colon 5-FU) or amino-acid (bladder cisplatin)
biomarker in patients, and the two later papers from the same group move to other tasks. This supports the novelty
claim, but it is not a true citation search. One cited-by search (Google Scholar, Web of Science, or OpenAlex
`filter=cites:` from a machine with access) should be run before submission.

## G16. Joint organoid + patient model: paired bootstrap and extended lab-weight grid

**What was done.** (1) For each of the 11 units (network-shortlist features, the main analysis of
`nar/03_organoid_plus_patient`), I took the saved held-out (cross-fitted) patient scores and recomputed the pooled
C-index (survival) or AUROC (response) exactly as the original "POOLED" rows (concordant pairs summed within cohorts).
All 31 recomputed pooled values match the original tables (largest difference 1e-16;
`tables/G16_check_pooled_C_vs_original.tsv`). Then 2,000 paired bootstrap resamples (patients resampled within each
cohort; the same resample for all models). (2) I re-ran the joint models (`run_joint.run_unit`, same data, folds and
seed) with the lab weight grid extended to w = 0, 0.05, 0.1, 0.25, 0.5, 0.75 (w = weight of the lab-model loss; w = 0
means patients only), plus a model with w fixed at 0.1. Code `code/g16_joint_extended_grid.py`,
`code/g16_bootstrap_cindex.py`. Tables `tables/G16_paired_bootstrap_cindex_diff.tsv`,
`tables/G16_extended_grid_chosen_weight.tsv`. Figure `figures/G16_joint_model_bootstrap`.

Counts over the 20 unit x outcome rows (difference with 95% paired bootstrap CI):

| Comparison | grid | CI above 0 | CI below 0 | median difference |
|---|---|---|---|---|
| joint minus patient-only | original (0.25/0.5/0.75) | 0 | 2 (BLCA cisplatin response -0.044, -0.084 to -0.006; PAAD gemcitabine survival -0.056, -0.093 to -0.024) | -0.018 |
| joint minus organoid-only | original | 6 | 0 | 0.036 |
| patient-only minus organoid-only | original | 6 | 0 | 0.055 |
| joint minus patient-only | extended (0-0.75) | 0 | 0 | -0.014 |
| joint (w = 0.1) minus patient-only | extended | 0 | 1 | -0.003 |

Examples (original grid): CRC 5-FU survival joint 0.532 vs patient-only 0.546, difference -0.014 (-0.036 to 0.008);
ovarian paclitaxel survival 0.553 vs 0.558, -0.006 (-0.031 to 0.018). Joint minus organoid-only is above 0 for bladder
gemcitabine (survival and response), CRC 5-FU survival, ovarian paclitaxel (both) and ovarian platinum survival.
With the extended grid the inner CV chose w = 0 in 18 of 61 fits, 0.05 in 14, 0.1 in 12, 0.25 in 12, 0.5 in 1 and 0.75
in 4: 44 of 61 fits (72%) chose a weight below the old minimum of 0.25. The extended-grid joint model is within
-0.08 to +0.02 of patient-only in every unit, and no interval excludes 0.
Note: the re-run gives slightly different organoid-only / patient-only C for CRC 5-FU (0.503 / 0.536 vs 0.511 / 0.546
in the saved run) and small changes for CRC irinotecan and oxaliplatin; the other units are identical. The likely cause
is a change in the 03 code or data after its saved run. Comparisons are paired within each run, so they stay valid.

What I think: no unit shows a measurable gain from adding organoid data to the patient model, and when the weight may
go to zero the models mostly use very little organoid information. The gain over organoid-only models comes from the
patient data. A fair statement for the paper: organoids act at best as a weak prior; with these sample sizes their added
value cannot be detected (interval widths typically 0.02-0.2 in C-index).

## G17. The 57-analysis survey with progression-free interval (PFI)

**What was done.** Same patients, frozen scores and median splits as A3; outcome = TCGA-CDR PFI (days / 30.44).
Code `code/g17_pfi_survey.py`. Tables `tables/G17_PFI_survey_57.tsv`, `tables/G17_PFI_vs_OS_counts.tsv`,
`tables/G17_PFI_events_power.tsv`. Figure `figures/G17_PFI_vs_OS`.

| Family | p < 0.05 (of 57) | expected | binomial p | BH q < 0.05 | smallest p |
|---|---|---|---|---|---|
| OS treated | 7 (4 expected direction) | 2.85 | 0.023 | 0 | 0.0011 |
| PFI treated | 8 (1 expected direction) | 2.85 | 0.007 | 0 | 0.0030 |
| OS untreated | 12 | 2.85 | 0.00002 | 0 | 0.0012 |
| PFI untreated | 4 | 2.85 | 0.32 | 0 | 0.0030 |

PFI gives more events (median 29 vs 17 per analysis; median ratio 1.6) and a smaller minimal detectable HR (median 2.8
vs 3.9). Only 2 analyses are below 0.05 with both endpoints: the paper's colon 5-FU on the 2019 files (OS p = 0.014;
PFI p = 0.015, 39 events, HR 2.21, expected direction) and Shi pancreas gemcitabine (opposite direction both times). The
other 6 PFI hits (colon irinotecan and oxaliplatin, sarcoma docetaxel and gemcitabine) are all in the opposite direction
(HR 0.32-0.52). Smallest PFI BH q = 0.17. Spearman correlation of p-values OS vs PFI 0.18; of continuous log HR 0.35;
same direction in 36 of 57.

What I think: the endpoint changes which analyses look "significant", as expected for noisy small tests. The one result
that holds with both endpoints in the expected direction is the paper's colon 5-FU result on the 2019 files. It is worth
reporting, read together with the current-data colon 5-FU result (OS p = 0.866) and the G18 result below.

## G18. Agreement between ssGSEA and the within-sample rank score

**What was done.** (a) Spearman correlation between the ssGSEA scores (all TCGA caches of the organoid-only agents plus
the GI result files) and the rank score of the training table (`robust/results/training_table/pathway_scores_tcga.tsv.gz`)
for every pathway in every cohort. (b) The same for the frozen rank-1 pathway of each of the 57 analyses. (c) All 57
treated tests re-run with the rank score: same pathway, same sign of the Ridge coefficient, same median-split rule, same
patients and survival. Code `code/g18_score_agreement.py`. Tables `tables/G18_pathway_correlation_ssgsea_vs_rank.tsv`,
`tables/G18_correlation_summary_by_cohort.tsv`, `tables/G18_rerun_57_with_rank_score.tsv`,
`tables/G18_hit_counts_ssgsea_vs_rank.tsv`. Figure `figures/G18_score_agreement`.

- 5,688 pathway x cohort pairs: median rho 0.98 (10th percentile 0.91); 98% above 0.7. Lowest for the paper-era 2019
  colon files (median 0.89; older expression processing) and COADREAD/READ (0.88).
- Frozen rank-1 pathways: median rho 0.97 (range 0.56-1.00); median share of treated patients in the same median-split
  half 0.93 (range 0.62-1.00).
- Treated p < 0.05: ssGSEA 7 of 57, rank score 3 of 57 (all 3 among the 7). Untreated: 12 and 12.

| ssGSEA hit | rho | same half | p ssGSEA | p rank score | direction (rank score) |
|---|---|---|---|---|---|
| melanoma PDX dacarbazine | 0.97 | 0.90 | 0.0011 | 0.011 | opposite |
| lung cell lines pemetrexed | 0.99 | 0.93 | 0.0048 | 0.013 | opposite |
| lung cell lines gemcitabine | 1.00 | 1.00 | 0.030 | 0.030 | expected |
| paper colon 5-FU (2019 files) | 0.86 | 0.84 | 0.014 | 0.32 | expected |
| Shi pancreas gemcitabine | 0.99 | 0.97 | 0.023 | 0.068 | opposite |
| paper bladder cisplatin | 0.96 | 0.94 | 0.044 | 0.053 | expected |
| PDX NSCLC erlotinib | 0.96 | 0.90 | 0.040 | 0.22 | opposite |

What I think: the two scoring methods agree closely, so the choice of score does not explain differences between the
organoid-only survey and the later rank-score analyses; the global picture (3-7 small treated p-values, 12 untreated) is
the same. Results near p = 0.05 move with the method: for the paper's colon 5-FU split, 16% of patients change half and p
goes from 0.014 to 0.32 (8 deaths). One more sign that this result is fragile.

---

## Figures

| Figure | Content |
|---|---|
| `figures/G9_permutation_null` | shuffle null of the count of p < 0.05 (treated, untreated) and of the smallest p |
| `figures/G10_meta_forest` | forest plot of 57 HRs per SD with pooled estimate and prediction interval; subset pooled HRs with equivalence bounds |
| `figures/G12_baseline_treated_vs_untreated` | age, stage III-IV and purity of treated vs no-drug-record patients per cohort |
| `figures/G16_joint_model_bootstrap` | paired bootstrap CIs of C-index differences (joint vs patient-only / organoid-only; extended grid) |
| `figures/G17_PFI_vs_OS` | events, p-values and HRs with PFI vs OS |
| `figures/G18_score_agreement` | ssGSEA vs rank score correlation, re-run p-values, split agreement |

## What could not be done

- G15: a true cited-by list. OpenAlex, Crossref and Semantic Scholar are blocked by the network proxy here, and PubMed
  offers only similarity links. Substitute: PubMed similar-articles plus targeted searches (above).
- Housekeeping: `_data/g16_stale_unpatched/` holds one CRC oxaliplatin output from a first launch in which the extended
  grid was not applied (joblib workers lost the patch). It is not used anywhere; it can be deleted.
