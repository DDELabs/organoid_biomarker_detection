# Analysis plan: registered patient test of the frozen DepMap hypotheses (ANCHOR Stages 1–2)

Written on 2026-11-17, **before any patient outcome was looked at** in this folder. The code in `code/` follows this
plan. Anything done later that is not in this plan is labelled "post hoc" in the report.

## 0. What is being tested

- Input: `nar/08_depmap/tables/frozen_hypotheses.tsv` (frozen 2026-11-17T13:03:40Z). Expected sha256
  `a0bdd620293e2559667cf7f4c218894ae3e70cf1eb404558a084d28aee878ab7`. The code stops if the hash differs. The file
  is read only; it is never changed.
- 7 drugs, 5 signed Reactome pathways each: fluorouracil (5-FU), oxaliplatin, irinotecan, gemcitabine, paclitaxel,
  pemetrexed, sorafenib.
- Sign "+" means "a higher pathway score means a more resistant tumour"; "−" means "more sensitive".

## 1. The per-patient hypothesis score

1. Pathway scores: the 673 Reactome rank scores in `robust/results/training_table/` (TCGA and GEO). GSE109211
   (STORM, not in the table) is scored with the same function (`robust/obd/scoring.rank_score`, Reactome sets from
   the repo MSigDB v6.1 file, genes mapped with `to_canonical`, sets with ≥ 5 genes).
2. Rows with a `possible_duplicate_of` value are dropped.
3. Each pathway is standardised (mean 0, SD 1) **within cohort** (each TCGA project is its own cohort).
4. Hypothesis score = mean over the 5 pathways of (sign × standardised score). It is then standardised again within
   cohort, so every hazard ratio (HR) or odds ratio (OR) is **per 1 SD of score**.
5. A higher score = predicted more resistant to the drug.

## 2. Stage 1: pooled TCGA test (one test per drug)

**Exposure.** A TCGA drug record of systemic therapy (the same filter as the existing pipeline: `THERAPY` types,
supportive drugs removed) naming the drug. Capecitabine is not counted as 5-FU in the primary analysis (sensitivity
analysis below). TCGA drug records are rebuilt with the existing loader `robust/obd/tcga_response.py` (sources:
GDC biotab clinical_drug tables via the kemplab mirror, TCGA-CDR), written to `_data/tcga_response/`.

**Cancers included.** TCGA projects with ≥ 10 exposed patients who have pathway scores and survival data.

**Comparator.** All other patients of the included projects. Exposure is **time-varying**: a patient is counted as
unexposed from diagnosis until the first start day of the drug, and exposed from then on. So exposed time is
measured from drug start, and patients who start the drug later contribute their earlier time to the comparator.
A record with no start day is counted as exposed from day 0.

**Primary endpoint.** Progression-free interval (PFI, TCGA-CDR, days). **Secondary endpoint.** Overall survival (OS).

**Primary model** (Cox, Breslow ties, stratified by TCGA project, own Newton–Raphson code checked against
statsmodels `PHReg`):

  log h(t) = h0_project(t) + b1·score + b2·exposed(t) + b3·score×exposed(t) + age/10 + stage (I–IV, numeric) + male
  + ABSOLUTE purity (+ missing-indicators for stage and purity; values centred within project, missing set to 0).

- Test statistic: b3. **Interaction HR = exp(b3)** per SD.
- **Pre-specified direction: interaction HR > 1** (a higher resistance score goes with a relatively worse PFI under
  the drug than without it).
- p-value: two-sided Wald. **Multiple testing: Holm across the 7 primary tests (7 drugs × PFI)**, family-wise
  α = 0.05. OS and every sensitivity analysis are secondary and are not part of the Holm family (raw p plus Holm
  within each secondary family of 7 are reported).
- Score HR among exposed = exp(b1 + b3); among unexposed = exp(b1).

**Gate G2 (per drug), all must hold:**
1. Holm-adjusted p < 0.05 and interaction HR > 1;
2. the no-purity (S1) and IPW (S2) estimates keep HR > 1 with p < 0.05;
3. both negative controls (below) for that drug have p ≥ 0.05;
4. the score HR in unexposed patients is closer to 1 than in exposed patients.

**Sensitivity analyses (secondary):**
- S1 no purity adjustment;
- S2 inverse-probability weighting (IPW): propensity of ever being exposed from age, stage, sex, purity, project
  (logistic); stabilised weights truncated at the 1st/99th percentile; weighted version of the primary model;
  95% CI and p from 200 bootstrap resamples of patients (propensity re-fitted in each);
- S3 active comparator: comparator restricted to patients with a record of some other systemic drug;
- S4 primary model plus a proliferation score as an extra covariate and its interaction with exposure;
- S5 5-FU only: exposure = fluorouracil or capecitabine.

**Random-signature null.** 1,000 random signed sets of 5 pathways (drawn from the 673, random signs; seed 20261117;
the same 1,000 sets for every drug). Each is scored and tested with the primary model. Empirical one-sided p =
(1 + #{z_random ≥ z_observed}) / 1,001, where z = b3 / SE (positive = pre-specified direction).

**Negative controls** (same primary model, the control score in place of the hypothesis score):
- proliferation score = mean of standardised REACTOME_CELL_CYCLE_MITOTIC, REACTOME_DNA_REPLICATION,
  REACTOME_MITOTIC_M_M_G1_PHASES, REACTOME_G2_M_CHECKPOINTS, REACTOME_S_PHASE (re-standardised within project);
- tumour purity (ABSOLUTE, standardised within project; purity removed from covariates for this model; patients
  without purity dropped for this model only).
- Also reported (secondary): a 7 × 7 specificity matrix — every drug's frozen score tested against every drug's
  exposure.

**Per-cancer estimates.** The primary model fitted in each included project (no strata) when the project has ≥ 5
events in exposed time and ≥ 5 in unexposed time. DerSimonian–Laird random-effects meta-analysis of log
interaction HRs; I² and τ² reported. Forest plot per drug.

**KM curves (descriptive).** Score split at the within-project median. Exposed patients: PFI from drug start
(patients whose first start is before their PFI event or censoring). Comparator: never-exposed patients of the
included projects, PFI from diagnosis. Log-rank p for high vs low within each group.

**Power (each test).** D1 = PFI events in exposed time, D0 = events in unexposed time. Expected SE of b3 ≈
√(1/D1 + 1/D0). Reported: minimal detectable interaction HR (MDHR, 80% power) at α = 0.05 and at the Holm
first-step α = 0.05/7; power for an interaction HR of 1.25 (= 1/0.8) and 1.43 (= 1/0.7); and the same from the
observed SE.

## 3. Stage 2: external cohorts (registered, one-sided)

Each drug's frozen score (section 1) in every curated cohort with a comparator for that drug. **One-sided p in the
pre-specified direction** (time-to-event: interaction HR > 1; binary response: interaction OR < 1, because a higher
resistance score should lower the chance of response under the drug). Interaction term = score × (arm receives
the drug). Covariates: age, sex and stage where recorded (time-to-event models stratified by stage group where
stage is recorded; binary models adjust for stage group). A model is fitted only when each arm has ≥ 5 events
(time-to-event) or ≥ 5 responders and ≥ 5 non-responders (binary).

| Drug | Cohort | Contrast (drug arm vs comparator) | Endpoint (primary; secondary) | Randomised |
|---|---|---|---|---|
| 5-FU | GSE39582 | any 5-FU-based adjuvant (5FU, FUFOL, FOLFOX, FOLFIRI) vs no chemotherapy | RFS; OS | no |
| 5-FU | GSE103479 | 5-FU adjuvant vs surgery only | PFS; OS | no |
| 5-FU | GSE14333 | adjuvant chemotherapy (type not recorded, 5-FU-based in this setting) vs none | DFS | no |
| oxaliplatin | GSE72970 | FOLFOX(±bevacizumab) vs FOLFIRI(±bevacizumab/cetuximab) | PFS; OS; response | no |
| oxaliplatin | GSE39582 | FOLFOX vs 5FU/FUFOL | RFS; OS | no |
| oxaliplatin | GSE104645 | oxaliplatin-based (FOLFOX/SOX/XELOX ± biologic) vs irinotecan-based (FOLFIRI/IRIS ± biologic) | response | no |
| oxaliplatin | GSE87211 | 5-FU + oxaliplatin + RT vs 5-FU + RT | DFS; response | no |
| irinotecan | GSE72970 | FOLFIRI(±biologic) vs FOLFOX(±bevacizumab) | PFS; OS; response | no |
| irinotecan | GSE104645 | irinotecan-based vs oxaliplatin-based | response | no |
| irinotecan | GSE39582 | FOLFIRI vs 5FU/FUFOL | RFS | no |
| paclitaxel | GSE20271 | T/FAC vs FAC (randomised arm) | pCR response | yes |
| paclitaxel | GSE41998 | AC→paclitaxel vs AC→ixabepilone (randomised arm) | response | yes |
| paclitaxel | ISPY2 | HER2+ arms containing paclitaxel vs T-DM1 + pertuzumab (no paclitaxel) | pCR response | yes |
| paclitaxel | GSE9891 | platinum + paclitaxel vs platinum | RFS; OS | no |
| paclitaxel | GSE63885 | taxane/platinum vs platinum/cyclophosphamide | OS | no |
| paclitaxel | GSE15622 | paclitaxel vs carboplatin | response | no |
| sorafenib | GSE109211 (STORM) | sorafenib vs placebo | responder (recurrence-based, as deposited) | yes |
| gemcitabine, pemetrexed | GSE42127 | adjuvant platinum doublet (partner not recorded) vs observation | OS | no — exploratory only, not counted for G3 |

Single-arm cohorts (no comparator) are reported as **supportive only** (score main effect, one-sided; not counted
for G3): 5-FU — GSE87211 (5-FU+RT arm, DFS); oxaliplatin — GSE28702, GSE19860 (FOLFOX response); irinotecan —
GSE62080 (FOLFIRI response); paclitaxel — BRIGHTNESS paclitaxel-alone arm, ISPY2 paclitaxel-alone control arm,
GSE32646, GSE51373 (response), GSE25066 (RFS). Gemcitabine and pemetrexed have no curated cohort with the drug recorded.

**Gate G3 (per drug):** one-sided p < 0.025 in at least one controlled cohort whose power (from its observed
interaction SE) to detect a log effect of 0.357 (HR 1.43 or OR 0.7 per SD) is ≥ 0.80, **and** the same direction
in every other controlled cohort for that drug. Secondary: Holm across all controlled Stage-2 primary tests;
inverse-variance fixed-effect and random-effects pooling per drug (time-to-event and binary pooled separately).

**Power (each cohort).** Events (or responders/non-responders) per arm; power for a log effect of 0.357 at one-sided
α = 0.025 from the expected SE (√(1/D1 + 1/D0) for time-to-event; √(1/(n1·p1(1−p1)) + 1/(n0·p0(1−p0))) for
binary) and from the observed SE; minimal detectable effect at 80% power.

## 4. Rules

- Seeds: 20261117. At most 2 CPU processes.
- Every number in the report comes from a table written by the code.
- Results are reported with HR/OR, 95% CI, p, Holm p, n and events, whatever the direction.
