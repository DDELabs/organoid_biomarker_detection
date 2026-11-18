# Registered patient test of the frozen DepMap hypotheses (ANCHOR Stages 1–2)

Folder `nar/09_registered_test/`. Code in `code/` (run in order `s00`–`s04`), tables in `tables/`, figures in
`figures/` (PNG 300 dpi + PDF + plain caption). Seed 20261117. At most 2 CPU processes were used.

## In short

- The 35 frozen pathways (7 drugs × 5) were tested **exactly as frozen**. The file's sha256 is
  `a0bdd620293e2559667cf7f4c218894ae3e70cf1eb404558a084d28aee878ab7`, the same as `nar/08_depmap/tables/frozen_hypotheses.sha256`.
  The code checks this hash at every run and stops if it differs. The file was never edited.
- The analysis plan (`ANALYSIS_PLAN.md`, sha256 `97b2a634…1ef1967`, written 2026-11-17T13:11:00Z) was written before
  any outcome was opened. Everything below follows it. Post-hoc remarks are marked as such. Wording note (audit I16):
  this is an **internally time-stamped analysis plan**, not an external registration. The proof of timing is the sha256
  hashes and file times only; the TCGA survival cache was built at 13:08, before the plan (13:11), and the plan was
  committed to git after the results (13:22). Before any further (Stage 3) test, the plan and its hash should be
  deposited publicly (OSF or Zenodo).
- **Stage 1 (TCGA, 7 primary tests):** no drug's score showed the predicted score × drug interaction. All 7
  Holm-adjusted p-values are 1.0. The closest was gemcitabine (HR 1.12, 95% CI 0.94–1.32, p = 0.199; random-signature
  p = 0.070). **Gate G2: 0 of 7 drugs pass.**
- **Stage 2 (external cohorts, one-sided):** 17 primary controlled tests in 13 cohorts. One had one-sided p < 0.025
  (oxaliplatin, GSE72970 FOLFOX vs FOLFIRI, PFS interaction HR 1.63, 95% CI 1.00–2.64, p = 0.024). It had only 30%
  power, and the other oxaliplatin cohorts did not go in the same direction. **Gate G3: 0 of 7 drugs pass.**
- Power was the main limit. Only paclitaxel, gemcitabine and 5-FU had ≥ 80% power in TCGA for an interaction HR of
  1.43 per SD. No external cohort had ≥ 80% power for that effect (the best was 54%).

What I think: with a clean, internally time-stamped design, these cell-line-derived pathway hypotheses do not show
drug-specific predictive value in patients. For paclitaxel, gemcitabine and 5-FU this is an informative estimate: the
confidence intervals exclude interaction HRs above about 1.13–1.32 per SD in TCGA. For irinotecan, sorafenib,
oxaliplatin and pemetrexed the patient data cannot tell either way.

## Data used (`tables/data_versions.tsv`)

| File | Bytes | sha256 (first 16) |
|---|---|---|
| nar/08_depmap/tables/frozen_hypotheses.tsv | 11,626 | a0bdd620293e2559 |
| robust/results/training_table/patients.tsv.gz | 345,036 | a68da5ad8a0aa0ae |
| robust/results/training_table/pathway_scores_tcga.tsv.gz | 11,150,748 | 0949e9fbe202dfd3 |
| robust/results/training_table/pathway_scores_geo.tsv.gz | 11,879,734 | ce3e32a577cd2999 |
| robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt (ABSOLUTE purity) | 901,812 | f430a975433d82e0 |
| _data/tcga_response/drug_response_long.tsv.gz (TCGA drug records, rebuilt) | 185,454 | 10695990742b1020 |
| _data/tcga_response/survival_cdr.tsv.gz (TCGA-CDR) | 175,611 | 4017d88e5ed0e509 |
| data/curated/trials/GSE109211/expression.tsv.gz (STORM) | 7,313,144 | 0f577406d247f2f3 |

The shared TCGA drug-record cache (`data/external/tcga_response/`) was missing. I rebuilt it with the existing loader
`robust/obd/tcga_response.py` (same sources: TCGA biotab drug tables and TCGA-CDR), writing only into `_data/`
(`code/s00_fetch_tcga_drug.py`). Rows flagged `possible_duplicate_of` were dropped. Pathway scores were standardised
within each cohort (each TCGA project is one cohort). STORM (GSE109211) is not in the training table, so I scored it with
the same rank-score function. The Cox code (`code/common.py`) gives the same coefficients and standard errors as
statsmodels `PHReg` to 8 decimals (`code/t_cox_check.py`).

## How the score works

For each patient: score = mean of (sign × standardised pathway score) over the drug's 5 frozen pathways, then
standardised again within cohort. A higher score means "predicted more resistant". All effects are per 1 standard
deviation (SD) of score.

## Stage 1: TCGA, pooled per drug (Figures 2–6)

Design (from the plan): Cox model stratified by cancer type (each cancer has its own baseline risk). Drug exposure is
time-varying: a patient counts as "on the drug" from the first drug start day. The model includes score, exposure,
score × exposure, age, stage, sex and tumour purity. **The interaction HR is the test.** HR > 1 = the predicted
direction. Primary endpoint: PFI (progression-free interval). Holm correction across the 7 drugs.

### Primary results (`tables/stage1_primary.tsv`, Figure 2a)

| Drug | Cancers | n | Patients exposed before progression | PFI events on / off drug | Interaction HR (95% CI) | p | Holm p | Random-set p (one-sided) |
|---|---|---|---|---|---|---|---|---|
| 5-FU | BRCA, COAD, ESCA, PAAD, READ, STAD | 2,427 | 406 | 121 / 500 | 0.93 (0.76–1.13) | 0.471 | 1.0 | 0.717 |
| Oxaliplatin | COAD, PAAD, READ, STAD | 1,163 | 159 | 48 / 341 | 0.93 (0.67–1.29) | 0.667 | 1.0 | 0.636 |
| Irinotecan | COAD, GBM, LGG, PAAD | 1,281 | 23 | 14 / 520 | 0.64 (0.32–1.26) | 0.197 | 1.0 | 0.888 |
| Gemcitabine | BLCA, LUAD, LUSC, OV, PAAD, SARC | 2,258 | 248 | 154 / 916 | 1.12 (0.94–1.32) | 0.199 | 1.0 | 0.070 |
| Paclitaxel | BLCA, BRCA, CESC, HNSC, LUAD, LUSC, OV, UCEC, UCS | 4,311 | 838 | 405 / 997 | 1.00 (0.89–1.13) | 0.988 | 1.0 | 0.507 |
| Pemetrexed | LUAD, MESO | 586 | 78 | 40 / 225 | 0.98 (0.70–1.37) | 0.899 | 1.0 | 0.544 |
| Sorafenib | KIRC, LIHC | 895 | 8 | 6 / 333 | 0.82 (0.45–1.50) | 0.520 | 1.0 | 0.681 |

What I think: none of the 7 tests point the predicted way with any strength. Paclitaxel is the best-powered test
(405 events on drug). Its estimate is exactly 1.00, with a tight interval (0.89–1.13). Gemcitabine leans the predicted
way (HR 1.12), and it beats 93% of random 5-pathway sets, but p = 0.199. Irinotecan and sorafenib in TCGA are mostly
given after progression. So only 23 and 8 patients had the drug before a PFI event, and their tests carry very little
information.

**Secondary endpoint OS** (Holm within the OS family; all Holm p = 1.0): 5-FU HR 0.95 (0.78–1.16, p = 0.609,
events 117/482); oxaliplatin 0.80 (0.59–1.11, p = 0.183, 44/327); irinotecan 1.01 (0.71–1.43, p = 0.958, 52/386);
gemcitabine 0.99 (0.86–1.15, p = 0.919, 223/800); paclitaxel 1.01 (0.89–1.15, p = 0.868, 358/1,042); pemetrexed
0.87 (0.62–1.22, p = 0.417, 42/211); sorafenib 0.75 (0.50–1.12, p = 0.162, 26/278).
What I think: OS agrees with PFI; nothing points the predicted way.

**Score effect in each group (PFI).** Score HR on drug vs off drug: 5-FU 1.03 vs 1.11; oxaliplatin 1.08 vs 1.16;
irinotecan 0.56 vs 0.88; gemcitabine 1.05 vs 0.94; paclitaxel 0.99 vs 0.99; pemetrexed 1.13 vs 1.16; sorafenib 0.78
vs 0.95. What I think: where the score relates to outcome at all (for example 5-FU and oxaliplatin, HR 1.11–1.16
off drug), it does so just as much without the drug. That points to a general (prognostic) signal, not a drug-specific one.

### Sensitivity analyses (`tables/stage1_sensitivity.tsv`, Figure 2b)

| Drug | No purity | IPW (bootstrap 200) | Active comparator | + proliferation | 5-FU or capecitabine |
|---|---|---|---|---|---|
| 5-FU | 0.93 (0.76–1.13), p 0.458 | 0.91 (0.72–1.16), p 0.417 | 0.95 (0.76–1.18), p 0.627 | 0.94 (0.77–1.15), p 0.532 | 0.98 (0.83–1.17), p 0.861 (159/659 ev) |
| Oxaliplatin | 0.93 (0.68–1.29), p 0.676 | 0.96 (0.62–1.41), p 0.845 | 0.95 (0.67–1.33), p 0.751 | 0.93 (0.67–1.28), p 0.643 | – |
| Irinotecan | 0.66 (0.33–1.31), p 0.232 | 0.51 (0.12–1.34), p 0.278 | 0.82 (0.42–1.62), p 0.571 | 0.69 (0.30–1.60), p 0.387 | – |
| Gemcitabine | 1.11 (0.94–1.31), p 0.211 | 1.19 (0.97–1.42), p 0.071 | 1.04 (0.87–1.25), p 0.682 | **1.49 (1.12–1.99), p 0.006** (Holm within family 0.043) | – |
| Paclitaxel | 1.00 (0.88–1.12), p 0.937 | 0.90 (0.74–1.12), p 0.313 | 0.89 (0.78–1.03), p 0.111 | 0.80 (0.63–1.02), p 0.071 | – |
| Pemetrexed | 0.94 (0.67–1.32), p 0.713 | 0.91 (0.58–1.29), p 0.637 | 0.82 (0.54–1.25), p 0.356 | 0.98 (0.69–1.38), p 0.886 | – |
| Sorafenib | 0.82 (0.45–1.50), p 0.515 | 0.82 (0.19–2.89), p NA* | 0.84 (0.43–1.65), p 0.619 | 1.40 (0.47–4.21), p 0.546 | – |

\*Sorafenib IPW: with 6 events on the drug, some bootstrap fits were unstable. That inflated the bootstrap SE
(z ≈ 0, which printed as p = 1.000). The percentile interval is shown and the p-value is reported as NA (audit I20);
`tables/stage1_sensitivity.tsv` still holds the raw bootstrap-SE value. IPW weights were well behaved: range 0.18–5.56, effective sample size
47–99% of n (`tables/stage1_ipw_weights.tsv`).

What I think: no adjustment changes the picture. Gemcitabine is the one to note. When proliferation and its own
interaction with gemcitabine are added, the gemcitabine score interaction rises to HR 1.49 (p = 0.006). This is a
secondary analysis, so it is not a G2 pass. It suggests the gemcitabine score carries a signal that is partly hidden
by proliferation. That is a hypothesis for a new pre-registered test, not a finding.

### Random-signature null (`tables/stage1_random_null_z.tsv.gz`, Figure 4)

1,000 random signed sets of 5 pathways per drug, same model. One-sided empirical p: 5-FU 0.717, oxaliplatin 0.636,
irinotecan 0.888, gemcitabine 0.070, paclitaxel 0.507, pemetrexed 0.544, sorafenib 0.681.
What I think: apart from gemcitabine (top 7%), the frozen scores sit in the middle of what random pathway sets give.

### Negative controls (`tables/stage1_negative_controls.tsv`, Figure 5a)

| Drug | Proliferation × drug HR (95% CI), p | Purity × drug HR (95% CI), p |
|---|---|---|
| 5-FU | 0.83 (0.67–1.02), 0.080 | 1.13 (0.92–1.39), 0.232 |
| Oxaliplatin | 1.04 (0.75–1.44), 0.824 | 1.14 (0.84–1.53), 0.408 |
| Irinotecan | 1.24 (0.72–2.14), 0.446 | 1.98 (1.01–3.90), 0.046 |
| Gemcitabine | 1.02 (0.86–1.21), 0.840 | 1.03 (0.85–1.26), 0.736 |
| Paclitaxel | 0.94 (0.84–1.06), 0.336 | 0.92 (0.81–1.03), 0.152 |
| Pemetrexed | 1.48 (0.99–2.22), 0.056 | 1.11 (0.77–1.60), 0.580 |
| Sorafenib | 1.99 (0.90–4.41), 0.089 | 0.95 (0.43–2.06), 0.888 |

What I think: the controls behave as expected (13 of 14 have p ≥ 0.05). The irinotecan purity control reached
p = 0.046 on 14 events, which fits chance (Holm within the family: 0.325).

**Specificity matrix** (`tables/stage1_specificity_matrix.tsv`, Figure 5b): each drug's score was tested against
every drug's exposure (49 tests). The diagonal (pre-registered pairs) z-values are −0.72, −0.43, −1.29, 1.29, 0.02,
−0.13, −0.64. The largest off-diagonal value is the gemcitabine score against 5-FU exposure (z = 2.35).
What I think: no score is specific to its own drug. With 42 off-diagonal tests, one |z| > 2.3 is about what chance gives.

### Per-cancer and meta-analysis (`tables/stage1_per_cancer.tsv`, `stage1_meta_analysis.tsv`, Figure 6)

Random-effects interaction HR (k cancers): 5-FU 0.91 (0.73–1.14, p 0.422, k 6, I² 8%); oxaliplatin 1.16
(0.80–1.67, p 0.436, k 2); irinotecan 1.02 (0.41–2.58, k 1, COAD only); gemcitabine 1.11 (0.93–1.33, p 0.265, k 6,
I² 0%); paclitaxel 0.92 (0.71–1.20, p 0.527, k 9, I² 53%); pemetrexed 0.99 (0.43–2.30, k 2, I² 78%); sorafenib 1.23
(0.58–2.58, k 1, KIRC only). Single-cancer values with p < 0.05: paclitaxel LUAD 2.44 (1.25–4.78, p 0.009, 16 events on
drug) and paclitaxel CESC 0.45 (0.22–0.89, p 0.021, 6 events). They point in opposite directions.
What I think: the cancers agree with the pooled result. Paclitaxel shows moderate spread between cancers (I² 53%),
driven by small strata that point in opposite directions. Bladder gemcitabine (1.35, 0.95–1.93, p 0.099, 31 events)
is the largest gemcitabine stratum leaning the predicted way.

### KM curves (`tables/stage1_km_logrank.tsv`, Figure 3, descriptive)

High vs low score (log-rank p), on drug from drug start / comparator from diagnosis: 5-FU 0.233 / 0.075; oxaliplatin
0.737 / 0.103; irinotecan 0.004 (n = 23, with high score doing *better*) / 0.210; gemcitabine 0.932 / 0.241;
paclitaxel 0.891 / 0.326; pemetrexed 0.655 / 0.034; sorafenib 0.406 / 0.137.
What I think: the curves show no consistent separation in the treated groups. The irinotecan split rests on 23
patients and goes against the hypothesis.

### Power (`tables/stage1_power.tsv`)

| Drug | Events on / off | Minimal detectable interaction HR (80% power, α 0.05) | at Holm α 0.05/7 | Power for HR 1.25 | Power for HR 1.43 |
|---|---|---|---|---|---|
| 5-FU | 121 / 500 | 1.33 | 1.43 | 0.60 | 0.94 |
| Oxaliplatin | 48 / 341 | 1.54 | 1.72 | 0.30 | 0.64 |
| Irinotecan | 14 / 520 | 2.14 | 2.60 | 0.13 | 0.26 |
| Gemcitabine | 154 / 916 | 1.28 | 1.36 | 0.73 | 0.98 |
| Paclitaxel | 405 / 997 | 1.18 | 1.23 | 0.97 | 1.00 |
| Pemetrexed | 40 / 225 | 1.62 | 1.83 | 0.26 | 0.55 |
| Sorafenib | 6 / 333 | 3.17 | 4.28 | 0.08 | 0.14 |

(Expected from events; the observed-SE versions are in the table.) What I think: paclitaxel, gemcitabine and 5-FU
were adequately powered for a moderate effect (HR 1.43), and paclitaxel even for a small one (HR 1.25). Irinotecan and
sorafenib were not testable in TCGA in any useful sense.

### Gate G2 (`tables/stage1_gate_G2.tsv`)

No drug meets condition 1 (Holm p < 0.05 with HR > 1), so **no drug passes G2**. Conditions 3 (negative controls null)
and 4 (weaker effect in comparators) hold for some drugs; the table lists each condition.

## Stage 2: external cohorts (Figures 7–8)

One-sided tests in the frozen direction: interaction HR > 1 for survival, and interaction OR < 1 for response.
Interaction = score × being in the drug arm. Covariates: age, sex and stage where recorded.

### Controlled cohorts, primary tests (`tables/stage2_controlled.tsv`)

| Drug | Cohort (contrast) | n | Events or R/NR (drug / comparator) | Effect (95% CI) | One-sided p | Holm (17 tests) | Power* |
|---|---|---|---|---|---|---|---|
| 5-FU | GSE39582 (5-FU-based vs none), RFS | 454 | 49 / 66 | HR 1.06 (0.74–1.52) | 0.367 | 1.0 | 0.49 |
| 5-FU | GSE103479 (5-FU vs surgery), PFS | 150 | 23 / 18 | HR 0.94 (0.50–1.78) | 0.572 | 1.0 | 0.19 |
| 5-FU | GSE14333 (chemo vs none), DFS | 226 | 27 / 23 | HR 0.84 (0.47–1.52) | 0.714 | 1.0 | 0.22 |
| Oxaliplatin | GSE72970 (FOLFOX vs FOLFIRI), PFS | 119 | 31 / 79 | **HR 1.63 (1.00–2.64)** | **0.024** | 0.386 | 0.30 |
| Oxaliplatin | GSE39582 (FOLFOX vs 5FU/FUFOL), RFS | 148 | 7 / 42 | HR 0.63 (0.24–1.64) | 0.830 | 1.0 | 0.11 |
| Oxaliplatin | GSE104645 (ox- vs iri-based), response | 167 | 66R/47NR vs 38R/16NR | OR 1.67 (0.82–3.43) | 0.920 | 1.0 | 0.16 |
| Oxaliplatin | GSE87211 (5-FU+ox+RT vs 5-FU+RT), DFS | 171 | 13 / 21 | HR 0.83 (0.39–1.79) | 0.682 | 1.0 | 0.15 |
| Irinotecan | GSE72970 (FOLFIRI vs FOLFOX), PFS | 119 | 79 / 31 | HR 1.26 (0.82–1.93) | 0.148 | 1.0 | 0.37 |
| Irinotecan | GSE104645 (iri- vs ox-based), response | 167 | 38R/16NR vs 66R/47NR | OR 1.06 (0.52–2.16) | 0.567 | 1.0 | 0.17 |
| Irinotecan | GSE39582 (FOLFIRI vs 5FU/FUFOL), RFS | 137 | 0 / 42 | not fitted (< 5 events in FOLFIRI arm) | – | – | – |
| Paclitaxel | GSE20271 (T/FAC vs FAC, randomised), pCR | 178 | 18R/66NR vs 8R/86NR | OR 0.95 (0.37–2.45) | 0.459 | 1.0 | 0.11 |
| Paclitaxel | GSE41998 (AC→T vs AC→ixabepilone, randomised) | 253 | 34R/87NR vs 35R/97NR | OR 1.01 (0.54–1.87) | 0.507 | 1.0 | 0.20 |
| Paclitaxel | I-SPY2 HER2+ (T-containing vs T-DM1+P, randomised) | 180 | 57R/71NR vs 30R/22NR | OR 1.44 (0.69–3.01) | 0.832 | 1.0 | 0.15 |
| Paclitaxel | GSE9891 (platinum+T vs platinum), RFS | 238 | 138 / 41 | HR 0.65 (0.42–1.00) | 0.976 | 1.0 | 0.37 |
| Paclitaxel | GSE63885 (T/platinum vs platinum/cyclo), OS | 75 | 36 / 30 | HR 1.25 (0.76–2.07) | 0.188 | 1.0 | 0.29 |
| Paclitaxel | GSE15622 (paclitaxel vs carboplatin), response | 35 | 13R/7NR vs 9R/6NR | OR 10.0 (1.27–78.5) | 0.986 | 1.0 | 0.05 |
| Sorafenib | GSE109211 STORM (sorafenib vs placebo, randomised) | 139 | 21R/46NR vs 21R/51NR | OR 1.20 (0.45–3.23) | 0.641 | 1.0 | 0.11 |

\*Power (from the observed SE) to detect HR 1.43 or OR 0.7 per SD at one-sided α 0.025. The smallest detectable
effects (80% power) range from 1.67 (GSE39582) to 19 (GSE15622) on the HR/1-over-OR scale (`min_detectable_obs`).

Secondary endpoints: 5-FU GSE39582 OS 1.07 (0.76–1.49, p 0.356); GSE103479 OS 1.03 (0.57–1.86, p 0.457); oxaliplatin
GSE72970 OS 1.28 (0.74–2.21, p 0.184) and response OR 0.49 (0.19–1.23, p 0.065); GSE39582 OS 0.35 (0.12–0.96,
p 0.979); GSE87211 response OR 1.31 (0.56–3.03, p 0.734); irinotecan GSE72970 OS 0.91 (0.56–1.49, p 0.642) and response
OR 0.58 (0.24–1.39, p 0.110); paclitaxel GSE9891 OS 0.61 (0.38–0.98, p 0.981). Exploratory GSE42127 (platinum doublet,
partner not recorded, vs observation, OS, 9/52 events): gemcitabine score HR 0.43 (0.15–1.19, p 0.948); pemetrexed
score 0.63 (0.33–1.23, p 0.911).

What I think: the oxaliplatin score in GSE72970 is the one positive-direction signal. FOLFOX patients with a high
score had relatively shorter PFS than FOLFIRI patients (HR 1.63), and OS and response lean the same way. But it does not
survive Holm across the 17 tests (0.386), and two other oxaliplatin comparisons (GSE39582, GSE87211) lean the opposite
way. GSE72970 is also the cohort behind the earlier "COAD oxaliplatin" result, so it is not independent of past work
on this project (the frozen pathways themselves were chosen without patient data). The paclitaxel ovarian cohort
GSE9891 leans clearly against the hypothesis (HR 0.65, upper limit 1.00). The tiny GSE15622 gives an unstable OR of 10.

**Pooled per drug** (`tables/stage2_pooled.tsv`, one-sided p, fixed effect): 5-FU survival HR 0.99 (0.75–1.30,
p 0.533, 3 cohorts); oxaliplatin survival 1.20 (0.82–1.74, p 0.174, 3 cohorts, I² 52%); paclitaxel response OR 1.22
(0.81–1.85, p 0.828, 4 cohorts); paclitaxel survival HR 0.86 (0.62–1.19, p 0.820, 2 cohorts, I² 74%). Single-cohort
pools for irinotecan and sorafenib are the same as their single cohorts above.

### Single-arm cohorts, supportive only (`tables/stage2_single_arm_supportive.tsv`)

| Drug | Cohort | n | Events or R/NR | Effect (95% CI) | One-sided p |
|---|---|---|---|---|---|
| 5-FU | GSE87211 5-FU+RT arm, DFS | 95 | 21 ev | HR 0.90 (0.61–1.32) | 0.702 |
| Oxaliplatin | GSE28702 FOLFOX, response | 83 | 42/41 | OR 0.85 (0.55–1.32) | 0.236 |
| Oxaliplatin | GSE19860 FOLFOX, response | 40 | 15/25 | OR 1.78 (0.92–3.45) | 0.956 |
| Irinotecan | GSE62080 FOLFIRI, response | 21 | 9/12 | OR 1.03 (0.42–2.48) | 0.521 |
| Paclitaxel | BrighTNess paclitaxel-alone arm, pCR | 123 | 40/83 | OR 0.67 (0.44–1.01) | 0.029 |
| Paclitaxel | I-SPY2 paclitaxel control arm, pCR | 179 | 31/148 | OR 0.62 (0.40–0.94) | 0.013 |
| Paclitaxel | GSE32646 (P→FEC), pCR | 115 | 27/88 | OR 0.44 (0.25–0.77) | 0.002 |
| Paclitaxel | GSE51373, response | 28 | 16/12 | OR 0.50 (0.21–1.19) | 0.058 |
| Paclitaxel | GSE25066 (T + anthracycline), RFS | 457 | 102 ev | HR 0.95 (0.78–1.15) | 0.698 |

What I think (post hoc): in three breast cohorts given paclitaxel, a high paclitaxel score goes with less pCR. That is
the frozen direction. But the same score shows no interaction in the randomised paclitaxel-vs-no-paclitaxel
comparisons (GSE20271, I-SPY2 HER2+). Two of the paclitaxel pathways are cell-cycle pathways with a "−" sign, so the
score is partly "low proliferation". Low proliferation predicts less pCR after any breast chemotherapy. These
single-arm associations are therefore most likely general chemosensitivity, not paclitaxel-specific.

### Gate G3 (`tables/stage2_gate_G3.tsv`)

| Drug | Controlled primary tests | Adequately powered (≥ 0.8) | One-sided p < 0.025 | Same direction | G3 |
|---|---|---|---|---|---|
| 5-FU | 3 | 0 | 0 | 1 of 3 | no |
| Oxaliplatin | 4 | 0 | 1 (GSE72970) | 1 of 4 | no |
| Irinotecan | 2 | 0 | 0 | 1 of 2 | no |
| Paclitaxel | 6 | 0 | 0 | 1 of 6 | no |
| Sorafenib | 1 | 0 | 0 | 0 of 1 | no |
| Gemcitabine, pemetrexed | 0 (only the exploratory GSE42127) | – | – | – | not testable |

What I think: no available cohort could have passed G3 on its own, because none reached 80% power for a moderate
effect. A real test for G3 needs pooled randomised trial data (for example several FOLFOX-vs-FOLFIRI trials for
oxaliplatin and irinotecan) or new data (Stage 3).

## Figures

| Figure | Content |
|---|---|
| `fig1_preregistration_flow` | flow of the pre-registration, from frozen hypotheses to gates |
| `fig2_stage1_pooled_interaction_forest` | per-drug pooled interaction HR (PFI, OS) and sensitivity analyses |
| `fig3_km_high_vs_low_score` | KM, high vs low score, on drug vs comparator, each drug |
| `fig4_random_signature_null` | 1,000 random-set histograms with the observed z |
| `fig5_negative_controls` | proliferation/purity controls; 7 × 7 specificity matrix |
| `fig6_per_cancer_forest_meta` | per-cancer forest plots with random-effects meta-analysis |
| `fig7_external_cohort_forest` | external controlled and single-arm cohorts |
| `fig8_summary_heatmap` | drug × cohort summary of z-values in the predicted direction |

## Limits and deviations

- TCGA treatment is not randomised. The interaction is still open to confounding by indication, even with purity
  adjustment, IPW and an active comparator.
- In TCGA, irinotecan and sorafenib are mostly given after progression. Under the pre-specified PFI design (exposure
  must come before the event), only 23 and 8 patients counted as exposed.
- Drug records without a start day are counted as exposed from day 0 (1–6% of patients per drug). Sensitivity
  (audit I19, `tables/stage1_sensitivity_no_missing_start.tsv`, made by `nar/15_fixes/code/f09_sensitivity.py`):
  dropping these patients changes the primary PFI interaction HRs only slightly: 5-FU 0.93 → 0.93 (7 dropped),
  oxaliplatin 0.93 → 0.91 (2), irinotecan unchanged (0), gemcitabine 1.12 (p 0.199) → 1.09 (p 0.326) (12),
  paclitaxel 1.00 → 1.01 (59), pemetrexed 0.98 → 0.95 (4), sorafenib 0.82 → 0.71 (1). No drug changes its conclusion.
- GSE39582 treated group (audit I01). The primary Stage 2 test (frozen) uses 5FU/FUFOL/FOLFOX/FOLFIRI at all stages
  (156 vs 298 with RFS). The project-wide shared definition (`nar/15_fixes/definitions/GSE39582_treatment.tsv`:
  5FU/FUFOL/FOLFOX, stage II–III, unknown regimen excluded) gives, as a sensitivity analysis
  (`tables/stage2_sensitivity_GSE39582_harmonised.tsv`): RFS interaction HR 1.17 (0.80–1.71), one-sided p = 0.21,
  137 vs 252 patients (47 vs 61 events); OS 1.20 (0.81–1.80), p = 0.18. The primary was 1.06 (0.74–1.52), p = 0.37.
  What I think: the harmonised definition moves the estimate a little towards the expected direction, but it stays
  far from p < 0.025.
- GSE14333 "adjuvant chemotherapy" has no regimen recorded. It was assumed to be 5-FU-based, as stated in the plan.
- The only change after the plan was cosmetic (figure layout), plus one numerical fix. The first Stage 1 run crashed
  in the sorafenib random-set null. By then it had already printed the primary PFI estimates, so I had seen them. I
  then added step-halving and NaN handling to the Cox fit and re-ran everything. The primary estimates were identical
  to the first run. All 1,000 null fits were valid for every drug (`random_null_n_valid` in `stage1_primary.tsv`).
