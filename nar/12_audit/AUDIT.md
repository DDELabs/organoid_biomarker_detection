# Statistics and consistency audit (before NAR Cancer submission)

Folder `nar/12_audit/`. Auditor: statistics and consistency agent. Date 2026-11-18.
Scope: `nar/02`, `03`, `05`, `07` (code and finished outputs only; still running), `08`, `09`, `10`.
Every issue is in `issues.tsv` (id, severity, evidence, exact fix, owner). No other agent's folder was edited.

## Summary

- **24 issues: 4 high, 13 medium, 7 low.**
- The core statistics are sound. Log-rank uses statsmodels `survdiff` and Cox uses statsmodels `PHReg` (02, 05).
  The custom Cox in 09 matches `PHReg`, with entry times and strata, to 8 decimals (I re-ran `t_cox_check.py`).
  The BH and Holm codes in 05 and 09 match `statsmodels.multipletests` on 200 random vectors. A3's q-values and
  binomial p-values recompute exactly from the table.
- The headline numbers agree across reports and source tables. Examples: colon 5-FU n = 114, p = 0.01435;
  untreated p = 0.1605; bladder p = 0.0435; A3 7/57 and 12/57; B1 196 tests; navitoclax rho; Stage 1 HRs.
  The differences that do exist come from different **patient definitions**, not from arithmetic.
- The biggest risks for referees are these:
  1. Inconsistent treated-group definitions (GSE39582; GSE72970 FOLFIRINOX patients).
  2. Wording in the 02 report text that contradicts its own numbers.
  3. Figures that are too large to print legibly.
  4. Different direction conventions (03 versus everything else).

## Top 5

1. **I01 (high), GSE39582 "5-FU-treated" is defined three ways.**
   - 03: 200 patients. This includes 66 "chemo type not given" and FOLFIRI patients; stage IV is removed.
   - 05 B1: 202 patients.
   - 09: 156 patients. Stage I-IV are included; unknown regimens are not.
   - Fix: choose one definition (5FU/FUFOL/FOLFOX, stage II-III, unknown regimen excluded). Re-run 03 and 05 with it.
     09's primary analysis is frozen, so keep it and add the harmonised version as a sensitivity analysis.
2. **I02 (high), GSE72970 FOLFIRINOX patients sit in both arms.** Four patients got FOLFIRINOX (both oxaliplatin and
   irinotecan).
   - 03 and 05 count them as FOLFOX in the oxaliplatin unit (40 vs 84) and as FOLFIRI in the irinotecan unit (88 vs 36).
   - 09 excludes them (36 vs 83). 09 is correct.
3. **I03 (high), 02 auto-text contradicts its numbers.** For COAD, COADREAD and Shi pancreas oxaliplatin, the best of
   30 settings beats 99.2%, 99.3% and 97.3% of shuffles. The text still says "picking the best setting is not enough".
   Those best settings point the opposite way (for example COAD k = 8, HR 0.31). This needs a direction-aware sentence.
4. **I04 (high), main figures are not journal-sized.** Examples: 02 colon panels are 17 × 31 in; BLCA_full is
   16 × 19.5 in; B1_Fig5 is 9.5 × 18 in. Some titles are truncated or overlap the panel letters.
5. **I05 / I06 (medium), 03 conventions and alignment.**
   - 03 uses "HR < 1 = expected". Every other folder uses "HR > 1 = expected".
   - 03's label-free alignment uses the held-out cohort's expression. This is transductive, but no outcomes are used.

## 1. Number consistency (gap G11 included)

`tables/MASTER_NUMBERS.tsv` gives one authoritative value per headline number. Each value is read by code
(`code/a02_master_numbers.py`, `code/a03_reconcile.py`) directly from the owning agent's table, and the source path is
given. The `reports_quoting` column is only a loose substring search, so treat it as a pointer. The rows starting
`G11_` reconcile the numbers that differ between agents:

| Number | Values seen | Correct value and reason |
|---|---|---|
| BRAF V600E vs cetuximab model, current data | 0.00072 (02) vs 0.00048 (10) | Both are correct but use different cohorts. 02 uses TCGA-COAD only (50 vs 377), the same design as the paper, so quote it. 10 uses COAD+READ from Xena calls (51 vs 522). |
| Colon 5-FU, current data | 0.866 (02) vs 0.62 (`docs/PROJECT_SUMMARY.md`) | Use 0.866 (n = 129, 25 deaths; the rank-1 pathway changes to purine biosynthesis). The 0.62 comes from the earlier `robust/` run, and no table row reproduces it. Correct the docs. |
| Broutier liver sorafenib | 0.139 (02, 05) vs 0.089 (10) | Use 0.139. 10 used the Ji gene-universe bundle, because B3 built no Broutier bundle. 05 rebuilt 0.139 exactly. |
| Bladder cisplatin n | 79 / 288 here vs 77 / 294 in the paper | The paper's patient list is not public. Report 79 and state the paper's 77. 03 uses 68 treated, with a different survival source. |
| TCGA-COAD cetuximab | 15/7 deaths (02 GI) vs 14/6 deaths (02 OTHER) | The two use different survival sources and parsers. Harmonise on TCGA-CDR. |
| GSE72970 arm sizes | 36/83 (09) vs 40/84 (03, 05) | 09 is correct (I02). |
| GSE39582 treated group | 156 / 200 / 202 | Harmonise (I01). |

Other consistency findings:

- **I10.** 08 says that, within colorectal lines, every navitoclax confidence interval crosses 0. This is false for
  PRISM: rho 0.405 (0.030-0.681).
- **I15.** "Bootstrap" means two different statistics:
  - the paper's top-rank statistic: colon 0.0074, bladder 0.0437;
  - the rank-1 p share: colon 0.0072, bladder 0.0663.

  02 quotes one of them and 05 quotes the other.

## 2. Method correctness

**Tests and implementations.**
- Log-rank and Cox are verified, as described in the summary.
- Median split: `pred > median` = predicted non-responder, as in the paper's code (`<= median` = responder). Ties are
  not an issue for continuous scores.
- HR direction: "HR > 1 = expected" is used consistently in 02, 05, 07 and 09. 03 uses the opposite convention (I05).
- Ties handling is `efron` in 02 GI and `breslow` in 02 OTHER (I14).

**Censoring and time units.**
- The paper does **not** censor at 5 years. It reports the 5-year KM value only (`kmf.predict(60)`, in months), and
  every agent does the same, so nothing needs to change.
- Time units are months in 02, 05 and 03 (GEO), and days from TCGA-CDR in 09. Rank-based tests and HRs are unaffected.
- 03 does pool KM curves across different endpoints (I13).

**Treated-patient definitions.** I01, I02 and I14 above. In 09, drug records without a start date count as exposed
from day 0. This affects 1-6% of patients per drug (I19).

**Duplicate patients.** There are none across TCGA-COAD and TCGA-READ in the training table (0 duplicated
barcodes). The COAD+READ pooled n equals the sum of the two (188 = 129 + 59).

**Leakage.**
- 03 uses no patient outcomes outside the training folds: the inner CV for w and lambda runs on training blocks
  only. Two smaller problems remain:
  - the alignment and CORAL step uses all of the unit's patient expression, including the held-out cohort (I06);
  - the learning-curve hyper-parameters come from the main run (I06).
- 07 B3: the Bolis gene selection is done inside the folds; patient outcomes are never used for selection. The DL grid
  is chosen on all lab models before the organoid out-of-fold predictions (I21).

**Pseudo-replication.**
- Bladder: the 9 organoids include 2 recurrence lines from the same patients, and the SCBO-3 response pools SMBO-3
  screens (I07).
- Ovarian: 3 patients give 2 organoids each (already stated).

**Multiple-testing families.** The A3 binomial p-values ignore correlation between tests (I08). Adjusted Cox models
run at EPV < 1-3 (I09).

**Proxies and substitutions.** These are documented in the reports, and the manuscript must keep them visible:
- Oxaliplatin and carboplatin use the cisplatin proximity. This is stated in 02, 08 and 03.
- Bladder uses CoderData AUC instead of IC50; higher AUC = more resistant (checked in `SOURCE.md`).
- The Broutier gene-set mismatch: I12. The fix is a Broutier bundle in 07 (owner 07/10); otherwise exclude Broutier
  explicitly.

**Pre-registration (09).**
- The frozen-file sha256 (`a0bdd620…`) and the plan sha256 (`97b2a634…`) both verify.
- The timing proof is file mtimes only (I16).

## 3. Reproducibility

Re-runs were done in a sandbox mirror, so no other folder was written. The mirror is in the session scratchpad, not
kept in the repo.

| Script | Result |
|---|---|
| 09 `s02_stage2.py` | all 4 tables identical to the stored ones |
| 09 `s01_stage1.py` (1,000-set nulls, 2 processes) | all 9 tables identical, including the random nulls |
| 09 `t_cox_check.py` | own Cox = PHReg to 8 decimals |
| 05 A3 rebuild (`_data/a3_check_*.tsv`, by 05) | 114/114 log-rank p identical to 02's stored values |
| 10 `c01` (independent re-derivation from B3 bundles) | 15/17 headline p identical. The other 2: Broutier (I12) and the Gao pancreas label lookup |

**Seeds** are fixed in every folder (02: 0 and 1; 03: 20240611; 05, 08, 09, 07: 20261117).

**Data versions.**
- sha256 tables exist for 02 OTHER, 03, 05, 08, 09 and 10.
- 02 GI had sizes only. I wrote `tables/data_versions_02GI_sha256.tsv` (24 files).
- 07 has no data-versions table yet.

## 4. Figures (21 spot-checked)

I viewed these 21 figures:
- 02: COAD authors' 5-FU panels; GI and OTHER cross-set forests; BLCA_CISPLATIN_full; regimen split.
- 03: main_colorectal; forest_cross_cancer; alignment_bladder.
- 05: A3 Fig1, Fig2, Fig5; B1 Fig5, Fig7.
- 08: fig3, fig4, fig5.
- 09: fig2, fig7.
- 10: Fig1c, SuppFig4.
- 07: SuppFig7 centrality.

What is good:
- The palettes are colour-blind safe throughout: Okabe-Ito, viridis, and a purple-orange diverging map.
- Axes are labelled with units and direction ("> 1 = expected").
- The captions I checked match their tables (for example colon p = 0.014, bladder 0.0435, Stage 1 HRs).

Problems:
- **Size (I04):** the multi-panel figures are far larger than a journal page.
- **Overlaps and truncation (I23):**
  - BLCA_full: truncated titles in panels c, f, i and l; panel letters over titles in e, h and j; the panel d legend
    covers points.
  - A3_Fig2b: the "interaction" tick label overlaps.
  - A3_Fig5b: drug labels overlap.
  - Regimen split, panels i and k: p-value text overlaps.
  - GI forest: the "shuffle" header overlaps the plot frame.
  - 03 main panels c, g and h: legend and title overlaps.
  - 08 fig3: the panel letter "e" is used twice.
- **Tiny fonts (I23):** B1_Fig7 and 10 SuppFig4 drop to about 4 pt at print size.
- **Other details (I23):** log axes show "10^-1" instead of 0.1; KM curves have no numbers-at-risk.
- **Pooled KM (I13):** 03's pooled KM mixes OS, RFS, DFS and PFS endpoints.

## 5. Statistical hierarchy (gap G6)

`tables/STATISTICAL_HIERARCHY.tsv` gives one row per analysis family: its tier, folder, number of tests, primary
statistic, correction, what it can claim, headline number and source. The tiers are:

1. **Descriptive:** exact reproduction (02 paper cases) and paper completeness (10, 07).
2. **Exploratory survey:**
   - A3, 57 analyses: BH/Holm and the chance count;
   - B1, 196 tests: BH;
   - 03 joint models: held-out C-index;
   - 07 benchmark.
3. **Mechanistic and outcome-free:** 08 Part A, and the Stage 0 frozen list.
4. **Pre-registered confirmatory:** 09 Stage 1 (Holm across 7) and Stage 2 (Holm across 17). These are the only
   confirmatory claims. Result: G2 0/7 and G3 0/7.

The scorecard figure the review asked for is not drawn here. The table holds its content. One figure would need the
07 results, which are still running.

## 6. Edits other agents need (owners)

- **02:** I03, I07, I14, I15, I17, I18, I24; resize figures (I04).
- **03:** I01, I02, I05, I06, I13; resize figures.
- **05:** I01, I02, I08, I09, I15; resize B1_Fig5; fix A3_Fig2b/5b labels.
- **07:** I12 (Broutier bundle), I21, and a data_versions.tsv.
- **08:** I10; duplicate panel letter in fig3; I24.
- **09:** I16, I19, I20; add the harmonised GSE39582 sensitivity (I01).
- **10:** I12 (re-run c01 after the 07 bundle exists); font sizes in SuppFig4.
- **Coordinator:** correct `docs/PROJECT_SUMMARY.md` "colon 5-FU p = 0.62" to 0.866 (I11).

## Files written

- `AUDIT.md` (this file) and `issues.tsv`.
- `tables/MASTER_NUMBERS.tsv`, `tables/STATISTICAL_HIERARCHY.tsv`, `tables/data_versions_02GI_sha256.tsv`.
- `code/a01_sha256_02GI.py`, `code/a02_master_numbers.py`, `code/a03_reconcile.py`. Re-run them in this order;
  each takes seconds.
