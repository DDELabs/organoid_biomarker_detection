# Fixes that belong to agents still running (not applied by the audit fixer)

Folders `nar/07_B3_models/`, `nar/13_gaps_bio/` and `nar/14_gaps_stats/` were running while the fixes in
`nar/15_fixes/` were applied, so nothing in them was edited. Each item below names the audit issue
(`nar/12_audit/issues.tsv`) and the exact change.

## nar/07_B3_models (B3 model benchmark)

1. **I12, Broutier bundle.** `build_bundles.py` has no ("broutier", "LIHC") bundle, so nar/10 re-derived Broutier
   sorafenib with the Ji gene-universe bundle (p = 0.089 instead of the authoritative 0.139). Add a bundle scored on the
   Broutier gene universe, then nar/10 must re-run `c01_verify_rank1.py`. If this is not done, state in the 07 and 10
   reports that Broutier is excluded from the cross-set analyses.
2. **I22, data versions.** 07 has no `tables/data_versions.tsv`. Write one (file, bytes, sha256) for every input, as in
   `nar/03_organoid_plus_patient/tables/data_versions.tsv`.
3. **I21, DL hyper-parameters.** The MOLI grid is chosen by CV on all lab models and then used for out-of-fold lab
   predictions (no nesting). Nest the grid inside the organoid folds, or state "organoid rho for DL is optimistic
   (grid not nested); patient tests unaffected".
4. **I01 / I02, patient definitions.** If 07 evaluates on GSE39582 or GSE72970, import
   `nar/15_fixes/code/definitions.py`: GSE39582 treated = 5FU|FUFOL|FOLFOX, stage II-III, unknown regimen excluded
   (`nar/15_fixes/definitions/GSE39582_treatment.tsv`); FOLFIRINOX (regex `DEF.BOTH_OX_IRI`) out of both arms of any
   oxaliplatin-vs-irinotecan comparison.
5. **I05, direction.** Report HR with "HR > 1 = expected direction" (risk score), as all other folders now do.

## nar/14_gaps_stats (statistics gaps)

1. **I08 / G9.** nar/05 REPORT now says the A3 binomial p (7/57 treated, 12/57 untreated) assumes independence and
   points to G9. G9 should give the dependence-aware version: a permutation over drugs/pathways keeping shared
   patients, or the count of unique (cohort, pathway) splits ("12 of 57, from k unique splits"). When done, replace the
   binomial p values quoted in `nar/05_A3_B1/REPORT.md` ("In short" and the "Expected by chance" table) with it.
2. **I01.** Any re-analysis of GSE39582 should use the shared definition above (139 treated stage II-III, 138 with RFS;
   264 untreated stage II-III, 258 with RFS).
3. **I09.** If G-analyses re-use A3 adjusted Cox models, use `nar/05_A3_B1/tables/A3_adjusted_cox_epv.tsv`
   (flag `low_epv_full`; reduced models with >= 3 deaths per term).
4. If 14 draws the statistical-hierarchy scorecard, take nar/03 HRs from the new `HR_risk` columns
   (`nar/03_organoid_plus_patient/tables/summary_main.tsv`), which follow "HR > 1 = expected".

## nar/13_gaps_bio (biology gaps)

1. **I07.** Bladder organoids: 9 lines from 7 patients (SCBO-3_2 and SCBO-11_2 are recurrence lines). Any new bladder
   analysis should state n patients and, if feasible, keep one line per patient
   (sensitivity: `nar/02_organoid_only/BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv`).
   SCBO-3 response pools SMBO-3 screens (`data/curated/preclinical/bladder_lee2018/response.tsv`, column screens);
   check in the curation whether SMBO-3 is the same patient.
2. **I24.** If 13 re-uses GSE81005 (HCT-8 5-FU isogenic), state which samples and time points enter the test
   (08 uses 3 vs 3 per time point: 0 h p = 0.18, 24 h p = 0.007, 48 h p = 0.019).
3. **I01 / I02 / I05** as for 07 above, if GSE39582, GSE72970 or HRs appear.

## nar/10_completeness (not running, but depends on 07)

- **I12:** re-run `c01_verify_rank1.py` after the 07 Broutier bundle exists (then Broutier p should be 0.139).
- **I23:** font size in SuppFig4 (~4 pt at print size) and the truncated row label in Fig1c were not changed here.
