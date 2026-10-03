# Cross-cancer summary (robust pipeline, 500 permutations each)

| Study | Pre-clinical model | n | Candidate pathways | Treated pts / deaths | Adj. HR (95% CI) | Permutation p (robust) | Paper-style best p → permutation p | Direction | Status |
|---|---|---|---|---|---|---|---|---|---|
| COAD 5-FU | organoids (van de Wetering) | 19 | 37 | 114 / 8 | 2.08 (0.86–5.07) | **0.048** | 0.076 → 0.48 | as expected | promising, not validated |
| LIHC sorafenib | organoids (LICOB) | 50 | 167 | 26 / 15 | 0.39 (0.18–0.83) | 0.95 | 0.003 → 0.25 | reversed | not supported |
| LUAD cisplatin | GDSC cell lines (proxy) | 65 | 69 | 83 / 27 | 0.81 (0.54–1.21) | 0.90 | 0.147 → 0.81 | reversed | not supported |
| BRCA doxorubicin | GDSC cell lines (proxy) | 37 | 8 | 334 / 29 | 1.01 (0.69–1.49) | 0.51 | 0.092 → 0.49 | inconsistent | not supported |
| BLCA gemcitabine | GDSC cell lines (proxy) | 17 | 35 | 88 / 34 | 0.96 (0.67–1.37) | 0.61 | 0.029 → 0.22 | reversed | not supported |
| PAAD gemcitabine | GDSC cell lines (proxy) | 15 | 35 | 103 / 50 | 0.98 (0.73–1.32) | 0.57 | 0.198 → 0.72 | reversed | not supported |
| STAD cisplatin | GDSC cell lines (proxy) | 16 | 69 | 46 / 16 | 0.93 (0.55–1.60) | 0.60 | 0.175 → 0.69 | reversed | not supported |

HR is per SD of the predicted-IC50 score in treated patients, adjusted for stage, age and sex
where estimable (`adj_covariates` in each `report.json`). Primary data version for the six new
studies: current GDC STAR FPKM-UQ with ssGSEA; for COAD: the paper's committed matrices.

## Reading the table

1. **The original readout produces false positives.** Picking the best of 3 models x 9 values
   of k gives nominal p < 0.1 in 4 of 7 studies (LIHC p = 0.003, BLCA p = 0.029). None survives
   a whole-pipeline permutation of the pre-clinical response (permutation p 0.22–0.81). The
   log-rank test also ignores direction: the LIHC "hit" separates patients the wrong way round.
2. **Only COAD/5-FU keeps a signal**, and only weakly. The robust 3-pathway signature (BH3-only
   protein activation, hyaluronan uptake/degradation, reversible hydration of CO2) passes
   5/5 scorecard checks and permutation p = 0.048, with direction preserved in all three data
   versions. But it rests on 8 deaths, the treatment-interaction p is 0.12, and the effect
   shrinks on current GDC data (HR about 1.2). It is a hypothesis for an external cohort, not
   a validated biomarker.
3. **Cell-line proxies did not transfer.** They are not organoids: no stroma or 3-D context, and
   GDSC 2012 has only 15–17 lines for bladder, pancreas and stomach. In-sample organoid fit is
   also unreliable at this n: leave-one-out Spearman is negative for BLCA (-0.75) and BRCA.
4. **TCGA treatment confounding** (sorafenib and gemcitabine are given in advanced disease) can
   reverse associations. The treatment-interaction test is there to catch this.

## Not run

- **Glioblastoma**: no reachable temozolomide response data (GDSC 2012 lacks TMZ; HGCC and GBO
  sets are on hosts blocked from this environment).
- **Bladder organoids (Lee 2018), PDAC organoids (Tiriac 2018)**: expression on GEO / dbGaP,
  blocked here. Add a loader with `cohorts.generic_organoids` once downloaded.
