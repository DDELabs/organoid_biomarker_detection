# Organoid / cancer drug-response biomarker project — results and conclusions (as of 2026-11-19)

Repository: github.com/DDELabs/organoid_biomarker_detection, branch `claude/zealous-darwin-zeiopw`, folder `nar/`.
Start point: Kong et al. 2020 Nat Commun (DOI 10.1038/s41467-020-19313-8): network-proximity pathways + organoid drug
response + ML → predict patient survival on a drug (colon 5-FU, bladder cisplatin).

## 1. Reproduction of the paper
- Colon 5-FU reproduces exactly on the authors' own files: rank-1 pathway "Activation of BH3-only proteins", 114 treated
  patients, 8 deaths, log-rank p = 0.0144; untreated p = 0.161; isogenic HCT-8 cells p = 0.0085; BRAF V600E p = 0.031.
- On current TCGA data the same analysis gives p = 0.866. Cause (decomposition): the patient expression scoring; putting
  the authors' 2019 patient scores back gives p = 0.005. 94 of the 114 "5-FU" patients actually received FOLFOX.
- Bladder cisplatin: paper's pathway is rank 2, not rank 1; p = 0.031 (paper 0.01); current data p = 0.067 with shuffle
  p = 0.036; ERCC2 mutation link p = 0.003 (specific: only 0.2% of pathways do as well); isogenic cells p = 0.0002.
- 77 paper items checked (checklist in nar/00_checklist): all reproduced or explained.

## 2. Does the organoid-network method generalise? (57 organoid-only analyses, 11 lab-model sets, 7+ cancers)
- 7/57 treated tests p<0.05 (2.9 expected by chance; permutation p = 0.03); none survives multiple-testing correction.
- Untreated patients show MORE small p-values (12/57) → the signal is mostly general prognosis, not drug response.
- Pooled effect across all 57: HR 0.96 per SD (0.90–1.03); an average HR above 1.02 is ruled out.
- Median 17 deaths per test → only very large effects (HR ~3.9) detectable. Main limit = patient data, not method.
- External cohorts (196 tests): 14 at p<0.05 after audit fixes; 1 survives correction (rectal GSE87211, HR 1.92, q 0.002).
- Pre-registered test of DepMap-derived hypotheses (7 drugs, pooled TCGA + trials): 0/7 pass; rules out moderate effects
  for 5-FU, gemcitabine, paclitaxel.
- Model benchmark (46 pairs): XGBoost/neural nets beat the paper's pathway in cell lines and held-out organoids, but all
  models are near chance on patients (C-index 0.46–0.53). Organoid + patient joint models: organoids add little.

## 3. POSITIVE RESULTS (rule: held-out p<0.05 in expected direction, or q<0.10) — nar/POSITIVE_RESULTS/
- BREAST (strongest): pathway model predicts chemotherapy response (pCR) in 13/13 held-out cohorts (3,511 patients,
  mean AUC 0.752, range 0.67–0.85); 8 of 15 new cohorts also pass (mean 0.682). Every model type gives 0.74–0.76
  (best linear SVM 0.764). Beats ER/HER2 alone (0.693) and proliferation alone (0.646). Holds in ER+/HER2− 0.67,
  HER2+ 0.73, triple-negative 0.72. Uses MANY pathways (~100 needed); top weights: P2Y receptors +0.175, mycobacterium
  immune pathway +0.154, HIV transcription elongation +0.143, E2F pre-replication −0.132, licensing factors +0.130;
  cell-cycle pathways only ~5% of weight (0.721 with proliferation removed). General chemo response, NOT drug-specific.
- LUNG: untreated prognosis score HR 1.36 per SD (1.29–1.44), 2,938 patients, 14/14 cohorts (prognostic, not predictive).
  Adjuvant chemo benefit (JBR.10 randomised) HR 1.02, p 0.95 — not predictive.
- BLADDER: ERCC2 concordance p 0.003 / 0.00055 (current data); isogenic T24 p 0.0002; cisplatin split nominal.
- COLORECTAL: cetuximab vs BRAF V600E p 0.00048 (stable 87–100% of leave-out refits).
- RECTAL: GSE87211 5-FU signature HR 1.91, q 0.002 (no comparator arm).
- OVARIAN: response GSE30161 AUC 0.728, GSE9891 0.624; joint platinum HR 0.85 (p 0.001) and paclitaxel HR 0.88
  (p 0.010) — prognostic in treated patients; olaparib signature GSE63885 p 0.004 (nominal, 6 lines); BRCA PFS HR 0.42.
- PAN-CANCER (DepMap/GDSC/CTRP): BH3-only pathway predicts BH3-mimetic (navitoclax, ABT-737) sensitivity, rank 1 of 672
  pathways, q ≤ 1e-6 → an apoptosis-priming marker, not a 5-FU marker.
- BORDERLINE (not positive): HER2 trials NOAH/CHER-LOB — score works in chemo arm but reverses in anti-HER2 arm
  (p 0.025, 0.039; direction not pre-specified); ICON7, APEX, bladder GSE13507 p≈0.06–0.09.

## 4. User-uploaded organoid data (nar/22_user_data_results)
- Lung GSE165611: drug "IC50" values are per CLUSTER (4), not per sample; 72/120 cells missing → per-sample IC50 needed.
- Ovarian GSE229119: 7 lines incl. MCF7 (breast); olaparib signature nominal in GSE63885 only.
- Liver GSE188541: 3 organoids, mostly non-cancer drugs, no expression in series matrix → no model possible.

## 5. Data assembled
- Patient training table v2: 149 cohorts, 32,226 tumour samples (v1: 73 cohorts, 20,766), 673 Reactome pathway scores.
- 11 lab-model drug-response sets (organoids, PDX, CTRP, PRISM, GDSC); DepMap/Project Score CRISPR.

## 6. Recommended paper directions (user to decide)
1. Rigorous re-evaluation of the organoid-network biomarker (exact reproduction + failure to generalise + power/
   prognostic explanation + pre-registered test) — mostly done.
2. Breast chemotherapy-response pathway model across 28 cohorts (strongest positive), with lung prognosis and
   mechanism results (ERCC2, BRAF, BH3-only priming) as supporting findings.
3. Collect better data (per-sample organoid IC50; paired organoid–patient cohorts on EGA) before claiming drug-specific
   biomarkers.

## 7. Still open
- Needs user: per-sample lung IC50; citation search (Google Scholar); optional Zenodo deposit of frozen hypotheses;
  controlled-access data applications.
- Small fixes: one liver verify script (0.089 vs 0.139), oversized organoid-only figures, bladder decomposition, stale
  myeloma summary flag.
- Not started (by user's choice): manuscript, final figures, reproducibility package.

Key files: nar/POSITIVE_RESULTS/INDEX.md, nar/11_review/STORY_ASSESSMENT.md, nar/12_audit/tables/MASTER_NUMBERS.tsv,
nar/15_fixes/CHANGELOG.md, nar/10_completeness/RESULTS_INDEX.md, each nar/<folder>/REPORT.md.
