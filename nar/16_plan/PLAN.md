# Plan: more cohorts, every model, and what drives the breast model

Folder `nar/16_plan/`. Written 2026-11-18 by the planning advisor. Light work only: I streamed the header of each
GEO series file (no full downloads), read the existing code and results, and wrote this plan.

Files here
- `PLAN.md` (this file).
- `COHORTS_TO_ADD.tsv`: 105 candidate cohorts. Columns: accession, cancer, platform, n, treatment, endpoint,
  verified, why_useful, priority (1 = do first, 4 = optional, `-` = do not use), note.
- `code/probe_geo_headers.py`: the checker. It reads only the top of each GEO series-matrix file and lists every
  sample field with its most common values. Re-run: `python3 -I code/probe_geo_headers.py _verify GSE14814 ...`
- `_verify/<GSE>.txt`: what the checker saw for each cohort (the evidence for the `verified` column).

What "verified" means
- `yes`: the series exists, tumour expression is there, and both a treatment field and an outcome field are in the
  public GEO record. The numbers in the table (n, events, arm sizes) are copied from that record.
- `yes (outcome); treatment no`: expression and outcome are there, but no treatment field. Use these only as
  untreated or unknown-treatment comparison groups.
- `partial`: the series exists, but the treatment or the outcome is not in the GEO record. It may be in the paper's
  supplement. Check before building.
- `no`: not checked from here (controlled access, for example POPLAR/OAK and IMvigor210 at EGA).
- "per paper" in the treatment column means the GEO record does not say it; the original paper does (from memory,
  not re-checked here).

The Seqout and PubMed tools were not needed: GEO FTP is open from this machine, and reading the real sample
fields is a stronger check than a search hit.

## Where we start (numbers from existing files)

- Breast, leave-one-cohort-out (train on all other breast cohorts, test on one never seen): L2 logistic regression
  on 673 Reactome pathway scores, AUC 0.67-0.85 in 13 of 13 cohorts
  (`robust/results/train_combined/part1_response_leave_one_cohort_out.tsv`).
- Lung: only TCGA-LUAD (n = 79, AUC 0.55) and TCGA-LUSC (n = 55, AUC 0.65) have response labels. GSE42127 is the only
  GEO lung cohort in the table, with survival only.
- Positive results so far: breast response model; bladder ERCC2 concordance (p = 0.003); colorectal cetuximab vs
  BRAF V600E (p = 0.031, 0.00072 on current data); BH3-only pathway vs navitoclax in about 800-960 cell lines
  (rho about -0.20); isogenic HCT-8 and T24 pairs (`nar/11_review/STORY_ASSESSMENT.md`, C2 and C6).
  The user also named ovarian. I found no ovarian result file that meets the "positive" rule below. W5 checks this
  and reports it either way.

## Verified new cohorts (count of `yes` rows with treatment and outcome; full detail in the TSV)

| Cancer | With treatment + outcome | Untreated/unknown-treatment controls | Partial (need supplement) |
|---|---|---|---|
| Lung | 14 (1 is a NanoString panel, GSE93157) | 10 | 5, plus CPTAC-3 on GDC |
| Breast | 26 (incl. 2 untreated sets) | 0 | 3 |
| Ovarian | 7 | 0 | 1 |
| Rectal (chemoradiotherapy) | 7 | 0 | 0 |
| Colorectal | 1 | 2 | 1 (GSE5851 cetuximab) |
| Bladder | 2 | 1 | 0 |
| Stomach / oesophagus | 1 / 1 | 1 | 3 |
| Melanoma | 3 (all small) | 0 | 0 |
| Myeloma (APEX, randomised) | 1 | 0 | 0 |
| Head and neck | 2 | 0 | 0 |
| Prostate | 1 | 0 | 1 |
| Pancreas | 0 | 3 | 0 |

Top lung picks: GSE14814 (JBR.10, randomised cisplatin-vinorelbine 71 vs observation 62, 60 deaths), GSE68465
(adjuvant chemo 89 vs none 341, 236 deaths), GSE161537 (immunotherapy, 82, RECIST response + OS),
GSE218989 (immunotherapy, 355, responder 168 vs 187; first check it is bulk RNA, not single-cell),
GSE33072 (BATTLE, erlotinib 25 / sorafenib 39), GSE37745 and GSE29013 (adjuvant yes/no).
Rows marked `-` (do not use) are not counted above. Not usable: GSE61676 (blood, not tumour), GSE43580 (no outcome), GSE28796 (2 responders).

Honest expectation for lung: almost all lung GEO cohorts have survival, not response. So the lung test is mainly
"does a pathway score predict who benefits from adjuvant chemo" (score x treatment interaction, with GSE14814 as the
only randomised test), plus a small immunotherapy response test (GSE161537, GSE218989, GSE135222, GSE126044).

---

## Shared rules for every workstream

- Follow `nar/AGENT_BRIEF.md` (plain English, real numbers, "What I think:", no git, seeds, 300 dpi PNG + PDF,
  caption files, at most 2 worker processes).
- New training table goes to a NEW folder `robust/results/training_table_v2/` so running agents keep the old one.
  The only edits allowed outside an agent's own folder: (a) add cohort ids to `COHORTS` in
  `robust/fetch_geo_matrices.py`, (b) add builder functions + names to `BUILDERS` in
  `robust/build_training_cohorts.py`, (c) add an output-folder argument and new `GROUP` keys in
  `robust/build_training_table.py`. The coordinator approves these three edits.
- Pathway scores: always the same 673 Reactome within-sample rank scores (`build_training_table.score`), then
  z-score inside each cohort (`obd.survlearn.zscore_within`). Drop `possible_duplicate_of` rows before pooling.

---

## W1. Lung cohort harvest and analysis

Folder: `nar/17_lung/` (code, tables, figures, REPORT.md).

Steps
1. Download: add a `"lung"` list to `COHORTS` in `robust/fetch_geo_matrices.py` with all priority 1-3 lung rows of
   `COHORTS_TO_ADD.tsv`; run it (matrix + supplementary listing only).
2. GSE218989: open its supplementary listing; keep it only if it holds a bulk gene-by-sample matrix. Record the
   decision in REPORT.md.
3. GSE232569, GSE162520, GSE81089: look in the supplementary listings for a clinical table. Use only if the outcome
   is there; otherwise list in `skipped`.
4. Write one builder per cohort in `robust/build_training_cohorts.py` (same pattern as `GSE48276()`: `array_expr` for
   arrays, symbol-level matrix for RNA-seq, `base()` for clinical). Field mapping:
   - `drugs`: GSE14814 ACT -> `CISPLATIN;VINORELBINE`, OBS -> `NONE`; GSE68465/GSE37745/GSE29013 adjuvant yes ->
     `PLATINUM` (regimen not given), no -> `NONE`, unknown -> empty; GSE33072 -> `ERLOTINIB` / `SORAFENIB`;
     immunotherapy sets -> `ICI_ANTI_PD1_PDL1`; untreated-control sets -> empty (not `NONE`, because not recorded).
   - `randomised=True` only for GSE14814.
   - `responder`: GSE161537 CR/PR = 1, SD/PD = 0; GSE218989 Responder = 1; GSE126044 responder = 1;
     GSE33072 8-week disease control = 1.
   - Survival in months (GSE14814 OS time is in years: multiply by 12; GSE41271 from dates; GSE3141 parse text).
   - Stage, age, sex, histology (adeno / squamous) into existing columns; add `histology` to `COLS`.
5. Add `("lung", "lung")` to `GROUP` in `build_training_table.py`; rebuild into `training_table_v2`.
6. Analyses (new script `nar/17_lung/code/run_lung.py`, importing `run_train_combined` functions):
   a. Response, leave-one-cohort-out, lung only: immunotherapy cohorts as one group; chemo response cohorts
      (TCGA-LUAD/LUSC) as another. Also apply the frozen breast model and a pan-cancer model to each lung cohort.
   b. Adjuvant-chemo benefit: `part3`-style Cox model `treated + score + treated x score` (+ age, stage), training
      on all other lung cohorts, testing on each held-out cohort with >= 15 treated, >= 15 untreated, >= 15 events.
      GSE14814 is the key test because it is randomised.
   c. Prognosis-only check: the same score in untreated patients alone (GSE31210, GSE50081, GSE30219, GSE72094...).
      If it predicts survival there too, it is prognostic, not drug-specific.
   d. Repeat (b) with the paper's network-proximity shortlist for cisplatin (existing STRING code) vs all pathways.
Inputs: GEO matrices under `data/external/geo_train/`, `robust/results/training_table_v2/`.
Outputs: `nar/17_lung/tables/lung_loco_response.tsv`, `lung_adjuvant_interaction.tsv`, `lung_untreated_prognosis.tsv`,
`lung_cohorts_built.tsv`; figures `lung_forest_interaction.png/.pdf`, `lung_auc_by_cohort.png/.pdf` + captions.
Acceptance
- every priority 1-2 lung row is either in `training_table_v2/cohorts.tsv` or in `skipped.tsv` with a reason;
- GSE14814 shows 71 ACT / 62 OBS and 60 deaths after building (matches `_verify/GSE14814.txt`);
- duplicate check run against GSE42127 and TCGA (same hospitals possible);
- each test row reports n, events, HR with 95% CI and p.
Runtime: download 30-60 min (about 1-1.5 GB), builders 30 min, scoring 20 min, analysis 10 min.

## W2. Other-cancer cohort harvest

Folder: `nar/18_other_cohorts/`.

Steps
1. Same download and builder routine as W1, for priority 1-3 non-lung rows. Order: breast (GSE16446, GSE50948,
   GSE34138, GSE17705, GSE66305, GSE42822, GSE6861, GSE32603, GSE6532, GSE58644, GSE7390, GSE11121), ovarian
   (GSE140082, GSE32062, GSE26712, GSE49997, GSE17260, GSE30161), bladder (GSE13507, GSE176307), myeloma (GSE9782),
   rectal chemoradiotherapy (GSE35452, GSE119409, GSE150082, GSE94104, GSE45404), oesophagus (GSE104958),
   stomach (GSE26253), head and neck (GSE65858), prostate (GSE116918).
2. Before building, de-duplicate: GSE173839 vs GSE194040 (I-SPY2), GSE4779 vs GSE6861, GSE42822 vs
   GSE20194/GSE25066, GSE25055 is part of GSE25066 (skip). Use the existing `DUP_CHECK` correlation (> 0.95).
3. For `partial` rows (GSE5851 cetuximab, GSE14210, GSE62254, GSE20181, GSE87455): read the supplementary listing;
   build only if a per-sample outcome table exists.
4. New `GROUP` keys: `("stomach","stomach")`, `("gastric","stomach")`, `("rectal","colorectal")` already there,
   `("myeloma","myeloma")`, `("head","head and neck")`, `("prostate","prostate")`, `("pancrea","pancreas")`.
5. Re-run `run_train_combined.py` parts 1-3 on `training_table_v2` (output to `nar/18_other_cohorts/tables/`).
   New tests it unlocks: ICON7 bevacizumab x score (randomised), NOAH trastuzumab x score (randomised), APEX
   bortezomib vs dexamethasone (randomised), CHER-LOB arms, tamoxifen survival cohorts, rectal chemoradiotherapy
   response leave-one-cohort-out (7 cohorts), ovarian platinum survival.
Outputs: `nar/18_other_cohorts/tables/cohorts_built.tsv`, `part1_v2.tsv`, `part2_v2.tsv`, `part3_v2.tsv`,
`new_randomised_tests.tsv`; REPORT.md with before/after counts per cancer.
Acceptance: arm sizes and event counts after building match `_verify/*.txt`; every attempted cohort is built or
skipped with a reason; the old breast 13-cohort AUCs are reproduced exactly from `training_table_v2` before adding
new breast cohorts (guards against breaking the builder).
Runtime: download 1-2 h, builders 1-2 h (about 30 cohorts), training 20 min.

## W3. All-model benchmark on breast and every positive result

Folder: `nar/19_model_benchmark/`. Reuse `nar/07_B3_models/code/learning.py`, `nets.py` (network models) and
`b3common.py` where they fit; do not rewrite them.

Models (all with the same folds):
L2, L1 and elastic-net logistic regression; random forest; XGBoost; linear and RBF SVM; small MLP (sklearn);
graph model with a network prior (GCN or network-smoothed logistic, from `nets.py`, CPU torch); and the network
proximity shortlist vs all 673 pathways, for each linear model. Baselines: ER+HER2 only (breast), random
pathway sets of equal size (200 draws), and a proliferation-only score (see W4).

Settings
- Breast response: leave-one-cohort-out over all breast response cohorts (13 now, more after W2).
- Positive results with patient labels: bladder ERCC2 (mutant vs wild type as the label, TCGA-BLCA), colorectal
  cetuximab BRAF V600E (organoid-trained model, TCGA-COAD), rectal chemoradiotherapy (after W2), lung
  immunotherapy (after W1). For single-cohort settings use 5-fold x 10 repeats instead.
- Hyperparameters chosen by inner cross-validation inside the training cohorts only (never the test cohort).
Steps
1. `bench_breast.py`: build X (z-scored within cohort), y, cohort ids; loop models x held-out cohort; save
   per-cohort AUC, AUC within ER/HER2 groups, and test-set predictions.
2. `bench_positives.py`: same model loop for the other settings.
3. Compare models: mean AUC across cohorts, paired DeLong or bootstrap difference vs L2, number of cohorts where
   each model beats L2.
Outputs: `tables/breast_loco_all_models.tsv` (model, test_cohort, n, responders, auc, auc_within_er_her2),
`tables/positives_all_models.tsv`, `tables/model_comparison.tsv`; figure `fig_model_heatmap.png/.pdf`
(models x cohorts) and `fig_model_mean_auc.png/.pdf`.
Acceptance: the L2 row reproduces `part1_response_leave_one_cohort_out.tsv` to 3 decimals for all 13 cohorts;
no test cohort used in tuning; seeds fixed; runtime log saved.
Runtime: linear models 5 min; RF/XGB/SVM about 30-40 min; MLP 15 min; graph model 30-60 min (2 workers).

## W4. Which pathways drive the breast prediction

Folder: `nar/20_breast_pathways/`. Start from `robust/run_breast_pathways.py` and
`robust/results/train_combined/breast_model_pathways.tsv` (weights, sign stability, per-cohort odds ratios exist).

Questions and how to answer each
1. One, few or many? Leave-one-cohort-out AUC for: the single best pathway (chosen in training cohorts only);
   the top k pathways for k = 1, 2, 3, 5, 10, 20, 50, 100, 200, 673 (ranked by training-fold weight or
   univariate p); plot AUC vs k. Report the smallest k whose mean AUC is within 0.01 of the 673-pathway model.
2. Weights and stability: L2 weights refit in each leave-one-out fold; L1 / elastic-net selection frequency over
   100 bootstrap fits (stability selection); list pathways selected in >= 80% of fits.
3. Minimal subset: greedy forward selection inside training folds only (nested), so the reported AUC stays honest.
4. Is it mostly proliferation? Define a proliferation score = mean of cell-cycle Reactome pathways (names containing
   CELL_CYCLE, MITOTIC, G1_S, G2_M, E2F, DNA_REPLICATION, PRE_REPLICATIVE, CHECKPOINT; list saved to a table) and,
   separately, an MKI67-like gene score from expression. Then: (a) AUC of the proliferation score alone per
   cohort; (b) AUC of the full model after removing all proliferation pathways; (c) partial test: does the full
   model add to proliferation + ER + HER2 in a logistic model (likelihood-ratio p per cohort, combined);
   (d) correlation between the model score and the proliferation score per cohort.
5. Pathway families: group the 673 pathways into Reactome top-level families (cell cycle, immune system,
   signal transduction, metabolism, ...) and report the summed absolute weight per family.
Outputs: `tables/auc_vs_k.tsv`, `tables/stability_selection.tsv`, `tables/minimal_subset.tsv`,
`tables/proliferation_tests.tsv`, `tables/family_weights.tsv`; figures `fig_auc_vs_k`, `fig_stability`,
`fig_proliferation` (+ captions); REPORT.md answers "single vs few vs many" and "mostly proliferation?" in two
short paragraphs with the numbers.
Acceptance: every subset choice is made inside training folds only; the 673-pathway row equals the W3 L2 row.
Runtime: about 30-45 min (bootstrap L1 fits dominate).

## W5. Positive-results folder

Folder: `nar/POSITIVE_RESULTS/`, one subfolder per cancer: `breast/`, `bladder/`, `colorectal/`, `ovarian/` (only if
something qualifies), `pan_cancer_cell_lines/` (BH3-only vs navitoclax), plus `lung/` etc. if W1/W2 produce one.

Rule for "positive" (written in `nar/POSITIVE_RESULTS/README.md`): a held-out or independent test with p < 0.05 in
the pre-stated direction, or BH q < 0.10 across its own family of tests, with n and events reported. Everything else
stays where it is.
Steps
1. Script `nar/POSITIVE_RESULTS/code/collect.py` scans result tables (paths listed in the script, starting from
   `nar/10_completeness/RESULTS_INDEX.md` and the W1-W4 outputs) and writes `POSITIVE_INDEX.tsv`
   (cancer, result, n, events, effect, CI, p, q, source file).
2. Copy (do not move) each qualifying figure, caption and table into its subfolder; write a short `README.md` per
   subfolder: what was tested, the numbers, "What I think:", and the source path.
3. Ovarian: report the best ovarian result found and whether it meets the rule.
Acceptance: every number in each README is copied from a file listed next to it; no original file moved or edited.
Runtime: 20-30 min.

## W6. Loader for the user's own organoid data

Folder: `nar/21_user_organoids/`. Code: `code/load_user_organoids.py`.

Input contract (one folder):
- `expression.*` (tsv/csv/xlsx, optionally .gz): genes x samples or samples x genes (detected), gene symbols or
  Ensembl ids (mapped to symbols), raw counts, TPM/FPKM or log values (detected; counts -> log2 CPM, TPM -> log2(x+1)).
- `drug_response.*`: long format `sample, drug, metric, value` (metric = AUC, IC50, viability, ...), or wide
  (samples x drugs). Drug names normalised with the existing name map used in `robust/obd/`.
- optional `samples.*`: sample, patient, cancer_type, organoid_line, replicate.
Steps the loader runs
1. Read, detect orientation and value type, map genes, report sample name mismatches between files.
2. Score the 673 Reactome pathways with the same `score()` function as the patient table.
3. Write `organoid_set.json` + `expression.tsv.gz` + `response.tsv` + `pathway_scores.tsv.gz` in the layout of
   `data/curated/preclinical/<set>/` so `obd.preclinical_sets` can load it like the existing sets.
4. One command runs the existing pipeline on it: organoid-only model -> TCGA/GEO patients (paper method), and the
   frozen breast/other models applied to the organoids.
Tests: a synthetic folder (made from `colorectal_vandewetering2015`, saved in three formats) must reproduce that
set's existing pathway scores exactly; unit tests for orientation, Ensembl ids, counts vs log values, wide vs long.
Outputs: loader, `tests/`, `README.md` (how the user drops files in). Runtime: 1-2 h agent time, seconds to run.

---

## Recommended order

1. W6 (independent, small) and W1 steps 1-4 (lung download + builders) in parallel.
2. W1 steps 5-6 (rebuild `training_table_v2`, lung analysis), then W2 (same machinery, more cohorts).
3. W3 and W4 on breast as soon as `training_table_v2` reproduces the old breast numbers (they can start on the
   current table and be re-run on v2 later; scripts take the table folder as an argument).
4. W3 on the other positives after W1/W2.
5. W5 last, because it collects everything.

## Token-efficient execution notes

- Scripted vs agent work: downloads, builders, table rebuild, benchmarks and the W5 collector are all scripts; the
  agent writes the script once, runs it in the background, and reads only the summary TSV (head/tail), never full
  logs. The checker `code/probe_geo_headers.py` is the cheap way to look at any cohort's fields.
- Run at most 2 Sonnet execution agents at a time (4 CPUs, other agents running). Pairs: (W6, W1) -> (W2, W3-breast)
  -> (W4, W3-positives) -> W5.
- Each agent reads only: `nar/AGENT_BRIEF.md`, this PLAN section, the named TSV rows, and the named code files.
- Reuse, do not rewrite: `robust/fetch_geo_matrices.py`, `robust/build_training_cohorts.py`,
  `robust/build_training_table.py`, `robust/run_train_combined.py`, `robust/run_breast_pathways.py`,
  `nar/07_B3_models/code/`.

Short prompts for the execution agents (paste as is):

- W1a: "Read nar/AGENT_BRIEF.md and nar/16_plan/PLAN.md section W1. Add the priority 1-3 lung accessions of
  nar/16_plan/COHORTS_TO_ADD.tsv to COHORTS in robust/fetch_geo_matrices.py and run it. Write one builder per cohort
  in robust/build_training_cohorts.py following GSE48276(), with the field mapping in W1 step 4. Check arm sizes and
  events against nar/16_plan/_verify/<GSE>.txt. Report built/skipped per cohort. Folder nar/17_lung/."
- W1b: "Rebuild the training table into robust/results/training_table_v2 (W1 step 5) and run W1 step 6 analyses
  a-d. Report n, events, HR (95% CI), p per test. Folder nar/17_lung/."
- W2: "Do PLAN.md W2 for the listed non-lung cohorts in the listed order. De-duplicate first. Then re-run
  run_train_combined parts 1-3 on training_table_v2. Folder nar/18_other_cohorts/."
- W3: "Do PLAN.md W3. First reproduce the L2 breast leave-one-cohort-out AUCs to 3 decimals, then add the other
  models with inner-CV tuning only. Reuse nar/07_B3_models/code. Folder nar/19_model_benchmark/."
- W4: "Do PLAN.md W4 questions 1-5 on breast, all subset choices inside training folds. Folder
  nar/20_breast_pathways/."
- W5: "Do PLAN.md W5 with the stated positive rule. Copy, never move. Folder nar/POSITIVE_RESULTS/."
- W6: "Do PLAN.md W6. Build the loader and its tests using a synthetic copy of colorectal_vandewetering2015.
  Folder nar/21_user_organoids/."

Total estimated wall time with 2 agents at a time: about 1.5-2 working days of compute and agent time, most of it
downloads (W1/W2) and the slow models (W3).
