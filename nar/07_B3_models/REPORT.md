# B3. How does the paper's network method compare with other models?

All numbers below are copied from the files in `tables/` (code in `code/`, seed fixed in `code/b3common.py`).
Figures are in `figures/` (each as PNG + PDF + a plain-language `.caption.md`).

## 1. What was compared
Every model uses the same lab-model data (organoids, PDX, or cell lines standing in for organoids), the same splits and
the same drug-treated TCGA patients. There are **46 pairs** (a lab-model set plus a drug with at least 15 treated patients).

| Group | Models |
|---|---|
| The paper's method | single network pathway; top-2 network pathways (ridge on pathway scores near the drug's targets) |
| Linear and tree models | ridge on all network pathways; elastic net; random forest; XGBoost (all pathways) |
| Other inputs | whole transcriptome; all pathways; first-neighbour genes; Bolis correlation genes; network-propagation genes |
| Deep learning | MOLI-style network (genes); fingerprint MLP trained on cell lines only (zero-shot) and fine-tuned on lab models |
| Cell-line-pretrained | the fingerprint MLP above (pretrained on GDSC + CTRPv2 cell lines) |
| Graph network | GCN on the STRING network |
| Controls | 100 random pathway sets of the same size; degree/betweenness/closeness-central pathways; lab-model-to-tumour alignment (nar/03) |

Three tasks: (1) predict drug response in lab models left out of training; (2) in drug-treated TCGA patients, does the
predicted resistance order survival (C-index, and hazard ratio per standard deviation; above 0.5 and above 1 = expected
direction); (3) cell-line cross-validation (GDSC + CTRPv2) where cells are new (cell-blind) or drugs are new (drug-blind).

## 2. Main results (table: `tables/benchmark_model_summary.tsv`, `tables/cellline_cv_summary.tsv`)
Mean over the 46 pairs, with 95% CI across pairs.

| Model | Lab models: held-out Spearman | Patients: C-index | Patients: HR per SD (geo-mean) | Cell-blind rho | Drug-blind rho |
|---|---|---|---|---|---|
| Paper: single network pathway | 0.072 (-0.015, 0.158) | 0.489 (0.457, 0.521) | 0.954 (0.847, 1.073) | 0.023 | n/a |
| Paper: top-2 network pathways | 0.093 (0.006, 0.181) | 0.493 (0.464, 0.523) | 0.999 (0.908, 1.099) | 0.062 | n/a |
| Ridge: all network pathways | 0.119 (0.029, 0.208) | 0.487 (0.460, 0.514) | 1.009 (0.918, 1.109) | 0.406 | n/a |
| Elastic net (42 pairs) | 0.159 (0.064, 0.254) | 0.492 (0.465, 0.519) | 1.047 (0.937, 1.170) | 0.352 | 0.461 |
| Random forest | 0.214 (0.125, 0.303) | 0.486 (0.451, 0.520) | 1.008 (0.895, 1.135) | 0.425 | 0.554 |
| XGBoost | 0.220 (0.129, 0.311) | 0.464 (0.434, 0.493) | 0.953 (0.859, 1.056) | 0.452 | 0.558 |
| Deep learning, MOLI-style | 0.151 (0.068, 0.235) | 0.495 (0.459, 0.531) | 0.981 (0.880, 1.092) | n/a | n/a |
| Fingerprint MLP, zero-shot (cell-line pretrained) | 0.223 (0.145, 0.302) | 0.494 (0.462, 0.525) | 1.054 (0.959, 1.159) | 0.440 | 0.535 |
| Fingerprint MLP, fine-tuned on lab models | 0.261 (0.174, 0.347) | 0.492 (0.456, 0.527) | 1.057 (0.951, 1.174) | 0.440 | 0.535 |
| GCN on STRING | -0.012 (-0.101, 0.076) | 0.528 (0.499, 0.558) | 1.096 (0.975, 1.232) | 0.350 (16 units) | n/a |
| Random pathway sets (median of 100) | 0.062 (0.031, 0.093) | 0.499 (0.488, 0.511) | 1.024 (0.993, 1.056) | n/a | n/a |

What I think: the clearest differences are in the lab-model and cell-line tasks. The fingerprint MLP (fine-tuned) has the
highest lab-model Spearman (0.261). On cell lines, XGBoost is highest for cell-blind (0.452) and for drug-blind (0.558),
with random forest (0.554) and the fingerprint MLP (0.535) close behind. The paper's single pathway is at 0.023 cell-blind,
so for new cell lines the models that use many features carry the signal.

## 3. Best model per task
| Task | Highest mean | Value | Comment |
|---|---|---|---|
| Held-out lab models (Spearman) | Fingerprint MLP, fine-tuned | 0.261 | next: Fingerprint MLP zero-shot 0.223, XGBoost 0.220, random forest 0.214 |
| Patients, C-index | GCN on STRING | 0.528 (0.499, 0.558) | all other models are between 0.464 and 0.507; random sets 0.499 |
| Patients, HR per SD | GCN on STRING | 1.096 (0.975, 1.232) | fine-tuned MLP 1.057, zero-shot MLP 1.054 |
| Cell-blind (cell lines) | XGBoost | 0.452 | fingerprint MLP 0.440, random forest 0.425 |
| Drug-blind (cell lines) | XGBoost | 0.558 | random forest 0.554, fingerprint MLP 0.535 |

What I think: no single model leads everywhere. The patient tasks sit close to 0.5 for all models, so differences there are
small compared with the width of the CI. Per pair, the model with the top C-index varies: ridge on all network pathways is top
in 6 of 46 pairs, network-propagation genes in 5, GCN, top-2 pathways, zero-shot MLP and MOLI-style in 4 each
(counts from `tables/benchmark_long.tsv`, random-set and alignment rows excluded).

## 4. The paper's two headline pairs (`tables/benchmark_long.tsv`)
**Colorectal organoids to 5-FU-treated TCGA-COAD (paper data; 19 organoids, 114 patients, 8 events).**
Paper single pathway: C = 0.727 (0.502, 0.860), HR per SD 1.519 (0.936, 2.465), p = 0.091, log-rank of the median split p = 0.014.
Highest C: GCN 0.732 (0.538, 0.918), HR 2.675 (0.939, 7.619), p = 0.065. Ridge on all network pathways: C = 0.698, HR 2.241 (0.971, 5.169), p = 0.058.
Elastic net 0.634, random forest 0.629, XGBoost 0.576, fine-tuned MLP 0.649, MOLI-style 0.676.
What I think: with only 8 events the intervals are wide (C from about 0.5 to 0.9). The paper's pathway and the GCN give about the same C-index,
and ridge on all network pathways has a higher held-out lab-model Spearman than the paper pathway (0.489 versus 0.239).

**Bladder organoids to cisplatin-treated TCGA-BLCA (paper data; 9 organoids, 79 patients, 20 events).**
Highest C: ridge on all network pathways 0.585 (0.472, 0.746), HR 1.644 (1.031, 2.620), p = 0.037. Paper single pathway: C = 0.558 (0.446, 0.727),
HR 1.289 (0.805, 2.066), p = 0.291, log-rank p = 0.044. Elastic net 0.490, random forest 0.473, XGBoost 0.519, GCN 0.478.
What I think: the paper's pathway and the top-2 pathways (0.558 each) sit between ridge-on-all-pathways and the other modern models here.

**Liver organoids (Broutier) to sorafenib-treated TCGA-LIHC (5 organoids, 30 patients, 15 events; new bundle for fix I12).**
Paper single pathway: log-rank p = 0.139 (this matches the authoritative value), C = 0.354, HR 0.544 (0.283, 1.045).
Random forest has the highest C (0.703, HR 2.268, 1.133 to 4.540, p = 0.021); GCN C = 0.662, HR 2.361 (1.169, 4.769).
What I think: five organoids is very few for training. Treat these numbers as a description of this one set; they are in the 46-pair averages above.

## 5. Other checks
- **Random pathway sets** (`tables/random_sets_vs_paper_pathway.tsv`): across the 46 pairs, on average 0.551 of random sets had a C-index at least as high as the paper pathway
  (median C of random sets 0.499; paper pathway mean C 0.489). For lab-model Spearman the share is 0.503. What I think: the chosen pathway is not set apart from random pathways on average, though single pairs differ (range 0.010 to 1.000).
- **Network centrality** (SuppFig 7; `tables/SuppFig7_centrality_profile.tsv`): drug-proximal pathway genes are more central than the average pathway gene. In COAD/5-FU the median degree of proximal genes is 96 versus 43 for all pathway genes (Mann-Whitney p = 2.7e-30); in BLCA/cisplatin 36 versus 43 (p = 0.011). The rank-1 centrality models have mean C = 0.494 (degree), 0.463 (betweenness), 0.507 (closeness).
- **Lab-model-to-tumour alignment** (`tables/alignment_with_vs_without.tsv`): mean C changes by +0.002 (paper pathway), -0.011 (ridge), -0.006 (elastic net), -0.027 (random forest), +0.015 (XGBoost); Wilcoxon p between 0.13 and 0.98. Held-out lab-model Spearman falls for all five models (for example 0.072 to -0.016 for the paper pathway, p = 0.041). What I think: alignment does not add a measurable gain in these pairs.
- **Learning curves** (`tables/learning_curve_lab_models.tsv`, 6 larger sets): held-out Spearman at 25% / 100% of the training data: XGBoost 0.165 / 0.347, random forest 0.175 / 0.323, elastic net 0.103 / 0.297, fine-tuned MLP 0.221 / 0.278, paper pathway 0.065 / 0.141. Patient C-index stays near 0.5 for every model at every size (0.46 to 0.53). What I think: the many-feature models gain from more lab models; the fingerprint MLP is already strong with a quarter of the data.

## 6. Fixes applied from `nar/15_fixes/PENDING_FOR_RUNNING_AGENTS.md`
- **I12.** Added the Broutier liver bundle (`LIHC_broutier`, Broutier gene universe) to `code/build_bundles.py`, `code/b3common.py`; benchmark pair run and in the tables and figures (`figures/fig3_replica/Fig3_LIHC_broutier_SORAFENIB.*`). Log-rank p for the paper pathway = 0.139. nar/10 should now re-run `c01_verify_rank1.py`.
  Note: the bundle builder re-saved `COAD_paper.pkl` as a side effect (112 bytes smaller; same inputs, ordering differences only); its new sha256 is in the data-versions table. The earlier COAD_paper benchmark results were not re-run.
- **I22.** `tables/data_versions.tsv` lists file, bytes and sha256 for every input (27 rows including the Broutier bundle).
- **I21.** The MOLI-style grid is chosen by 3-fold CV on all lab models of a pair and then used for the out-of-fold lab-model predictions (grid not nested). So the lab-model Spearman for the MOLI-style model (0.151) is optimistic; patient tests are unaffected.
- **I01 / I02.** B3 evaluates on TCGA patients only (GSE39582 and GSE72970 are not used), so the shared GSE39582 definitions do not apply.
- **I05.** HR is reported as a risk score (HR above 1 = expected direction) throughout.

## 7. Paper checklist items covered
- **Fig3b-f**: `figures/fig3_replica/` (one replica per pair: 46 PNG + PDF + caption), the benchmark heatmaps and `figures/B3_main_figure.*`.
- **SuppFig6**: `figures/SuppFig6_replica_COAD_paper_FLUOROURACIL.*` and `tables/SuppFig6_7_centrality_selection.tsv`.
- **SuppFig7**: `figures/SuppFig7_replica_BLCA_paper_CISPLATIN.*`, `figures/SuppFig7_network_centrality_profile.*`, `tables/SuppFig7_centrality_profile.tsv`.
- **SuppTable4**: `tables/SuppTable4_DL_hyperparameters.md`, `tables/SuppTable4_DL_selected.tsv`, `tables/SuppTable4_DL_grid_all_pairs.tsv`.

## 8. Figures produced
`B3_main_figure`, `B3_benchmark_heatmap`, `B3_supp_heatmap_Cindex`, `B3_supp_heatmap_LOO`, `B3_task_bars`, `B3_learning_curves`,
`B3_random_pathway_sets`, `B3_alignment_with_without`, `SuppFig6_replica_*`, `SuppFig7_replica_*`, `SuppFig7_network_centrality_profile`, and the Fig3 replicas.
Each has PNG, PDF and a `.caption.md`.

## 9. Limits to keep in mind
Pairs share some patients (for example the paper-data and current-data bladder sets), so CIs across pairs are approximate. Several sets have 5 to 19 lab models.
`learning.py` was not re-run (its output and `_data/learning.done` were already complete).
