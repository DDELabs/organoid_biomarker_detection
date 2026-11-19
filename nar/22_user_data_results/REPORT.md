# Results on the user's own lung, ovarian and liver data (nar/22_user_data_results)

All numbers come from files in `tables/` written by `code/`. Seeds fixed. One CPU process. Settings of the paper method: 500 bootstraps, 300 shuffles (engine of nar/02, unchanged). Figures: `figures/*.png|pdf` with `.caption.md`.
HR in the paper-method tables = predicted-less-sensitive half vs predicted-more-sensitive half (HR above 1 is the expected direction). In the joint (nar/03) tables HR is per 1 SD of predicted benefit (below 1 is expected).

## 1. What the uploaded files really contain (read this first)
| Set | What I found | What I think |
|---|---|---|
| Lung GSE165611 | 102 patient-derived cell lines (pleural fluid etc.; mixed histology), RNA-seq, 30 drugs. `median_response.txt` and `supp_table_S2b.txt` carry the same numbers (checked row by row). The column called IC50 holds **one value per cluster and drug** (the same number for every line in a cluster, repeated 8 times). Only 48 of 120 cluster x drug cells are listed; the other cells are blank. Your `luad drug.xlsx` writes 0 in blank cells (matches exactly). The `.ods` is an older version: 13 of 120 cells differ (e.g. vemurafenib C4 -0.30 vs 0.40) (`luad_ods_vs_S2b_cluster_values.tsv`). The expression file `geneID_expression_median.txt` equals the GEO series matrix (99.995% identical values; `input_checks.tsv`). | There is no line-level IC50 here. A model trained on it can only learn "which cluster is this line in". Drugs with the same cluster pattern (gemcitabine, paclitaxel, etoposide: C4 higher, others 0) give the same model. I filled blanks with 0 (as your Excel does) as the main version; the "listed only" version is the sensitivity check. If you have the line-level IC50s, they would be much better. |
| Ovary GSE229119 | 7 lines (Hey, OVCA432, SKOV3, 3 organoids 2414/2445/3008, plus MCF7 which is a **breast** line); 28 samples are time points. IC50 is identical across time points of a line (56 rows collapse to 14). Baseline = DMSO day 3 expression. Olaparib and niraparib IC50 agree across lines: Spearman 0.81 (p 0.050, 6 ovarian lines). | Main analysis uses 6 ovarian lines; 7 lines is a check. |
| Liver GSE188541 | The series matrix has **no expression rows**. GEO lists 6 samples of ONE iPSC-derived liver organoid line on chips (control, APAP, FIAU, tenofovir, tenofovir+inarigivir, plain organoid). Your IC50 table is for 3 other organoids (72.3, 2E, CC3), 12 mostly non-cancer drugs. | No organoid has both expression and IC50. A model cannot be trained. |
Loader (`nar/21_user_organoids/code/import_organoids.py`): used unchanged on reshaped inputs (`code/prepare_inputs.py`). Issues recorded: (1) it would take the natural log of any column named IC50, so the lung log2FC values are imported as metric "log2FC" with `--direction lower`; (2) it needs a plain symbol column, so the uniprot column was dropped; (3) the raw series matrices use Ensembl_symbol ids, not accepted, so the symbol file was used; (4) the liver set cannot be imported (no expression). Import reports: `_data/loader_out/imports/*/IMPORT_REPORT.md` (lung 19,505 of 19,861 genes scored, 673 pathways; ovary 671 pathways).
TCGA drug tables used are your uploaded files (LUAD, OV); LUSC from the repo. Treated counts (biotab): LUAD cisplatin 87, carboplatin 77, paclitaxel 44, pemetrexed 60, gemcitabine 25, etoposide 15, docetaxel 22. The combined training table v2 (used by `apply` and joint runs) records fewer drug-specific patients (e.g. LUAD gemcitabine 10), so counts differ between method tables.
Drugs without a patient group: navitoclax (no lung patient got it), and 25 more cluster-table drugs. Proximity: paper table for gemcitabine, paclitaxel, docetaxel, cisplatin (carboplatin uses cisplatin), camptothecin; etoposide, olaparib, niraparib recomputed on the paper's STRING network (`u22.proximity`, targets TOP2A/B, PARP1/2).

## 2. Lung
### 2a. Paper method (organoid-only), TCGA, `tables/paper_method_lung.tsv` (all rows, incl. bootstrap p, 27-test minimum, leave-one-line-out)
| Organoid drug -> patient drug | Cohort | n treated / events | log-rank p | HR (95% CI) | Untreated n, p, HR | Shuffle p |
|---|---|---|---|---|---|---|
| gemcitabine -> gemcitabine | LUAD | 25 / 16 | 0.60 | 1.31 (0.48-3.53) | 331, 0.021, 0.64 | 0.88 |
| gemcitabine -> gemcitabine | LUSC | 35 / 20 | 0.62 | 0.80 (0.33-1.93) | 354, 0.40, 1.14 | 0.62 |
| paclitaxel -> paclitaxel | LUAD | 44 / 21 | 0.047 | 0.40 (0.16-1.01) | 331, 0.52, 0.89 | 0.080 |
| paclitaxel -> paclitaxel | LUSC | 28 / 10 | 0.51 | 0.65 (0.18-2.33) | 354, 0.036, 0.72 | 0.48 |
| paclitaxel -> docetaxel | LUAD | 22 / 13 | 0.006 | 0.14 (0.03-0.69) | 331, 0.52, 0.89 | 0.063 |
| paclitaxel -> docetaxel | LUSC | 31 / 14 | 0.58 | 0.73 (0.23-2.29) | 354, 0.036, 0.72 | 0.69 |
| etoposide -> etoposide | LUAD | 15 / 4 | 0.51 | 2.20 (0.20-24.3) | 331, 0.030, 1.51 | 0.40 |
| etoposide -> etoposide | LUSC | 14 / 6 | 0.94 | 1.07 (0.17-6.61) | 354, 0.18, 0.81 | 0.99 |
| gemcitabine -> cisplatin (proxy) | LUAD / LUSC | 84 / 28; 72 / 23 | 0.91; 0.85 | 0.96 (0.46-2.02); 0.93 (0.41-2.10) | | 0.91; 0.93 |
| gemcitabine -> carboplatin (proxy) | LUAD / LUSC | 77 / 33; 57 / 23 | 0.72; 0.93 | 1.14 (0.56-2.30); 0.96 (0.41-2.27) | | 0.77; 0.93 |
| etoposide -> cisplatin (proxy) | LUAD / LUSC | 84 / 28; 72 / 23 | 0.87; 0.44 | 1.07 (0.51-2.25); 1.40 (0.60-3.27) | | 0.84; 0.21 |
| etoposide -> carboplatin (proxy) | LUAD / LUSC | 77 / 33; 57 / 23 | 0.24; 0.061 | 0.66 (0.33-1.33); 0.44 (0.18-1.07) | | 0.11; 0.12 |
What I think: the one nominal hit (paclitaxel model, docetaxel LUAD, p 0.006) points the opposite way to the method's expectation (HR 0.14) and its shuffle p is 0.063. 16 tests are shown, so one p below 0.05 is within what chance gives. Gemcitabine, etoposide and platinum proxies give HRs near 1. Leave-one-line-out Spearman of the organoid model is 0.52-0.81, but this mostly shows the model recovers the cluster.

### 2b. Same models on the lung cohorts of training table v2 (rank-scored pathways; `tables/apply_lung.tsv`, `apply_lung_listed.tsv`; figure 2)
Platinum proxy (top pathway, model drug gemcitabine; paclitaxel and etoposide give identical numbers because they share the same cluster pattern): 
| Cohort | n treated / events | log-rank p | HR (95% CI) | Shuffle p | Untreated/observation n, p |
|---|---|---|---|---|---|
| TCGA-LUAD | 124 / 43 | 0.57 | 0.84 (0.46-1.55) | 0.69 | 369, 0.88 |
| TCGA-LUSC | 103 / 30 | 0.89 | 0.95 (0.46-1.96) | 0.88 | 386, 0.023 |
| GSE14814 (JBR.10, cisplatin/vinorelbine vs observation) | 71 / 30 | 0.56 | 1.24 (0.60-2.53) | 0.70 | 62, 0.53 |
| GSE29013 | 34 / 12 | 0.89 | 0.92 (0.29-2.88) | 0.93 | 21, 0.85 |
| GSE68465 | 89 / 64 | 0.95 | 1.02 (0.62-1.67) | 0.96 | 340, 0.21 |
| GSE37745 | 29 / 23 | 0.56 | 0.78 (0.34-1.79) | 0.44 | 71, 0.52 |
| GSE42127 | 49 / 9 | 0.73 | 0.79 (0.21-2.97) | 0.81 | 124, 0.69 |
Camptothecin model (platinum proxy): GSE37745 p 0.047 (HR 0.44, 0.19-1.01, shuffle 0.048); GSE42127 p 0.051; GSE14814 p 0.13 (HR 1.74); others p above 0.2. Interaction (treated vs observation) in GSE14814: p 0.86 (gemcitabine model), 0.048 (camptothecin model). GSE207422 (anti-PD-1 + chemo, 24 patients, 9 responders): AUROC for predicting non-response 0.24, p 0.039, shuffle p 0.07-0.10 (0.55 with camptothecin). Anti-PD-1-only cohorts (GSE161537, GSE218989) were not used: the drug models do not apply to them. The listed-only version (gemcitabine, 70 lines with entries) gives the same picture (GSE14814 HR 1.77, p 0.12; shuffle 0.07).
What I think: direction and size vary from cohort to cohort and no result survives the shuffle check across cohorts; in the randomised JBR.10 trial the model drug does not separate benefit (interaction p 0.86). With cluster-level drug values this is the expected limit.

### 2c. Joint organoid + patient model (nar/03 code; `tables/joint_lung_pooled.tsv`, figure 3)
Leave-one-cohort-out, 7 cohorts, 499 platinum-treated patients, 233 events, shortlist = paper's cisplatin pathways (47):
| Lab drug | organoid-only HR/SD (p) | patient-only | joint |
|---|---|---|---|
| gemcitabine | 0.945 (0.41) | 0.881 (0.051) | 0.908 (0.16) |
| etoposide / paclitaxel | 0.941 (0.37) | 0.881 (0.051) | 0.924 (0.24) |
| camptothecin | 1.117 (0.094) | 0.881 (0.051) | 0.958 (0.53) |
What I think: adding the organoid information does not move the patient-only result (0.88). Matched-drug joint units (gemcitabine, paclitaxel, etoposide in TCGA only) did not run: fewer than two trainable cohorts (IndexError, `_data/joint/lung/failed.tsv`; earlier failures there were a table being rewritten by another job and were rerun). No shuffle test was run for the joint models (CPU limit); the paper-method shuffles above are the shuffle evidence.

### 2d. Clusters C1-C4 (`tables/lung_cluster_*.tsv`, figures 1 and 4)
- Sizes C1 13, C2 45, C3 19, C4 25; histology table in `lung_cluster_histology.tsv`. Pathway scores differ strongly by cluster (top: FGFR ligand-binding pathways, Kruskal p 7e-16, BH q 3e-13).
- A classifier on pathway scores recovers the cluster: balanced accuracy 0.89 (5-fold, 5 repeats) vs 0.25 chance and 0.25 for shuffled labels (shuffle p 0.01, 100 shuffles).
- Drug vs cluster: because values are defined by cluster, there is nothing left to test line by line; the cluster x drug values are in `luad_cluster_by_drug_value_listed.tsv` (e.g. gemcitabine C2 -0.2, C4 0.5; paclitaxel C4 0.3; etoposide C4 0.1; navitoclax C2 0.3, C3 -0.2, C4 -0.1).
- Predicted clusters in patients are similar everywhere (C2 about 37%, C3 28%, C4 20%, C1 13%).
- Expected drug value from predicted cluster, treated patients (OS): TCGA-LUAD paclitaxel 31/14 events p 0.56 HR 0.72 (0.23-2.23); LUAD platinum proxy with etoposide/paclitaxel values p 0.13 HR 0.62 (0.33-1.16); GSE14814 camptothecin values interaction p 0.048; GSE29013 and GSE42127 HR 0.2-0.6 with p 0.05-0.3 (shuffle p 0.07-0.2). Full table with 500 shuffles: `lung_cluster_expected_value_survival.tsv`.
What I think: the cluster is readable from expression, but I see no consistent link between the predicted cluster and survival after chemotherapy.

## 3. Ovary (`tables/ovary_*.tsv`, `apply_ovary*.tsv`, figures 5-7)
With 6 ovarian lines (7 with MCF7) a model can only be exploratory. With 6 lines the strongest correlation of any pathway with olaparib ln IC50 reaches rho = 1 by chance in 69% of the 720 possible orderings (family-wise p 0.69; niraparib 0.88; 7 lines 0.21 and 0.34). So no single pathway is reliably tied to IC50. DNA-repair pathways vs olaparib IC50: base excision repair rho -0.71, DNA repair -0.37 (`ovary_dna_repair_pathways_vs_ic50.tsv`).
Paper-method (PARP1/2 proximity shortlist; 6 lines) in TCGA-OV (`paper_method_ovary.tsv`): olaparib -> carboplatin 348 / 218 events, p 0.020, HR 0.73 (0.56-0.95), untreated 37 p 0.68, shuffle p 0.096; niraparib -> carboplatin p 0.88, HR 0.98, shuffle 1.0; paclitaxel p 0.20 and 0.48. On the v2 cohorts (olaparib model, platinum): TCGA-OV OS p 0.051, HR 0.77 (0.59-1.00), shuffle 0.13; GSE9891 p 0.005, HR 0.58 (0.39-0.85), shuffle 0.008; GSE63885 p 0.99. Both HR below 1 mean predicted-less-sensitive did better, the opposite of the expected direction. Niraparib: all p above 0.3. Response cohorts GSE15622 / GSE51373: AUROC 0.57 / 0.37 (p 0.69 / 0.26) for the signature.
Correlation signature (top 10 positive minus top 10 negative pathways; exact shuffle over all orders) in platinum-treated patients: TCGA-OV OS p 0.65 (HR 1.06, 0.82-1.38), PFS p 0.60; GSE9891 p 0.54; GSE63885 olaparib p 0.004, HR 2.09 (1.24-3.50), exact permutation p 0.014 (this one points in the expected direction; niraparib p 0.084); other cohorts p above 0.28. 
BRCA1/2 somatic mutation (17 of 288 TCGA-OV patients with mutation data): signature score vs mutation p 0.13-0.63 (6 lines); mutated patients live longer on platinum: OS HR 0.47 (0.22-1.02, p 0.049), PFS HR 0.42 (0.21-0.82, p 0.009), 16 mutated of 253 (germline BRCA is not public).
What I think: two cohorts give p below 0.05 for the olaparib model and signature but in opposite directions across cohorts, so I would not read them as evidence yet. Olaparib was given to no TCGA patient, so platinum is only a stand-in.

## 4. Liver (`tables/liver_*.tsv`, figure 8)
Valid analyses only: the 3 organoids differ in drug profile (rank agreement between organoid pairs across 12 drugs: 2E vs 72.3 rho -0.10, 2E vs CC3 0.25, 72.3 vs CC3 -0.34; all p above 0.28), drug spread `liver_drug_spread.tsv` (thioguanine ranges 2.9 to 37). Of the 12 drugs only tamoxifen has a TCGA-LIHC patient (1). No model, no survival test. What I think: these data describe liver toxicity screening; linking them to hepatocellular cancer outcome needs line-matched expression and cancer drugs (e.g. sorafenib).

## 5. Limits
Cluster-level drug values (not per line); few treated patients for drug-specific lung runs (14-44); proxies for platinum; 500 bootstraps and 300 shuffles (not 10,000/1,000); no shuffle for joint models; many tests, nominal p values. Code: `code/` (run order in `run_all_chain.sh`, `run_rerun.sh`, `run_rerun2.sh`). Inputs: `data/user_upload/*` (sizes in `tables/input_checks.tsv`); soft files in `nar/_cache/user_soft/`. No file above 50 MB created (`_data` 41 MB).
