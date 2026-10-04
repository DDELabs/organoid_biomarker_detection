# GDSC / CTRPv2 cell-line pre-training data

Built 2026-10-04 by `scratchpad/build_gdsc.py` (not in repo). Loader: `robust.obd.preclinical_sets.load_gdsc`.

## Files

| file | content | sha256 |
|---|---|---|
| `expression.tsv.gz` | GDSC RMA basal expression (Affymetrix HG-U219), 17,419 gene symbols x 1,014 cell lines; columns are COSMIC IDs; first column `gene`; probes averaged per symbol, duplicated array columns (`<id>.1`) averaged; rounded to 3 dp | `2db4f146983aadfee5dc402b9d4035b877ba3a7b7403274fedafeafdad1810a9` |
| `response.tsv.gz` | long table, 744,789 rows: `cosmic_id, cell_line, sanger_model_id, depmap_id, tissue, tissue_sub, tcga_label, gdsc_cancer_type, depmap_lineage, drug (generic upper-case via reference.common_drug_name), drug_name_gdsc, drug_id, ln_ic50, auc, max_conc_uM, dataset` (gzipped to keep the repo small; the loader reads either `.tsv.gz` or `.tsv`) | `484156f0036b98a503381cd4d81f96caf98e33257f5f9af7d39075aa3604ac92` |
| `cell_line_map.tsv` | 1,048 lines: COSMIC ID <-> DepMap ACH <-> Sanger SIDM <-> name <-> GDSC tissue / TCGA label / DepMap lineage, plus has_expression / has_gdsc1 / has_gdsc2 / has_ctrpv2 flags | `baccfd215578a7cbd8b136633b77026fb575172f8c3590cc676c02f3e2a01007` |

Datasets in `response.tsv.gz`: GDSC1 (333,161 rows, 378 drug names, 970 lines), GDSC2 (242,036 rows, 286 drug names, 969 lines), CTRPv2 (184,415 rows, 481 compounds, AUC only; 471 of 664 CTRP lines mapped to a COSMIC ID, 454 of them with responses).

## Raw sources (copies in `scratchpad/pre_raw/`)

| raw file | URL | sha256 |
|---|---|---|
| GDSC2_fitted_dose_response_27Oct23.xlsx | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/old/sensitivity/GDSC2_fitted_dose_response_27Oct23.xlsx | `0e99120a14a003dcbd3c27c020f1cfd35e38ee091067a2406f66c0d4dbc1c899` |
| GDSC1_fitted_dose_response_27Oct23.xlsx | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/old/sensitivity/GDSC1_fitted_dose_response_27Oct23.xlsx | `6dfad242263d93b6814cc51146b15346d3c69abe1aba834a0977ba9e3fd3b6a5` |
| Cell_Lines_Details.xlsx (GDSC) | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/old/Cell_Lines_Details.xlsx | `5520e779d25469f864de9cc5d42b99bf9632bd346607c0ee6a9ebad7ece17ac3` |
| Cell_line_RMA_proc_basalExp.txt (306 MB) | https://media.githubusercontent.com/media/hwr9912/pRRophetic/master/data/Cell_line_RMA_proc_basalExp.txt | `379c226502533e7f54eee84886ef2bad5f296164c1444cce7b511a140f4dda2d` |
| DepMap Model.csv (2026 release, has COSMICID/SangerModelID) | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/annotation/Model.csv | `ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8` |
| DepMap sample_info.csv (22Q1/PRISM 220228) | https://raw.githubusercontent.com/lasseignelab/Cancer_Signature_Reversion/main/data/DEPMAP_PRISM_220228/sample_info.csv | `cfc85cd9b2cb51f879dfe6f6ddaa6122eb3ba60758bc9a3ba6e76b5ca98fb93a` |
| CTRPv2.2 v22.data.auc_sensitivities.txt | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/DepMap/sensitivity/CTRPv2.2_2015_pub_CancerDisc_5_1210/v22.data.auc_sensitivities.txt | `c4174ac2deb7c16e64f3a898cf3234f94cef48cf82f71c57453e433b986f2976` |
| CTRPv2.2 v22.meta.per_cell_line.txt | same folder | `45787280209203b2d5eeb3c96a12077b86484790fea3d6a244e5a0c0b8d83c98` |
| CTRPv2.2 v22.meta.per_compound.txt | same folder | `666ca98fb7f955e15dee48c0200fb44c60986447c7491ae26f13c10c06862814` |
| DepMap "Drug sensitivity AUC (Sanger GDSC2)" (ACH IDs; join test only, not merged) | https://media.githubusercontent.com/media/zuhany/cell_lines_pv_26/main/data/raw/GDSC/sensitivity/Drug_sensitivity_AUC_%28Sanger_GDSC2%29_subsetted.csv | `a9bd955a1b0f9ba937f932f63fb9e353f83799899eaa41a268a34fb946331729` |

zuhany/cell_lines_pv_26 commit `301ddcb05294f380499690a86b43a098be7d9a44`. The 27Oct23 GDSC release has no COSMIC_ID column. COSMIC IDs come from SANGER_MODEL_ID via DepMap Model.csv / sample_info, with a fallback on the normalised cell-line name via Cell_Lines_Details.xlsx. All GDSC1/GDSC2 rows were mapped.

## Join stats (lines with GDSC RMA expression + response)

| drug | GDSC2 | GDSC1 | CTRPv2 (AUC) |
|---|---|---|---|
| FLUOROURACIL | 940 | 885 | 435 |
| CISPLATIN | 757 | 924 | — |
| GEMCITABINE | 933 | 933 | 390 |
| SORAFENIB | 930 | 385 | 425 |
| OXALIPLATIN | 939 | — | — |
| IRINOTECAN | 937 (SN-38: 930) | — (SN-38: 924) | — (SN-38: 390) |
| TEMOZOLOMIDE | 938 | 881 | 428 |
| DOXORUBICIN | — | 934 | 435 |
| PACLITAXEL | 930 | 383 | 424 |

Temozolomide in glioma lines (`tissues=['glioma']`, GDSC2 then GDSC1): 50 lines with expression. CTRPv2: 31. With `['glioma','GBM','LGG','CNS']` and DepMap CNS lineage: up to 54.
DepMap ACH GDSC2 AUC file: 946 of 947 ACH rows map to a COSMIC ID through `cell_line_map.tsv`, and 919 have RMA expression. Its TEMOZOLOMIDE AUC correlates r = 0.74 with the GDSC2 `AUC` here (n = 930). DepMap re-fits the curves, so do not mix the two.

## Direction check (lower = more sensitive)
`ln_ic50` is the natural log of IC50 in µM. `auc` is the area under the fitted viability curve (fraction, ≤1 for GDSC). A lower value means more sensitive for both. Checked on GDSC2: SLFN11 expression vs irinotecan ln_ic50 gives Spearman −0.49 (SLFN11-high lines are sensitive, as expected). MGMT vs temozolomide ln_ic50 gives +0.16 (MGMT-high lines are resistant, as expected).

## Caveats
- Some drugs have several DRUG_IDs within a dataset: GDSC1 cisplatin 1005/1496, gemcitabine 135/1393, doxorubicin 133/1386, SN-38 1490/1494; GDSC2 oxaliplatin 1089/1806. The loader takes the median.
- GDSC1 and GDSC2 use different assays and concentration ranges. `load_gdsc(datasets=('GDSC2','GDSC1'))` fills from GDSC1 only where GDSC2 has no value. For strict analyses, use one dataset.
- CTRPv2 AUC is an unnormalised area over a 16-point curve (median 13.7, max 27) and is not comparable to GDSC AUC. Not all CTRP lines have a COSMIC ID, so 193 CTRP lines are dropped.
- Expression is microarray RMA (log2), not RNA-seq. Gene symbols are mapped to canonical names by the loader.
- Tissue labels: `tissue` and `tissue_sub` are the GDSC descriptors, e.g. `nervous_system` / `glioma`, `digestive_system` / `stomach`. `tcga_label` is the TCGA code, `depmap_lineage` is the Oncotree lineage.
