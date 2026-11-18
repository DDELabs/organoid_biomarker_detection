# Results index: every figure and table under nar/

354 figures (each also has a PDF and a .caption.md next to it) and 910 tables. Paths are relative to `nar/`. Descriptions come from the caption's first line (figures) or the column names (tables). Built by `10_completeness/code/c09_index.py`; the same list is in `10_completeness/tables/results_index.tsv`. Items being made by the B3 agent (`07_B3_models/`) will appear after a re-run.

## Part A. By paper section

### Fig. 1 (method, schematics, data map)

**colon / rectum**

- `08_depmap/figures/fig1_workflow_bh3_vs_5fu.png` — **Figure 1. Testing the organoid 5-FU biomarker in colorectal cancer cell lines.**

**all sets / cross-cancer**

- `08_depmap/figures/S1_data_matching_flow.png` — **Supplementary Figure S1. How many cell lines survive each matching step.** From all lines in each screen, to
- `08_depmap/tables/S1_data_matching_flow.tsv` — table: screen, all_lines_in_screen, lines_with_5FU, 5FU_and_RNAseq, colorectal_5FU_RNAseq, colorectal_plus_CRISPR, colorectal_plus_mutation_data
- `09_registered_test/figures/fig1_preregistration_flow.png` — **Figure 1. How the test was pre-registered.** Each box is one step, top to bottom. The 35 pathway hypotheses were
- `10_completeness/figures/Fig1a_workflow_schematic.png` — Fig. 1a (redrawn): the whole method in one picture
- `10_completeness/figures/Fig1c_data_map_all_sets.png` — Fig. 1c (new): map of every analysis
- `10_completeness/tables/Fig1c_data_map.tsv` — table: set, cohort, drug, n_models, treated_n

### Fig. 2 (network pathway biomarker, survival)

**PDX (Gao) / breast / melanoma**

- `02_organoid_only/COAD_isella_pdx/figures/Fig2d_main_composite_COAD_isella_pdx.png` — Main figure: Colorectal PDX (Isella 2017, cetuximab)
- `02_organoid_only/COAD_isella_pdx/tables/cutoff_sensitivity_COAD_CETUXIMAB.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PDX_novartis_gao2015/MAIN_PDX_novartis_gao2015_forest.png` — Novartis PDX (Gao 2015), several tumour types: all drugs (rank-1 Ridge pathway)
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_PACLITAXEL_full.png` — Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, paclitaxel
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_PACLITAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA paclitaxel
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_TAMOXIFEN_full.png` — Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, tamoxifen
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_TAMOXIFEN_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA tamoxifen
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_TRASTUZUMAB_full.png` — Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA, trastuzumab
- `02_organoid_only/PDX_novartis_gao2015/figures/BRCA_TRASTUZUMAB_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (BRCA, Gao 2015) -> TCGA-BRCA trastuzumab
- `02_organoid_only/PDX_novartis_gao2015/figures/COAD_CETUXIMAB_full.png` — Novartis PDX (CRC, Gao 2015) -> TCGA-COAD, cetuximab
- `02_organoid_only/PDX_novartis_gao2015/figures/COAD_CETUXIMAB_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (CRC, Gao 2015) -> TCGA-COAD cetuximab
- `02_organoid_only/PDX_novartis_gao2015/figures/COAD_FLUOROURACIL_full.png` — Novartis PDX (CRC, Gao 2015) -> TCGA-COAD, fluorouracil
- `02_organoid_only/PDX_novartis_gao2015/figures/COAD_FLUOROURACIL_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (CRC, Gao 2015) -> TCGA-COAD fluorouracil
- `02_organoid_only/PDX_novartis_gao2015/figures/LUAD_ERLOTINIB_full.png` — Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD, erlotinib
- `02_organoid_only/PDX_novartis_gao2015/figures/LUAD_ERLOTINIB_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD erlotinib
- `02_organoid_only/PDX_novartis_gao2015/figures/LUAD_PACLITAXEL_full.png` — Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD, paclitaxel
- `02_organoid_only/PDX_novartis_gao2015/figures/LUAD_PACLITAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (NSCLC, Gao 2015) -> TCGA-LUAD paclitaxel
- `02_organoid_only/PDX_novartis_gao2015/figures/PAAD_GEMCITABINE_full.png` — Novartis PDX (PDAC, Gao 2015) -> TCGA-PAAD, gemcitabine
- `02_organoid_only/PDX_novartis_gao2015/figures/PAAD_GEMCITABINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (PDAC, Gao 2015) -> TCGA-PAAD gemcitabine
- `02_organoid_only/PDX_novartis_gao2015/figures/SKCM_DACARBAZINE_full.png` — Novartis PDX (CM, Gao 2015) -> TCGA-SKCM, dacarbazine
- `02_organoid_only/PDX_novartis_gao2015/figures/SKCM_DACARBAZINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Novartis PDX (CM, Gao 2015) -> TCGA-SKCM dacarbazine
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/PDX_novartis_gao2015/tables/PDX_summary.tsv` — table: set, set_label, project, drug, model_drug, n_models, proximity_source
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/figures/Fig2d_main_composite_COAD_vandewetering.png` — Main figure: Colorectal organoids (van de Wetering 2015)
- `02_organoid_only/COAD_vandewetering/figures/Fig2d_supp_composite_READ_and_pooled.png` — Main figure: Colorectal organoids (van de Wetering 2015)
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COADREAD_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COADREAD_IRINOTECAN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COADREAD_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_CETUXIMAB.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_IRINOTECAN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_authors_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_authors_IRINOTECAN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_COAD_authors_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_READ_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/COAD_vandewetering/tables/cutoff_sensitivity_READ_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_shi/tables/cutoff_sensitivity_PAAD_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_tiriac/tables/cutoff_sensitivity_PAAD_FLUOROURACIL.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `03_organoid_plus_patient/figures/main_colorectal.png` — **Colorectal, 5-FU-based chemotherapy: a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/CRC_5FU__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `05_A3_B1/figures/B1_Fig1_forest_colorectal.png` — **External validation, colorectal.** Each row = one frozen organoid-trained signature (right label; [rank1] = the
- `08_depmap/tables/A_e_5fu_rank_summary.tsv` — table: screen, drugs, 5FU_colon_rho, 5FU_colon_rank_of, 5FU_colon_percentile, 5FU_pancancer_rho, 5FU_pancancer_rank_of
- `15_fixes/backup/03/figures/main_colorectal.png` — **Colorectal, 5-FU-based chemotherapy: a model that learns from lab models and patients together.**
- `15_fixes/backup/03/tables/CRC_5FU__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/CRC_5FU__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**bladder**

- `02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/figures/cross_set_forest.png` — All sets and drugs in one plot
- `02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv` — table: folder, set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/MAIN_BLCA_lee2018_current_forest.png` — Bladder organoids (Lee 2018), current GDC data: all drugs (rank-1 Ridge pathway)
- `02_organoid_only/BLCA_lee2018/MAIN_BLCA_lee2018_paper_forest.png` — Bladder organoids (Lee 2018), paper source-data NES: all drugs (rank-1 Ridge pathway)
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_CARBOPLATIN_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, carboplatin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_CARBOPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA carboplatin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_CISPLATIN_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, cisplatin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA cisplatin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_DOXORUBICIN_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, doxorubicin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_DOXORUBICIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA doxorubicin
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_GEMCITABINE_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, gemcitabine
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_GEMCITABINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA gemcitabine
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_METHOTREXATE_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, methotrexate
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_METHOTREXATE_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA methotrexate
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_PACLITAXEL_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, paclitaxel
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_PACLITAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA paclitaxel
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_VINBLASTINE_full.png` — Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA, vinblastine
- `02_organoid_only/BLCA_lee2018/current/figures/BLCA_VINBLASTINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), current GDC data -> TCGA-BLCA vinblastine
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_paper_biomarker_KM.png` — The pathway the paper reported for bladder cisplatin
- `02_organoid_only/BLCA_lee2018/extras/tables/paper_biomarker_amino_acid_KM.tsv` — table: version, group, ridge_rank_of_AA, coef, n, events, n_nonresp
- `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` — Bladder organoids (Lee 2018), paper source-data NES -> TCGA-BLCA, cisplatin
- `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Bladder organoids (Lee 2018), paper source-data NES -> TCGA-BLCA cisplatin
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/figures/BLCA_CISPLATIN_full.png` — all 9 lines (7 patients), as main run -> TCGA-BLCA, cisplatin
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/figures/BLCA_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, all 9 lines (7 patients), as main run -> TCGA-BLCA cisplatin
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/figures/BLCA_CISPLATIN_full.png` — 7 lines, one per patient (recurrence lines SCBO-3_2, SCBO-11_2 dropped) -> TCGA-BLCA, cisplatin
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/figures/BLCA_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, 7 lines, one per patient (recurrence lines SCBO-3_2, SCBO-11_2 dropped) -> TCGA-BLCA cisplatin
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/BLCA_lee2018/tables/BLCA_summary.tsv` — table: set, set_label, project, drug, model_drug, n_models, proximity_source
- `03_organoid_plus_patient/figures/main_bladder.png` — **Bladder, cisplatin-based chemotherapy: a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `05_A3_B1/figures/B1_Fig2_forest_bladder.png` — **External validation, bladder.** Each row = one frozen organoid-trained signature (right label; [rank1] = the
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**pancreas**

- `02_organoid_only/PAAD_shi/figures/Fig2d_main_composite_PAAD_shi.png` — Main figure: Pancreatic organoids (Shi 2022)
- `02_organoid_only/PAAD_shi/tables/cutoff_sensitivity_PAAD_GEMCITABINE.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_shi/tables/cutoff_sensitivity_PAAD_IRINOTECAN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_shi/tables/cutoff_sensitivity_PAAD_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_tiriac/figures/Fig2d_main_composite_PAAD_tiriac.png` — Main figure: Pancreatic organoids (Tiriac 2018)
- `02_organoid_only/PAAD_tiriac/tables/cutoff_sensitivity_PAAD_GEMCITABINE.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_tiriac/tables/cutoff_sensitivity_PAAD_OXALIPLATIN.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/PAAD_tiriac/tables/cutoff_sensitivity_PAAD_SN-38.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/SARC_alshihabi2024/MAIN_SARC_alshihabi2024_forest.png` — Sarcoma organoids (Al Shihabi 2024): all drugs (rank-1 Ridge pathway)
- `02_organoid_only/SARC_alshihabi2024/figures/SARC_DOCETAXEL_full.png` — Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC, docetaxel
- `02_organoid_only/SARC_alshihabi2024/figures/SARC_DOCETAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC docetaxel
- `02_organoid_only/SARC_alshihabi2024/figures/SARC_GEMCITABINE_full.png` — Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC, gemcitabine
- `02_organoid_only/SARC_alshihabi2024/figures/SARC_GEMCITABINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Sarcoma organoids (Al Shihabi 2024) -> TCGA-SARC gemcitabine
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_summary.tsv` — table: set, set_label, project, drug, model_drug, n_models, proximity_source
- `03_organoid_plus_patient/figures/main_pancreas.png` — **Pancreas, gemcitabine: a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**liver**

- `02_organoid_only/LIHC_broutier/figures/Fig2d_main_composite_LIHC_broutier.png` — Main figure: Liver cancer organoids (Broutier 2017)
- `02_organoid_only/LIHC_broutier/tables/cutoff_sensitivity_LIHC_SORAFENIB.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/LIHC_ji/figures/Fig2d_main_composite_LIHC_ji.png` — Main figure: Liver cancer organoids (Ji 2023)
- `02_organoid_only/LIHC_ji/tables/cutoff_sensitivity_LIHC_SORAFENIB.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `03_organoid_plus_patient/figures/main_liver.png` — **Liver, sorafenib: a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `05_A3_B1/figures/B1_Fig4_forest_liver_and_lung.png` — **External validation, liver and lung.** Each row = one frozen organoid-trained signature (right label; [rank1] = the
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**ovarian**

- `02_organoid_only/OV_vias2023/MAIN_OV_vias2023_forest.png` — Ovarian HGSOC organoids (Vias 2023): all drugs (rank-1 Ridge pathway)
- `02_organoid_only/OV_vias2023/figures/OV_CARBOPLATIN_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, carboplatin
- `02_organoid_only/OV_vias2023/figures/OV_CARBOPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV carboplatin
- `02_organoid_only/OV_vias2023/figures/OV_CISPLATIN_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, cisplatin
- `02_organoid_only/OV_vias2023/figures/OV_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV cisplatin
- `02_organoid_only/OV_vias2023/figures/OV_DOCETAXEL_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, docetaxel
- `02_organoid_only/OV_vias2023/figures/OV_DOCETAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV docetaxel
- `02_organoid_only/OV_vias2023/figures/OV_DOXORUBICIN_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, doxorubicin
- `02_organoid_only/OV_vias2023/figures/OV_DOXORUBICIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV doxorubicin
- `02_organoid_only/OV_vias2023/figures/OV_GEMCITABINE_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, gemcitabine
- `02_organoid_only/OV_vias2023/figures/OV_GEMCITABINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV gemcitabine
- `02_organoid_only/OV_vias2023/figures/OV_PACLITAXEL_full.png` — Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV, paclitaxel
- `02_organoid_only/OV_vias2023/figures/OV_PACLITAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Ovarian HGSOC organoids (Vias 2023) -> TCGA-OV paclitaxel
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/OV_vias2023/tables/OV_summary.tsv` — table: set, set_label, project, drug, model_drug, n_models, proximity_source
- `03_organoid_plus_patient/figures/main_ovarian.png` — **Ovarian, paclitaxel (with platinum): a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/OV_PLATINUM__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/OV_PLATINUM__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `05_A3_B1/figures/B1_Fig3_forest_ovarian.png` — **External validation, ovarian.** Each row = one frozen organoid-trained signature (right label; [rank1] = the
- `15_fixes/backup/03/tables/OV_PACLITAXEL__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/OV_PACLITAXEL__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/OV_PLATINUM__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/OV_PLATINUM__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**lung**

- `02_organoid_only/LUAD_GSE165611_substitute/MAIN_LUAD_GSE165611_substitute_forest.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611): all drugs (rank-1 Ridge pathway)
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_CARBOPLATIN_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, carboplatin
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_CARBOPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD carboplatin
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_CISPLATIN_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, cisplatin
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_CISPLATIN_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD cisplatin
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_DOCETAXEL_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, docetaxel
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_DOCETAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD docetaxel
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_ERLOTINIB_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, erlotinib
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_ERLOTINIB_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD erlotinib
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_ETOPOSIDE_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, etoposide
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_ETOPOSIDE_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD etoposide
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_GEMCITABINE_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, gemcitabine
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_GEMCITABINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD gemcitabine
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_PACLITAXEL_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, paclitaxel
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_PACLITAXEL_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD paclitaxel
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_PEMETREXED_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, pemetrexed
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_PEMETREXED_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD pemetrexed
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_VINORELBINE_full.png` — Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD, vinorelbine
- `02_organoid_only/LUAD_GSE165611_substitute/figures/LUAD_VINORELBINE_supp_zcutoff_KM.png` — Cut-off sensitivity, Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611) -> TCGA-LUAD vinorelbine
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_rank1_tests.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_summary.json` — JSON summary: set, set_label, project, drug, model_drug, n_models
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_zcutoff_sensitivity.tsv` — table: cutoff, n_proximal, pathway, treated_p, treated_hr, treated_hr_lo, treated_hr_hi
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_summary.tsv` — table: set, set_label, project, drug, model_drug, n_models, proximity_source
- `03_organoid_plus_patient/figures/main_lung.png` — **Lung, platinum (cell-line proxy): a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/LUNG_PLATINUM__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/LUNG_PLATINUM__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**sarcoma**

- `03_organoid_plus_patient/figures/main_sarcoma.png` — **Sarcoma, gemcitabine (+ docetaxel): a model that learns from lab models and patients together.**
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist_foldalign__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__all__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__shortlist__interaction.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison

**all sets / cross-cancer**

- `02_organoid_only/figures/cross_set_cutoff_sensitivity.png` — Cut-off sensitivity across all sets
- `02_organoid_only/figures/cross_set_forest.png` — Cross-set forest plot
- `02_organoid_only/tables/cross_set_best_direction.tsv` — table: folder, cohort, model_drug, best_setting, best_p, best_HR, best_HR_lo
- `02_organoid_only/tables/cross_set_cutoff_sensitivity.tsv` — table: cutoff, n_pathways, rank1_pathway, treated_p, treated_HR, treated_HR_lo, treated_HR_hi
- `02_organoid_only/tables/cross_set_summary.tsv` — table: set, folder, cohort, model_drug, patient_drug, proximity_from, metric
- `03_organoid_plus_patient/figures/forest_cross_cancer.png` — **Cross-cancer summary (network-proximal pathways).**
- `03_organoid_plus_patient/tables/interaction_all.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `03_organoid_plus_patient/tables/summary_all.tsv` — table: cohort, model, outcome, n, events, C, logrank_p
- `03_organoid_plus_patient/tables/summary_main.tsv` — table: cancer, unit, drug, features, n_features, model, outcome
- `05_A3_B1/figures/A3_Fig1_master_forest.png` — **Master forest plot of all 57 organoid-only analyses.** Each row is one lab-model set x drug x TCGA cohort
- `05_A3_B1/figures/A3_Fig4_interaction_forest.png` — **Is the score linked to survival more strongly in drug-treated than in untreated patients?** For each analysis a
- `05_A3_B1/figures/B1_Fig5_KM_main_cohorts.png` — **Kaplan-Meier curves in the main external cohorts.** Left column = patients who got the drug; right column = the
- `05_A3_B1/figures/B1_Fig5b_KM_main_cohorts.png` — (Second part of Fig. 5: cohorts 4-6, panels continue the lettering.)
- `05_A3_B1/tables/A3_adjusted_cox_epv_summary.tsv` — table: set, n, same_direction
- `05_A3_B1/tables/B1_summary_counts.tsv` — table: group, tests, treated_p<0.05, treated_p<0.05_expected_dir, treated_expected_dir, with_comparator, interaction_p<0.05
- `07_B3_models/figures/B3_main_figure.png` — **Main figure: the paper's network method next to its comparison methods and modern models, on identical splits.**
- `07_B3_models/tables/cellline_cv_summary.tsv` — table: task, model, n_units, mean_rho, lo, hi, share_units_rho_pos
- `08_depmap/figures/fig5_anchor_stage0_frozen_hypotheses.png` — **Figure 5. ANCHOR Stage 0: pre-registered, outcome-free mechanism hypotheses.**
- `08_depmap/tables/B_stage0_summary.tsv` — table: drug, screens, pathways_tested, S1_pass, S2_pass, S3_pass, two_plus_sources
- `09_registered_test/figures/fig2_stage1_pooled_interaction_forest.png` — **Figure 2. Stage 1: does the frozen score predict a worse outcome specifically under the drug?** Each square is a
- `09_registered_test/figures/fig3_km_high_vs_low_score.png` — **Figure 3. Kaplan–Meier curves of progression-free interval (PFI) by frozen score, per drug (TCGA, pooled
- `09_registered_test/figures/fig6_per_cancer_forest_meta.png` — **Figure 6. Stage 1 per cancer type, with a random-effects meta-analysis.** For each drug, the score x drug
- `09_registered_test/figures/fig7_external_cohort_forest.png` — **Figure 7. Stage 2: the frozen scores in outside patient cohorts.** **a**, cohorts with a comparison arm. Each line
- `09_registered_test/figures/fig8_summary_heatmap.png` — **Figure 8. Summary: every test of every frozen score, drug by cohort.** Each cell is a z-value (estimate divided by
- `09_registered_test/tables/stage1_gate_G2.tsv` — table: drug, G2_1_holm_and_direction, G2_2_purity_ipw_robust, G2_3_negative_controls_null, G2_4_weaker_in_comparator, G2_pass
- `09_registered_test/tables/stage1_ipw_weights.tsv` — table: drug, w_min, w_max, w_mean, ess, n
- `09_registered_test/tables/stage1_km_logrank.tsv` — table: drug, group, n, n_high, events, logrank_chi2, logrank_p
- `09_registered_test/tables/stage1_meta_analysis.tsv` — table: drug, k, HR_re, lo_re, hi_re, p_re, HR_fe
- `09_registered_test/tables/stage1_negative_controls.tsv` — table: drug, control, n, HR_int, lo, hi, z
- `09_registered_test/tables/stage1_per_cancer.tsv` — table: drug, project, n, n_exposed, events_exposed, events_unexposed, fitted
- `09_registered_test/tables/stage1_power.tsv` — table: drug, endpoint, test, n, events_exposed, events_unexposed, expected_se
- `09_registered_test/tables/stage1_primary.tsv` — table: drug, projects, n_projects, PFI_HR_int, PFI_lo, PFI_hi, PFI_z
- `09_registered_test/tables/stage1_sensitivity.tsv` — table: drug, analysis, n, events_exposed, events_unexposed, HR_int, lo
- `09_registered_test/tables/stage1_sensitivity_no_missing_start.tsv` — table: drug, analysis, n, n_exposed, n_dropped, events_exposed, events_unexposed
- `09_registered_test/tables/stage1_specificity_matrix.tsv` — table: exposure, score, HR_int, lo, hi, z, p
- `09_registered_test/tables/stage2_controlled.tsv` — table: drug, cohort, drug_arm, comparator, endpoint, kind, role
- `09_registered_test/tables/stage2_gate_G3.tsv` — table: drug, n_controlled_primary, n_adequately_powered, n_one_sided_p_lt_0.025, n_same_direction, G3_powered_hit, G3_all_same_direction
- `09_registered_test/tables/stage2_pooled.tsv` — table: drug, kind, k, cohorts, FE, FE_lo, FE_hi
- `09_registered_test/tables/stage2_sensitivity_GSE39582_harmonised.tsv` — table: definition, endpoint, n, n_drug_arm, n_comparator, events_drug_arm, events_comparator
- `09_registered_test/tables/stage2_single_arm_supportive.tsv` — table: drug, cohort, arm, endpoint, kind, n, events
- `09_registered_test/tables/summary_drug_by_cohort.tsv` — table: drug, col, z, p
- `10_completeness/tables/verify_headline_biomarkers.tsv` — table: key, label, bundle, drug, pathway, ridge_coef, treated_n
- `14_gaps_stats/figures/G10_meta_forest.png` — **Figure G10. Pooled effect of the organoid scores across the 57 analyses.**
- `14_gaps_stats/tables/G18_correlation_summary_by_cohort.tsv` — table: project, count, mean, std, min, 10%, 25%
- `15_fixes/backup/03/figures/forest_cross_cancer.png` — **Cross-cancer summary (network-proximal pathways).**
- `15_fixes/backup/03/tables/interaction_all.tsv` — table: cohort, model, outcome, n_treated, n_comparison, events_treated, events_comparison
- `15_fixes/backup/03/tables/summary_all.tsv` — table: cohort, model, outcome, n, events, C, logrank_p
- `15_fixes/backup/03/tables/summary_main.tsv` — table: cancer, unit, drug, features, n_features, model, outcome
- `15_fixes/backup/05/tables/B1_summary_counts.tsv` — table: group, tests, treated_p<0.05, treated_p<0.05_expected_dir, treated_expected_dir, with_comparator, interaction_p<0.05
- `15_fixes/backup/09/tables/stage1_gate_G2.tsv` — table: drug, G2_1_holm_and_direction, G2_2_purity_ipw_robust, G2_3_negative_controls_null, G2_4_weaker_in_comparator, G2_pass
- `15_fixes/backup/09/tables/stage1_ipw_weights.tsv` — table: drug, w_min, w_max, w_mean, ess, n
- `15_fixes/backup/09/tables/stage1_km_logrank.tsv` — table: drug, group, n, n_high, events, logrank_chi2, logrank_p
- `15_fixes/backup/09/tables/stage1_meta_analysis.tsv` — table: drug, k, HR_re, lo_re, hi_re, p_re, HR_fe
- `15_fixes/backup/09/tables/stage1_negative_controls.tsv` — table: drug, control, n, HR_int, lo, hi, z
- `15_fixes/backup/09/tables/stage1_per_cancer.tsv` — table: drug, project, n, n_exposed, events_exposed, events_unexposed, fitted
- `15_fixes/backup/09/tables/stage1_power.tsv` — table: drug, endpoint, test, n, events_exposed, events_unexposed, expected_se
- `15_fixes/backup/09/tables/stage1_primary.tsv` — table: drug, projects, n_projects, PFI_HR_int, PFI_lo, PFI_hi, PFI_z
- `15_fixes/backup/09/tables/stage1_sensitivity.tsv` — table: drug, analysis, n, events_exposed, events_unexposed, HR_int, lo
- `15_fixes/backup/09/tables/stage1_specificity_matrix.tsv` — table: exposure, score, HR_int, lo, hi, z, p
- `15_fixes/backup/09/tables/stage2_controlled.tsv` — table: drug, cohort, drug_arm, comparator, endpoint, kind, role
- `15_fixes/backup/09/tables/stage2_gate_G3.tsv` — table: drug, n_controlled_primary, n_adequately_powered, n_one_sided_p_lt_0.025, n_same_direction, G3_powered_hit, G3_all_same_direction
- `15_fixes/backup/09/tables/stage2_pooled.tsv` — table: drug, kind, k, cohorts, FE, FE_lo, FE_hi
- `15_fixes/backup/09/tables/stage2_single_arm_supportive.tsv` — table: drug, cohort, arm, endpoint, kind, n, events
- `15_fixes/backup/09/tables/summary_drug_by_cohort.tsv` — table: drug, col, z, p

### Fig. 3 / Supp 6-9 (comparison feature sets, centrality, deep learning)

**PDX (Gao) / breast / melanoma**

- `07_B3_models/figures/fig3_replica/Fig3_PDX_BRCA_PACLITAXEL.png` — **Paper Fig. 3 (b-f) replica for BRCA PDX | paclitaxel.** Each panel trains one model on the 37 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_BRCA_TAMOXIFEN.png` — **Paper Fig. 3 (b-f) replica for BRCA PDX | tamoxifen.** Each panel trains one model on the 37 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_BRCA_TRASTUZUMAB.png` — **Paper Fig. 3 (b-f) replica for BRCA PDX | trastuzumab.** Each panel trains one model on the 36 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_CM_DACARBAZINE.png` — **Paper Fig. 3 (b-f) replica for melanoma PDX | dacarbazine.** Each panel trains one model on the 31 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_CRC_FLUOROURACIL.png` — **Paper Fig. 3 (b-f) replica for CRC PDX (Gao) | fluorouracil.** Each panel trains one model on the 41 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_NSCLC_ERLOTINIB.png` — **Paper Fig. 3 (b-f) replica for NSCLC PDX | erlotinib.** Each panel trains one model on the 25 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_NSCLC_PACLITAXEL.png` — **Paper Fig. 3 (b-f) replica for NSCLC PDX | paclitaxel.** Each panel trains one model on the 24 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PDX_PDAC_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for PDAC PDX | gemcitabine.** Each panel trains one model on the 33 lab models and applies it to the

**colon / rectum**

- `07_B3_models/figures/SuppFig6_replica_COAD_paper_FLUOROURACIL.png` — **Paper Supp Fig 6 replica (COAD organoids, fluorouracil; paper source data).** Instead of the
- `07_B3_models/figures/fig3_replica/Fig3_COADREAD_vdW_CETUXIMAB.png` — **Paper Fig. 3 (b-f) replica for CRC organoids | cetuximab.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COADREAD_vdW_FLUOROURACIL.png` — **Paper Fig. 3 (b-f) replica for CRC organoids | fluorouracil.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COADREAD_vdW_IRINOTECAN.png` — **Paper Fig. 3 (b-f) replica for CRC organoids | irinotecan.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COADREAD_vdW_OXALIPLATIN.png` — **Paper Fig. 3 (b-f) replica for CRC organoids | oxaliplatin.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COAD_isella_CETUXIMAB.png` — **Paper Fig. 3 (b-f) replica for CRC PDX (Isella) | cetuximab.** Each panel trains one model on the 192 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COAD_paper_FLUOROURACIL.png` — **Paper Fig. 3 (b-f) replica for COAD organoids (paper data) | fluorouracil.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COAD_paper_IRINOTECAN.png` — **Paper Fig. 3 (b-f) replica for COAD organoids (paper data) | irinotecan.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_COAD_paper_OXALIPLATIN.png` — **Paper Fig. 3 (b-f) replica for COAD organoids (paper data) | oxaliplatin.** Each panel trains one model on the 19 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PAAD_shi_FLUOROURACIL.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Shi) | fluorouracil.** Each panel trains one model on the 39 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PAAD_tiriac_FLUOROURACIL.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Tiriac) | fluorouracil.** Each panel trains one model on the 36 lab models and applies it to the

**bladder**

- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png` — Other ways of choosing features, cisplatin-treated TCGA-BLCA (paper Fig. 3b-d, Supp. Fig. 4b, 8b, 9)
- `02_organoid_only/BLCA_lee2018/extras/tables/comparison_feature_sets_treated_KM.tsv` — table: model, n_features, survival, n, events, n_nonresp, n_resp
- `07_B3_models/figures/SuppFig7_replica_BLCA_paper_CISPLATIN.png` — **Paper Supp Fig 7 replica (BLCA organoids, cisplatin; paper source data).** Instead of the
- `07_B3_models/figures/fig3_replica/Fig3_BLCA_current_CARBOPLATIN.png` — **Paper Fig. 3 (b-f) replica for BLCA organoids | carboplatin.** Each panel trains one model on the 9 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_BLCA_current_CISPLATIN.png` — **Paper Fig. 3 (b-f) replica for BLCA organoids | cisplatin.** Each panel trains one model on the 9 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_BLCA_current_DOXORUBICIN.png` — **Paper Fig. 3 (b-f) replica for BLCA organoids | doxorubicin.** Each panel trains one model on the 9 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_BLCA_current_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for BLCA organoids | gemcitabine.** Each panel trains one model on the 11 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_BLCA_paper_CISPLATIN.png` — **Paper Fig. 3 (b-f) replica for BLCA organoids (paper data) | cisplatin.** Each panel trains one model on the 9 lab models and applies it to the

**pancreas**

- `07_B3_models/figures/fig3_replica/Fig3_PAAD_shi_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Shi) | gemcitabine.** Each panel trains one model on the 39 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PAAD_shi_OXALIPLATIN.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Shi) | oxaliplatin.** Each panel trains one model on the 39 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PAAD_tiriac_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Tiriac) | gemcitabine.** Each panel trains one model on the 36 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_PAAD_tiriac_OXALIPLATIN.png` — **Paper Fig. 3 (b-f) replica for PAAD org. (Tiriac) | oxaliplatin.** Each panel trains one model on the 36 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_SARC_alshihabi_DOCETAXEL.png` — **Paper Fig. 3 (b-f) replica for SARC organoids | docetaxel.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_SARC_alshihabi_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for SARC organoids | gemcitabine.** Each panel trains one model on the 14 lab models and applies it to the

**liver**

- `07_B3_models/figures/fig3_replica/Fig3_LIHC_broutier_SORAFENIB.png` — **Paper Fig. 3 (b-f) replica for LIHC org. (Broutier) | sorafenib.** Each panel trains one model on the 5 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LIHC_ji_SORAFENIB.png` — **Paper Fig. 3 (b-f) replica for LIHC org. (Ji) | sorafenib.** Each panel trains one model on the 61 lab models and applies it to the

**ovarian**

- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_CARBOPLATIN.png` — **Paper Fig. 3 (b-f) replica for OV organoids | carboplatin.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_CISPLATIN.png` — **Paper Fig. 3 (b-f) replica for OV organoids | cisplatin.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_DOCETAXEL.png` — **Paper Fig. 3 (b-f) replica for OV organoids | docetaxel.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_DOXORUBICIN.png` — **Paper Fig. 3 (b-f) replica for OV organoids | doxorubicin.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for OV organoids | gemcitabine.** Each panel trains one model on the 14 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_OV_vias_PACLITAXEL.png` — **Paper Fig. 3 (b-f) replica for OV organoids | paclitaxel.** Each panel trains one model on the 14 lab models and applies it to the

**lung**

- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_CARBOPLATIN.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | carboplatin.** Each panel trains one model on the 43 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_CISPLATIN.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | cisplatin.** Each panel trains one model on the 40 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_DOCETAXEL.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | docetaxel.** Each panel trains one model on the 46 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_ERLOTINIB.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | erlotinib.** Each panel trains one model on the 46 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_ETOPOSIDE.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | etoposide.** Each panel trains one model on the 46 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_GEMCITABINE.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | gemcitabine.** Each panel trains one model on the 46 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_PACLITAXEL.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | paclitaxel.** Each panel trains one model on the 46 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_PEMETREXED.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | pemetrexed.** Each panel trains one model on the 45 lab models and applies it to the
- `07_B3_models/figures/fig3_replica/Fig3_LUAD_cells_VINORELBINE.png` — **Paper Fig. 3 (b-f) replica for LUAD cell lines | vinorelbine.** Each panel trains one model on the 46 lab models and applies it to the

**all sets / cross-cancer**

- `07_B3_models/figures/B3_benchmark_heatmap.png` — **Benchmark heatmap: every model (rows) on every task (columns).** Numbers are averages over the 46 lab-model set x drug pairs
- `07_B3_models/figures/SuppFig7_network_centrality_profile.png` — **Network centrality of the drug targets and of the selected pathways (STRING>700, the paper's network).** Rows: colorectal
- `07_B3_models/tables/SuppFig6_7_centrality_selection.tsv` — table: pair, centrality, n_selected, rank, pathway, coef, n
- `07_B3_models/tables/SuppFig7_centrality_profile.tsv` — table: pair, centrality, n_targets, targets, median_targets, median_proximal_genes, median_all_pathway_genes
- `07_B3_models/tables/benchmark_heatmap_values.tsv` — table: Unnamed: 0, Lab models: held-out Spearman, Patients: C-index - 0.5, Patients: log HR per SD, Cell lines: cell-blind Spearman, Cell lines: dr
- `07_B3_models/tables/benchmark_long.tsv` — table: bundle, drug, pair, n_models, cv, n_proximal, model
- `07_B3_models/tables/benchmark_model_summary.tsv` — table: model, family, n_pairs, loo_rho_mean, loo_rho_lo, loo_rho_hi, C_mean
- `08_depmap/tables/B_drug_targets_used.tsv` — table: drug, target_source, targets
- `10_completeness/figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png` — Fig. 2c and 2g (and the same picture for every other set): where the biomarker pathway sits in the protein network
- `10_completeness/figures/Fig3a_comparison_feature_sets_schematic.png` — Fig. 3a (redrawn): the other ways to build a biomarker that the network method is compared with
- `10_completeness/figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png` — Supp. Figs. 8 and 9 (extended to every set): simpler network features instead of the proximal pathway
- `10_completeness/tables/Fig2c_2g_network_neighbourhoods.tsv` — table: key, label, drug, pathway, targets, n_targets_in_network, n_pathway_genes
- `10_completeness/tables/SuppFig8_target_genes_alone.tsv` — table: key, label, drug, gene, lab_rho, lab_rho_p, n
- `10_completeness/tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv` — table: key, label, drug, hops, n_genes, n, events
- `15_fixes/tables/03_foldalign_comparison.tsv` — table: unit, model, outcome, n, C_transductive, C_fold_aligned, C_diff

### Fig. 4 (isogenic cells, bootstrap)

**PDX (Gao) / breast / melanoma**

- `02_organoid_only/COAD_isella_pdx/tables/bootstrap_summary_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_permutation_null.tsv.gz` — table: p, hr, pathway

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COADREAD_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COADREAD_IRINOTECAN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COADREAD_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_IRINOTECAN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_authors_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_authors_IRINOTECAN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_COAD_authors_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_READ_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/COAD_vandewetering/tables/bootstrap_summary_READ_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_shi/tables/bootstrap_summary_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_tiriac/tables/bootstrap_summary_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p

**bladder**

- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_isogenic_T24.png` — Isogenic cisplatin-sensitive vs -resistant bladder cells (paper Fig. 4c, Supp. Table 3)
- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_paper_bootstrap_recount.png` — Paper's bootstrap (Fig. 4f) recounted from its source data
- `02_organoid_only/BLCA_lee2018/extras/tables/isogenic_T24_tests.tsv` — table: label, pathway, proximity_z, mean_sensitive, mean_resistant, t_sensitive_minus_resistant, p
- `02_organoid_only/BLCA_lee2018/extras/tables/paper_bootstrap_Fig4f_recount.tsv` — table: iterations, not_found, rank1, rank2, empirical_p_rank<=1, empirical_p_rank<=2
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway

**pancreas**

- `02_organoid_only/PAAD_shi/tables/bootstrap_summary_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_shi/tables/bootstrap_summary_PAAD_IRINOTECAN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_shi/tables/bootstrap_summary_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_tiriac/tables/bootstrap_summary_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_tiriac/tables/bootstrap_summary_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/PAAD_tiriac/tables/bootstrap_summary_PAAD_SN-38.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_permutation_null.tsv.gz` — table: p, hr, pathway

**liver**

- `02_organoid_only/LIHC_broutier/tables/bootstrap_summary_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p
- `02_organoid_only/LIHC_ji/tables/bootstrap_summary_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, bootstrap_random_pathways_p

**ovarian**

- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_permutation_null.tsv.gz` — table: p, hr, pathway

**lung**

- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_permutation_null.tsv.gz` — table: p, hr, pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_bootstrap_random_pathways.tsv.gz` — table: iteration, top_rank, feature, rank1_p, rank1_pathway
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_permutation_null.tsv.gz` — table: p, hr, pathway

**all sets / cross-cancer**

- `07_B3_models/figures/B3_random_pathway_sets.png` — **Paper network pathway (blue line) versus 100 random pathway sets of the same size (grey histogram).** For each random set the
- `07_B3_models/tables/random_sets_vs_paper_pathway.tsv` — table: pair, paper_C, paper_lr_p, paper_loo_rho, random_C_median, share_random_C_ge_paper, share_random_lr_p_le_paper
- `08_depmap/figures/S4_permutation_and_random_controls.png` — **Supplementary Figure S4. Is the BH3-only result different from chance?** Red line = observed Spearman rho of the
- `08_depmap/figures/S7_isogenic_pairs.png` — **Supplementary Figure S7. Matched 5-FU-resistant and parental colorectal lines (GEO).** For each pair we compare
- `08_depmap/tables/S4_permutation_random_controls.tsv` — table: screen, n, rho_BH3, label_permutation_p_two_sided, label_permutation_p_one_sided_negative, reactome_pathways, BH3_percentile_among_all_pathw
- `08_depmap/tables/S7_isogenic_pairs.tsv` — table: dataset, gene_set, n_resistant, n_parental, mean_resistant, mean_parental, difference_res_minus_par
- `09_registered_test/figures/fig4_random_signature_null.png` — **Figure 4. Is the frozen score better than a random score?** For each drug, 1,000 random scores were built, each
- `09_registered_test/tables/stage1_random_null_z.tsv.gz` — table: FLUOROURACIL, OXALIPLATIN, IRINOTECAN, GEMCITABINE, PACLITAXEL, PEMETREXED, SORAFENIB
- `10_completeness/figures/Fig4a_isogenic_validation_schematic.png` — Fig. 4a (redrawn): the isogenic cell-line check
- `10_completeness/figures/Fig4bc_isogenic_other_cancers.png` — Fig. 4b/c equivalents for other cancers: biomarker pathways in drug-resistant cell lines
- `10_completeness/figures/Fig4d_bootstrap_schematic.png` — Fig. 4d (redrawn): the random-pathway bootstrap
- `10_completeness/tables/Fig4bc_isogenic_other_cancers.tsv` — table: dataset, drug, biomarker, pathway, organoid_coef, expected, n_sensitive
- `14_gaps_stats/figures/G16_joint_model_bootstrap.png` — **Figure G16. Does adding organoid data to a patient model improve prediction in unseen patients?**
- `14_gaps_stats/figures/G9_permutation_null.png` — **Figure G9. How many small p-values would the 57-analysis survey give by chance, taking shared patients into account?**
- `14_gaps_stats/tables/G10_random_effects_meta.tsv` — table: effect, subset, k, HR, lo_HKSJ, hi_HKSJ, p_HKSJ
- `14_gaps_stats/tables/G16_paired_bootstrap_cindex_diff.tsv` — table: unit, outcome, n_cohorts, n_patients, cohorts, C_organoid-only, C_patient-only
- `14_gaps_stats/tables/G9_per_analysis_permutation.tsv` — table: analysis_id, cluster, treated_p_obs, untreated_p_obs, treated_perm_p_single, familywise_minP_adjusted_p, null_share_p_below_0.05
- `14_gaps_stats/tables/G9_permutation_null_summary.tsv` — table: statistic, observed, null_mean, null_median, null_2.5pct, null_97.5pct, null_95pct
- `15_fixes/backup/09/tables/stage1_random_null_z.tsv.gz` — table: FLUOROURACIL, OXALIPLATIN, IRINOTECAN, GEMCITABINE, PACLITAXEL, PEMETREXED, SORAFENIB

### Fig. 5 / Supp 12, 13, 18 (mutation concordance)

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/figures/Fig5a_cetuximab_BRAF_V600E.png` — Fig. 5a: cetuximab model vs BRAF V600E in TCGA colon tumours
- `02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv` — table: label, n_models, pathway, coef, n_mut, n_wt, p_one_sided

**bladder**

- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_ERCC2_concordance.png` — ERCC2 mutation and predicted cisplatin resistance (paper Fig. 5b and Supp. Fig. 13)
- `02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv` — table: feature_set, n_features, n_mut, n_wt, median_mut, median_wt, p_one_sided

**all sets / cross-cancer**

- `08_depmap/figures/S5_mutation_adjustment.png` — **Supplementary Figure S5. Adjusting for colorectal driver mutations.** Linear regression in colorectal lines of
- `08_depmap/tables/S5_mutation_adjustment.tsv` — table: screen, model, n, beta_BH3_per_SD, ci_low, ci_high, p
- `10_completeness/figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png` — Fig. 5 equivalents for other cancers, Supp. Fig. 12 and Supp. Fig. 18
- `10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv` — table: model, pathway, mutation, expected_in_mutants, n_mut, n_wt, median_mut
- `10_completeness/tables/SuppFig12_cetuximab_BRAF_feature_sets.tsv` — table: feature_set, n_mut, n_wt, median_mut, median_wt, p_one_sided
- `10_completeness/tables/SuppFig18_leave_k_out_cetuximab_BRAF.tsv` — table: k, rep, pathway, n_mut, n_wt, median_mut, median_wt
- `10_completeness/tables/SuppFig18_leave_k_out_summary.tsv` — table: k, share_p_below_005, median_p, share_same_rank1
- `13_gaps_bio/figures/G8_mutation_concordance_null_and_adjustment.png` — Are the mutation checks special to the organoid pathways, and do they survive adjustment?
- `13_gaps_bio/tables/G8_BRAF_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G8_BRAF_null_and_adjustment.tsv` — table: test, version, model, pathway, coef, n_proximal, n_mut
- `13_gaps_bio/tables/G8_ERCC2_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G8_ERCC2_null_and_adjustment.tsv` — table: test, version, model, pathway, coef, ridge_rank, n_proximal

### Supp. Figs 1, 14 (batch / PCA)

**PDX (Gao) / breast / melanoma**

- `02_organoid_only/COAD_isella_pdx/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling
- `03_organoid_plus_patient/figures/alignment_colorectal.png` — **Alignment of lab models with patient tumours (colorectal).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for colorec

**bladder**

- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_PCA_batch.png` — Organoids vs patients before and after scaling (paper Supp. Fig. 1c,d,g,h)
- `02_organoid_only/BLCA_lee2018/extras/tables/SuppFig1_PCA_variance.tsv` — table: level, scaled, pc1, pc2
- `03_organoid_plus_patient/figures/alignment_bladder.png` — **Alignment of lab models with patient tumours (bladder).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for bladder; c

**pancreas**

- `02_organoid_only/PAAD_shi/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling
- `02_organoid_only/PAAD_tiriac/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling
- `03_organoid_plus_patient/figures/alignment_pancreas.png` — **Alignment of lab models with patient tumours (pancreas).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for pancreas;

**liver**

- `02_organoid_only/LIHC_broutier/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling
- `02_organoid_only/LIHC_ji/figures/SuppFig1_PCA_models_vs_TCGA.png` — Supp. Fig. 1: lab models vs patients, before and after scaling
- `03_organoid_plus_patient/figures/alignment_liver.png` — **Alignment of lab models with patient tumours (liver).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for liver; colou

**ovarian**

- `03_organoid_plus_patient/figures/alignment_ovarian.png` — **Alignment of lab models with patient tumours (ovarian).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for ovarian; c

**lung**

- `02_organoid_only/LUAD_GSE165611_substitute/figures/GSE165611_PDC_vs_TCGA_PCA.png` — Lung patient-derived cells (GSE165611) next to TCGA-LUAD tumours
- `02_organoid_only/LUAD_GSE165611_substitute/tables/GSE165611_PCA_variance.tsv` — table: scaled, pc1, pc2
- `03_organoid_plus_patient/figures/alignment_lung.png` — **Alignment of lab models with patient tumours (lung).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for lung; coloure

**sarcoma**

- `03_organoid_plus_patient/figures/alignment_sarcoma.png` — **Alignment of lab models with patient tumours (sarcoma).** Each dot is one sample. Grey dots are patient tumours from all cohorts used for sarcoma; c

**all sets / cross-cancer**

- `03_organoid_plus_patient/tables/alignment_metrics.tsv` — table: cancer, stage, n_tumours, n_models, model_vs_tumour_AUROC, tumour_share_of_model_neighbours, tumour_share_expected
- `07_B3_models/figures/B3_alignment_with_without.png` — **With versus without the lab-model-to-tumour alignment (B3 item h).** Each dot is one set x drug pair. x = model trained on the
- `07_B3_models/tables/alignment_with_vs_without.tsv` — table: model, n_pairs, mean_C_without, mean_C_with, mean_delta_C, wilcoxon_p_C, mean_loo_without
- `10_completeness/figures/SuppFig14_cell_lines_vs_TCGA_PCA.png` — Supp. Fig. 14 (with the explained-variance numbers; extended to 8 cancers)
- `10_completeness/tables/SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv` — table: cancer, tcga, level, scaled, n_cell_lines, n_patients, n_features
- `10_completeness/tables/SuppFig1_PCA_variance_all_sets.tsv` — table: bundle, label, level, scaled, n_lab, n_patients, n_features
- `15_fixes/backup/03/tables/alignment_metrics.tsv` — table: cancer, stage, n_tumours, n_models, model_vs_tumour_AUROC, tumour_share_of_model_neighbours, tumour_share_expected

### Supp. Figs 2-5, 10, 11, 15-17 (models, multi-pathway, size, cell lines, proximity, inner validation)

**PDX (Gao) / breast / melanoma**

- `02_organoid_only/COAD_isella_pdx/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_CETUXIMAB_SuppFig11.png` — Supplementary figure: Colorectal PDX (Isella 2017, cetuximab) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_isella_pdx/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_CETUXIMAB.png` — Model coefficients: Colorectal PDX (Isella 2017, cetuximab) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_isella_pdx/figures/SuppFig10_size_vs_p_COAD_CETUXIMAB.png` — Supp. Fig. 10: pathway size vs survival p (cetuximab, COAD (current GDC))
- `02_organoid_only/COAD_isella_pdx/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_CETUXIMAB.png` — Kaplan-Meier grid: Colorectal PDX (Isella 2017, cetuximab) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_isella_pdx/figures/supp_proximity_heatmap.png` — Proximity heatmap: Colorectal PDX (Isella 2017, cetuximab)
- `02_organoid_only/COAD_isella_pdx/tables/SuppFig10_size_vs_p_COAD_CETUXIMAB.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_isella_pdx/tables/coefficients_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_isella_pdx/tables/loo_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_isella_pdx/tables/proximity_COAD_CETUXIMAB.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_isella_pdx/tables/proximity_all_drugs.tsv` — table: Pathway, CETUXIMAB
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_PACLITAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TAMOXIFEN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/BRCA_TRASTUZUMAB_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_CETUXIMAB_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/COAD_FLUOROURACIL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_ERLOTINIB_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/LUAD_PACLITAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/PAAD_GEMCITABINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/PDX_novartis_gao2015/tables/SKCM_DACARBAZINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COADREAD_FLUOROURACIL.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COADREAD_IRINOTECAN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COADREAD_OXALIPLATIN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_CETUXIMAB_SuppFig11.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_FLUOROURACIL.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_IRINOTECAN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_OXALIPLATIN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_IRINOTECAN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_OXALIPLATIN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_READ_FLUOROURACIL.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_READ_OXALIPLATIN.png` — Supplementary figure: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COADREAD_FLUOROURACIL.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COADREAD_IRINOTECAN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COADREAD_OXALIPLATIN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_CETUXIMAB.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_FLUOROURACIL.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_IRINOTECAN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_OXALIPLATIN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_IRINOTECAN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_OXALIPLATIN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_READ_FLUOROURACIL.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_READ_OXALIPLATIN.png` — Model coefficients: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/Fig2d_regimen_split_5FU_vs_FOLFOX.png` — 5-FU-treated colorectal patients split by regimen
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COADREAD_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, COAD+READ pooled (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COADREAD_IRINOTECAN.png` — Supp. Fig. 10: pathway size vs survival p (irinotecan, COAD+READ pooled (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COADREAD_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, COAD+READ pooled (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_CETUXIMAB.png` — Supp. Fig. 10: pathway size vs survival p (cetuximab, COAD (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, COAD (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_IRINOTECAN.png` — Supp. Fig. 10: pathway size vs survival p (irinotecan, COAD (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, COAD (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_authors_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, COAD (authors' 2019 files))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_authors_IRINOTECAN.png` — Supp. Fig. 10: pathway size vs survival p (irinotecan, COAD (authors' 2019 files))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_authors_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, COAD (authors' 2019 files))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_READ_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, READ (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_READ_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, READ (current GDC))
- `02_organoid_only/COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png` — Supp. Fig. 15 (colon): 5-FU models trained on cancer cell lines instead of organoids
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COADREAD_FLUOROURACIL.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COADREAD_IRINOTECAN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COADREAD_OXALIPLATIN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD+READ pooled (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_CETUXIMAB.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | CETUXIMAB
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_FLUOROURACIL.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_IRINOTECAN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_OXALIPLATIN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_FLUOROURACIL.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_IRINOTECAN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | IRINOTECAN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_OXALIPLATIN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | COAD (authors' 2019 files) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_READ_FLUOROURACIL.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | FLUOROURACIL
- `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_READ_OXALIPLATIN.png` — Kaplan-Meier grid: Colorectal organoids (van de Wetering 2015) | READ (current GDC) | OXALIPLATIN
- `02_organoid_only/COAD_vandewetering/figures/supp_proximity_heatmap.png` — Proximity heatmap: Colorectal organoids (van de Wetering 2015)
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COADREAD_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COADREAD_IRINOTECAN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COADREAD_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_CETUXIMAB.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_IRINOTECAN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_authors_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_authors_IRINOTECAN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_COAD_authors_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_READ_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig10_size_vs_p_READ_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/COAD_vandewetering/tables/SuppFig15_GDSC_cell_line_models_5FU.tsv` — table: training, n_lines, patients, rank, pathway, n, events
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COADREAD_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COADREAD_IRINOTECAN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COADREAD_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_IRINOTECAN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_authors_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_authors_IRINOTECAN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_COAD_authors_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_READ_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/coefficients_READ_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/COAD_vandewetering/tables/loo_COADREAD_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COADREAD_IRINOTECAN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COADREAD_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_CETUXIMAB.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_IRINOTECAN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_authors_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_authors_IRINOTECAN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_COAD_authors_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_READ_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/loo_READ_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/COAD_vandewetering/tables/proximity_COADREAD_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COADREAD_IRINOTECAN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COADREAD_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_CETUXIMAB.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_IRINOTECAN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_authors_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_authors_IRINOTECAN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_COAD_authors_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_READ_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_READ_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/COAD_vandewetering/tables/proximity_all_drugs.tsv` — table: Pathway, FLUOROURACIL, IRINOTECAN, OXALIPLATIN, CETUXIMAB
- `02_organoid_only/COAD_vandewetering/tables/regimen_split_5FU_FOLFOX.tsv` — table: cohort, subgroup, n, events, n_resp, n_nonresp, events_resp
- `02_organoid_only/PAAD_shi/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_FLUOROURACIL.png` — Supplementary figure: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_shi/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_FLUOROURACIL.png` — Model coefficients: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_shi/figures/SuppFig10_size_vs_p_PAAD_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, PAAD (current GDC))
- `02_organoid_only/PAAD_shi/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_FLUOROURACIL.png` — Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_shi/tables/SuppFig10_size_vs_p_PAAD_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_shi/tables/coefficients_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_shi/tables/loo_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_shi/tables/proximity_PAAD_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_tiriac/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_FLUOROURACIL.png` — Supplementary figure: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_tiriac/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_FLUOROURACIL.png` — Model coefficients: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_tiriac/figures/SuppFig10_size_vs_p_PAAD_FLUOROURACIL.png` — Supp. Fig. 10: pathway size vs survival p (fluorouracil, PAAD (current GDC))
- `02_organoid_only/PAAD_tiriac/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_FLUOROURACIL.png` — Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | FLUOROURACIL
- `02_organoid_only/PAAD_tiriac/tables/SuppFig10_size_vs_p_PAAD_FLUOROURACIL.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_tiriac/tables/coefficients_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_tiriac/tables/loo_PAAD_FLUOROURACIL.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_tiriac/tables/proximity_PAAD_FLUOROURACIL.tsv` — table: Pathway, proximity_z, proximal
- `03_organoid_plus_patient/tables/CRC_5FU__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/CRC_5FU__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**bladder**

- `02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/proximity_recomputation_check_cisplatin.tsv` — table: pathways, pearson_r, spearman_r, proximal_paper, proximal_ours, proximal_both
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CARBOPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_DOXORUBICIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_GEMCITABINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_METHOTREXATE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_PACLITAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/current/tables/BLCA_VINBLASTINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_GDSC_trained.png` — Models trained on cancer cell lines instead of organoids (paper Supp. Fig. 15c,d)
- `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_multipathway_vs_paper.png` — Multi-pathway models (k = 2..10) next to the paper's source data (MOESM3 'FigureS5', Supp. Fig. 5b)
- `02_organoid_only/BLCA_lee2018/extras/tables/SuppFig15_GDSC_trained.tsv` — table: training, n_lines, rank, pathway, coef, n, events
- `02_organoid_only/BLCA_lee2018/extras/tables/multipathway_k_vs_paper_sourcedata.tsv` — table: model, k, p_paper, p_ours_paper_data, p_ours_current_data
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/paper/tables/BLCA_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/all9/tables/BLCA_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/onepp/tables/BLCA_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**pancreas**

- `02_organoid_only/PAAD_shi/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_GEMCITABINE.png` — Supplementary figure: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_shi/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_IRINOTECAN.png` — Supplementary figure: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | IRINOTECAN
- `02_organoid_only/PAAD_shi/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_OXALIPLATIN.png` — Supplementary figure: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_shi/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_GEMCITABINE.png` — Model coefficients: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_shi/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_IRINOTECAN.png` — Model coefficients: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | IRINOTECAN
- `02_organoid_only/PAAD_shi/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_OXALIPLATIN.png` — Model coefficients: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_shi/figures/SuppFig10_size_vs_p_PAAD_GEMCITABINE.png` — Supp. Fig. 10: pathway size vs survival p (gemcitabine, PAAD (current GDC))
- `02_organoid_only/PAAD_shi/figures/SuppFig10_size_vs_p_PAAD_IRINOTECAN.png` — Supp. Fig. 10: pathway size vs survival p (irinotecan, PAAD (current GDC))
- `02_organoid_only/PAAD_shi/figures/SuppFig10_size_vs_p_PAAD_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, PAAD (current GDC))
- `02_organoid_only/PAAD_shi/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_GEMCITABINE.png` — Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_shi/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_IRINOTECAN.png` — Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | IRINOTECAN
- `02_organoid_only/PAAD_shi/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_OXALIPLATIN.png` — Kaplan-Meier grid: Pancreatic organoids (Shi 2022) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_shi/figures/supp_proximity_heatmap.png` — Proximity heatmap: Pancreatic organoids (Shi 2022)
- `02_organoid_only/PAAD_shi/tables/SuppFig10_size_vs_p_PAAD_GEMCITABINE.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_shi/tables/SuppFig10_size_vs_p_PAAD_IRINOTECAN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_shi/tables/SuppFig10_size_vs_p_PAAD_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_shi/tables/coefficients_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_shi/tables/coefficients_PAAD_IRINOTECAN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_shi/tables/coefficients_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_shi/tables/loo_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_shi/tables/loo_PAAD_IRINOTECAN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_shi/tables/loo_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_shi/tables/proximity_PAAD_GEMCITABINE.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_shi/tables/proximity_PAAD_IRINOTECAN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_shi/tables/proximity_PAAD_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_shi/tables/proximity_all_drugs.tsv` — table: Pathway, FLUOROURACIL, GEMCITABINE, IRINOTECAN, OXALIPLATIN
- `02_organoid_only/PAAD_tiriac/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_GEMCITABINE.png` — Supplementary figure: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_tiriac/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_OXALIPLATIN.png` — Supplementary figure: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_tiriac/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_PAAD_SN-38.png` — Supplementary figure: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | SN-38 (models) / IRINOTECAN (patients)
- `02_organoid_only/PAAD_tiriac/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_GEMCITABINE.png` — Model coefficients: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_tiriac/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_OXALIPLATIN.png` — Model coefficients: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_tiriac/figures/Fig2b_SuppFig2_SuppFig3_coefficients_PAAD_SN-38.png` — Model coefficients: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | SN-38 (models) / IRINOTECAN (patients)
- `02_organoid_only/PAAD_tiriac/figures/SuppFig10_size_vs_p_PAAD_GEMCITABINE.png` — Supp. Fig. 10: pathway size vs survival p (gemcitabine, PAAD (current GDC))
- `02_organoid_only/PAAD_tiriac/figures/SuppFig10_size_vs_p_PAAD_OXALIPLATIN.png` — Supp. Fig. 10: pathway size vs survival p (oxaliplatin, PAAD (current GDC))
- `02_organoid_only/PAAD_tiriac/figures/SuppFig10_size_vs_p_PAAD_SN-38.png` — Supp. Fig. 10: pathway size vs survival p (sn-38, PAAD (current GDC))
- `02_organoid_only/PAAD_tiriac/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_GEMCITABINE.png` — Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | GEMCITABINE
- `02_organoid_only/PAAD_tiriac/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_OXALIPLATIN.png` — Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | OXALIPLATIN
- `02_organoid_only/PAAD_tiriac/figures/SuppFig5_SourceData-multipathway_KMgrid_PAAD_SN-38.png` — Kaplan-Meier grid: Pancreatic organoids (Tiriac 2018) | PAAD (current GDC) | SN-38 (models) / IRINOTECAN (patients)
- `02_organoid_only/PAAD_tiriac/figures/supp_proximity_heatmap.png` — Proximity heatmap: Pancreatic organoids (Tiriac 2018)
- `02_organoid_only/PAAD_tiriac/tables/SuppFig10_size_vs_p_PAAD_GEMCITABINE.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_tiriac/tables/SuppFig10_size_vs_p_PAAD_OXALIPLATIN.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_tiriac/tables/SuppFig10_size_vs_p_PAAD_SN-38.tsv` — table: pathway, size, p_alone
- `02_organoid_only/PAAD_tiriac/tables/coefficients_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_tiriac/tables/coefficients_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_tiriac/tables/coefficients_PAAD_SN-38.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/PAAD_tiriac/tables/loo_PAAD_GEMCITABINE.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_tiriac/tables/loo_PAAD_OXALIPLATIN.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_tiriac/tables/loo_PAAD_SN-38.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/PAAD_tiriac/tables/proximity_PAAD_GEMCITABINE.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_tiriac/tables/proximity_PAAD_OXALIPLATIN.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_tiriac/tables/proximity_PAAD_SN-38.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/PAAD_tiriac/tables/proximity_all_drugs.tsv` — table: Pathway, FLUOROURACIL, GEMCITABINE, OXALIPLATIN, IRINOTECAN
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_DOCETAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/SARC_alshihabi2024/tables/SARC_GEMCITABINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**liver**

- `02_organoid_only/LIHC_broutier/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_LIHC_SORAFENIB.png` — Supplementary figure: Liver cancer organoids (Broutier 2017) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_broutier/figures/Fig2b_SuppFig2_SuppFig3_coefficients_LIHC_SORAFENIB.png` — Model coefficients: Liver cancer organoids (Broutier 2017) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_broutier/figures/SuppFig10_size_vs_p_LIHC_SORAFENIB.png` — Supp. Fig. 10: pathway size vs survival p (sorafenib, LIHC (current GDC))
- `02_organoid_only/LIHC_broutier/figures/SuppFig5_SourceData-multipathway_KMgrid_LIHC_SORAFENIB.png` — Kaplan-Meier grid: Liver cancer organoids (Broutier 2017) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_broutier/figures/supp_proximity_heatmap.png` — Proximity heatmap: Liver cancer organoids (Broutier 2017)
- `02_organoid_only/LIHC_broutier/tables/SuppFig10_size_vs_p_LIHC_SORAFENIB.tsv` — table: pathway, size, p_alone
- `02_organoid_only/LIHC_broutier/tables/coefficients_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/LIHC_broutier/tables/loo_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/LIHC_broutier/tables/proximity_LIHC_SORAFENIB.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/LIHC_broutier/tables/proximity_all_drugs.tsv` — table: Pathway, SORAFENIB
- `02_organoid_only/LIHC_ji/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_LIHC_SORAFENIB.png` — Supplementary figure: Liver cancer organoids (Ji 2023) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_ji/figures/Fig2b_SuppFig2_SuppFig3_coefficients_LIHC_SORAFENIB.png` — Model coefficients: Liver cancer organoids (Ji 2023) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_ji/figures/SuppFig10_size_vs_p_LIHC_SORAFENIB.png` — Supp. Fig. 10: pathway size vs survival p (sorafenib, LIHC (current GDC))
- `02_organoid_only/LIHC_ji/figures/SuppFig5_SourceData-multipathway_KMgrid_LIHC_SORAFENIB.png` — Kaplan-Meier grid: Liver cancer organoids (Ji 2023) | LIHC (current GDC) | SORAFENIB
- `02_organoid_only/LIHC_ji/figures/supp_proximity_heatmap.png` — Proximity heatmap: Liver cancer organoids (Ji 2023)
- `02_organoid_only/LIHC_ji/tables/SuppFig10_size_vs_p_LIHC_SORAFENIB.tsv` — table: pathway, size, p_alone
- `02_organoid_only/LIHC_ji/tables/coefficients_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, z
- `02_organoid_only/LIHC_ji/tables/loo_LIHC_SORAFENIB.tsv` — table: Unnamed: 0, observed, loo_Ridge, loo_SVR, loo_LinearRegression
- `02_organoid_only/LIHC_ji/tables/proximity_LIHC_SORAFENIB.tsv` — table: Pathway, proximity_z, proximal
- `02_organoid_only/LIHC_ji/tables/proximity_all_drugs.tsv` — table: Pathway, SORAFENIB
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**ovarian**

- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_CARBOPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_DOCETAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_DOXORUBICIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_GEMCITABINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/OV_vias2023/tables/OV_PACLITAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `03_organoid_plus_patient/tables/OV_PLATINUM__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/OV_PACLITAXEL__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/OV_PLATINUM__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**lung**

- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CARBOPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_CISPLATIN_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_DOCETAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ERLOTINIB_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_ETOPOSIDE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_GEMCITABINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PACLITAXEL_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_PEMETREXED_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_coefficients.tsv` — table: pathway, coef_Ridge, coef_SVR, coef_LinearRegression, proximity_z
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_loo_predictions.tsv` — table: Unnamed: 0, Ridge, SVR, LinearRegression, observed
- `02_organoid_only/LUAD_GSE165611_substitute/tables/LUAD_VINORELBINE_multipathway_k.tsv` — table: n, events, n_nonresp, n_resp, events_nonresp, events_resp, p
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/LUNG_PLATINUM__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**sarcoma**

- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__learning_curve.tsv` — table: unit, axis, fraction, n, rep, model, concordance

**all sets / cross-cancer**

- `03_organoid_plus_patient/tables/learning_curves.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `03_organoid_plus_patient/tables/learning_curves_mean.tsv` — table: unit, axis, model, 0.1, 0.25, 0.5, 0.75
- `07_B3_models/figures/B3_learning_curves.png` — **Learning curves: how performance changes with the amount of training data.** (a, b) The larger lab-model sets
- `07_B3_models/tables/learning_curve_cell_lines.tsv` — table: frac, rep, n_cells, model, mean_rho
- `07_B3_models/tables/learning_curve_lab_models.tsv` — table: bundle, drug, fold, frac, rep, n_train, model
- `08_depmap/tables/B_network_proximity_z.tsv` — table: Unnamed: 0, FLUOROURACIL, OXALIPLATIN, IRINOTECAN, CISPLATIN, CARBOPLATIN, GEMCITABINE
- `10_completeness/figures/Fig1b_network_proximity_schematic.png` — Fig. 1b (redrawn): how "network proximity" is measured
- `10_completeness/figures/SuppFig10_size_vs_p_all_sets.png` — Supp. Fig. 10 (extended to every set): does pathway size drive the survival split?
- `10_completeness/figures/SuppFig16_proximity_all_drugs.png` — Supp. Fig. 16 (paper's 44 drug columns + drugs added for the other sets)
- `10_completeness/figures/SuppFig17_split_603010_all_sets.png` — Supp. Fig. 17 (paper's 60/10/30 split, every headline set)
- `10_completeness/figures/SuppFig4_cross_test_matrix_all_sets.png` — Supp. Fig. 4 (extended to every set): is each biomarker specific to its own cancer and drug?
- `10_completeness/tables/SuppFig10_size_vs_p_summary_all_sets.tsv` — table: key, label, n_pathways, R2, slope_p, share_p_below_005, rank1_pathway_p_any_direction
- `10_completeness/tables/SuppFig16_paper_proximity_44_columns.tsv.gz` — table: Pathway, AICA RIBONUCLEOTIDE, AKT INHIBITOR VIII, AXITINIB, BI-2536, BMS-536924, BORTEZOMIB
- `10_completeness/tables/SuppFig16_proximity_all_drugs.tsv.gz` — table: Unnamed: 0, AICA RIBONUCLEOTIDE, AKT INHIBITOR VIII, AXITINIB, BI-2536, BMS-536924, BORTEZOMIB
- `10_completeness/tables/SuppFig17_split_603010_summary.tsv` — table: key, label, n_models, n_test, median_spearman, share_positive
- `10_completeness/tables/SuppFig4_cross_test_matrix.tsv` — table: biomarker, cohort, pathway, own, n, events, p
- `10_completeness/tables/SuppFig4_cross_test_summary.json` — JSON summary: n_own, own_p_below_005, n_off, off_p_below_005, off_share, paper_SuppFig4a_COAD_with_BLCA_amino
- `13_gaps_bio/figures/G14_SuppFig15_cell_line_models.png` — Supp. Fig. 15 for the other cancers: models trained on cancer cell lines instead of organoids
- `13_gaps_bio/tables/G14_SuppFig15_cell_line_models.tsv` — table: setting, project, patient_drugs, line_drug, proximity_column, training, n_lines
- `15_fixes/backup/03/tables/learning_curves.tsv` — table: unit, axis, fraction, n, rep, model, concordance
- `15_fixes/backup/03/tables/learning_curves_mean.tsv` — table: unit, axis, model, 0.1, 0.25, 0.5, 0.75
- `15_fixes/backup/10/SuppFig4_cross_test_matrix_all_sets.png` — figure (no caption file)

### Supp. Tables 1-4 (clinical balance, isogenic, bootstrap)

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/figures/Fig4b_SuppTable3_isogenic_HCT8_5FU.png` — Fig. 4b (and Supp. Table 3): isogenic 5-FU-sensitive vs -resistant colorectal cells
- `02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv` — table: cohort, nes_version, pathway, ridge_rank, ridge_coef, proximity_z, mean_sensitive

**bladder**

- `02_organoid_only/BLCA_lee2018/extras/tables/SuppTable2_clinical_balance.tsv` — table: version, biomarker, feature, level, responders, non_responders, chi2_p

**all sets / cross-cancer**

- `02_organoid_only/tables/SuppTable1_clinical_features_by_predicted_group.tsv` — table: folder, cohort, drug, feature, level, responder, non_responder
- `07_B3_models/tables/SuppTable4_DL_grid_all_pairs.tsv` — table: pair, batch_size, nodes, lr, dropout, weight_decay, gamma
- `07_B3_models/tables/SuppTable4_DL_selected.tsv` — table: pair, batch_size, nodes, lr, dropout, weight_decay, gamma
- `10_completeness/tables/SuppTable1_2_clinical_balance_all_sets.tsv` — table: key, label, feature, level, responder, non_responder, chi2_p

### DepMap / registered test / other

**PDX (Gao) / breast / melanoma**

- `02_organoid_only/COAD_isella_pdx/tables/all_tests_COAD_CETUXIMAB.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp

**colon / rectum**

- `02_organoid_only/COAD_vandewetering/tables/all_tests_COADREAD_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COADREAD_IRINOTECAN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COADREAD_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_CETUXIMAB.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_IRINOTECAN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_authors_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_authors_IRINOTECAN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_authors_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_READ_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/COAD_vandewetering/tables/all_tests_READ_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_shi/tables/all_tests_PAAD_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_tiriac/tables/all_tests_PAAD_FLUOROURACIL.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `03_organoid_plus_patient/figures/supp/unit_CRC_5FU.png` — **Supplementary figure for CRC_5FU (colorectal, 5-FU-based chemotherapy).**
- `03_organoid_plus_patient/figures/supp/unit_CRC_IRINOTECAN.png` — **Supplementary figure for CRC_IRINOTECAN (colorectal, irinotecan (FOLFIRI)).**
- `03_organoid_plus_patient/figures/supp/unit_CRC_OXALIPLATIN.png` — **Supplementary figure for CRC_OXALIPLATIN (colorectal, oxaliplatin (FOLFOX)).**
- `03_organoid_plus_patient/tables/CRC_5FU__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_5FU__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_5FU__all__meta.tsv` — table: unit, CRC_5FU
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist__meta.tsv` — table: unit, CRC_5FU
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_5FU__shortlist_foldalign__meta.tsv` — table: unit, CRC_5FU
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__all__meta.tsv` — table: unit, CRC_IRINOTECAN
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_IRINOTECAN__shortlist__meta.tsv` — table: unit, CRC_IRINOTECAN
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__all__meta.tsv` — table: unit, CRC_OXALIPLATIN
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/CRC_OXALIPLATIN__shortlist__meta.tsv` — table: unit, CRC_OXALIPLATIN
- `08_depmap/figures/S3_colon_vs_pancancer.png` — **Supplementary Figure S3. Colorectal lines compared with other cancer types.** Spearman rho between the BH3-only
- `08_depmap/tables/A_b_bh3_vs_5fu_colon.tsv` — table: screen, metric, lines, n, spearman_rho, ci_low, ci_high
- `08_depmap/tables/A_c_frozen_model_colon_lines.tsv` — table: screen, signature, k_pathways, n, spearman_rho, p, null_q025
- `08_depmap/tables/A_h_crispr_vs_5fu.tsv` — table: screen, gene, colon_rho, colon_p, colon_n, pancancer_partial_rho, pancancer_p
- `08_depmap/tables/S3_colon_vs_pancancer.tsv` — table: screen, group, n, rho, p
- `13_gaps_bio/figures/G2_decomposition_colon_5FU.png` — What changed between the 2019 colon 5-FU result and today's data
- `13_gaps_bio/tables/G2_decomposition_colon_5FU.tsv` — table: frozen_sign, n_proximal, bh3_proximal, rank1, n, deaths, p_rank1
- `15_fixes/backup/03/tables/CRC_5FU__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_5FU__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_5FU__all__meta.tsv` — table: unit, CRC_5FU
- `15_fixes/backup/03/tables/CRC_5FU__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_5FU__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_5FU__shortlist__meta.tsv` — table: unit, CRC_5FU
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__all__meta.tsv` — table: unit, CRC_IRINOTECAN
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_IRINOTECAN__shortlist__meta.tsv` — table: unit, CRC_IRINOTECAN
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__all__meta.tsv` — table: unit, CRC_OXALIPLATIN
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/CRC_OXALIPLATIN__shortlist__meta.tsv` — table: unit, CRC_OXALIPLATIN

**bladder**

- `02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/data_versions.tsv` — table: file, bytes, sha256
- `02_organoid_only/BLCA_lee2018/extras/tables/extras_log.json` — JSON summary: named_biomarker, paper_bootstrap, isogenic, ercc2, comparisons, clinical
- `02_organoid_only/BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv` — table: analysis, lines, n_models, n_proximal, biomarker, treated_n, treated_events
- `03_organoid_plus_patient/figures/supp/unit_BLCA_CISPLATIN.png` — **Supplementary figure for BLCA_CISPLATIN (bladder, cisplatin-based chemotherapy).**
- `03_organoid_plus_patient/figures/supp/unit_BLCA_GEMCITABINE.png` — **Supplementary figure for BLCA_GEMCITABINE (bladder, gemcitabine).**
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__all__meta.tsv` — table: unit, BLCA_CISPLATIN
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist__meta.tsv` — table: unit, BLCA_CISPLATIN
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist_foldalign__meta.tsv` — table: unit, BLCA_CISPLATIN
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__all__meta.tsv` — table: unit, BLCA_GEMCITABINE
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/BLCA_GEMCITABINE__shortlist__meta.tsv` — table: unit, BLCA_GEMCITABINE
- `08_depmap/figures/fig3_drug_specificity_and_bladder.png` — **Figure 3. Is the BH3-only signal specific to 5-FU, and does the bladder organoid model carry over?**
- `08_depmap/tables/A_f_bladder_cisplatin.tsv` — table: screen, drug, organoids_trained, proximal_pathways, signature_pathways, n_bladder_lines, spearman_rho
- `08_depmap/tables/A_f_bladder_organoid_model_table.tsv` — table: Unnamed: 0, selection_freq, sign_agreement, mean_coef
- `08_depmap/tables/A_f_bladder_organoid_signature.tsv` — table: Unnamed: 0, weight
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__all__meta.tsv` — table: unit, BLCA_CISPLATIN
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/BLCA_CISPLATIN__shortlist__meta.tsv` — table: unit, BLCA_CISPLATIN
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__all__meta.tsv` — table: unit, BLCA_GEMCITABINE
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/BLCA_GEMCITABINE__shortlist__meta.tsv` — table: unit, BLCA_GEMCITABINE

**pancreas**

- `02_organoid_only/PAAD_shi/tables/all_tests_PAAD_GEMCITABINE.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_shi/tables/all_tests_PAAD_IRINOTECAN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_shi/tables/all_tests_PAAD_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_tiriac/tables/all_tests_PAAD_GEMCITABINE.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_tiriac/tables/all_tests_PAAD_OXALIPLATIN.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/PAAD_tiriac/tables/all_tests_PAAD_SN-38.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `03_organoid_plus_patient/figures/supp/unit_PAAD_GEMCITABINE.png` — **Supplementary figure for PAAD_GEMCITABINE (pancreas, gemcitabine).**
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__all__meta.tsv` — table: unit, PAAD_GEMCITABINE
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist__meta.tsv` — table: unit, PAAD_GEMCITABINE
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/PAAD_GEMCITABINE__shortlist_foldalign__meta.tsv` — table: unit, PAAD_GEMCITABINE
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__all__meta.tsv` — table: unit, PAAD_GEMCITABINE
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/PAAD_GEMCITABINE__shortlist__meta.tsv` — table: unit, PAAD_GEMCITABINE

**liver**

- `02_organoid_only/LIHC_broutier/tables/all_tests_LIHC_SORAFENIB.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `02_organoid_only/LIHC_ji/tables/all_tests_LIHC_SORAFENIB.tsv` — table: model, k, pathways, coefs, treated_n, treated_events, treated_n_resp
- `03_organoid_plus_patient/figures/supp/unit_LIHC_SORAFENIB.png` — **Supplementary figure for LIHC_SORAFENIB (liver, sorafenib).**
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__all__meta.tsv` — table: unit, LIHC_SORAFENIB
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist__meta.tsv` — table: unit, LIHC_SORAFENIB
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LIHC_SORAFENIB__shortlist_foldalign__meta.tsv` — table: unit, LIHC_SORAFENIB
- `05_A3_B1/figures/B1_Fig6_STORM_sorafenib.png` — **Liver: STORM trial (GSE109211, sorafenib vs placebo, randomised).** Share of responders (as defined in the
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__all__meta.tsv` — table: unit, LIHC_SORAFENIB
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/LIHC_SORAFENIB__shortlist__meta.tsv` — table: unit, LIHC_SORAFENIB

**ovarian**

- `02_organoid_only/OV_GSE229119_olaparib/figures/OV_GSE229119_olaparib_timecourse.png` — Olaparib time course in ovarian organoids (GSE229119) and what it says about TCGA-OV
- `02_organoid_only/OV_GSE229119_olaparib/tables/day7_change_similarity_spearman.tsv` — table: Unnamed: 0, Hey|D7, MCF7|D7, ORG2414|D7, ORG2445|D7, ORG3008|D7, SKOV3|D7
- `02_organoid_only/OV_GSE229119_olaparib/tables/olaparib_signature_pathways.json` — JSON summary: up, down
- `02_organoid_only/OV_GSE229119_olaparib/tables/organoid_olaparib_pathway_changes.tsv` — table: pathway, mean_change_d7, same_sign_all3, t_p, proximity_z, mean_change_d3, mean_change_recovery
- `02_organoid_only/OV_GSE229119_olaparib/tables/pathway_change_vs_DMSO.tsv` — table: Term, Hey|D3, Hey|D7, Hey|R, MCF7|D3, MCF7|D7, MCF7|R
- `02_organoid_only/OV_GSE229119_olaparib/tables/tcga_ov_association.json` — JSON summary: n_organoids, n_models, proximal_pathways_olaparib, signature_up, signature_down, noise_sd_d7
- `03_organoid_plus_patient/figures/supp/unit_OV_PACLITAXEL.png` — **Supplementary figure for OV_PACLITAXEL (ovarian, paclitaxel (with platinum)).**
- `03_organoid_plus_patient/figures/supp/unit_OV_PLATINUM.png` — **Supplementary figure for OV_PLATINUM (ovarian, carboplatin / platinum (cell-line proxy)).**
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__all__meta.tsv` — table: unit, OV_PACLITAXEL
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist__meta.tsv` — table: unit, OV_PACLITAXEL
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/OV_PACLITAXEL__shortlist_foldalign__meta.tsv` — table: unit, OV_PACLITAXEL
- `03_organoid_plus_patient/tables/OV_PLATINUM__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/OV_PLATINUM__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/OV_PLATINUM__all__meta.tsv` — table: unit, OV_PLATINUM
- `03_organoid_plus_patient/tables/OV_PLATINUM__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/OV_PLATINUM__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/OV_PLATINUM__shortlist__meta.tsv` — table: unit, OV_PLATINUM
- `15_fixes/backup/03/tables/OV_PACLITAXEL__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/OV_PACLITAXEL__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/OV_PACLITAXEL__all__meta.tsv` — table: unit, OV_PACLITAXEL
- `15_fixes/backup/03/tables/OV_PACLITAXEL__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/OV_PACLITAXEL__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/OV_PACLITAXEL__shortlist__meta.tsv` — table: unit, OV_PACLITAXEL
- `15_fixes/backup/03/tables/OV_PLATINUM__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/OV_PLATINUM__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/OV_PLATINUM__all__meta.tsv` — table: unit, OV_PLATINUM
- `15_fixes/backup/03/tables/OV_PLATINUM__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/OV_PLATINUM__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/OV_PLATINUM__shortlist__meta.tsv` — table: unit, OV_PLATINUM

**lung**

- `02_organoid_only/LUAD_GSE165611_substitute/tables/GSE165611_sample_annotation.tsv` — table: Unnamed: 0, tissue, histopathology
- `03_organoid_plus_patient/figures/supp/unit_LUNG_PLATINUM.png` — **Supplementary figure for LUNG_PLATINUM (lung, platinum (cell-line proxy)).**
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__all__meta.tsv` — table: unit, LUNG_PLATINUM
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist__meta.tsv` — table: unit, LUNG_PLATINUM
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/LUNG_PLATINUM__shortlist_foldalign__meta.tsv` — table: unit, LUNG_PLATINUM
- `15_fixes/backup/03/tables/LUNG_PLATINUM__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/LUNG_PLATINUM__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/LUNG_PLATINUM__all__meta.tsv` — table: unit, LUNG_PLATINUM
- `15_fixes/backup/03/tables/LUNG_PLATINUM__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/LUNG_PLATINUM__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/LUNG_PLATINUM__shortlist__meta.tsv` — table: unit, LUNG_PLATINUM

**sarcoma**

- `03_organoid_plus_patient/figures/supp/unit_SARC_GEMCITABINE.png` — **Supplementary figure for SARC_GEMCITABINE (sarcoma, gemcitabine (+ docetaxel)).**
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__all__meta.tsv` — table: unit, SARC_GEMCITABINE
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist__meta.tsv` — table: unit, SARC_GEMCITABINE
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist_foldalign__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist_foldalign__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `03_organoid_plus_patient/tables/SARC_GEMCITABINE__shortlist_foldalign__meta.tsv` — table: unit, SARC_GEMCITABINE
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__all__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__all__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__all__meta.tsv` — table: unit, SARC_GEMCITABINE
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__shortlist__eval.tsv` — table: cohort, model, outcome, n, events, C, num
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__shortlist__folds.tsv` — table: heldout, model, train_cohorts, w, lam, unit, features
- `15_fixes/backup/03/tables/SARC_GEMCITABINE__shortlist__meta.tsv` — table: unit, SARC_GEMCITABINE

**all sets / cross-cancer**

- `00_checklist/checklist.tsv` — table: ID, type, cancer_drug, what_it_shows, reported_numbers, data_needed, how_to_generalise
- `00_checklist/tables/source_data_facts.tsv` — table: fact, value
- `03_organoid_plus_patient/tables/data_versions.tsv` — table: file, bytes, sha256
- `03_organoid_plus_patient/tables/joint_weight_chosen.tsv` — table: unit, w
- `03_organoid_plus_patient/tables/storm_auroc_by_arm.tsv` — table: model, arm, n, responders, AUROC, p
- `05_A3_B1/figures/A3_Fig2_QQ_expected_by_chance.png` — **How many small p-values would chance alone give?** (a) QQ plot: each dot is one of the 57 analyses, sorted by
- `05_A3_B1/figures/A3_Fig3_adjusted_vs_unadjusted.png` — **Does adjusting for clinical factors change the hazard ratios?** Each dot is one analysis (dot size = number of
- `05_A3_B1/figures/A3_Fig5_power_vs_events.png` — **How much could each analysis detect?** (a) Power of the median-split log-rank test (Schoenfeld formula, equal
- `05_A3_B1/figures/B1_Fig7_heatmap_signature_by_cohort.png` — **Summary heatmap: every frozen signature in every matching cohort.** Rows = signatures (organoid set, drug,
- `05_A3_B1/tables/A3_adjusted_cox_all_terms.tsv` — table: analysis_id, term, coef, se, HR, lo, hi
- `05_A3_B1/tables/A3_adjusted_cox_epv.tsv` — table: analysis_id, deaths, n_terms_full, deaths_per_term_full, low_epv_full, reduced_covariates, n_terms_reduced
- `05_A3_B1/tables/A3_expected_by_chance.tsv` — table: family, n_tests, observed_p_below_0.05, expected_by_chance, binomial_p_at_least_observed, observed_p_below_0.01, expected_0.01
- `05_A3_B1/tables/A3_master_results.tsv` — table: analysis_id, set_label, project, drug, model_drug, n_models, pathway
- `05_A3_B1/tables/A3_winners_curse.tsv` — table: analysis_id, treated_n, events, treated_p, observed_HR, observed_HR_lo, observed_HR_hi
- `05_A3_B1/tables/B1_external_results.tsv` — table: n, events, effect, logrank_p_median, est, lo, hi
- `05_A3_B1/tables/B1_frozen_signatures.tsv` — table: group, analysis_id, kind, sig_id, pathways, coefs
- `05_A3_B1/tables/data_versions.tsv` — table: file, bytes, sha256
- `07_B3_models/figures/B3_supp_heatmap_Cindex.png` — **C-index in treated patients for every model (rows) and every lab-model set x drug pair (columns).** Printed numbers are the raw values;
- `07_B3_models/figures/B3_supp_heatmap_LOO.png` — **Held-out spearman in lab models for every model (rows) and every lab-model set x drug pair (columns).** Printed numbers are the raw values;
- `07_B3_models/figures/B3_task_bars.png` — **Each task separately, with 95% confidence intervals.** (a) Held-out Spearman correlation in lab models, mean over the
- `07_B3_models/tables/bolis_selected_genes.tsv` — table: pair, n_genes, genes
- `07_B3_models/tables/cellline_cv_per_unit.tsv` — table: task, fold, model, screen, drug, n, rho
- `07_B3_models/tables/cellline_cv_per_unit_gcn.tsv` — table: task, fold, model, screen, drug, n, rho
- `07_B3_models/tables/data_versions.tsv` — table: file, bytes, sha256
- `08_depmap/figures/S2_gdsc1_vs_gdsc2.png` — **Supplementary Figure S2. Do the screens agree on 5-FU?** (a) 5-FU ln IC50 of the same cell lines in GDSC1 and
- `08_depmap/figures/S6_single_bh3_genes.png` — **Supplementary Figure S6. Each gene on its own.** Spearman rho between the expression (log2 TPM+1) of each gene and
- `08_depmap/figures/S8_metric_and_geneset_robustness.png` — **Supplementary Figure S8. Does the answer depend on the response measure or the gene-set version?** Spearman rho
- `08_depmap/figures/fig2_frozen_model_and_sign_agreement.png` — **Figure 2. The organoid model applied unchanged to colorectal cell lines.**
- `08_depmap/figures/fig4_bh3_mimetics_and_crispr.png` — **Figure 4. Mechanism checks: BH3 mimetics and CRISPR knockouts.**
- `08_depmap/tables/A_c_predictions.tsv` — table: Unnamed: 0, Unnamed: 1, predicted_resistance, observed
- `08_depmap/tables/A_d_sign_agreement_37.tsv` — table: screen, pathways_tested, sign_agree, fraction, binomial_p_one_sided
- `08_depmap/tables/A_d_sign_agreement_per_pathway.tsv` — table: screen, pathway, organoid_coef, cell_line_rho, p, n, sign_agrees
- `08_depmap/tables/A_e_bh3_vs_all_drugs.tsv` — table: screen, drug, colon_rho, colon_p, colon_n, pancancer_partial_rho, pancancer_p
- `08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv` — table: screen, drug, scope, n, rho, ci_low, ci_high
- `08_depmap/tables/B_candidates_two_or_more_sources.tsv` — table: drug, pathway, expr_meta_z, expr_p, expr_q, expr_sources, screens_used
- `08_depmap/tables/B_source1_expression_association.tsv.gz` — table: drug, pathway, expr_meta_z, expr_p, expr_q, expr_sources, screens_used
- `08_depmap/tables/B_source2_crispr_codependency.tsv.gz` — table: pathway, crispr_z, mean_gene_z, genes_in_screen, crispr_p, crispr_q, drug
- `08_depmap/tables/S2_gdsc1_vs_gdsc2_agreement.tsv` — table: pair, n_all, rho_all, p_all, n_colorectal, rho_colorectal, p_colorectal
- `08_depmap/tables/S6_single_genes.tsv` — table: screen, gene, core_BH3_only, n, rho, p
- `08_depmap/tables/S8_metric_and_geneset_robustness.tsv` — table: screen, metric, gene_set, colon_n, colon_rho, colon_p, pancancer_n
- `08_depmap/tables/data_versions.tsv` — table: file, bytes, sha256
- `08_depmap/tables/frozen_hypotheses.tsv` — table: frozen_at_utc, drug, rank, pathway, sign, direction, n_sources
- `09_registered_test/figures/fig5_negative_controls.png` — **Figure 5. Negative controls and drug specificity (TCGA, PFI).** **a**, the frozen score (blue) next to two control
- `09_registered_test/tables/data_versions.tsv` — table: file, bytes, sha256
- `10_completeness/completeness.tsv` — table: ID, type, cancer_drug, paper_before, paper_after, paper_paths, paper_files_found
- `10_completeness/tables/completeness_counts.json` — JSON summary: paper_before, paper_after, generalisation_before, generalisation_after
- `10_completeness/tables/data_versions.tsv` — table: file, bytes, sha256
- `10_completeness/tables/results_index.tsv` — table: section, cancer, path, description
- `11_review/GAP_LIST.tsv` — table: rank, gap_id, priority, what, why, data_available, data_path
- `12_audit/issues.tsv` — table: id, severity, area, folder, issue, evidence, exact_fix
- `12_audit/tables/MASTER_NUMBERS.tsv` — table: id, description, authoritative_value, source_table, expected_text, reports_quoting, note
- `12_audit/tables/STATISTICAL_HIERARCHY.tsv` — table: tier, analysis_family, folder, n_tests, primary_statistic, multiplicity_correction, what_it_can_claim
- `12_audit/tables/data_versions_02GI_sha256.tsv` — table: file, bytes, sha256
- `13_gaps_bio/figures/G13_odd_pathway_composition.png` — What is inside the "HIV" pathways among the frozen hypotheses?
- `13_gaps_bio/figures/G3_bh3_vs_bh3_mimetics_organoids.png` — BH3-only pathway score and drug response in the 19 colorectal organoids (van de Wetering 2015)
- `13_gaps_bio/figures/G5_folfox_question.png` — Is the paper's 5-FU marker really a 5-FU marker, or an oxaliplatin / FOLFOX marker?
- `13_gaps_bio/figures/G7_bh3_beyond_simple_markers.png` — Does the BH3-only pathway tell us more than BCL2-family genes and proliferation?
- `13_gaps_bio/tables/G13_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G13_frozen_pathway_composition.tsv` — table: pathway, n_genes, odd, n_proteasome (PSM*), share_proteasome (PSM*), n_nuclear pore (NUP*, RAN*, KPNB1, XPO1, SEH1L, SEC13, TPR, AAAS, RAE1)
- `13_gaps_bio/tables/G13_most_overlapping_pathways.tsv` — table: pathway, other, jaccard, overlap, n_other
- `13_gaps_bio/tables/G13_purity_immune_correlation.tsv` — table: pathway, odd, vs, n_cancer_types, median_rho, q10, q90
- `13_gaps_bio/tables/G14_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G2_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G2_gene_universe_sizes.tsv` — table: Unnamed: 0, 0
- `13_gaps_bio/tables/G3_bh3_vs_all_drugs_organoids_AUC.tsv` — table: drug, n, rho, p, rank_most_negative
- `13_gaps_bio/tables/G3_bh3_vs_drugs_organoids.tsv` — table: score, metric, drug, n, rho, lo, hi
- `13_gaps_bio/tables/G3_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G3_response_range_organoids.tsv` — table: drug, n, max_conc_uM, auc_min, auc_median, auc_max, auc_sd
- `13_gaps_bio/tables/G3_single_genes_vs_bh3_mimetics_organoids.tsv` — table: gene, drug, n, rho, lo, hi, p
- `13_gaps_bio/tables/G5_cell_lines_bh3_vs_oxaliplatin.tsv` — table: screen, drug, scope, n, rho, lo, hi
- `13_gaps_bio/tables/G5_data_versions.tsv` — table: file, bytes, sha256
- `13_gaps_bio/tables/G5_folfox_patients.tsv` — table: frozen_sign, subset, n, events, cohort, setup, endpoint
- `13_gaps_bio/tables/G7_partial_correlations.tsv` — table: screen, drug, n, adjusted_for, var, r, lo
- `13_gaps_bio/tables/G7_regression_added_value.tsv` — table: screen, drug, n, R2_markers, R2_markers_plus_BH3, dR2, F
- `13_gaps_bio/tables/check_ssgsea_fast_vs_gseapy.tsv` — table: set, r, fast0, gseapy0
- `14_gaps_stats/figures/G12_baseline_treated_vs_untreated.png` — **Figure G12. Treated patients and patients with no drug record are different groups.**
- `14_gaps_stats/figures/G17_PFI_vs_OS.png` — **Figure G17. The 57-analysis survey with progression-free interval (PFI) instead of overall survival (OS).**
- `14_gaps_stats/figures/G18_score_agreement.png` — **Figure G18. Do the two pathway scores agree in TCGA?**
- `14_gaps_stats/tables/G10_equivalence_counts.tsv` — table: bound, n_analyses, n_upper_bound_below, n_two_sided_equivalent
- `14_gaps_stats/tables/G10_per_analysis_equivalence.tsv` — table: analysis_id, drug, drug_class, project, treated_n, treated_events, cont_HR
- `14_gaps_stats/tables/G12_baseline_readable.tsv` — table: cohort, treated with, n treated / untreated, age treated, age untreated, p age, % stage III-IV treated / untreated
- `14_gaps_stats/tables/G12_baseline_treated_vs_untreated.tsv` — table: cohort, treated_with, treated_n, treated_age_median, treated_age_q1, treated_age_q3, treated_pct_male
- `14_gaps_stats/tables/G15_literature_related_to_Kong2020.tsv` — table: pmid, doi, year, journal, title, found_by, relation_to_Kong2020
- `14_gaps_stats/tables/G16_check_pooled_C_vs_original.tsv` — table: unit, outcome, model, C_original_table, C_recomputed, abs_diff
- `14_gaps_stats/tables/G16_extended_grid_chosen_weight.tsv` — table: unit, n_fits, chosen_w_0.0, chosen_w_0.05, chosen_w_0.1, chosen_w_0.25, chosen_w_0.5
- `14_gaps_stats/tables/G17_PFI_events_power.tsv` — table: median_events_OS, median_events_PFI, median_ratio_PFI_over_OS, median_MDHR_OS, median_MDHR_PFI, n_power80_HR1.5_PFI, n_same_side_p0.05_OS_an
- `14_gaps_stats/tables/G17_PFI_survey_57.tsv` — table: analysis_id, drug, project, os_treated_n, os_treated_events, os_treated_p, os_cont_HR
- `14_gaps_stats/tables/G17_PFI_vs_OS_counts.tsv` — table: family, n_tests, p_below_0.05, expected, binomial_p, BH_q_below_0.05, min_p
- `14_gaps_stats/tables/G18_hit_counts_ssgsea_vs_rank.tsv` — table: family, ssgsea, rank_score, n_rank, both
- `14_gaps_stats/tables/G18_pathway_correlation_ssgsea_vs_rank.tsv` — table: project, pathway, spearman, n, n_sources
- `14_gaps_stats/tables/G18_rerun_57_with_rank_score.tsv` — table: analysis_id, pathway, ssgsea_treated_p, direction_ssgsea, n_with_rank_score, n_total, spearman_ssgsea_vs_rank_patients
- `14_gaps_stats/tables/G9_null_count_distribution.tsv` — table: k, perm_null_prob, binomial_57_0.05_prob, perm_null_untreated_prob
- `14_gaps_stats/tables/data_versions.tsv` — table: file, bytes, sha256
- `15_fixes/backup/03/tables/data_versions.tsv` — table: file, bytes, sha256
- `15_fixes/backup/03/tables/joint_weight_chosen.tsv` — table: unit, w
- `15_fixes/backup/03/tables/storm_auroc_by_arm.tsv` — table: model, arm, n, responders, AUROC, p
- `15_fixes/backup/05/b1_patient_scores.tsv.gz` — table: cohort, comparison, endpoint, sig_id, sample, arm1, score
- `15_fixes/backup/05/tables/A3_adjusted_cox_all_terms.tsv` — table: analysis_id, term, coef, se, HR, lo, hi
- `15_fixes/backup/05/tables/A3_expected_by_chance.tsv` — table: family, n_tests, observed_p_below_0.05, expected_by_chance, binomial_p_at_least_observed, observed_p_below_0.01, expected_0.01
- `15_fixes/backup/05/tables/A3_master_results.tsv` — table: analysis_id, set_label, project, drug, model_drug, n_models, pathway
- `15_fixes/backup/05/tables/A3_winners_curse.tsv` — table: analysis_id, treated_n, events, treated_p, observed_HR, observed_HR_lo, observed_HR_hi
- `15_fixes/backup/05/tables/B1_external_results.tsv` — table: n, events, effect, logrank_p_median, est, lo, hi
- `15_fixes/backup/05/tables/B1_frozen_signatures.tsv` — table: group, analysis_id, kind, sig_id, pathways, coefs
- `15_fixes/backup/05/tables/data_versions.tsv` — table: file, bytes, sha256
- `15_fixes/backup/09/tables/data_versions.tsv` — table: file, bytes, sha256
- `15_fixes/backup/MASTER_NUMBERS_before.tsv` — table: id, description, authoritative_value, source_table, expected_text, reports_quoting, note
- `15_fixes/definitions/GSE39582_treatment.tsv` — table: subset, rule_plain_words, drugs_regex, arms_included, n_all_stages, n_stage_II_III, n_stage_II_III_with_RFS
- `15_fixes/tables/MASTER_NUMBERS.tsv` — table: id, description, authoritative_value, source_table, expected_text, reports_quoting, note

## Part B. By cancer (counts; full list in Part A)

| Cancer | Figures | Tables | Folders |
|---|---|---|---|
| PDX (Gao) / breast / melanoma | 34 | 81 | 02_organoid_only/COAD_isella_pdx, 02_organoid_only/PDX_novartis_gao2015, 07_B3_models/figures |
| colon / rectum | 85 | 167 | 02_organoid_only/COAD_vandewetering, 02_organoid_only/PAAD_shi, 02_organoid_only/PAAD_tiriac, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 05_A3_B1/figures, 07_B3_models/figures, 08_depmap/figures, 08_depmap/tables, 13_gaps_bio/figures, 13_gaps_bio/tables, 15_fixes/backup |
| bladder | 43 | 138 | 02_organoid_only/BLCA_OV_LUAD_OTHER, 02_organoid_only/BLCA_lee2018, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 05_A3_B1/figures, 07_B3_models/figures, 08_depmap/figures, 08_depmap/tables, 15_fixes/backup |
| pancreas | 44 | 83 | 02_organoid_only/PAAD_shi, 02_organoid_only/PAAD_tiriac, 02_organoid_only/SARC_alshihabi2024, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 07_B3_models/figures, 15_fixes/backup |
| liver | 21 | 38 | 02_organoid_only/LIHC_broutier, 02_organoid_only/LIHC_ji, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 05_A3_B1/figures, 07_B3_models/figures, 15_fixes/backup |
| ovarian | 25 | 94 | 02_organoid_only/OV_GSE229119_olaparib, 02_organoid_only/OV_vias2023, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 05_A3_B1/figures, 07_B3_models/figures, 15_fixes/backup |
| lung | 32 | 97 | 02_organoid_only/LUAD_GSE165611_substitute, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 07_B3_models/figures, 15_fixes/backup |
| sarcoma | 3 | 22 | 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 15_fixes/backup |
| all sets / cross-cancer | 67 | 190 | 00_checklist/checklist.tsv, 00_checklist/tables, 02_organoid_only/figures, 02_organoid_only/tables, 03_organoid_plus_patient/figures, 03_organoid_plus_patient/tables, 05_A3_B1/figures, 05_A3_B1/tables, 07_B3_models/figures, 07_B3_models/tables, 08_depmap/figures, 08_depmap/tables, 09_registered_test |
