"""Build completeness.tsv: status of every checklist item for (a) the paper's own cancer/drug and (b) the other sets,
before and after this agent's fills. Every path is checked for existence (glob)."""
from common import *
import glob
ck = pd.read_csv(NAR / "00_checklist/checklist.tsv", sep="\t")
G = "02_organoid_only/"; C = "10_completeness/"; B = "07_B3_models/ (in progress by B3 agent)"
GI_P = G + "COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png"
GI_ALL = G + "*/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_*.png"; OT_ALL = G + "*/figures/*_full.png"
BL = G + "BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png"; BX = G + "BLCA_lee2018/extras/"
COEF = G + "*/figures/Fig2b_SuppFig2_SuppFig3_coefficients_*.png"
SETS_ALL = "colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao)"
D, P, M, IP, NA = "DONE", "PARTIAL", "MISSING", "IN_PROGRESS_B3", "N/A"
# ID: (paper_before, paper_after, paper_paths, paper_check, gen_before, gen_after, gen_paths, note)
S = {
 "Fig1a": (M, D, C + "figures/Fig1a_workflow_schematic.png", "drawing", M, D, C + "figures/Fig1c_data_map_all_sets.png", "Fig1c maps all 57 set x cohort x drug analyses"),
 "Fig1b": (M, D, C + "figures/Fig1b_network_proximity_schematic.png", "drawing", M, D, C + "figures/Fig1b_network_proximity_schematic.png", "one drawing serves all sets"),
 "Fig2a": (D, D, GI_P, "37 proximal pathways, rank-1 BH3-only", D, D, GI_ALL + ";" + OT_ALL, SETS_ALL),
 "Fig2b": (D, D, G + "COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png", "", D, D, COEF + ";" + OT_ALL, SETS_ALL),
 "Fig2c": (M, D, C + "figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png", "TYMS + 15 BH3-only genes in network, mean d 2.00", M, D, C + "figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png", "18 panels, every headline drug"),
 "Fig2d-left": (D, D, GI_P, "n=114, p=0.014 (paper 0.014) verified in tables/verify_headline_biomarkers.tsv", D, D, GI_ALL + ";" + OT_ALL + ";" + G + "figures/cross_set_forest.png", SETS_ALL),
 "Fig2d-right": (D, D, GI_P, "n=298, p=0.161 (paper 0.16)", D, D, GI_ALL + ";" + OT_ALL, SETS_ALL),
 "Fig2e": (D, D, BL, "47 -> 30 proximal pathways (paper 30)", D, D, OT_ALL, "as Fig2a"),
 "Fig2f": (D, D, BL, "9 organoids", D, D, OT_ALL, "as Fig2b"),
 "Fig2g": (M, D, C + "figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png", "4 targets, amino acid pathway 16 genes in network", M, D, C + "figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png", ""),
 "Fig2h-left": (D, D, BL + ";" + BX + "figures/S_BLCA_paper_biomarker_KM.png", "rank-1 p=0.044, paper's pathway p=0.031 (paper 0.01); AUC used, IC50 not public", D, D, OT_ALL, ""),
 "Fig2h-right": (D, D, BL, "n=288, p=0.93", D, D, OT_ALL, ""),
 "Fig3a": (M, D, C + "figures/Fig3a_comparison_feature_sets_schematic.png", "drawing", M, D, C + "figures/Fig3a_comparison_feature_sets_schematic.png", ""),
 "Fig3b-COAD": (M, IP, B, "", M, IP, B, "B3 bench.py: whole transcriptome for every set"),
 "Fig3b-BLCA": (D, D, BX + "figures/S_BLCA_comparison_feature_sets.png", "whole transcriptome p=0.71", M, IP, B, ""),
 "Fig3c-COAD": (M, IP, B, "", M, IP, B, "B3 bench.py: all pathways for every set"),
 "Fig3c-BLCA": (D, D, BX + "figures/S_BLCA_comparison_feature_sets.png", "whole pathway p=0.59 / 0.86", M, IP, B, ""),
 "Fig3d-COAD": (M, D, C + "figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png", "first neighbours of TYMS p=0.719", M, D, C + "tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv", "all 16 headline sets (also in B3 bench)"),
 "Fig3d-BLCA": (D, D, BX + "figures/S_BLCA_comparison_feature_sets.png", "neighbours 1 step p=0.70", M, D, C + "tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv", ""),
 "Fig3e-COAD": (M, IP, B, "", M, IP, B, ""), "Fig3e-BLCA": (M, IP, B, "", M, IP, B, ""),
 "Fig3f-COAD": (M, IP, B, "", M, IP, B, ""), "Fig3f-BLCA": (M, IP, B, "", M, IP, B, ""),
 "Fig4a": (M, D, C + "figures/Fig4a_isogenic_validation_schematic.png", "drawing", M, D, C + "figures/Fig4a_isogenic_validation_schematic.png", ""),
 "Fig4b": (D, D, G + "COAD_vandewetering/figures/Fig4b_SuppTable3_isogenic_HCT8_5FU.png;08_depmap/figures/S7_isogenic_pairs.png", "p=0.0085 (paper 0.0085)", M, P, C + "figures/Fig4bc_isogenic_other_cancers.png", "ovarian A2780 cisplatin (GSE15709) and liver HepG2/Huh7 sorafenib (GSE322742) added; no public pair found for pancreas gemcitabine, sarcoma, melanoma, PDX drugs"),
 "Fig4c": (D, D, BX + "figures/S_BLCA_isogenic_T24.png", "amino acid p=0.00022 (paper's NES; raw T24 data not public)", M, P, C + "figures/Fig4bc_isogenic_other_cancers.png", "as Fig4b"),
 "Fig4d": (M, D, C + "figures/Fig4d_bootstrap_schematic.png", "drawing", M, D, C + "figures/Fig4d_bootstrap_schematic.png", ""),
 "Fig4e": (D, D, GI_P, "share 0.0074 (paper 0.0012)", D, D, GI_ALL + ";" + OT_ALL, SETS_ALL),
 "Fig4f": (D, D, BX + "figures/S_BLCA_paper_bootstrap_recount.png", "recount 0.014 (paper 0.014)", D, D, OT_ALL, ""),
 "Fig5a": (D, D, G + "COAD_vandewetering/figures/Fig5a_cetuximab_BRAF_V600E.png", "p=0.031 (paper 0.037); recomputed here 0.00048 on current data", M, P, C + "figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png", "KRAS (colon), EGFR (lung), BRCA1/2 (ovary) added; no accepted mutation marker for gemcitabine, sorafenib, dacarbazine"),
 "Fig5b": (D, D, BX + "figures/S_BLCA_ERCC2_concordance.png", "p=0.0030 (paper 0.002)", M, P, C + "figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png", "as Fig5a"),
 "SuppFig1a-d": (D, D, G + "COAD_vandewetering/figures/SuppFig1_PCA_models_vs_TCGA.png;" + BX + "figures/S_BLCA_PCA_batch.png", "BLCA pathways PC1/PC2 27.1/17.8 (paper 27.0/17.8)", P, D, C + "tables/SuppFig1_PCA_variance_all_sets.tsv", "variance numbers for all 15 bundles"),
 "SuppFig1e-h": (D, D, BX + "tables/SuppFig1_PCA_variance.tsv", "scaled 20.7/18.8 (paper 20.7/18.8)", P, D, C + "tables/SuppFig1_PCA_variance_all_sets.tsv", ""),
 "SuppFig2a": (D, D, G + "COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png", "", D, D, COEF + ";" + OT_ALL, ""),
 "SuppFig2b": (D, D, G + "COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png", "", D, D, COEF, ""),
 "SuppFig2c": (D, D, BL, "", D, D, OT_ALL, ""), "SuppFig2d": (D, D, BL, "", D, D, OT_ALL, ""),
 "SuppFig3a-c": (D, D, G + "COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png", "SVR and LR p=0.014 too", D, D, COEF + ";" + OT_ALL, ""),
 "SuppFig3d-f": (D, D, BL, "same rank 1 in all 3 models", D, D, OT_ALL, ""),
 "SuppFig4a": (M, D, C + "figures/SuppFig4_cross_test_matrix_all_sets.png", "p=0.164 (paper 0.16)", M, D, C + "tables/SuppFig4_cross_test_matrix.tsv", "17 biomarkers x 16 cohorts"),
 "SuppFig4b": (D, D, BX + "figures/S_BLCA_comparison_feature_sets.png;" + C + "tables/SuppFig4_cross_test_matrix.tsv", "p=0.98 (paper 0.86)", M, D, C + "figures/SuppFig4_cross_test_matrix_all_sets.png", ""),
 "SuppFig5a": (D, D, G + "COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_FLUOROURACIL.png", "", D, D, G + "*/figures/SuppFig5_*.png;" + OT_ALL, "OTHER sets: p vs k panel g"),
 "SuppFig5b": (D, D, BX + "figures/S_BLCA_multipathway_vs_paper.png", "Ridge k=2 p=0.060 (paper 0.022)", D, D, OT_ALL, ""),
 "SourceData-multipathway": (D, D, BX + "tables/multipathway_k_vs_paper_sourcedata.tsv", "", D, D, G + "*/tables/*multipathway*", ""),
 "SuppFig6a": (M, IP, B, "", M, IP, B, "B3 centrality.py"), "SuppFig6b": (M, IP, B, "", M, IP, B, ""), "SuppFig6c": (M, IP, B, "", M, IP, B, ""),
 "SuppFig7a": (M, IP, B, "", M, IP, B, ""), "SuppFig7b": (M, IP, B, "", M, IP, B, ""), "SuppFig7c": (M, IP, B, "", M, IP, B, ""),
 "SuppFig8a": (M, D, C + "tables/SuppFig8_target_genes_alone.tsv", "TYMS alone", M, D, C + "figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png", "all 16 headline sets"),
 "SuppFig8b": (D, D, BX + "tables/comparison_feature_sets_treated_KM.tsv", "A2M 0.74, ATOX1 0.37, MPG 0.68, TF 0.65", M, D, C + "tables/SuppFig8_target_genes_alone.tsv", ""),
 "SuppFig9a": (P, D, BX + "tables/comparison_feature_sets_treated_KM.tsv;" + C + "tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv", "2 steps: BLCA 0.94, COAD 0.660", M, D, C + "figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png", ""),
 "SuppFig9b": (P, D, C + "tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv", "3 steps: BLCA 0.99, COAD 0.774", M, D, C + "tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv", ""),
 "SuppFig10a": (D, D, G + "COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_authors_FLUOROURACIL.png;" + C + "figures/SuppFig10_size_vs_p_all_sets.png", "R2 0.012 (GI) / 0.003 (here); paper 0.0026", P, D, C + "tables/SuppFig10_size_vs_p_summary_all_sets.tsv", ""),
 "SuppFig10b": (M, D, C + "figures/SuppFig10_size_vs_p_all_sets.png", "R2 0.000", P, D, C + "tables/SuppFig10_size_vs_p_summary_all_sets.tsv", ""),
 "SuppFig11": (D, D, G + "COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_CETUXIMAB_SuppFig11.png", "", D, D, G + "COAD_isella_pdx/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_*.png", "cetuximab only exists in colorectal sets"),
 "SuppFig12a-e": (P, D, C + "figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png", "transcriptome 0.0030, all pathways 0.94, neighbours 0.9991/0.047/0.0069", NA, NA, "", "cetuximab-specific; mutation tests for other drugs under Fig5a"),
 "SuppFig13a-e": (D, D, BX + "tables/ERCC2_mutation_concordance.tsv", "transcriptome 0.065 (paper 0.04)", NA, NA, "", "cisplatin-ERCC2 specific"),
 "SuppFig14a-h": (M, D, C + "figures/SuppFig14_cell_lines_vs_TCGA_PCA.png", "COAD genes raw 38.7/6.8, scaled 11.1/8.5", M, D, C + "tables/SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv", "8 cancers"),
 "SuppFig14i-p": (M, D, C + "figures/SuppFig14_cell_lines_vs_TCGA_PCA.png", "BLCA pathways raw 64.1/12.5, scaled 22.4/17.1 (mean-expression score, not ssGSEA)", M, D, C + "tables/SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv", ""),
 "SuppFig15a": (D, D, G + "COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png", "p=0.583 authors' patients (paper 0.58)", P, P, G + "LUAD_GSE165611_substitute/figures/*_full.png", "lung set is itself GDSC-trained; not run for pancreas, liver, ovary, sarcoma"),
 "SuppFig15b": (D, D, G + "COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png", "p=0.506 (paper 0.51)", P, P, "", "as SuppFig15a"),
 "SuppFig15c": (D, D, BX + "figures/S_BLCA_GDSC_trained.png", "0.82/0.012/0.011", P, P, "", "as SuppFig15a"),
 "SuppFig15d": (D, D, BX + "figures/S_BLCA_GDSC_trained.png", "0.21/0.69/0.24", P, P, "", "as SuppFig15a"),
 "SuppFig16": (P, D, C + "figures/SuppFig16_proximity_all_drugs.png", "44 paper columns + 8 added drugs", P, D, C + "tables/SuppFig16_proximity_all_drugs.tsv.gz;" + G + "*/figures/supp_proximity_heatmap.png", ""),
 "SuppFig17a": (P, D, GI_P + ";" + C + "figures/SuppFig17_split_603010_all_sets.png", "60/10/30 median rho 0.37 (LOO 0.44)", P, D, C + "tables/SuppFig17_split_603010_summary.tsv", ""),
 "SuppFig17b": (P, D, BL + ";" + C + "figures/SuppFig17_split_603010_all_sets.png", "60/10/30 median rho 0.50 (paper 0.89 with IC50)", P, D, C + "tables/SuppFig17_split_603010_summary.tsv", ""),
 "SuppFig18": (M, D, C + "tables/SuppFig18_leave_k_out_summary.tsv", "share p<0.05: k=1 1.00, k=2 0.93, k=3 0.94, k=4 0.87", NA, NA, "", "cetuximab-specific"),
 "SuppTable1": (D, D, G + "tables/SuppTable1_clinical_features_by_predicted_group.tsv", "age 0.260, sex 0.188, T 0.435, N 0.440 (paper 0.26/0.19/0.44/0.44)", P, D, C + "tables/SuppTable1_2_clinical_balance_all_sets.tsv", "GI sets by 02; all headline sets here"),
 "SuppTable2": (D, D, BX + "tables/SuppTable2_clinical_balance.tsv", "chi2 p 0.11-1.0", P, D, C + "tables/SuppTable1_2_clinical_balance_all_sets.tsv", ""),
 "SuppTable3": (D, D, G + "COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv;" + BX + "tables/isogenic_T24_tests.tsv", "", M, P, C + "tables/Fig4bc_isogenic_other_cancers.tsv", "as Fig4b"),
 "SuppTable4": (P, IP, B, "bootstrap draws exist in 02 tables; formatted table by B3", P, IP, B, ""),
 "SourceData-STRING": (D, D, G + "BLCA_OV_LUAD_OTHER/_data/moesm3/Figure_2a_f.tsv", "network re-used by all agents", D, D, G + "BLCA_OV_LUAD_OTHER/_data/moesm3/Figure_2a_f.tsv", ""),
 "SourceData-DLhyper": (M, IP, B, "", M, IP, B, ""),
 "Text-method-params": (D, D, G + "REPORT_GI.md;" + G + "REPORT_BLCA_OV_LUAD_OTHER.md", "cut-off -1.2816, RidgeCV 0.1-0.9, 10,000 bootstraps", D, D, G + "tables/cross_set_cutoff_sensitivity.tsv", ""),
 "Text-claim-specificity": (P, D, C + "tables/SuppFig4_cross_test_summary.json", "own 6/17 p<0.05; off-diagonal 9/245 (3.7%)", P, D, C + "tables/SuppFig4_cross_test_matrix.tsv", ""),
}
rows = []
for r in ck.itertuples():
    pb, pa, pp, pc, gb, ga, gp, note = S[r.ID]
    def exists(paths):
        if not paths or paths.startswith("07_B3"):
            return ""
        return ";".join(str(len(glob.glob(str(NAR / p.strip())))) for p in paths.split(";") if p.strip())
    rows.append(dict(ID=r.ID, type=r.type, cancer_drug=r.cancer_drug, paper_before=pb, paper_after=pa, paper_paths=pp,
                     paper_files_found=exists(pp), paper_number_check=pc, generalisation_before=gb, generalisation_after=ga,
                     generalisation_paths=gp, generalisation_files_found=exists(gp), note=note))
T = pd.DataFrame(rows); T.to_csv(OUT / "completeness.tsv", sep="\t", index=False)
bad = T[T.paper_files_found.str.contains(r"(^|;)0(;|$)") | T.generalisation_files_found.str.contains(r"(^|;)0(;|$)")]
print("missing files:\n", bad[["ID", "paper_paths", "paper_files_found", "generalisation_paths", "generalisation_files_found"]].to_string())
cnt = {c: T[c].value_counts().to_dict() for c in ("paper_before", "paper_after", "generalisation_before", "generalisation_after")}
json.dump(cnt, open(TAB / "completeness_counts.json", "w"), indent=1); print(cnt, len(T))
