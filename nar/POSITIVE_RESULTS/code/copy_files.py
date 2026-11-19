"""Copies (never moves) figures, captions and tables into nar/POSITIVE_RESULTS/<CANCER>/{figures,tables}.
Figure entries are base paths (no extension): PNG, PDF and .caption.md are copied; a missing caption is written
from the source report text and marked 'caption written by W5, not by the source agent'.
Writes COPY_MANIFEST.tsv (source, destination, sha256). Source files are opened read-only."""
import os, shutil, hashlib, csv
ROOT = '/home/user/organoid_biomarker_detection'; NAR = f'{ROOT}/nar'; OUT = f'{NAR}/POSITIVE_RESULTS'
def sha(p):
    h = hashlib.sha256(); h.update(open(p, 'rb').read()); return h.hexdigest()[:16]
FIG = {  # cancer -> list of figure base paths relative to ROOT
 'BREAST': ['nar/19_model_benchmark/figures/fig_model_mean_auc', 'nar/19_model_benchmark/figures/fig_model_heatmap', 'nar/19_model_benchmark/figures/fig_roc_top3',
            'nar/18_other_cohorts/figures/fig1_breast_loco_forest',
            'nar/20_breast_pathways/figures/fig_auc_vs_k', 'nar/20_breast_pathways/figures/fig_top30_weights', 'nar/20_breast_pathways/figures/fig_stability',
            'nar/20_breast_pathways/figures/fig_themes', 'nar/20_breast_pathways/figures/fig_proliferation', 'nar/20_breast_pathways/figures/fig_subgroup_forest'],
 'LUNG': ['nar/17_lung/figures/lung_untreated_prognosis', 'nar/17_lung/figures/lung_auc_by_cohort', 'nar/17_lung/figures/lung_forest_interaction', 'nar/17_lung/figures/lung_gse14814_km',
          'nar/19_model_benchmark/figures/fig_bar_lung_untreated_prognosis', 'nar/05_A3_B1/figures/A3_Fig1_master_forest', 'nar/18_other_cohorts/figures/fig2_other_cancers_loco_forest'],
 'BLADDER': ['nar/02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full', 'nar/02_organoid_only/BLCA_lee2018/current/figures/BLCA_CISPLATIN_full',
             'nar/02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_ERCC2_concordance', 'nar/02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_isogenic_T24',
             'nar/13_gaps_bio/figures/G8_mutation_concordance_null_and_adjustment', 'nar/03_organoid_plus_patient/figures/main_bladder',
             'nar/19_model_benchmark/figures/fig_bar_bladder_ERCC2_currentdata', 'nar/19_model_benchmark/figures/fig_bar_bladder_ERCC2_paperdata',
             'nar/19_model_benchmark/figures/fig_bar_bladder_logrank_paperdata', 'nar/19_model_benchmark/figures/fig_bar_bladder_logrank_currentdata',
             'nar/05_A3_B1/figures/B1_Fig2_forest_bladder', 'nar/05_A3_B1/figures/A3_Fig1_master_forest'],
 'COLORECTAL': ['nar/02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL',
                'nar/02_organoid_only/COAD_vandewetering/figures/Fig4b_SuppTable3_isogenic_HCT8_5FU', 'nar/02_organoid_only/COAD_vandewetering/figures/Fig5a_cetuximab_BRAF_V600E',
                'nar/10_completeness/figures/Fig5_SuppFig12_SuppFig18_mutation_concordance', 'nar/13_gaps_bio/figures/G2_decomposition_colon_5FU',
                'nar/13_gaps_bio/figures/G8_mutation_concordance_null_and_adjustment', 'nar/19_model_benchmark/figures/fig_bar_colon_cetuximab_BRAF',
                'nar/05_A3_B1/figures/B1_Fig1_forest_colorectal', 'nar/05_A3_B1/figures/A3_Fig1_master_forest', 'nar/18_other_cohorts/figures/fig2_other_cancers_loco_forest'],
 'RECTAL': ['nar/05_A3_B1/figures/B1_Fig1_forest_colorectal', 'nar/19_model_benchmark/figures/fig_bar_rectal_GSE87211_5FU', 'nar/05_A3_B1/figures/B1_Fig7_heatmap_signature_by_cohort'],
 'OVARIAN': ['nar/18_other_cohorts/figures/fig2_other_cancers_loco_forest', 'nar/19_model_benchmark/figures/fig_bar_ovarian_olaparib_GSE63885', 'nar/19_model_benchmark/figures/fig_bar_ovarian_olaparib_GSE9891',
             'nar/19_model_benchmark/figures/fig_bar_ovarian_olaparib_TCGAOV', 'nar/19_model_benchmark/figures/fig_bar_ovarian_platinum_joint_survival',
             'nar/19_model_benchmark/figures/fig_bar_ovarian_paclitaxel_joint_survival', 'nar/19_model_benchmark/figures/fig_bar_ovarian_platinum_joint_response',
             'nar/03_organoid_plus_patient/figures/main_ovarian', 'nar/10_completeness/figures/Fig4bc_isogenic_other_cancers',
             'nar/22_user_data_results/figures/fig5_ovary_apply_forest', 'nar/22_user_data_results/figures/fig7_ovary_ic50_signature', 'nar/05_A3_B1/figures/B1_Fig3_forest_ovarian'],
 'PAN_CANCER_DEPMAP': ['nar/08_depmap/figures/fig3_drug_specificity_and_bladder', 'nar/08_depmap/figures/fig4_bh3_mimetics_and_crispr', 'nar/13_gaps_bio/figures/G7_bh3_beyond_simple_markers',
                       'nar/08_depmap/figures/S4_permutation_and_random_controls', 'nar/08_depmap/figures/S2_gdsc1_vs_gdsc2'],
}
TAB = {
 'BREAST': ['robust/results/train_combined/part1_response_leave_one_cohort_out.tsv', 'robust/results/train_combined/breast_model_pathways.tsv',
            'nar/19_model_benchmark/tables/breast_loco_all_models.tsv', 'nar/19_model_benchmark/tables/model_comparison.tsv', 'nar/19_model_benchmark/tables/random_pathway_sets_auc.tsv',
            'nar/19_model_benchmark/tables/l2_reproduction_check.tsv', 'nar/18_other_cohorts/tables/part1_v2.tsv', 'nar/18_other_cohorts/tables/part2_v2.tsv',
            'nar/20_breast_pathways/tables/auc_vs_k.tsv', 'nar/20_breast_pathways/tables/top30_pathways.tsv', 'nar/20_breast_pathways/tables/stability_selection.tsv',
            'nar/20_breast_pathways/tables/proliferation_tests.tsv', 'nar/20_breast_pathways/tables/subgroup_auc_pooled.tsv', 'nar/20_breast_pathways/tables/minimal_subset.tsv',
            'nar/20_breast_pathways/tables/family_weights.tsv'],
 'LUNG': ['nar/17_lung/tables/lung_untreated_prognosis.tsv', 'nar/17_lung/tables/lung_loco_response.tsv', 'nar/17_lung/tables/lung_adjuvant_interaction.tsv', 'nar/17_lung/tables/lung_cohorts_built.tsv',
          'nar/17_lung/tables/lung_arm_counts.tsv', 'nar/19_model_benchmark/tables/positives_all_models_lung.tsv', 'nar/05_A3_B1/tables/A3_master_results.tsv'],
 'BLADDER': ['nar/02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv', 'nar/02_organoid_only/BLCA_lee2018/extras/tables/isogenic_T24_tests.tsv',
             'nar/02_organoid_only/BLCA_lee2018/extras/tables/paper_biomarker_amino_acid_KM.tsv', 'nar/02_organoid_only/BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv',
             'nar/13_gaps_bio/tables/G8_ERCC2_null_and_adjustment.tsv', 'nar/03_organoid_plus_patient/tables/summary_main.tsv', 'nar/03_organoid_plus_patient/tables/interaction_all.tsv',
             'nar/03_organoid_plus_patient/tables/BLCA_CISPLATIN__shortlist__eval.tsv', 'nar/03_organoid_plus_patient/tables/BLCA_GEMCITABINE__shortlist__eval.tsv',
             'nar/05_A3_B1/tables/A3_master_results.tsv', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'nar/19_model_benchmark/tables/positives_all_models_rest.tsv'],
 'COLORECTAL': ['nar/02_organoid_only/tables/cross_set_summary.tsv', 'nar/02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv',
                'nar/02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv', 'nar/10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv',
                'nar/10_completeness/tables/SuppFig18_leave_k_out_summary.tsv', 'nar/13_gaps_bio/tables/G8_BRAF_null_and_adjustment.tsv', 'nar/13_gaps_bio/tables/G2_decomposition_colon_5FU.tsv',
                'nar/18_other_cohorts/tables/part1_v2.tsv', 'nar/05_A3_B1/tables/A3_master_results.tsv', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'nar/19_model_benchmark/tables/positives_all_models_rest.tsv'],
 'RECTAL': ['nar/05_A3_B1/tables/B1_external_results.tsv', 'nar/05_A3_B1/tables/B1_frozen_signatures.tsv', 'nar/19_model_benchmark/tables/positives_all_models_rest.tsv',
            'nar/03_organoid_plus_patient/tables/interaction_all.tsv', 'nar/18_other_cohorts/tables/part1_v2.tsv'],
 'OVARIAN': ['nar/18_other_cohorts/tables/part1_v2.tsv', 'nar/22_user_data_results/tables/ovary_signature_survival.tsv', 'nar/22_user_data_results/tables/paper_method_ovary.tsv',
             'nar/22_user_data_results/tables/ovary_BRCA_somatic_survival_platinum.tsv', 'nar/10_completeness/tables/Fig4bc_isogenic_other_cancers.tsv',
             'nar/03_organoid_plus_patient/tables/summary_main.tsv', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'nar/19_model_benchmark/tables/positives_all_models_rest.tsv'],
 'PAN_CANCER_DEPMAP': ['nar/08_depmap/tables/A_e_bh3_vs_all_drugs.tsv', 'nar/08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv', 'nar/13_gaps_bio/tables/G7_regression_added_value.tsv',
                       'nar/13_gaps_bio/tables/G7_partial_correlations.tsv', 'nar/08_depmap/tables/S4_permutation_random_controls.tsv', 'nar/13_gaps_bio/tables/G3_bh3_vs_drugs_organoids.tsv'],
}
REPORTS = {'BREAST': ['nar/19_model_benchmark/REPORT_breast.md', 'nar/20_breast_pathways/REPORT.md'], 'LUNG': ['nar/17_lung/REPORT.md'],
           'BLADDER': ['nar/03_organoid_plus_patient/REPORT.md'], 'PAN_CANCER_DEPMAP': ['nar/08_depmap/REPORT.md'], 'OVARIAN': ['nar/22_user_data_results/REPORT.md']}
man = []; miss = []
def cp(src, dst):
    if not os.path.exists(f'{ROOT}/{src}'): miss.append(src); return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(f'{ROOT}/{src}', dst); man.append((src, dst.replace(OUT + '/', ''), sha(dst)))
def prefix(src):  # unique name: "<folder number>__<file>" when two sources share a file name
    parts = src.split('/'); mid = [p for p in parts[2:-1] if p not in ('figures', 'tables')]
    return '__'.join([parts[1].split('_')[0]] + mid + [parts[-1]]) if mid else f'{parts[1].split("_")[0]}__{parts[-1]}'
for c, lst in FIG.items():
    for b in lst:
        pre = prefix(b + '.png')
        stem = pre[:-4]
        for ext in ['.png', '.pdf', '.caption.md']:
            if os.path.exists(f'{ROOT}/{b}{ext}'): cp(b + ext, f'{OUT}/{c}/figures/{stem}{ext}')
            elif ext == '.caption.md' and os.path.exists(f'{ROOT}/{b}.png'):
                os.makedirs(f'{OUT}/{c}/figures', exist_ok=True)
                open(f'{OUT}/{c}/figures/{stem}.caption.md', 'w').write(f'Caption written by W5 (the source agent saved this figure without a caption file): {b}.png. See the source report listed in the README of this folder for what it shows.\n')
                man.append((b + '.png (no source caption)', f'{c}/figures/{stem}.caption.md', 'written'))
            elif ext != '.caption.md': miss.append(b + ext)
for c, lst in TAB.items():
    for t in lst:
        pre = prefix(t)
        cp(t, f'{OUT}/{c}/tables/{pre}')
for c, lst in REPORTS.items():
    for t in lst: cp(t, f'{OUT}/{c}/source_reports/{prefix(t)}')
with open(f'{OUT}/COPY_MANIFEST.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['source', 'copied_to', 'sha256_16']); w.writerows(man)
print(len(man), 'files copied;', 'missing:', miss)
