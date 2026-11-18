# Completeness check against the Kong et al. 2020 checklist (77 items)

Folder `nar/10_completeness/`. Status is given twice for every checklist item: (a) the paper's own cancer and drug (exact
reproduction: colon 5-FU, bladder cisplatin, colon cetuximab) and (b) the same item for the other sets (colon current data,
bladder current data, pancreas Tiriac and Shi, liver Ji and Broutier, ovarian Vias, lung cell-line substitute, sarcoma,
PDX: Isella colorectal and Gao breast/colorectal/lung/melanoma/pancreas). Full table: `completeness.tsv` (one row per item,
with the files that cover it, how many files the path matched, and the number that was checked).

Status words: DONE = figure/table exists and its number was checked; PARTIAL = done for some sets only, or a close
substitute; MISSING = not produced; IN_PROGRESS_B3 = being produced now by the B3 agent (`nar/07_B3_models/`: Fig 3b-f
for every set, Supp Figs 6-7 centrality, Supp Table 4, deep-learning grid); N/A = the item only exists for one drug.

## Counts before and after this agent's work

| | DONE | PARTIAL | MISSING | IN_PROGRESS_B3 | N/A |
|---|---|---|---|---|---|
| (a) paper's own case, before | 42 | 8 | 27 | 0 | 0 |
| (a) paper's own case, after | 63 | 0 | 0 | 14 | 0 |
| (b) other sets, before | 22 | 15 | 37 | 0 | 3 |
| (b) other sets, after | 49 | 9 | 0 | 16 | 3 |

Before = what the earlier agents (02, 03, 08, 09) had produced. In the 'before' rows, items now marked IN_PROGRESS_B3 are
counted as MISSING (14 for (a), 16 for (b)).

## How the earlier numbers were checked

- Every path in `completeness.tsv` was checked with a file search; all paths exist (column `*_files_found`, all > 0).
- Headline survival tests re-derived from the shared data bundles and the agents' coefficient tables (`tables/verify_headline_biomarkers.tsv`): 15 of 17 reproduce the reported treated-patient p exactly (to 6 decimals). The two others: Gao pancreas PDX (reported 0.14 in the report; recomputed 0.137, the table lookup used a
different project label) and Broutier liver (0.089 here vs 0.139 reported, because the shared bundle holds TCGA-LIHC
scores computed for the Ji gene universe; Broutier was left out of the new cross-set analyses for this reason).
- Spot checks of report numbers against tables (all matched): colon 5-FU p = 0.014 (n = 114); bladder p = 0.0435;
Supp Table 1 age 0.260, sex 0.188, T 0.435, N 0.440; ERCC2 p = 0.0030; T24 isogenic p = 0.00022; HCT-8 p = 0.0085 (paper) /
0.039 (recomputed); BRAF V600E p = 0.031; bladder PCA 27.1/17.8 and 20.7/18.8; GDSC 5-FU 0.708 / 0.583; organoid+patient
CRC_5FU joint HR 0.89 (0.80-1.00), C 0.532; registered test gemcitabine HR 1.12 (0.94-1.32), p = 0.199; DepMap GDSC1
rho 0.054 (n = 44).

## What was added here (all in `nar/10_completeness/figures` and `tables`)

1. **Schematics** (Fig 1a workflow, 1b network proximity, 3a comparison feature sets, 4a isogenic check, 4d bootstrap) and
   a new Fig 1c map of all 57 set x cohort x drug analyses.
2. **Fig 2c / 2g network pictures** for the paper's two drugs and every other headline drug (18 panels).
3. **Supp Fig 4 cross-test matrix** (17 biomarkers x 16 treated cohorts). Paper's two cross-tests: colon 5-FU patients with the bladder amino acid pathway p = 0.164 (paper 0.16); bladder cisplatin patients with the colon BH3-only pathway p = 0.98 (paper 0.86). Own-cohort tests p < 0.05: 6 of 17; other-cohort tests p < 0.05: 9 of 245 (3.7%, about the 5% expected by chance).
   What I think: the biomarkers are not general survival markers across cancers; this supports the paper's specificity claim.
4. **Supp Fig 8 (target genes alone) and Supp Fig 9 (1-, 2-, 3-step network neighbours)** for colon (paper data: TYMS
   neighbours 1/2/3 steps p = 0.719 / 0.660 / 0.774) and all 16 headline sets. Neighbour models with p < 0.05: 1 of 48.
5. **Supp Fig 10 size vs p** for all 16 sets: R2 from 0.000 to 0.012 (paper 0.0026 and 0.0008).
   What I think: pathway size does not drive the survival splits in any set.
6. **Supp Fig 14 with explained-variance numbers**, cell lines vs TCGA for 8 cancers. Colon genes raw 38.7/6.8%, scaled 11.1/8.5%. Supp Fig 1
   variance numbers for all 15 lab-model bundles (`tables/SuppFig1_PCA_variance_all_sets.tsv`).
7. **Supp Fig 16** proximity heatmap: the paper's 44 columns plus 8 drugs added for the other sets.
8. **Supp Fig 17 as the paper's 60/10/30 split** (100 repeats), every set; median held-out Spearman rho: BLCA_cur_CIS 1.00; BLCA_paper_CIS 0.50; COADREAD_5FU 0.43; COAD_isella_CET 0.45; COAD_paper_5FU 0.37; LIHC_ji_SOR -0.10; LUAD_cells_CIS 0.00; OV_vias_CARBO 0.40; PAAD_shi_GEM 0.00; PAAD_tiriac_GEM 0.21; PDX_BRCA_PTX 0.10; PDX_CM_DTIC 0.08; PDX_CRC_5FU 0.16; PDX_NSCLC_ERL 0.05; PDX_PDAC_GEM 0.09; SARC_GEM 0.00.
9. **Fig 5 for other cancers** (one-sided Mann-Whitney, all TCGA tumours with mutation calls): Colorectal organoids (vdW), current BRAF V600E p = 0.00048 (51 vs 522); Colorectal organoids (vdW), current KRAS G12/G13/Q61/A146 p = 0.88 (229 vs 344); Colorectal PDX (Isella) KRAS G12/G13/Q61/A146 p = 0.97 (172 vs 252); Colorectal PDX (Isella) BRAF V600E p = 0.27 (50 vs 374); Lung PDX (Gao) EGFR activating p = 1 (52 vs 442); Lung adenocarcinoma cell lines EGFR activating p = 0.66 (52 vs 442); Ovarian organoids (Vias) BRCA1/2 somatic p = 0.83 (17 vs 271); Ovarian organoids (Vias) BRCA1/2 somatic p = 0.83 (17 vs 271).
   What I think: the colon organoid cetuximab model recovers BRAF V600E again (p = 0.00048) but does not pick up KRAS
   mutations; the erlotinib and platinum models do not separate EGFR or BRCA1/2 mutants in the expected direction.
   No accepted mutation marker exists for gemcitabine, sorafenib or dacarbazine, so no such test was possible there.
10. **Supp Fig 12** (colon cetuximab, other feature sets vs BRAF V600E) and **Supp Fig 18** (leave-k-out organoids): share of
    refits still giving p < 0.05: k=1: 1.00, k=2: 0.93, k=3: 0.94, k=4: 0.87.
11. **Fig 4b/c for other cancers**: GSE15709 ovarian A2780 | Bladder organoids, paper data (rank 1): p = 0.22, predicted direction; GSE15709 ovarian A2780 | Bladder, paper's amino acid pathway: p = 0.053, opposite direction; GSE15709 ovarian A2780 | Bladder organoids, current data (rank 1): p = 0.14, predicted direction; GSE15709 ovarian A2780 | Ovarian organoids, carboplatin (rank 1): p = 0.017, predicted direction; GSE15709 ovarian A2780 | Ovarian organoids, cisplatin (rank 1): p = 0.017, predicted direction; GSE15709 ovarian A2780 | Lung cell lines, cisplatin (rank 1): p = 0.029, predicted direction; GSE322742 liver HepG2 parental | Liver organoids (Ji) (rank 1): p = 0.12, predicted direction; GSE322742 liver HepG2 parental | Liver organoids (Broutier) (rank 1): p = 0.77, predicted direction; GSE322742 liver Huh7 parental | Liver organoids (Ji) (rank 1): p = 0.55, opposite direction; GSE322742 liver Huh7 parental | Liver organoids (Broutier) (rank 1): p = 0.21, opposite direction.
    What I think: the ovarian and lung platinum pathways differ between A2780 and its cisplatin-resistant copy in the
    predicted direction (p 0.017-0.029, 5 vs 5); the liver pairs have 2 vs 2 samples and cannot decide much.
12. **Supp Tables 1-2 for every set**: clinical balance of predicted groups; chi-square p < 0.05 in 2 of 45 tests (LUAD_cells_CIS stage 0.037, SARC_GEM age 0.010).

## What remains, in one line each

- Fig 3b-f, Supp Figs 6-7, Supp Table 4, deep-learning grid: being produced by the B3 agent (`nar/07_B3_models/`).
- Isogenic pairs for pancreas gemcitabine, sarcoma, melanoma dacarbazine and the PDX drugs: no clean public parental vs
  resistant expression pair was found in this work (GEO search, 2026-11-18).
- Mutation concordance for gemcitabine, sorafenib, dacarbazine: no accepted mutation marker exists.
- Supp Fig 15 (cell-line-trained models) for pancreas, liver, ovary and sarcoma was not run (colon, bladder, and lung,
  whose main analysis is itself cell-line-trained, are done).
- The paper's bladder IC50 values and T24 raw expression are not public (AUC and the paper's NES are used).

## Item by item

| ID | (a) paper case | (b) other sets | Where (a) | Number checked | Note |
|---|---|---|---|---|---|
| Fig1a | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig1a_workflow_schematic.png` | drawing | Fig1c maps all 57 set x cohort x drug analyses |
| Fig1b | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig1b_network_proximity_schematic.png` | drawing | one drawing serves all sets |
| Fig2a | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` | 37 proximal pathways, rank-1 BH3-only | colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao) |
| Fig2b | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png` |  | colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao) |
| Fig2c | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png` | TYMS + 15 BH3-only genes in network, mean d 2.00 | 18 panels, every headline drug |
| Fig2d-left | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` | n=114, p=0.014 (paper 0.014) verified in tables/verify_headline_biomarkers.tsv | colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao) |
| Fig2d-right | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` | n=298, p=0.161 (paper 0.16) | colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao) |
| Fig2e | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | 47 -> 30 proximal pathways (paper 30) | as Fig2a |
| Fig2f | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | 9 organoids | as Fig2b |
| Fig2g | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png` | 4 targets, amino acid pathway 16 genes in network |  |
| Fig2h-left | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | rank-1 p=0.044, paper's pathway p=0.031 (paper 0.01); AUC used, IC50 not public |  |
| Fig2h-right | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | n=288, p=0.93 |  |
| Fig3a | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig3a_comparison_feature_sets_schematic.png` | drawing |  |
| Fig3b-COAD | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  | B3 bench.py: whole transcriptome for every set |
| Fig3b-BLCA | DONE | IN_PROGRESS_B3 (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png` | whole transcriptome p=0.71 |  |
| Fig3c-COAD | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  | B3 bench.py: all pathways for every set |
| Fig3c-BLCA | DONE | IN_PROGRESS_B3 (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png` | whole pathway p=0.59 / 0.86 |  |
| Fig3d-COAD | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png` | first neighbours of TYMS p=0.719 | all 16 headline sets (also in B3 bench) |
| Fig3d-BLCA | DONE | DONE (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png` | neighbours 1 step p=0.70 |  |
| Fig3e-COAD | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| Fig3e-BLCA | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| Fig3f-COAD | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| Fig3f-BLCA | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| Fig4a | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig4a_isogenic_validation_schematic.png` | drawing |  |
| Fig4b | DONE | PARTIAL (was MISSING) | `02_organoid_only/COAD_vandewetering/figures/Fig4b_SuppTable3_isogenic_HCT8_5FU.png` | p=0.0085 (paper 0.0085) | ovarian A2780 cisplatin (GSE15709) and liver HepG2/Huh7 sorafenib (GSE322742) added; no public pair found for pancreas gemcitabine, sarcoma, melanoma, PDX drugs |
| Fig4c | DONE | PARTIAL (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_isogenic_T24.png` | amino acid p=0.00022 (paper's NES; raw T24 data not public) | as Fig4b |
| Fig4d | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/Fig4d_bootstrap_schematic.png` | drawing |  |
| Fig4e | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` | share 0.0074 (paper 0.0012) | colon (current), bladder (current), pancreas x2, liver x2, ovarian, lung substitute, sarcoma, PDX (Isella + Gao) |
| Fig4f | DONE | DONE | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_paper_bootstrap_recount.png` | recount 0.014 (paper 0.014) |  |
| Fig5a | DONE | PARTIAL (was MISSING) | `02_organoid_only/COAD_vandewetering/figures/Fig5a_cetuximab_BRAF_V600E.png` | p=0.031 (paper 0.037); recomputed here 0.00048 on current data | KRAS (colon), EGFR (lung), BRCA1/2 (ovary) added; no accepted mutation marker for gemcitabine, sorafenib, dacarbazine |
| Fig5b | DONE | PARTIAL (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_ERCC2_concordance.png` | p=0.0030 (paper 0.002) | as Fig5a |
| SuppFig1a-d | DONE | DONE (was PARTIAL) | `02_organoid_only/COAD_vandewetering/figures/SuppFig1_PCA_models_vs_TCGA.png` | BLCA pathways PC1/PC2 27.1/17.8 (paper 27.0/17.8) | variance numbers for all 15 bundles |
| SuppFig1e-h | DONE | DONE (was PARTIAL) | `02_organoid_only/BLCA_lee2018/extras/tables/SuppFig1_PCA_variance.tsv` | scaled 20.7/18.8 (paper 20.7/18.8) |  |
| SuppFig2a | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png` |  |  |
| SuppFig2b | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png` |  |  |
| SuppFig2c | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` |  |  |
| SuppFig2d | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` |  |  |
| SuppFig3a-c | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png` | SVR and LR p=0.014 too |  |
| SuppFig3d-f | DONE | DONE | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | same rank 1 in all 3 models |  |
| SuppFig4a | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/SuppFig4_cross_test_matrix_all_sets.png` | p=0.164 (paper 0.16) | 17 biomarkers x 16 cohorts |
| SuppFig4b | DONE | DONE (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png` | p=0.98 (paper 0.86) |  |
| SuppFig5a | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_FLUOROURACIL.png` |  | OTHER sets: p vs k panel g |
| SuppFig5b | DONE | DONE | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_multipathway_vs_paper.png` | Ridge k=2 p=0.060 (paper 0.022) |  |
| SourceData-multipathway | DONE | DONE | `02_organoid_only/BLCA_lee2018/extras/tables/multipathway_k_vs_paper_sourcedata.tsv` |  |  |
| SuppFig6a | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  | B3 centrality.py |
| SuppFig6b | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| SuppFig6c | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| SuppFig7a | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| SuppFig7b | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| SuppFig7c | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| SuppFig8a | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/tables/SuppFig8_target_genes_alone.tsv` | TYMS alone | all 16 headline sets |
| SuppFig8b | DONE | DONE (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/tables/comparison_feature_sets_treated_KM.tsv` | A2M 0.74, ATOX1 0.37, MPG 0.68, TF 0.65 |  |
| SuppFig9a | DONE (was PARTIAL) | DONE (was MISSING) | `02_organoid_only/BLCA_lee2018/extras/tables/comparison_feature_sets_treated_KM.tsv` | 2 steps: BLCA 0.94, COAD 0.660 |  |
| SuppFig9b | DONE (was PARTIAL) | DONE (was MISSING) | `10_completeness/tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv` | 3 steps: BLCA 0.99, COAD 0.774 |  |
| SuppFig10a | DONE | DONE (was PARTIAL) | `02_organoid_only/COAD_vandewetering/figures/SuppFig10_size_vs_p_COAD_authors_FLUOROURACIL.png` | R2 0.012 (GI) / 0.003 (here); paper 0.0026 |  |
| SuppFig10b | DONE (was MISSING) | DONE (was PARTIAL) | `10_completeness/figures/SuppFig10_size_vs_p_all_sets.png` | R2 0.000 |  |
| SuppFig11 | DONE | DONE | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_CETUXIMAB_SuppFig11.png` |  | cetuximab only exists in colorectal sets |
| SuppFig12a-e | DONE (was PARTIAL) |  | `10_completeness/figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png` | transcriptome 0.0030, all pathways 0.94, neighbours 0.9991/0.047/0.0069 | cetuximab-specific; mutation tests for other drugs under Fig5a |
| SuppFig13a-e | DONE |  | `02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv` | transcriptome 0.065 (paper 0.04) | cisplatin-ERCC2 specific |
| SuppFig14a-h | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/SuppFig14_cell_lines_vs_TCGA_PCA.png` | COAD genes raw 38.7/6.8, scaled 11.1/8.5 | 8 cancers |
| SuppFig14i-p | DONE (was MISSING) | DONE (was MISSING) | `10_completeness/figures/SuppFig14_cell_lines_vs_TCGA_PCA.png` | BLCA pathways raw 64.1/12.5, scaled 22.4/17.1 (mean-expression score, not ssGSEA) |  |
| SuppFig15a | DONE | PARTIAL | `02_organoid_only/COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png` | p=0.583 authors' patients (paper 0.58) | lung set is itself GDSC-trained; not run for pancreas, liver, ovary, sarcoma |
| SuppFig15b | DONE | PARTIAL | `02_organoid_only/COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png` | p=0.506 (paper 0.51) | as SuppFig15a |
| SuppFig15c | DONE | PARTIAL | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_GDSC_trained.png` | 0.82/0.012/0.011 | as SuppFig15a |
| SuppFig15d | DONE | PARTIAL | `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_GDSC_trained.png` | 0.21/0.69/0.24 | as SuppFig15a |
| SuppFig16 | DONE (was PARTIAL) | DONE (was PARTIAL) | `10_completeness/figures/SuppFig16_proximity_all_drugs.png` | 44 paper columns + 8 added drugs |  |
| SuppFig17a | DONE (was PARTIAL) | DONE (was PARTIAL) | `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png` | 60/10/30 median rho 0.37 (LOO 0.44) |  |
| SuppFig17b | DONE (was PARTIAL) | DONE (was PARTIAL) | `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` | 60/10/30 median rho 0.50 (paper 0.89 with IC50) |  |
| SuppFig18 | DONE (was MISSING) |  | `10_completeness/tables/SuppFig18_leave_k_out_summary.tsv` | share p<0.05: k=1 1.00, k=2 0.93, k=3 0.94, k=4 0.87 | cetuximab-specific |
| SuppTable1 | DONE | DONE (was PARTIAL) | `02_organoid_only/tables/SuppTable1_clinical_features_by_predicted_group.tsv` | age 0.260, sex 0.188, T 0.435, N 0.440 (paper 0.26/0.19/0.44/0.44) | GI sets by 02; all headline sets here |
| SuppTable2 | DONE | DONE (was PARTIAL) | `02_organoid_only/BLCA_lee2018/extras/tables/SuppTable2_clinical_balance.tsv` | chi2 p 0.11-1.0 |  |
| SuppTable3 | DONE | PARTIAL (was MISSING) | `02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv` |  | as Fig4b |
| SuppTable4 | IN_PROGRESS_B3 (was PARTIAL) | IN_PROGRESS_B3 (was PARTIAL) | `07_B3_models/ (in progress by B3 agent)` | bootstrap draws exist in 02 tables; formatted table by B3 |  |
| SourceData-STRING | DONE | DONE | `02_organoid_only/BLCA_OV_LUAD_OTHER/_data/moesm3/Figure_2a_f.tsv` | network re-used by all agents |  |
| SourceData-DLhyper | IN_PROGRESS_B3 (was MISSING) | IN_PROGRESS_B3 (was MISSING) | `07_B3_models/ (in progress by B3 agent)` |  |  |
| Text-method-params | DONE | DONE | `02_organoid_only/REPORT_GI.md` | cut-off -1.2816, RidgeCV 0.1-0.9, 10,000 bootstraps |  |
| Text-claim-specificity | DONE (was PARTIAL) | DONE (was PARTIAL) | `10_completeness/tables/SuppFig4_cross_test_summary.json` | own 6/17 p<0.05; off-diagonal 9/245 (3.7%) |  |

## Code and data

Code `code/` (run in order): c01_verify_rank1, c02_schematics, c03_networks, c04_feature_tests,
c05_mutations, c06_pca_clinical_split, c07_isogenic_other, c08_completeness, c09_index, c10_completeness_md. Seed 20261117.
Inputs read from other agents' folders (never modified): B3 bundles, 02 coefficient and summary tables, the paper's STRING
network, the DepMap cache. New downloads in `_data/raw/` (Xena GDC somatic mutations COAD/READ/LUAD; GEO GSE15709, GSE322742,
GPL570 table). Sizes and sha256: `tables/data_versions.tsv`.
