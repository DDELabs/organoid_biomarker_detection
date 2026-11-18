# Story assessment for NAR Cancer (senior editor's view)

Folder `nar/11_review/`. Written 2026-11-18 by the reviewer / senior-editor agent. I only read files; I ran no new
analysis. Every number below is copied from a report or table written by an earlier agent, and the path is given.
Paths are relative to `nar/` unless they start with `docs/` or `robust/`.

Companion files: `REFEREE_REPORTS.md` (two mock referee reports) and `GAP_LIST.tsv` (ranked list of extra work).

---

## 1. Short answer

**Is there a coherent story? Yes.** It is a reproducibility and method-evaluation story, which NAR Cancer does publish
(see section 8). It is not a "new biomarker" story.

**Is it sufficient today? Almost.** The core is complete and unusually thorough: exact reproduction, 57 organoid-only
analyses, 196 external-cohort tests, a power analysis, a cell-line mechanism study and a pre-registered test. Four
things should be finished before submission (details in `GAP_LIST.tsv`):

1. finish the model benchmark (`07_B3_models/`, still running: on 2026-11-18, 11 set x drug result files exist in
   `07_B3_models/_data/bench/`, but the summary tables still cover only bladder cisplatin);
2. a "what changed between 2019 and today" analysis for colon 5-FU (data version, patient set, follow-up, and the
   paper's own pathway frozen, not re-selected);
3. test the BH3-only score against navitoclax **inside the 19 colorectal organoids** (the data are already in
   `data/curated/preclinical/colorectal_vandewetering2015/response.tsv`; navitoclax is screened in all 19);
4. make the pre-registration checkable by outsiders (public timestamped deposit of the frozen file and plan) and
   deposit code + derived data.

What I think: with these four items the paper is a solid NAR Cancer submission. Without item 1 the "other methods do
worse" part of the paper's claim stays unanswered, and a referee will ask for it. Without item 4 the word
"pre-registered" will be challenged.

---

## 2. Proposed title (pick one; the first is my choice)

1. **"Organoid-derived network biomarkers of chemotherapy response: exact reproduction, generalisation across 57
   analyses, and a pre-registered patient test"**
2. "What carries over from organoids to patients? A systematic re-evaluation of network-based drug-response biomarkers"
3. "Reproducing and stress-testing organoid-to-patient drug biomarkers across seven cancer types"

Short running title: "Re-evaluating organoid drug biomarkers".

---

## 3. Central claims the numbers support

Each claim has its key numbers and the file they come from.

**C1. The original headline result reproduces exactly on the authors' own inputs.**
- Colon 5-FU, authors' files: n = 114, 8 deaths (0 vs 8), log-rank p = 0.014 (paper 0.014); untreated n = 298,
  p = 0.161 (paper 0.16). (`02_organoid_only/REPORT_GI.md`, verified in `10_completeness/tables/verify_headline_biomarkers.tsv`)
- Bladder cisplatin, paper's source data: 30 proximal pathways (paper 30), 9 organoids; paper's amino-acid pathway
  treated p = 0.031, HR 2.66 (1.06-6.69) (paper p = 0.01, which used IC50; we used AUC because IC50 is not public).
  (`02_organoid_only/REPORT_BLCA_OV_LUAD_OTHER.md` section 3.1)
- 63 of 77 checklist items are DONE for the paper's own case; the 14 others are being made by the B3 agent.
  (`10_completeness/COMPLETENESS.md`)

**C2. The mechanistic checks of the paper also reproduce.**
- ERCC2-mutant bladder tumours score as more cisplatin-sensitive: one-sided Mann-Whitney p = 0.0030 (paper 0.002);
  p = 0.00055 on current data. (`02_organoid_only/REPORT_BLCA_OV_LUAD_OTHER.md` 3.4)
- Cetuximab model separates BRAF V600E colon tumours: p = 0.031 (paper 0.037); stronger on current data
  (p = 0.00072 in `REPORT_GI.md`; p = 0.00048 in `10_completeness/COMPLETENESS.md` item 9 — two numbers for
  the same test; must be reconciled, see gap G11).
- Isogenic cells: HCT-8 BH3-only p = 0.0085 (paper's NES), 0.039 recomputed; T24 amino-acid pathway p = 0.00022.
- Specificity of the two biomarkers to their own cancer: 9 of 245 off-diagonal cross-tests have p < 0.05 (3.7%,
  about chance). (`10_completeness/COMPLETENESS.md` item 3)

**C3. The single-pathway organoid-to-patient result does not generalise across cancers, drugs and data versions,
and the study design explains why.**
- Same analysis on current GDC data: colon 5-FU n = 129, 25 deaths, p = 0.866 (rank-1 pathway changes); bladder
  cisplatin p = 0.067, HR 2.03 (0.94-4.38). (`02_organoid_only/REPORT_GI.md`; `REPORT_BLCA_OV_LUAD_OTHER.md` 3.2)
- 57 organoid-only analyses (7 cancer types, 11 lab-model sets): 7 have treated p < 0.05 vs 2.85 expected
  (binomial p = 0.023), 4 of 7 in the expected direction, none survives Benjamini-Hochberg (smallest q = 0.063) or
  Holm. Untreated patients give **more** small p-values (12 of 57, binomial p = 0.00002). (`05_A3_B1/REPORT.md`)
- Power: median 17 deaths per analysis, median minimal detectable HR 3.9 at 80% power; only 2 of 57 had 80% power
  for HR 1.5. For the paper's colon 5-FU (8 deaths): power 0.09 at HR 1.5, expected exaggeration 4.2x.
  (`05_A3_B1/REPORT.md`, `docs/revamp/AUDIT_AND_POWER.md`)
- Frozen organoid signatures in 21 external cohort settings: 196 tests, 17 treated-arm p < 0.05 (9.8 expected),
  1 survives BH. The paper's BH3-only pathway in GSE39582: RFS HR 1.01 (0.81-1.27), p = 0.93. (`05_A3_B1/REPORT.md` B1)

**C4. Prognostic vs predictive: most "signal" is general prognosis.**
- 12 of 57 untreated splits vs 7 of 57 treated splits with p < 0.05; score x treatment interaction 5 of 57 (2.85
  expected), smallest BH q = 0.36. (`05_A3_B1/REPORT.md`)
- STORM sorafenib: AUROC 0.947-0.965 in the sorafenib arm and 0.894-0.957 in the placebo arm; interaction OR 1.34,
  p = 0.80. (`03_organoid_plus_patient/REPORT.md`)

**C5. Regimen confounding: the paper's "5-FU" cohort is mostly FOLFOX.**
- 94 of 114 "5-FU" patients also had oxaliplatin; FOLFOX-type subset p = 0.036 (6 deaths, all in predicted
  non-responders); only 11 patients had no oxaliplatin or irinotecan record (0 deaths). (`02_organoid_only/REPORT_GI.md`)
- Across TCGA, a drug's response label equals its partner drug's label in 99% of shared regimens (median);
  oxaliplatin is in a multi-drug regimen in 89% of records. (`docs/revamp/AUDIT_AND_POWER.md` 1.2)

**C6. The BH3-only pathway is a real readout of apoptotic priming, not of 5-FU response.**
- Pan-cancer, lineage-adjusted: navitoclax rho -0.206 (GDSC1, n = 926), -0.202 (GDSC2, n = 961), -0.211 (CTRPv2,
  n = 772); ABT-737 -0.246 (GDSC2) and -0.238 (CTRPv2), q < 0.001. 5-FU ranks 203/378, 204/286, 409/460, 1179/1442.
  (`08_depmap/REPORT.md` e, g)
- In colorectal lines the BH3-only score does not track 5-FU (GDSC1 rho 0.054, n = 44, p = 0.73; GDSC2 0.016;
  CTRPv2 0.257; PRISM 0.405, opposite direction); the frozen organoid model sits inside the random-signature band
  in all 4 screens. (`08_depmap/REPORT.md` b, c)
- One isogenic pair supports the organoid direction after 5-FU exposure (HCT-8, 24 h p = 0.007, 48 h p = 0.019).

**C7. A pre-registered test of cell-line-derived hypotheses gives informative estimates.**
- 35 frozen pathways for 7 drugs (sha256 a0bdd620...). Stage 1 (TCGA, pan-cancer, time-varying exposure, PFI):
  0 of 7 drugs pass; best gemcitabine interaction HR 1.12 (0.94-1.32), p = 0.199; paclitaxel 1.00 (0.89-1.13).
  The intervals exclude interaction HRs above about 1.13-1.32 per SD for paclitaxel, gemcitabine and 5-FU.
- Stage 2 (17 controlled external tests): 1 one-sided p < 0.025 (oxaliplatin GSE72970, HR 1.63, 1.00-2.64,
  Holm 0.386). (`09_registered_test/REPORT.md`)

**C8. Joint organoid + patient models: patients carry most of the information; organoids act as a light prior.**
- Organoid-only held-out C-index 0.31-0.52; patient-only usually best; joint pooled HR < 1 with CI excluding 1 in 5
  of 11 units (e.g. CRC_5FU 0.89, 0.80-1.00, C = 0.532 vs patient-only 0.546). Inner CV chose the smallest lab weight
  (0.25) in 55 of 61 fits. 3 of 33 interaction tests p < 0.05 (1.7 expected). (`03_organoid_plus_patient/REPORT.md`)

---

## 4. Which results carry the paper (my ranking)

| Rank | Result | Why it carries weight | Strength |
|---|---|---|---|
| 1 | 57-analysis generalisation + power + prognostic-vs-predictive (C3, C4) | The main new knowledge. Large, systematic, with multiplicity control, untreated control, interaction and power. Few papers do this. | Strong |
| 2 | Exact reproduction + regimen confounding (C1, C5) | Shows the original analysis is computationally correct, then shows what it actually measured (FOLFOX patients, 8 deaths). Fair to the original authors. | Strong |
| 3 | BH3-only = apoptotic priming (C6) | The one clearly positive biological finding, replicated in 3 independent screens with ~800-960 lines each. | Strong in cell lines; untested in organoids (gap G3) |
| 4 | Pre-registered test (C7) | Methodologically the most modern part. Gives confidence intervals that rule out moderate effects for 3 drugs. | Strong if the pre-registration is made public (gap G4) |
| 5 | ERCC2 / BRAF / isogenic reproduce (C2) | Shows the organoid models do learn real biology. Balances the paper. | Moderate (need specificity controls, gap G8) |
| 6 | Joint organoid + patient models (C8) | Useful and novel framing, but effects are small. | Moderate |
| 7 | External cohorts B1 (part of C3) | Large, but heterogeneous endpoints; 1 BH survivor without a comparator arm. | Supporting |
| 8 | Model benchmark (`07_B3_models/`) | Needed to answer "maybe a better model would work". | Pending |

What I think: the paper should lead with C1 + C3 + C4 (the "reproduces, then does not generalise, and here is why"
arc), and give C6 and C7 as the constructive half (what the organoid signal really is, and how to test such
hypotheses properly). C8 and the benchmark go at the end as "can better modelling rescue it?".

---

## 5. Abstract draft (about 250 words)

> Patient-derived organoids are widely used to discover drug-response biomarkers, but few such biomarkers have been
> tested outside the study that found them. We re-evaluated a published network-based machine-learning framework
> that links organoid drug response to patient survival through drug-proximal Reactome pathways. Using the authors'
> inputs, the colorectal 5-fluorouracil result reproduced exactly (n = 114, log-rank P = 0.014), as did the bladder
> cisplatin pathway (P = 0.031) and its ERCC2-mutation concordance (P = 0.003). We then applied the same pipeline to
> 57 lab-model–drug–cohort analyses across seven cancer types (11 organoid, xenograft and cell-line sets) and to 196
> tests in 21 external cohort settings. Seven of 57 analyses reached P < 0.05 (2.85 expected) and none passed
> false-discovery control; the same splits separated untreated patients more often (12 of 57), indicating prognostic
> rather than drug-specific signal. The median analysis had 17 deaths and could detect only hazard ratios above 3.9.
> In TCGA, 94 of 114 "5-fluorouracil" patients received FOLFOX. In about 800–960 cancer cell lines from three
> independent screens, the organoid biomarker "activation of BH3-only proteins" tracked sensitivity to the
> BCL-2/BCL-XL inhibitors navitoclax and ABT-737 (ρ ≈ −0.20 to −0.25, q < 0.001) but not 5-fluorouracil, identifying
> it as a marker of apoptotic priming. A pre-registered, pan-cancer test of 35 cell-line-derived pathway hypotheses
> for seven drugs found no drug-specific interaction (e.g. paclitaxel HR 1.00, 95% CI 0.89–1.13). Models that learn
> from organoids and patients jointly performed close to patient-only models. Organoid-derived pathway biomarkers
> capture real biology, but testing their clinical value requires drug-specific comparator arms and several hundred
> events; we provide code, frozen hypotheses and power tables to support such tests.

(About 260 words. Numbers from C1-C8. "BH3-only line" is the user's call: it says what the data show without
claiming a clinical marker.)

---

## 6. Main figures (6 + 1 optional) with exact source paths

Paths below are existing files. Several main figures must be **assembled** from panels of existing figures (a
layout job, no new analysis), except where "NEW" is marked (these need the gap analyses).

**Figure 1. Study design and data map.**
- a: `10_completeness/figures/Fig1a_workflow_schematic.png` (the method)
- b: `10_completeness/figures/Fig1b_network_proximity_schematic.png`
- c: `10_completeness/figures/Fig1c_data_map_all_sets.png` (all 57 analyses)
- d: `09_registered_test/figures/fig1_preregistration_flow.png` (the pre-registered arm)

**Figure 2. Exact reproduction, and what the original cohort contained.**
- a: colon 5-FU, authors' inputs: `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_authors_FLUOROURACIL.png`
- b: bladder cisplatin, paper source data: `02_organoid_only/BLCA_lee2018/paper/figures/BLCA_CISPLATIN_full.png` (KM panels) and `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_paper_biomarker_KM.png`
- c: FOLFOX split: `02_organoid_only/COAD_vandewetering/figures/Fig2d_regimen_split_5FU_vs_FOLFOX.png`
- d: same pipeline on current data: `02_organoid_only/COAD_vandewetering/figures/Fig2a_Fig2d_Fig4e_SuppFig17_panels_COAD_FLUOROURACIL.png`
- e (NEW, gap G2): decomposition of the 2019 -> current change (expression pipeline / patients / follow-up / frozen pathway).

**Figure 3. Generalisation across 57 analyses, with power.**
- a: `05_A3_B1/figures/A3_Fig1_master_forest.png`
- b: `05_A3_B1/figures/A3_Fig2_QQ_expected_by_chance.png` (treated vs untreated)
- c: `05_A3_B1/figures/A3_Fig4_interaction_forest.png`
- d: `05_A3_B1/figures/A3_Fig5_power_vs_events.png`

**Figure 4. Frozen signatures in independent cohorts.**
- a: `05_A3_B1/figures/B1_Fig7_heatmap_signature_by_cohort.png`
- b: `05_A3_B1/figures/B1_Fig1_forest_colorectal.png`
- c: `05_A3_B1/figures/B1_Fig5_KM_main_cohorts.png` (GSE39582 panel)
- d: `05_A3_B1/figures/B1_Fig6_STORM_sorafenib.png` (prognostic, not predictive)

**Figure 5. What the organoid biomarkers measure: mechanism and mutation concordance.**
- a: `08_depmap/figures/fig1_workflow_bh3_vs_5fu.png` (panel b: BH3-only vs 5-FU in colorectal lines)
- b: `08_depmap/figures/fig3_drug_specificity_and_bladder.png` (panel e: BH3-only vs every drug)
- c: `08_depmap/figures/fig4_bh3_mimetics_and_crispr.png` (panel g: BH3 mimetics)
- d (NEW, gap G3): BH3-only score vs navitoclax in the 19 colorectal organoids
- e: `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_ERCC2_concordance.png`
- f: `02_organoid_only/COAD_vandewetering/figures/Fig5a_cetuximab_BRAF_V600E.png`
- g: `08_depmap/figures/S7_isogenic_pairs.png` (HCT-8 after 5-FU)

**Figure 6. Pre-registered test of cell-line-derived hypotheses.**
- a: `08_depmap/figures/fig5_anchor_stage0_frozen_hypotheses.png`
- b: `09_registered_test/figures/fig2_stage1_pooled_interaction_forest.png`
- c: `09_registered_test/figures/fig4_random_signature_null.png`
- d: `09_registered_test/figures/fig7_external_cohort_forest.png`
- e: `09_registered_test/figures/fig8_summary_heatmap.png`

**Figure 7 (optional, or fold into Fig. 6). Can better models rescue transfer?**
- a: `03_organoid_plus_patient/figures/alignment_colorectal.png` (alignment works)
- b: `03_organoid_plus_patient/figures/forest_cross_cancer.png` (organoid-only vs patient-only vs joint)
- c: `03_organoid_plus_patient/figures/main_colorectal.png` (learning-curve panel)
- d (PENDING): model benchmark across all sets from `07_B3_models/figures/` (only `fig3_replica/Fig3_BLCA_paper_CISPLATIN.png` exists now)

What I think: six main figures is the right size for NAR Cancer. If the editor wants fewer, merge Fig. 4 into Fig. 3
and put Fig. 7 into the supplement.

---

## 7. Supplementary material

**Supplementary figures (suggested order).**
- S1 Batch / PCA, organoids vs TCGA: `02_organoid_only/COAD_vandewetering/figures/SuppFig1_PCA_models_vs_TCGA.png`, `10_completeness/figures/SuppFig14_cell_lines_vs_TCGA_PCA.png`
- S2 Network neighbourhoods: `10_completeness/figures/Fig2c_Fig2g_network_neighbourhoods_all_sets.png`
- S3 Coefficients across Ridge/SVR/LR: `02_organoid_only/COAD_vandewetering/figures/Fig2b_SuppFig2_SuppFig3_coefficients_COAD_authors_FLUOROURACIL.png`
- S4 Cross-test specificity matrix: `10_completeness/figures/SuppFig4_cross_test_matrix_all_sets.png`
- S5 Multi-pathway models: `02_organoid_only/COAD_vandewetering/figures/SuppFig5_SourceData-multipathway_KMgrid_COAD_authors_FLUOROURACIL.png`, `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_multipathway_vs_paper.png`
- S6-S7 Centrality and DL replicas: `07_B3_models/figures/SuppFig6_replica_COAD_paper_FLUOROURACIL.png`, `07_B3_models/figures/SuppFig7_network_centrality_profile.png`, `07_B3_models/figures/SuppFig7_replica_BLCA_paper_CISPLATIN.png`
- S8-S9 Target genes and 1/2/3-step neighbours: `10_completeness/figures/SuppFig8_SuppFig9_targets_and_neighbourhoods_all_sets.png`
- S10 Pathway size vs p: `10_completeness/figures/SuppFig10_size_vs_p_all_sets.png`
- S11 Comparison feature sets, bladder: `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_comparison_feature_sets.png`
- S12 Isogenic pairs, other cancers: `10_completeness/figures/Fig4bc_isogenic_other_cancers.png`, `02_organoid_only/COAD_vandewetering/figures/Fig4b_SuppTable3_isogenic_HCT8_5FU.png`, `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_isogenic_T24.png`
- S13 Mutation concordance, all sets: `10_completeness/figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png`
- S14 Bootstraps: `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_paper_bootstrap_recount.png`
- S15 Cell-line-trained models: `02_organoid_only/COAD_vandewetering/figures/SuppFig15_GDSC_cell_line_models_5FU.png`, `02_organoid_only/BLCA_lee2018/extras/figures/S_BLCA_GDSC_trained.png`
- S16 Proximity heatmap: `10_completeness/figures/SuppFig16_proximity_all_drugs.png`
- S17 60/10/30 inner validation: `10_completeness/figures/SuppFig17_split_603010_all_sets.png`
- S18 Cut-off sensitivity: `02_organoid_only/figures/cross_set_cutoff_sensitivity.png`
- S19 Per-set forests: `02_organoid_only/figures/cross_set_forest.png`, `02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/figures/cross_set_forest.png`
- S20 Adjusted vs unadjusted Cox: `05_A3_B1/figures/A3_Fig3_adjusted_vs_unadjusted.png`
- S21 External forests, bladder / ovarian / liver+lung: `05_A3_B1/figures/B1_Fig2_forest_bladder.png`, `B1_Fig3_forest_ovarian.png`, `B1_Fig4_forest_liver_and_lung.png`
- S22 DepMap controls: `08_depmap/figures/S1_data_matching_flow.png` to `S8_metric_and_geneset_robustness.png`, and `08_depmap/figures/fig2_frozen_model_and_sign_agreement.png`
- S23 Registered test extras: `09_registered_test/figures/fig3_km_high_vs_low_score.png`, `fig5_negative_controls.png`, `fig6_per_cancer_forest_meta.png`
- S24 Joint models per cancer: `03_organoid_plus_patient/figures/main_{bladder,pancreas,liver,ovarian,lung,sarcoma}.png`, `alignment_*.png`
- S25 Olaparib time course and lung PDC substitution (user's own attempts): `02_organoid_only/OV_GSE229119_olaparib/figures/OV_GSE229119_olaparib_timecourse.png`, `02_organoid_only/LUAD_GSE165611_substitute/figures/GSE165611_PDC_vs_TCGA_PCA.png`

**Supplementary tables.**
- T1 Checklist of all 77 paper items with status: `10_completeness/completeness.tsv` (+ `00_checklist/checklist.tsv`)
- T2 Master table of 57 analyses: `05_A3_B1/tables/A3_master_results.tsv`
- T3 External tests: `05_A3_B1/tables/B1_external_results.tsv`, `B1_frozen_signatures.tsv`
- T4 Clinical balance: `02_organoid_only/tables/SuppTable1_clinical_features_by_predicted_group.tsv`, `02_organoid_only/BLCA_lee2018/extras/tables/SuppTable2_clinical_balance.tsv`
- T5 Isogenic: `02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv`
- T6 DL hyper-parameters: `07_B3_models/tables/SuppTable4_DL_grid_all_pairs.tsv` (pending complete)
- T7 Frozen hypotheses + sha256: `08_depmap/tables/frozen_hypotheses.tsv`, `frozen_hypotheses.sha256`; plan `09_registered_test/ANALYSIS_PLAN.md`
- T8 Registered results: `09_registered_test/tables/stage1_primary.tsv`, `stage1_sensitivity.tsv`, `stage1_power.tsv`, `stage2_controlled.tsv`
- T9 Joint models: `03_organoid_plus_patient/tables/summary_main.tsv`, `interaction_all.tsv`, `alignment_metrics.tsv`
- T10 Power targets and label audit: `robust/revamp_audit/out/power_targets.tsv`, `regimen_confounding.tsv`
- T11 Data versions (every agent's `tables/data_versions.tsv`; merge into one).

---

## 8. Literature context (searched PubMed 2020-2026, 2026-11-18)

How I searched: PubMed MCP keyword searches. PubMed has no "cited by" function, so I could **not** list papers that
cite Kong et al. 2020 (PMID 33127883). The "similar articles" list for it is mostly organoid-biobank and
drug-screen papers, not re-evaluations. A Google Scholar or Web of Science citation search should be done by the
user before submission (gap G15).

What exists (all found in PubMed in this session):
- Kong et al. 2020, Nat Commun, PMID 33127883, [10.1038/s41467-020-19313-8](https://doi.org/10.1038/s41467-020-19313-8). A later review from the same group: PMID 36284440, BMB Rep 2023, [10.5483/BMBRep.2022-0155](https://doi.org/10.5483/BMBRep.2022-0155).
- Organoid-as-avatar evidence (functional, not transcriptomic, biomarkers): Vlachogiannis 2018 Science, PMID 29472484, [10.1126/science.aao2774](https://doi.org/10.1126/science.aao2774); Ooft 2019 Sci Transl Med, PMID 31597751, [10.1126/scitranslmed.aay2574](https://doi.org/10.1126/scitranslmed.aay2574) (prediction specific to irinotecan, not 5-FU/oxaliplatin); Tiriac 2018 Cancer Discov, PMID 29853643, [10.1158/2159-8290.CD-18-0349](https://doi.org/10.1158/2159-8290.CD-18-0349); Wensink 2021 npj Precis Oncol review of 17 PDO studies, PMID 33846504, [10.1038/s41698-021-00168-1](https://doi.org/10.1038/s41698-021-00168-1); 2024 mCRC PDO screen optimisation, PMID 38414064, [10.1186/s13046-024-02980-6](https://doi.org/10.1186/s13046-024-02980-6).
- 2026 meta-analyses: PDO vs PDX avatars, 411 patient-model pairs, 70% concordance, PMID 42497819, [10.1016/j.ctrv.2026.103190](https://doi.org/10.1016/j.ctrv.2026.103190); PDTO meta-analysis of 54 studies, PMID 43041741, [10.1002/btm2.70188](https://doi.org/10.1002/btm2.70188).
- Lab-to-patient transfer methods: Celligner, PMID 33397959, [10.1038/s41467-020-20294-x](https://doi.org/10.1038/s41467-020-20294-x); few-shot transfer (TCRP), PMID 34223192, [10.1038/s43018-020-00169-2](https://doi.org/10.1038/s43018-020-00169-2); CellHit web server (NAR 2025), PMID 40377071, [10.1093/nar/gkaf414](https://doi.org/10.1093/nar/gkaf414).
- Statistics context: random signatures associate with outcome (Venet 2011), PMID 22028643, [10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240); type S / type M errors (Gelman & Carlin 2014), PMID 26186114, [10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642); Reproducibility Project: Cancer Biology, PMID 34874005, [10.7554/eLife.71601](https://doi.org/10.7554/eLife.71601) and PMID 34874008, [10.7554/eLife.67995](https://doi.org/10.7554/eLife.67995).
- NAR Cancer scope (from PubMed records of the journal): it publishes computational drug-response papers, e.g. NeurixAI (PMID 40918644, [10.1093/narcan/zcaf029](https://doi.org/10.1093/narcan/zcaf029)), DVFGCDR (PMID 41268024, [10.1093/narcan/zcaf043](https://doi.org/10.1093/narcan/zcaf043)), alt-EJ transcriptional predictor of treatment (PMID 40061566, [10.1093/narcan/zcaf007](https://doi.org/10.1093/narcan/zcaf007)), and organoid papers such as bladder chemoresistance m6A (PMID 38023731, [10.1093/narcan/zcad054](https://doi.org/10.1093/narcan/zcad054)). I could not open the journal's author-guidelines page from here (oup.com not reachable); the user should check the current article types (standard article vs "Survey and summary") before choosing.

What I did **not** find: a published, systematic re-evaluation of an organoid-derived *transcriptomic* drug biomarker
across many cancers with a pre-registered patient test. The closest work evaluates *functional* organoid response
(the meta-analyses above) or proposes new transfer methods (Celligner, TCRP). What I think: this gap is the paper's
novelty claim, and it is defensible, but the citation search in G15 must confirm it.

---

## 9. Honest weak points (for the user to decide on)

- **The paper rests on many small analyses.** That is the point (power), but a referee may say "of course it does
  not replicate with 8-25 deaths". The answer is already in the data: the paper's own claim was made with 8 deaths,
  and the pre-registered test (C7) has 121-405 events on drug for 5-FU, gemcitabine and paclitaxel. Put this
  contrast in the abstract (done above).
- **TCGA "untreated" is not a control arm.** Stated in every report. The randomised comparisons (STORM, GSE20271,
  GSE41998, I-SPY2, GSE15622) are few and small. Say so plainly.
- **Several numbers differ slightly between agents** (BRAF current p 0.00072 vs 0.00048; colon 5-FU current
  p 0.866 in nar/02 vs 0.62 in `docs/PROJECT_SUMMARY.md`; Broutier 0.139 vs 0.089). These must be reconciled to one
  source before writing (gap G11).
- **Bladder rank-1 vs rank-2.** The paper's bladder biomarker was rank 2 in its own model; our rank-1 reproduces at
  p = 0.044 and the paper's pathway at p = 0.031, but the permutation p is 0.30. This should be described neutrally:
  "the reported pathway was chosen from the top 2 because it was consistent across algorithms".
- **The "pre-registration" is internal.** The plan and the frozen file have hashes and timestamps written by the
  project itself. Partial outside evidence exists: git commit `cf0e585` (2026-11-17 13:05:11 UTC, pushed to
  `origin/claude/zealous-darwin-zeiopw`) contains `frozen_hypotheses.tsv` and its sha256, 17 minutes before the
  registered-test commit `44d7246` (13:22:00 UTC). The analysis plan, however, was first committed together with
  the results, so it has only an internal timestamp. Honest wording: "hypotheses frozen, hashed and committed
  before any outcome was opened; analysis plan written before analysis (internal timestamp); both deposited at
  submission". See G4.
- **The B3 benchmark is unfinished**, so the claim "no model family rescues transfer" is not yet supported.

What I think: none of these stops the paper. All are fixable or are honest limits that make the paper more useful.
