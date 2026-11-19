# W4: Which pathways drive the breast chemotherapy-response prediction

## What was done (plain words)
- Data: `robust/results/training_table/` (v1, untouched). 3,511 treated breast patients, 13 cohorts, 673 Reactome pathway scores, each z-scored inside its cohort. Same selection and same L2 logistic model (penalty lambda = 100) as W3/part1.
- Every result is "leave one cohort out": train on 12 cohorts, test on the 13th. Every choice (which pathways, which direction, which slope) is made on the 12 training cohorts only. The test cohort is never used. Seeds fixed (20261119).
- Check: the 673-pathway model here gives mean AUC 0.752, the same as `L2_fixed` in W3. Full-data weights equal `robust/results/train_combined/breast_model_pathways.tsv` exactly (max difference 0.0).
- Mean AUC = plain average over the 13 cohorts; its 95% CI uses the summed per-cohort variances (Hanley-McNeil), as in W3.
- Not done: the "MKI67-like gene score from expression" in PLAN W4 (the training table holds pathway scores only, no gene expression; you did not list it in this request).
- Code: `code/` (a_auc_vs_k, b_stability, c_greedy, d_proliferation, e_tables, f_subgroups, g_figures, themes). Run in that order.

## Themes (plain-language groups, by keywords in the pathway name)
Table `tables/family_weights.tsv`, full list `tables/pathway_themes.tsv` (please review; keyword grouping is crude). Figure `fig_themes`.

| theme | pathways | total absolute weight | share of weight | share if equal weights | in top 30 |
|---|---|---|---|---|---|
| other | 264 | 9.07 | 39.0% | 39.2% | 10 |
| signalling | 210 | 7.10 | 30.5% | 31.2% | 10 |
| metabolism | 119 | 4.47 | 19.2% | 17.7% | 5 |
| immune | 53 | 1.51 | 6.5% | 7.9% | 1 |
| proliferation | 27 | 1.11 | 4.8% | 4.0% | 4 |

What I think: the weight is spread across all themes roughly in proportion to how many pathways each has. No theme takes over. Proliferation and metabolism are slightly above their equal share, immune slightly below.

## Q1. One, a few or many pathways? (`tables/auc_vs_k.tsv`, `fig_auc_vs_k`)
Mean AUC (95% CI) for the top k pathways. "Weight" = ranked by training-fold L2 weight; "alone" = ranked by training-fold one-at-a-time strength.

| k | by weight | by one-at-a-time | random sets |
|---|---|---|---|
| 1 | 0.540 | 0.680 | 0.561 |
| 2 | 0.587 | 0.684 | 0.592 |
| 3 | 0.634 | 0.688 | 0.615 |
| 5 | 0.666 | 0.689 | 0.638 |
| 10 | 0.684 | 0.699 | 0.670 |
| 20 | 0.722 | 0.710 | 0.696 |
| 50 | 0.737 | 0.725 | 0.720 |
| 100 | 0.744 | 0.734 | 0.730 |
| 200 | 0.750 | 0.740 | 0.743 |
| 673 | 0.752 (0.726-0.777) | 0.752 | |

- The single best pathway, chosen on training cohorts, was EARLY_PHASE_OF_HIV_LIFE_CYCLE in all 13 folds (AUC 0.680, CI 0.652-0.708). Ranked by model weight the single top pathway was P2Y_RECEPTORS in 12 folds (AUC 0.540).
- Smallest k within 0.01 of the full model (0.742 or higher): k = 100 by weight (0.744); k = 200 for the one-at-a-time ranking (0.740 is just outside, so in practice all 673). Random sets reached it at 200 (0.743).
- What I think: it is "many". One pathway gives 0.68, a handful reach about 0.72, but the last 0.03 needs about 100 or more pathways. Random sets of 100 to 200 pathways do almost as well as chosen ones (0.730 and 0.743), which says the information is spread over many correlated pathways.

## Q2. Weights and stability (`tables/top30_pathways.tsv`, `tables/stability_selection.tsv`, `fig_top30_weights`, `fig_stability`)
- Stability = fraction of 100 bootstrap fits (resampling patients inside each cohort) in which an L1-penalised logistic model (C = 0.1, the value W3 chose by inner cross-validation) kept the pathway. Also 20 bootstraps per leave-one-out fold (`tables/stability_by_fold.tsv`).
- An L1 fit keeps about 234 of 673 pathways. 58 pathways were kept in at least 80% of fits, 184 in at least 50%.
- Of the top 30 by weight, 28 are kept in at least 80% of fits (median 0.91). All 30 keep the same weight sign in all 13 leave-one-out refits.

Top 30 (weight from the L2 model on all patients; OR = raw odds ratio per 1 SD of pathway score, pooled over cohorts, `raw_odds_ratio_per_sd`; stability = L1 selection frequency):

| # | pathway | theme | weight | OR per SD | stability |
|---|---|---|---|---|---|
| 1 | p2y_receptors | signalling | +0.175 | 1.22 | 1.00 |
| 2 | latent_infection_of_homo_sapiens_with_mycobacterium_tuberculosis | other | +0.154 | 1.33 | 0.97 |
| 3 | abortive_elongation_of_hiv1_transcript_in_the_absence_of_tat | other | +0.143 | 1.18 | 0.98 |
| 4 | e2f_enabled_inhibition_of_pre_replication_complex_formation | proliferation | -0.132 | 1.41 | 0.94 |
| 5 | association_of_licensing_factors_with_the_pre_replicative_complex | proliferation | +0.130 | 1.60 | 0.99 |
| 6 | transport_to_the_golgi_and_subsequent_modification | other | +0.119 | 1.23 | 0.89 |
| 7 | abacavir_transport_and_metabolism | metabolism | +0.119 | 1.05 | 0.99 |
| 8 | g_alpha_s_signalling_events | signalling | -0.117 | 0.80 | 0.87 |
| 9 | chondroitin_sulfate_biosynthesis | metabolism | -0.115 | 1.03 | 0.92 |
| 10 | shc_related_events | signalling | -0.115 | 1.04 | 0.90 |
| 11 | insulin_receptor_recycling | signalling | -0.113 | 0.95 | 0.98 |
| 12 | initial_triggering_of_complement | immune | +0.111 | 1.28 | 0.91 |
| 13 | nrage_signals_death_through_jnk | signalling | +0.111 | 1.28 | 0.88 |
| 14 | early_phase_of_hiv_life_cycle | other | +0.110 | 1.84 | 0.96 |
| 15 | mtorc1_mediated_signalling | signalling | +0.109 | 1.03 | 0.88 |
| 16 | mitotic_prometaphase | proliferation | +0.106 | 1.57 | 0.96 |
| 17 | signal_regulatory_protein_sirp_family_interactions | signalling | -0.104 | 1.46 | 0.82 |
| 18 | homologous_recombination_repair_of_replication_independent_double_strand_breaks | other | -0.103 | 1.17 | 0.91 |
| 19 | integrin_cell_surface_interactions | signalling | +0.103 | 1.19 | 0.86 |
| 20 | signal_transduction_by_l1 | signalling | -0.101 | 1.01 | 0.98 |
| 21 | branched_chain_amino_acid_catabolism | metabolism | -0.099 | 0.67 | 0.92 |
| 22 | nucleotide_like_purinergic_receptors | metabolism | -0.098 | 1.21 | 0.91 |
| 23 | other_semaphorin_interactions | signalling | +0.097 | 1.00 | 0.90 |
| 24 | fatty_acyl_coa_biosynthesis | metabolism | +0.094 | 1.11 | 0.77 |
| 25 | g1_s_specific_transcription | proliferation | +0.094 | 1.66 | 0.82 |
| 26 | sema3a_pak_dependent_axon_repulsion | other | -0.093 | 1.10 | 0.87 |
| 27 | abc_family_proteins_mediated_transport | other | +0.093 | 0.76 | 0.85 |
| 28 | proteolytic_cleavage_of_snare_complex_proteins | other | +0.093 | 1.02 | 0.91 |
| 29 | rna_pol_i_transcription_initiation | other | -0.091 | 1.04 | 0.72 |
| 30 | apoptotic_execution_phase | other | -0.090 | 1.13 | 0.91 |

What I think: the top pathways are stable, but they are not an obvious biology story. Several (latent infection with Mycobacterium, HIV transcript elongation, early HIV life cycle, abacavir transport) are virus or drug pathways in the Reactome catalogue whose genes mostly are general cell machinery (for example ribosome, transcription, proteasome genes). I read them as markers of a general high-activity state, not as infection biology. Also, 11 of the top 30 have a weight with the opposite sign to their one-at-a-time odds ratio (for example E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION: weight -0.132 but OR 1.41). That happens when pathways are correlated and the model uses one to correct another, so the weight of a single pathway should not be read alone.

## Q3. Smallest subset (`tables/minimal_subset.tsv`, `tables/minimal_subset_honest_stop.tsv`, `fig_auc_vs_k` b)
Greedy forward selection inside each training set (candidate pool: top 60 by training-only strength; inner AUC from 4 inner folds that keep whole cohorts together; up to 15 steps).

| pathways | 1 | 2 | 3 | 5 | 8 | 10 | 15 |
|---|---|---|---|---|---|---|---|
| mean test AUC | 0.680 | 0.691 | 0.694 | 0.717 | 0.719 | 0.720 | 0.719 |

- Plateau at about 0.72 after 5 pathways, 0.031 to 0.035 below the full model. No step up to 15 reaches within 0.01 of 0.752.
- Honest stopping rule (smallest step within 0.01 of the best inner AUC, chosen in training only): steps chosen 4 to 6, mean AUC 0.718 (0.691-0.744).
- The most often picked first five: EARLY_PHASE_OF_HIV_LIFE_CYCLE (13 of 13 folds), GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATION (12), G1_S_SPECIFIC_TRANSCRIPTION (10), UNWINDING_OF_DNA (9).
- What I think: five pathways give a model that is about 0.03 AUC worse. Within 0.01 of the full model took about 100 pathways in the ranking approach (Q1). The greedy search was limited to a 60-pathway pool and 15 steps, so a larger pool might do slightly better; I did not test that.

## Q4. Is it proliferation? (`tables/proliferation_tests.tsv`, `tables/proliferation_per_cohort.tsv`, `fig_proliferation`)
Proliferation = the 27 pathways whose names contain CELL_CYCLE, MITOTIC, G1_S, G2_M, E2F, DNA_REPLICATION, PRE_REPLICATIVE or CHECKPOINT (`tables/proliferation_pathways.tsv`). Score = their mean; direction set on training cohorts.

| model | mean AUC (95% CI) |
|---|---|
| all 673 pathways | 0.752 (0.726-0.777) |
| proliferation score alone | 0.646 (0.618-0.674) |
| all 673 with proliferation regressed out of every pathway (slopes from training) | 0.721 (0.694-0.747) |
| 673 minus the 27 proliferation pathways | 0.748 (0.723-0.774) |

- Share of weight: 27 proliferation pathways hold 4.8% of the total absolute weight (range over leave-one-out refits 4.3-5.2%). With a broader cell-cycle list (adds CDK, CDC6/20/25, ORC, MCM, kinetochore, spindle, DNA repair and similar, `broad_cellcycle` in `pathway_themes.tsv`) the share is 9.0% (8.2-9.5%). Immune pathways: 6.5% (6.3-6.9%).
- Does the model add to proliferation plus ER and HER2? Likelihood-ratio test per cohort (adjusted for proliferation + ER + HER2 in 9 cohorts, + ER only in 2, proliferation only in BRIGHTNESS and TCGA-BRCA, which have no usable receptor variation): p < 0.05 in 11 of 13 cohorts (not significant in GSE32646, p = 0.137, and TCGA-BRCA, p = 0.068). Combined signed z = 14.0 (Stouffer). For the proliferation-regressed model: 11 of 13, z = 14.1.
- Correlation of model score with proliferation score: mean Spearman 0.41 (0.15 to 0.55 across cohorts); for the regressed-out model 0.05.
- What I think: the model is not mostly proliferation. Proliferation alone gives 0.646; removing it from the model costs 0.004 (dropping the 27 pathways) or 0.031 (regressing it out of everything). So proliferation explains some of the signal but most of the prediction remains once it is removed, and it adds clearly to proliferation plus receptors. The regressing-out loss of 0.031 shows proliferation is a real, shared part of the signal.

## Q5. Subgroups (`tables/subgroup_auc_pooled.tsv`, `tables/subgroup_auc_per_cohort.tsv`, `tables/subgroup_top30_odds_ratios.tsv`, `fig_subgroup_forest`)
Subgroups from ER/HER2 flags; cohort-subgroups with at least 5 responders and 5 non-responders. GSE22093, GSE23988 and TCGA-BRCA have no receptor data and are left out. BRIGHTNESS contains only triple-negative patients.

| subgroup | cohorts | patients | all 673 | top 100 | regressed-out | proliferation only |
|---|---|---|---|---|---|---|
| all patients | 13 | 3,511 | 0.752 | 0.744 | 0.721 | 0.646 |
| ER+/HER2- | 8 | 1,111 | 0.670 (0.599-0.741) | 0.692 | 0.627 | 0.627 |
| HER2+ | 8 | 520 | 0.728 (0.669-0.788) | 0.744 | 0.702 | 0.585 |
| triple-negative | 10 | 1,470 | 0.720 (0.678-0.761) | 0.689 | 0.695 | 0.606 |

- The all-pathway model is above 0.5 in 7 of 8 ER+/HER2- cohorts and in all HER2+ and TNBC cohorts.
- Among the top 30 pathways, the one-at-a-time odds ratio points the same way as overall in 25 of 30 (ER+/HER2-), 24 of 30 (HER2+) and 25 of 30 (TNBC).
- What I think: yes, it holds in all three groups. It is weakest in ER+/HER2- (0.670, wide interval; in this group the proliferation score alone gives 0.627, so more of the signal there is proliferation: regressing it out takes the model down to 0.627). It is strongest in HER2+, where proliferation alone is weak (0.585). Confidence intervals overlap, so I would not claim the groups differ.

## Files
Tables: `tables/` (auc_vs_k, auc_vs_k_per_cohort, top30_pathways, stability_selection, stability_by_fold, l2_weights_per_fold, minimal_subset, minimal_subset_honest_stop, greedy_path_per_fold, proliferation_tests, proliferation_per_cohort, proliferation_pathways, family_weights, pathway_themes, theme_share_per_fold, subgroup_*, loco_scores.tsv.gz). Figures (PNG + PDF + `.caption.md`): `figures/` fig_auc_vs_k, fig_top30_weights, fig_themes, fig_stability, fig_subgroup_forest, fig_proliferation.
Input files: `robust/results/training_table/{patients.tsv.gz, pathway_scores_tcga.tsv.gz, pathway_scores_geo.tsv.gz}` (v1); weights file `robust/results/train_combined/breast_model_pathways.tsv`.
Limits: L1 stability uses one fixed C (0.1). Greedy search used a 60-pathway pool. Top-k refits use lambda = 100 unchanged. Theme grouping is by keywords, not by the official Reactome hierarchy (not available offline in this folder).
