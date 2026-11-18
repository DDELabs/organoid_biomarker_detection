# Mock referee reports (NAR Cancer)

Folder `nar/11_review/`. These are written as two referees would write them for the manuscript outlined in
`STORY_ASSESSMENT.md`. After each concern: **Answered?** says whether an existing result already answers it (with the
path), or what extra work is needed (with the gap ID from `GAP_LIST.tsv`). Paths are relative to `nar/` unless they
start with `docs/` or `robust/`.

---

## Editor's summary (senior editor)

The manuscript re-evaluates a published network-based organoid biomarker framework (Kong et al. 2020). It reproduces
the original results exactly, then tests the framework on 57 lab-model analyses, 196 external tests, cell-line
mechanism data and a pre-registered pan-cancer test. Both referees find the work careful and of value to the field.
Both ask for (i) a cleaner statistical hierarchy (which test is primary), (ii) a direct explanation of why the 2019
result changes on current data, (iii) the model benchmark to be completed, and (iv) the BH3-only finding to be shown
in organoids, not only cell lines. Referee 2 also asks for an outside timestamp for the pre-registration.
**Decision if submitted as is: major revision.** With the gaps G1-G6 closed I would expect minor revision.

---

## Referee 1 — computational cancer biology

**Summary.** The authors take a well-cited organoid ML framework, reproduce it to the decimal on the original
inputs, and then show that the single-pathway organoid-to-patient transfer does not hold across cancers, drugs or
data versions. They link the main colorectal pathway (activation of BH3-only proteins) to BH3-mimetic sensitivity in
cell lines. The study is thorough and the writing is clear. My concerns are mostly about biology and about whether
the comparisons are fair to the original method.

### Major concerns

**R1-M1. Why does the colon 5-FU result change between 2019 and current data?** The authors show p = 0.014 on the
authors' files and p = 0.866 on current GDC data, with a different rank-1 pathway. This could come from the
expression pipeline (HTSeq vs STAR, FPKM-UQ), the patient set (114 vs 129), longer follow-up (8 vs 25 deaths), or
re-selection of the pathway. Please separate these. Most importantly: what happens if the *paper's* BH3-only pathway
is frozen and applied to the current patients?
- **Answered?** Partly. Both endpoints exist (`02_organoid_only/REPORT_GI.md`; `robust/results/paper_single_pathway/summary.tsv`).
  The BH3-only pathway appears in the current-data Ridge k=5 model (`02_organoid_only/COAD_vandewetering/tables/all_tests_COAD_FLUOROURACIL.tsv`),
  but no table gives the frozen BH3-only pathway alone on current patients, and no one-factor-at-a-time swap was
  done. For bladder this was done (amino-acid pathway on current data: p = 0.063; `02_organoid_only/BLCA_lee2018/extras/tables/paper_biomarker_amino_acid_KM.tsv`).
  **Needed: G2.**

**R1-M2. Regimen confounding.** 94 of 114 "5-FU" patients received FOLFOX. Is the BH3-only pathway then an
oxaliplatin marker? The organoids were also screened with oxaliplatin.
- **Answered?** Partly. The split is shown (`02_organoid_only/COAD_vandewetering/figures/Fig2d_regimen_split_5FU_vs_FOLFOX.png`,
  `tables/regimen_split_5FU_FOLFOX.tsv`) and the regimen-level label problem is quantified (`docs/revamp/AUDIT_AND_POWER.md` 1.2).
  GSE39582 subsets were tested (FOLFOX n = 12, HR 0.57, p = 0.40; `05_A3_B1/REPORT.md` B1). Not done: BH3-only score
  vs oxaliplatin response in the 19 organoids and in DepMap. **Needed: G5** (small).

**R1-M3. The BH3-only = apoptotic priming claim needs organoid evidence and a comparison with simpler markers.** The
cell-line correlations with navitoclax and ABT-737 are convincing (rho about -0.20 to -0.25, n = 772-961,
`08_depmap/REPORT.md` e, g). But (a) the paper's claim was about organoids, and (b) the BH3-only Reactome set
contains E2F1 and 14-3-3 genes; is the effect driven by BCL2L1/MCL1 expression, by proliferation, or by the set
itself?
- **Answered?** (b) partly: single-gene results exist (`08_depmap/tables/S6_single_genes.tsv`), mutation and
  lineage adjustment exist (`S5_mutation_adjustment.tsv`, `S3_colon_vs_pancancer.tsv`). Not done: navitoclax in the
  19 van de Wetering organoids (data present: `data/curated/preclinical/colorectal_vandewetering2015/response.tsv`,
  NAVITOCLAX and OBATOCLAX, 19 organoids each); a "does the pathway add beyond BCL2/BCL2L1/MCL1 expression and a
  proliferation score" model; a random 16-gene-set null for the navitoclax association. **Needed: G3, G7.**

**R1-M4. Fairness to the method: would a better model or better alignment work?** The original paper compared
against deep learning, whole transcriptome, etc. The authors should show the same comparison for every set, plus
modern models (elastic net, random forest, XGBoost, GNN) and Celligner-type alignment.
- **Answered?** No — pending. `07_B3_models/` has code and 11 set x drug result files (`_data/bench/`), but the
  summary tables (`07_B3_models/tables/benchmark_model_summary.tsv`) cover only bladder cisplatin. The joint models in
  `03_organoid_plus_patient/REPORT.md` already show that alignment brings lab models and tumours together
  (neighbour share 0.00 -> 0.87 in colorectal) yet organoid-only held-out C-index stays 0.31-0.52. **Needed: G1.**

**R1-M5. Choice of patient cohorts and the comparator.** TCGA "untreated" = "no drug record". Many such patients are
stage I or had treatment not recorded. The prognostic-vs-predictive argument depends on this group.
- **Answered?** Mostly. Adjusted Cox and interactions with age, sex, stage, purity, MSI and regimen are in
  `05_A3_B1/tables/A3_master_results.tsv` and `A3_adjusted_cox_all_terms.tsv` (direction kept in 42 of 49). The
  registered test uses time-varying exposure, IPW and an active comparator (`09_registered_test/tables/stage1_sensitivity.tsv`).
  Add one sentence and a supplementary table of stage distribution in treated vs untreated per cohort. **Needed: G12 (small).**

**R1-M6. Mutation concordance needs a null.** ERCC2 (p = 0.0030) and BRAF V600E (p = 0.031) reproduce, but ERCC2
mutants have high tumour mutation burden and BRAF V600E tumours are mostly MSI/right-sided. How many random
pathways separate these mutants equally well?
- **Answered?** No. Feature-set comparisons exist (`02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv`;
  `10_completeness/figures/Fig5_SuppFig12_SuppFig18_mutation_concordance.png`), but no random-pathway null and no
  adjustment for TMB/MSI/sidedness. **Needed: G8.**

### Minor concerns

- **R1-m1.** Several numbers differ between sections (BRAF current p = 0.00072 in `02_organoid_only/REPORT_GI.md` vs
  0.00048 in `10_completeness/COMPLETENESS.md`; colon 5-FU current p = 0.866 vs 0.62 in `docs/PROJECT_SUMMARY.md`).
  **Needed: G11.**
- **R1-m2.** The lung "organoid" analysis uses GDSC cell lines because the GSE165611 drug table is not public. Label
  it clearly as a cell-line substitute in every figure. **Answered:** labelled in `02_organoid_only/REPORT_BLCA_OV_LUAD_OTHER.md` section 6; carry through to figure titles.
- **R1-m3.** Ovarian organoids were screened with oxaliplatin, not carboplatin; state "platinum proxy" everywhere.
  **Answered:** stated in `REPORT_BLCA_OV_LUAD_OTHER.md` section 2.
- **R1-m4.** Odd frozen pathways ("APOBEC3G-mediated resistance to HIV-1") should be checked for composition. **Answered
  partly:** noted in `08_depmap/REPORT.md` Part B; purity is a covariate in Stage 1. A one-line gene-overlap check is
  enough (G13).
- **R1-m5.** Supp Fig 15 (cell-line-trained models) is missing for pancreas, liver, ovary and sarcoma
  (`10_completeness/COMPLETENESS.md`). Either add or say why not needed. **Needed: G14 (low).**
- **R1-m6.** The olaparib time course (GSE229119) is not a response panel; keep it in the supplement only.
  **Answered:** `REPORT_BLCA_OV_LUAD_OTHER.md` section 5.
- **R1-m7.** Discuss the functional-avatar literature (Vlachogiannis 2018, Ooft 2019, the 2026 PDO/PDX meta-analysis
  with 70% concordance) to make clear the paper evaluates *transcriptomic* biomarkers learned from organoids, not
  organoid response testing itself. **Needed: text only (G15).**

---

## Referee 2 — biostatistics

**Summary.** The authors carry out a large replication with careful multiplicity control, power calculations and a
pre-registered test. I appreciate the honest reporting of every test. My concerns are about the statistical
hierarchy, the independence of the 57 tests, the pre-registration, and how "no evidence" is distinguished from
"evidence of no effect".

### Major concerns

**R2-M1. Define one primary question and one primary family.** There are many families (57 treated log-rank tests,
57 untreated, 57 interactions, 196 external, 7 + 17 registered, 33 joint-model interactions). Readers need to know
which numbers are confirmatory. I suggest: (1) the exact reproduction is descriptive; (2) the registered test is the
only confirmatory analysis; (3) the 57 + 196 analyses are a systematic exploratory survey summarised by a
pre-specified meta-statistic.
- **Answered?** The pieces exist (`05_A3_B1/tables/A3_expected_by_chance.tsv`; `09_registered_test/ANALYSIS_PLAN.md`),
  but no single "statistical hierarchy" table. **Needed: G6** (writing plus one summary table).

**R2-M2. The 57 tests are not independent.** Colon, rectum and colon+rectum share patients; Ridge/SVR/LR share
features. The binomial "7 vs 2.85" comparison overstates precision.
- **Answered?** Acknowledged in `05_A3_B1/REPORT.md` ("the binomial p is a guide"). Not done: a dependence-aware
  null, e.g. permuting drug labels / survival within cohort and recomputing all 57 tests jointly to get the null
  distribution of "number of p < 0.05", or restricting to one analysis per independent patient set.
  **Needed: G9.**

**R2-M3. Absence of evidence vs evidence of absence.** Most analyses are underpowered (median minimal detectable HR
3.9). The authors should not conclude "does not generalise" from underpowered tests alone. Please give equivalence
or Bayesian evidence where power allows, e.g. a random-effects pooled estimate across the 57 analyses in the expected
direction, and replication Bayes factors relative to the original effect.
- **Answered?** Partly. The registered test gives tight intervals (paclitaxel 1.00, 0.89-1.13; gemcitabine 1.12,
  0.94-1.32; `09_registered_test/tables/stage1_primary.tsv`) and the power/winner's-curse tables exist
  (`05_A3_B1/tables/A3_winners_curse.tsv`, `docs/revamp/AUDIT_AND_POWER.md`). Not done: a pooled random-effects
  estimate over the 57 organoid-only analyses (signed log HR with direction) and a "probability of replication
  given the original estimate" calculation. **Needed: G10.**

**R2-M4. Pre-registration must be verifiable.** A hash written by the same pipeline proves the file did not change,
but not *when* it was fixed relative to outcome access. Also, the report states that one Stage 1 run crashed after
printing primary estimates, and the code was then patched.
- **Answered?** The deviation is honestly disclosed (`09_registered_test/REPORT.md`, "Limits and deviations"); the
  sha256 and timestamps are recorded (`08_depmap/tables/frozen_hypotheses.sha256`, `09_registered_test/ANALYSIS_PLAN.md`).
  Git history of the branch can provide an outside timestamp if the frozen file was committed before the Stage 1
  run. **Needed: G4** (check git log; deposit on Zenodo/OSF; reword as "frozen and hashed before outcome analysis").

**R2-M5. Median split and log-rank as the main test.** The original method dichotomises at the median and uses
log-rank. With 8 deaths, and 0 deaths in one arm, the HR is not estimable. Please give continuous-score Cox
estimates as the main effect measure throughout, and treat the median split as the replication of the original
analysis only.
- **Answered?** Yes for A3: continuous-score HR given (`05_A3_B1/REPORT.md`: paper colon 5-FU continuous HR 1.52 per
  SD, p = 0.091; adjusted 1.97, p = 0.15). B1, the registered test and the joint models already use continuous
  per-SD effects. Make this explicit in the Methods. **Text only.**

**R2-M6. Selection of the rank-1 pathway is itself a selection step.** The response-shuffle test includes it
(paper colon 5-FU: rank-1 shuffle share 0.034, best-of-30 share 0.164; bladder permutation p = 0.30). This is good
and should be in the main text, not only the supplement.
- **Answered?** Yes: `02_organoid_only/REPORT_GI.md`, `REPORT_BLCA_OV_LUAD_OTHER.md` 3.1, `05_A3_B1/tables/A3_master_results.tsv`
  (perm_p column). Put it in Fig. 2 or 3. **Layout only.**

**R2-M7. Joint models: is the improvement over patient-only real?** Pooled C-indices differ by 0.01-0.03. Please give
paired bootstrap CIs for the C-index difference (joint minus patient-only) per unit, and try lab weights below 0.25
because inner CV almost always chose the smallest weight.
- **Answered?** No. `03_organoid_plus_patient/tables/joint_weight_chosen.tsv` shows w = 0.25 chosen in 55 of 61 fits;
  the report itself lists "lighter weights were not tried". **Needed: G16.**

### Minor concerns

- **R2-m1.** Report the number of events, not only n, in every KM figure. **Answered:** given in tables
  (`A3_master_results.tsv`); check figure annotations.
- **R2-m2.** OS vs PFI: the organoid-only survey uses OS, the registered test PFI. Give the A3 survey with PFI as a
  sensitivity analysis. **Needed: G17.**
- **R2-m3.** ssGSEA (organoid-only) vs rank score (B1, registered, joint) — show that scores agree (correlation per
  pathway in TCGA) so the switch does not explain differences. **Needed: G18.**
- **R2-m4.** The bootstrap share for colon 5-FU (0.0074 here vs 0.0012 in the paper) should be explained (random
  seed, set size distribution). **Answered partly** in `02_organoid_only/REPORT_GI.md`; one sentence on the likely
  cause (different random sets; the paper's draws are in MOESM3, recounted for bladder to 0.014).
- **R2-m5.** STORM "response" is recurrence-based and PC1 alone gives AUROC 0.96; good that this is flagged
  (`03_organoid_plus_patient/REPORT.md`, `docs/revamp/AUDIT_AND_POWER.md`). Keep the interaction, not the AUROC, as
  the headline.
- **R2-m6.** Data and code availability: a single re-run script and a Zenodo DOI for derived tables are expected.
  **Needed: G19.**

---

## Concern -> status summary

| Concern | Status | Gap |
|---|---|---|
| R1-M1 2019 vs current decomposition | partly answered | G2 |
| R1-M2 FOLFOX / oxaliplatin | partly answered | G5 |
| R1-M3 BH3-only in organoids, beyond BCL2-family genes | partly answered | G3, G7 |
| R1-M4 model benchmark | pending | G1 |
| R1-M5 untreated comparator | mostly answered | G12 |
| R1-M6 mutation-concordance null | not answered | G8 |
| R2-M1 statistical hierarchy | partly | G6 |
| R2-M2 dependence among 57 tests | acknowledged only | G9 |
| R2-M3 evidence of absence / pooled estimate | partly | G10 |
| R2-M4 verifiable pre-registration | partly | G4 |
| R2-M5 continuous Cox as main measure | answered | text |
| R2-M6 shuffle test in main text | answered | layout |
| R2-M7 joint vs patient-only CIs, lower weights | not answered | G16 |
