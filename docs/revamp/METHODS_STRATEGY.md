# Revamped methods strategy: label-free hypotheses, causal testing in patients

Status: research and design document. Nothing here has been run. It builds on
`robust/README.md`, `robust/ATLAS.md`, `robust/BRIDGE.md`,
`robust/results/survival_validation/FINAL_BIOMARKER_SURVIVAL.md` and
`docs/PIPELINE_STRATEGY.md`. The earlier survey is in `docs/BRIDGE_strategy_research.md`.

**Recommendation in one line.** Stop *learning* biomarkers from response labels. Instead:
1. Derive a short list of pre-specified pathway hypotheses for each drug, using only
   label-free mechanistic evidence: drug-induced transcription and genetic dependency.
2. Test those hypotheses in patients as **treatment-effect modifiers on survival**, using
   active-comparator, cross-fitted, doubly-robust estimators.
3. Pool the estimates across cancers and cohorts with random-effects meta-analysis.
4. Confirm them in pre-registered randomised or comparator cohorts.

The integrated plan, **ANCHOR**, is in Section 4.

---

## 0. Diagnosis: why the current pipeline stalls

| # | Problem (repo evidence) | Root cause | What the fix has to do |
|---|---|---|---|
| P1 | Supervised TCGA RECIST models reach held-out AUROC of about 0.56. The signal is mostly a shared, prognostic programme (`ATLAS.md` §1, §3). | Labels are few (about 2,800 patient × drug pairs over 21 drugs) and noisy (free-text RECIST, mixed settings). Searching 673 pathways × 21 drugs × 3 layers spends all the power on selection. | Generate hypotheses **without** patient labels, so patient data are used only to *test* a handful of hypotheses. |
| P2 | Signals are prognostic rather than predictive. Score × treatment interactions are rarely significant. | The models learn from treated patients only, so "good outcome under drug" ≈ "good outcome". | Make the **treatment-effect contrast** (treated vs comparator) the estimand, not response in treated patients. |
| P3 | Organoid/PDX → patient transfer fails: pathway effects are uncorrelated (Spearman \|ρ\| ≤ 0.17, `BRIDGE.md`). | Learned *weights* transfer badly: culture lacks immune cells and stroma, IC50 is dominated by proliferation, and there is one mouse per arm. | Transfer **mechanism** (which pathways the drug engages or depends on) rather than fitted weights. Mechanism is a property of the drug and target, so it should be more conserved across systems than a regression vector. |
| P4 | Survival validation: 20 of 41 rows are not testable and only colon oxaliplatin replicated externally (GSE72970). | Per-cancer × drug strata have 2–40 events. Detecting an interaction needs about 4× the events of a main effect (Brookes 2004). | Pool across cancers (meta-analysis), use more endpoints (PFI, TTNT), add cohorts, and **set power targets before testing**. |
| P5 | Selection inflates significance: the paper's best p = 0.076 has a whole-pipeline permutation p of 0.48. | Too many choices are made after seeing outcomes. | Pre-registration, a frozen hypothesis list, and an α budget. |
| P6 | Composition artefacts: bevacizumab → GABA/neuronal pathways; B-cell receptor (BCR) → paclitaxel may be infiltrate. | Bulk pathway scores mix tumour, immune and stroma signal. | Deconvolve. Score pathways in the malignant compartment and adjust for purity and immune/stroma fractions. |

A useful hint from the repo: the **only** externally validated row, colon oxaliplatin in
GSE72970, came from an **active-comparator** contrast (FOLFOX vs FOLFIRI) in which the biomarker
modified the effect of the drug. The strategy below generalises that design.

---

## 1. Evaluation of the candidate ideas

Scores: **Evidence** = published support that the approach yields patient-relevant predictive
biomarkers. **Feasibility** = on 4 CPUs, 15 GB RAM and about 20 GB free disk, with the data
already in the repo or downloadable. **Gain** = expected improvement over the current pipeline in
the number of *replicated, predictive* pathway–drug biomarkers. The ratings are judgements, not
measured values.

| # | Idea | Evidence | Feasibility (4 CPU) | Expected gain | Role in ANCHOR |
|---|---|---|---|---|---|
| 1 | Mechanism-first (L1000/CMap/Tahoe, organoid pre/post) | Medium | High (L1000 landmark subset) | **High** for hypothesis quality; indirect for AUROC | Hypothesis source A |
| 2 | Genetic-dependency anchoring (DepMap/PRISM, SL/SR) | Medium–high for targeted drugs; low–medium for chemo | High | Medium–high | Hypothesis source B |
| 3 | Predictive design / heterogeneous treatment effect (HTE) as primary estimand | High (methods); medium (observational oncology) | High | **High**, because it fixes P2 | Core test engine |
| 4 | Time-on-treatment / TTNT / multi-endpoint | Medium | Medium (TCGA dates are sparse) | Medium: more events, so more power | Endpoint layer |
| 5 | Cross-fitting over all patients + random-effects meta-analysis | High (statistics) | High | **High** for power and replication | Inference layer |
| 6 | Deconvolution / cell-state features | High that confounding exists; medium that removing it rescues biomarkers | Medium (InstaPrism) to high (purity covariates) | Medium; removes false leads | Feature layer |
| 7 | Foundation models / knowledge graphs | Low for bulk survival HTE; KGs are useful as priors | Low (foundation models need a GPU); high (KG lookups) | Low | Optional prior only |
| 8 | Pre-registered staged validation with power targets | High (REMARK, Simon 2009) | High | **High** for credibility | Governance |
| 9 | *(added)* Active-comparator new-user designs | High (pharmacoepi standard) | High | **High** | Primary contrast |
| 10 | *(added)* Randomised multi-arm platforms (I-SPY2, BrighTNess, GSE72970-type) as HTE ground truth | High | High (already curated) | High | Stage-2 validation |
| 11 | *(added)* Negative-control and specificity tests | High | High | Medium (prevents false claims) | Gate criteria |
| 12 | *(added)* More open cohorts with treatment lines (MMRF CoMMpass, POG570, CTR-DB; Hartwig by data request) | Medium | Medium | Medium–high (events) | Data expansion |
| 13 | *(added)* Paired pre/on-treatment patient biopsies (window-of-opportunity, neoadjuvant) | Medium | Medium | Medium; in-patient mechanism | Hypothesis source C |

### 1.1 Idea 1: Mechanism-first biomarkers from drug-induced transcription

**Rationale.** A drug's transcriptional *footprint* shows which pathways it engages. Two kinds
of baseline biomarker follow from this:

- **(a) Dependency / engagement state.** If the drug acts by shutting down pathway P, tumours
  whose growth depends on high baseline P activity should benefit. The template is the MEK
  functional-activation signature: genes *repressed by MEK inhibition* measured at baseline
  predicted selumetinib sensitivity better than RAS/RAF mutations did (Dry 2010).
- **(b) Adaptive-resistance programmes.** Pathways the drug *induces* as a survival response
  (NF-κB, interferon, NRF2, EMT, p53-independent stress responses) indicate pre-existing
  resistance when they are already high at baseline. The repo's credible hits fit this pattern:
  IFN/RIG-I → cisplatin resistance (cf. the IFN-related DNA-damage resistance signature,
  Weichselbaum 2008) and NF-κB → poor outcome on rituximab regimens.

Footprint scores (genes *downstream* of a pathway) measure pathway activity better than
membership scores (genes *in* the pathway), as shown for PROGENy (Schubert 2018) and DoRothEA
(Garcia-Alonso 2019). The repo currently scores Reactome by gene membership.

**Pipeline**

1. **Drug signatures.**
   - LINCS L1000 Level 5 consensus signatures (Subramanian 2017; GEO GSE92742/GSE70138) for each
     repo drug: cisplatin, carboplatin, oxaliplatin, 5-FU, gemcitabine, paclitaxel, docetaxel,
     doxorubicin, epirubicin, etoposide, irinotecan, pemetrexed, temozolomide, sorafenib and
     others. Use landmark genes plus inferred genes, at the dose and time point with the highest
     transcriptional activity score (TAS), per cell line.
   - Tahoe-100M (Zhang 2026, *Cell*): 50 cancer lines × 1,100 drug-dose conditions. Use
     pseudobulk only.
   - Organoid pre/post-treatment RNA-seq where it exists in GEO. Catalogue it with the Seqout
     search tools (search terms: "organoid", drug name, "treated vs vehicle").
2. **Separate mechanism from cytotoxicity.** Regress each signature on the generic
   cell-death/proliferation axis, which dominates L1000 signatures (Szalai 2019). Keep the
   residual "drug-specific" signature. Without this step every drug looks like "cell cycle
   down", which reproduces the shared, prognostic layer from ATLAS.
3. **Pathway-level footprint.** Run decoupler (Badia-i-Mompel 2022) with PROGENy (14
   pathways), DoRothEA/CollecTRI transcription factors, and Reactome/Hallmark sets as a
   secondary layer. Output: a drug × pathway engagement matrix *E*, with sign and a consistency
   score across cell lines and lineages.
4. **Selectivity filter.** A pathway is "engaged by drug d" if (i) \|E\| is in the top 5% for d
   in at least 50% of lines of the relevant lineage and (ii) d is more extreme than its
   mechanism-class siblings. This gives the *specific* engagement needed for predictive
   biomarkers.
5. **Pre-clinical sanity check (no patient data).** Within lineage, does baseline activity of
   the engaged pathway predict GDSC/PRISM/CTRP sensitivity (`data/curated/preclinical`) after
   adjusting for doubling time and PC1 of all drugs (the existing NIT deconfounding)? Keep only
   pathways that pass.
6. **Output:** at most 5 pre-specified pathway hypotheses per drug, with the expected direction.

**How it addresses the current problems.**
- P1/P5: patient labels are not used for selection, so the hypothesis space shrinks from about
  673 pathways to about 5. For 10 drugs that is roughly a 100-fold smaller multiplicity burden.
- P3: only the *identity and sign* of the mechanism is transferred, not regression weights.
- P2: biomarkers of specific engagement are predictive by construction, not prognostic.

**Risks.**
- L1000 is measured in cell lines at 6–24 h. Classic chemotherapies such as 5-FU and
  gemcitabine give weak, mostly p53/DNA-damage signatures.
- Mitigations: use TAS-filtered signatures; use Tahoe-100M for compounds L1000 lacks; accept
  that some drugs will have no specific footprint and fall back to Idea 2.

**Feasibility.**
- Landmark matrix: 978 genes × about 30k signatures for the relevant compounds. Read with h5py
  slicing from the GCTX file, which is about 2 GB for landmark-only data, without loading the
  full file.
- decoupler-py runs in minutes.
- Tahoe-100M: use only pseudobulk summaries. The full h5ad is hundreds of GB and is not
  feasible here; check what the Hugging Face dataset card exposes.

**Plan (3 weeks).**
- Week 1: `obd/mechanism.py`, which fetches L1000 signatures for the 21 ATLAS drugs, removes
  the proliferation/death axis and runs decoupler footprints.
- Week 2: selectivity filter, plus the GDSC/PRISM within-lineage baseline test
  (`run_mechanism_screen.py`).
- Week 3: freeze `results/revamp/hypotheses_v1.tsv` (drug, pathway, direction, evidence).

### 1.2 Idea 2: Genetic-dependency anchoring (DepMap, PRISM, SL/SR networks)

**Rationale.**
- When a drug's effect across cell lines correlates with the CRISPR knockout effect of its
  target or pathway, that pathway is a *functional* mediator. This drug–gene co-dependency
  recovers mechanism of action and target engagement (Gonçalves 2020, *Mol Syst Biol*;
  Tsherniak 2017; Corsello 2020; Pacini 2024).
- Synthetic-lethal / synthetic-rescue (SL/SR) partners mined from tumour data turn a drug target
  into a transcriptomic biomarker without response labels. SELECT predicted response in about
  80% of 35 trials (Lee 2021). ENLIGHT extends this (Dinstag 2023). The repo already holds the
  ISLE and SynLethDB networks (`data/external/gi_networks`) and the ENLIGHT cohorts.

**Pipeline**

1. Inputs: DepMap Chronos gene effects (Dempster 2021) plus PRISM secondary screen and GDSC2
   AUC.
2. For each drug, compute within-lineage partial correlations between drug AUC and gene effect,
   controlling for growth rate and the drug-sensitivity PC1. Aggregate to pathways with a
   rank-based set test. Result: a drug × pathway **co-dependency matrix C**.
3. Find **dependency markers.** For each co-dependent pathway, which baseline expression
   footprint predicts its essentiality in cell lines? Use cross-validated ridge on PROGENy and
   TF activities. This converts a dependency into a feature that can be measured in tumours.
   Translate to tumours with Celligner alignment (Warren 2021), which the BRIDGE alignment code
   already approximates.
4. SL/SR scores:
   - Targeted drugs: SELECT-style. A patient's tumour is predicted to respond when SL partners
     of the target are under-expressed and SR partners are not over-expressed. Restrict to
     edges with tumour-data support.
   - Chemotherapies: use the DNA-repair and replication-stress gene sets as "targets" (for
     example SLFN11 and the Fanconi/HR pathway for platinum and topoisomerase inhibitors).
5. **Triangulate.** A hypothesis is promoted if it is supported by both C (Idea 2) and the
   engagement matrix E (Idea 1), or by one of them plus the literature or knowledge graph.

**How it addresses the current problems.**
- P1: no patient labels are used.
- P3: dependency is measured with genetic, not pharmacological, perturbation. It is
  orthogonal to IC50 noise, and co-dependency is more conserved than fitted weights.
- P2: SL/SR scores are drug-specific by construction.

**Risks.**
- Chemotherapies with DNA targets have weak co-dependency. Some PRISM compounds are not used
  clinically.
- Cell-line lineages are poorly matched to some cancers.
- Commercial ENLIGHT is not reproducible. Re-implement only the public SELECT logic.

**Feasibility.** The DepMap release CSVs are about 1–3 GB, and the correlations are cheap.

**Plan (2 weeks).**
- Week 1: `obd/dependency.py`, computing the co-dependency matrix C and pathway aggregation.
- Week 2: dependency-marker models and SELECT-style scores, then merge into
  `hypotheses_v1.tsv`.

### 1.3 Idea 3: Predictive (not prognostic) design as primary — heterogeneous treatment effects on survival

**Rationale.** The quantity a clinician needs is the *conditional average treatment effect*
(CATE) of drug d versus its alternative, given baseline pathway state. Treatment-interaction
methods estimate this directly:
- modified covariate method (Tian 2014);
- causal forests and generalised random forests (Wager & Athey 2018; Athey 2019);
- R-learner (Nie & Wager 2021) and DR-learner (Kennedy 2023);
- causal survival forests for censored outcomes (Cui 2023);
- frameworks for subgroup identification (Chen 2017; Ballarini 2018).

Target-trial emulation (Hernán & Robins 2016) defines time zero and eligibility so that
observational contrasts mimic a trial, which avoids immortal-time bias (Suissa 2008).

**Pipeline for each (cancer, drug) target trial**

1. **Protocol.**
   - Eligibility: stage and setting, untreated with the drug class before.
   - Arms: drug d vs an **active comparator** (Idea 9) or vs no systemic therapy.
   - Time zero: treatment start, i.e. `days_to_drug_therapy_start` in the BCR biotab. For
     comparators without treatment, use a landmark at the median start day.
   - Outcome: OS, PFI or TTNT (Idea 4).
2. **Nuisance models, cross-fitted.**
   - Propensity e(x) from clinical covariates (age, stage, grade, histology, year, tissue
     source site) plus purity and immune/stroma fractions (Idea 6).
   - Survival outcome models μ₀ and μ₁ (penalised Cox or random survival forest) on clinical
     variables plus the PCs of pathway scores.
   - Censoring model for IPCW.
3. **Effect modification test for each pre-specified pathway (confirmatory).**
   - Doubly-robust pseudo-outcome (restricted-mean survival time difference at τ = 3 years, or
     the log-HR via the modified covariate method).
   - Regress the pseudo-outcome on the standardised pathway score. This gives the best linear
     projection (BLP) slope with a robust SE (Chernozhukov 2018). The slope is the interaction
     estimate for that stratum.
4. **Exploratory layer.** Causal survival forest over all footprint pathways. Report RATE/AUTOC
   (rank-weighted average treatment effect; Yadlowsky 2025) to test whether *any* HTE exists.
   Use this only to generate hypotheses for the next cycle.
5. **Overlap and positivity diagnostics.** Trim to 0.05 < e(x) < 0.95. Report effective sample
   size.

**How it addresses the current problems.**
- P2 directly: the estimand is the *difference* in survival between arms as a function of the
  pathway, so a purely prognostic pathway has an expected slope of 0.
- P1: survival time is a far more abundant label than RECIST. Every treated or comparator
  patient with follow-up contributes, not just the roughly 2,800 with RECIST.

**Risks.**
- **Confounding by indication.** Sicker patients get chemotherapy, or get carboplatin rather
  than cisplatin. Mitigations: active comparators, negative controls (Idea 11) and
  E-value-style sensitivity analysis.
- **Unreliable "untreated" status in TCGA.** A missing drug form does not mean no treatment.
  Restrict untreated comparators to cancers and stages where the clinical "no adjuvant therapy"
  field is explicit, or use active comparators only.
- **Limited power.** See Section 3.

**Feasibility.** The `grf` package (R) or `econml`/`scikit-survival` (Python) on n ≤ 1,000 per
stratum takes minutes on 4 CPUs. The repo already has Cox, interaction and bootstrap code
(`obd/survival.py`, `obd/study.py`).

**Plan (3 weeks).**
- Week 1: `obd/target_trial.py`, building protocol objects (eligibility, arms, time zero,
  endpoints) from the biotab files for about 15 contrasts.
- Week 2: `obd/hte.py`, with cross-fitted DR pseudo-outcomes, BLP slope, RATE and overlap
  diagnostics. Unit tests on simulated data with known HTE.
- Week 3: run the contrasts on hypotheses_v1, using an α plan written before any outcome is
  looked at.

### 1.4 Idea 4: Time-on-treatment, TTNT and pooled multi-endpoint labels

**Rationale.**
- Time to treatment discontinuation (TTD) correlates with PFS/OS across chemotherapy, targeted
  therapy and immunotherapy in NSCLC trials (Blumenthal 2019).
- Real-world genomic cohorts use TTD and TTNT as their main effectiveness endpoints (Singal
  2019; Pleasance 2020 POG570; Jee 2024 MSK-CHORD).
- Pooling endpoints raises the number of events, and events determine power.

**Pipeline**

1. Derive endpoints per patient and regimen:
   - TCGA biotab: treatment start and end, regimen number, `new_tumor_event` days. Gives TTD
     (metastatic/advanced only), TTNT (start of the next regimen or death) and PFI/OS (TCGA-CDR,
     Liu 2018).
   - External: PFS/OS/pCR as available.
2. Multi-state or multi-endpoint model (Putter 2007): treatment start → progression or next
   treatment → death. Alternatively, a joint frailty model sharing the patient's pathway effect
   across endpoints.
3. Simpler and robust option: **stack** (endpoint × patient) rows in the HTE estimator, with
   endpoint-specific baselines and a cluster-robust SE by patient. Report each endpoint
   separately as well.
4. Adjuvant fixed-course regimens: treat TTD as *uninformative*, because it is planned, and use
   PFI/TTNT only.

**How it addresses the current problems.**
- P4: PFI has about 1.5–2× the events of OS in most TCGA cancers.
- TTNT adds events in metastatic cohorts (POG570, Hartwig).

**Risks.**
- TCGA treatment dates are missing for a sizeable fraction of patients.
- TTD reflects toxicity as well as efficacy.
- Heterogeneous censoring.

**Feasibility.** High. This is table processing only.

**Plan (2 weeks).**
- Week 1: `obd/endpoints.py`, which audits date completeness per cancer and drug and builds
  TTD, TTNT and PFI.
- Week 2: multi-endpoint stacking in `hte.py`, plus a simulation check of type-I error with
  cluster-robust SEs.

### 1.5 Idea 5: Cross-fitting across all TCGA patients plus random-effects meta-analysis

**Rationale.** The current split wastes data (FINAL table: held-out n = 2–35 for many rows).
Use everyone for both stages:
- **Cross-fitting** gives out-of-fold nuisance and CATE estimates for every patient
  (Chernozhukov 2018).
- **Pre-validation** (Tibshirani & Efron 2002) and cross-validated survival risk evaluation
  (Simon 2011) give honest performance estimates without a held-out set.
- The cross-validated adaptive signature design (Freidlin 2010) gives honest tests of a
  sensitive subgroup within one dataset.
- Use nested cross-validation whenever anything is tuned (Bates 2024 shows naive CV
  intervals under-cover).

**Pipeline**

1. **Confirmatory hypotheses**, which are pre-specified: no fitting is needed for the biomarker,
   so cross-fitting is used only for nuisance models. Every eligible patient contributes to
   the test.
2. **Exploratory CATE models:**
   - K = 5 cross-fitting repeated 20 times;
   - out-of-fold CATE score for each patient;
   - test with a permutation null that **re-runs the whole fit**. The repo already does this
     for IC50.
3. **Meta-analysis per pathway × drug (or class).**
   - Estimate the interaction log-HR or BLP slope in each cancer, and in each external cohort.
   - Pool with REML random effects (DerSimonian & Laird 1986), with a Hartung–Knapp correction
     and a 95% prediction interval (Riley 2011).
   - Hierarchy: cancer within mechanism class. This maps onto the ATLAS shared / class / drug
     layers, but now applied to *interaction* effects instead of response.
4. Decision rule: the pooled effect is significant at the α budget, I² < 50% or the prediction
   interval excludes the opposite direction, and the sign agrees in at least 70% of strata.

**How it addresses the current problems.**
- P4: pooling 8 cancers with about 20 events each reaches about 160 events, enough for an
  interaction HR of 0.6 at α = 0.05 (Section 3).
- P5: a confirmatory analysis requires only out-of-fold quantities.

**Risks.**
- True heterogeneity across cancers, where the biomarker works in one tissue only. Mitigation:
  report per-cancer estimates and the prediction interval, and allow class-level pooling.
- Pooled strata are not independent when patients receive several drugs in a regimen.
  Analyse at regimen level.

**Plan (1–2 weeks).** `obd/meta.py` (REML + Hartung–Knapp, prediction interval, forest plots),
which feeds into the gates.

### 1.6 Idea 6: Cell-state and deconvolution features

**Rationale.**
- Tumour purity confounds expression–outcome associations (Aran 2015; Yoshihara 2013).
- BayesPrism (Chu 2022; fast implementation InstaPrism, Hu 2024) and CIBERSORTx (Newman 2019)
  give cell fractions and **malignant-cell-specific expression**.
- EcoTyper (Luca 2021) and TME subtypes (Bagaev 2021) give cell states. Pan-cancer malignant
  meta-programs (Gavish 2023) give tumour-intrinsic states.
- The B-cell → paclitaxel signal could reflect B cells or infiltrate generally. A deconvolved
  analysis tells these apart, and both readings are useful but different.

**Pipeline**

1. Fast, always: purity (ESTIMATE, or ABSOLUTE where TCGA provides it) plus immune and stroma
   scores as covariates in every propensity and outcome model.
2. Medium:
   - InstaPrism per cancer with a 3CA / CELLxGENE single-cell reference;
   - malignant-compartment expression gives malignant-only pathway footprints;
   - immune/stromal fractions are modelled as *separate* candidate modifiers (for example
     B-cell fraction × taxane).
3. Gavish meta-program scores replace some Reactome sets in the tumour-intrinsic layer, which
   also cuts redundancy among 673 pathways.
4. Sensitivity analysis: report each interaction raw, purity-adjusted and malignant-only. Claim
   a tumour-intrinsic effect only if it survives the malignant-only analysis.

**How it addresses the current problems.** P6, and partly P3: organoids lack TME, so
comparing malignant-only patient footprints with organoid footprints is a fairer contrast.

**Risks.**
- References are imperfect for rare cancers.
- Deconvolution error adds noise.
- InstaPrism with 20+ cancers × 300–1,100 samples takes hours but is feasible on 4 CPUs.

**Plan (2 weeks).**
- Week 1: purity/ESTIMATE covariates for all cancers (quick).
- Week 2: InstaPrism for the 6 priority cancers (COAD, LUAD, BRCA, BLCA, HNSC, OV) and
  malignant-only footprints.

### 1.7 Idea 7: Foundation-model and knowledge-graph features

**Evidence.**
- Single-cell foundation models (scGPT, Cui 2024; Geneformer, Theodoris 2023) underperform
  simple baselines in zero-shot settings (Kedzierska 2025, *Genome Biol*).
- Deep perturbation predictors do not yet beat linear baselines (Ahlmann-Eltze 2025,
  *Nat Methods*).
- No bulk-tumour foundation model has shown predictive (not prognostic) survival biomarkers.
- Knowledge graphs (Hetionet, Himmelstein 2017; PrimeKG, Chandak 2023; TxGNN, Huang 2024)
  encode drug–target–pathway links. They are cheap and interpretable as *priors* for ranking
  hypotheses.

**Recommendation.** Do not use foundation-model embeddings as features: they need a GPU, are
hard to interpret, and the evidence is negative. Use PrimeKG path counts and distances
(drug → target → pathway) only as a third tie-breaker when ranking Idea 1/2 hypotheses.

**Plan.** At most 3 days, inside the Idea-2 week.

### 1.8 Idea 8: Pre-registered, staged validation with power targets

**Rationale.**
- REMARK (McShane 2005) and TRIPOD+AI (Collins 2024) set reporting standards.
- The prospective–retrospective design (Simon, Paik & Hayes 2009) shows how archived trial
  specimens can give level-I evidence with a pre-specified marker and analysis.
- Predictive and prognostic markers need different tests (Ballman 2015).

**Pipeline.**
1. Freeze `hypotheses_vN.tsv`, the scoring code (git tag plus SHA-256 of the scoring outputs)
   and the analysis plan.
2. Register on OSF *before* any patient outcome is touched.
3. Split α across stages (Section 4). Stage-2 tests are one-sided in the pre-specified
   direction.

**How it addresses the current problems.** P5 directly. It also prevents the "best of 27"
pattern from the original paper.

**Plan.** Ongoing: 2 days to write the analysis plan template
(`docs/revamp/PREREG_TEMPLATE.md`, created later), then 1 day per freeze.

### 1.9 Added ideas

- **(9) Active-comparator, new-user designs.** Compare drug A with a clinically interchangeable
  drug B in the same setting, so confounding by indication and the missing-untreated problem
  largely cancel. Candidate contrasts available in TCGA:

  | Cancer | Contrast |
  |---|---|
  | COAD/READ | FOLFOX/CAPOX vs 5-FU/capecitabine alone |
  | LUAD/LUSC | cisplatin vs carboplatin doublets |
  | BRCA | anthracycline-containing vs non-anthracycline; taxane-containing vs not |
  | HNSC | cisplatin vs cetuximab with RT |
  | GBM/LGG | temozolomide + RT vs RT alone |
  | BLCA | gemcitabine–cisplatin vs MVAC |
  | OV | carboplatin + paclitaxel dose-dense vs 3-weekly (if recorded) |

  Externally: GSE72970 (FOLFOX vs FOLFIRI), GSE41998 (paclitaxel vs ixabepilone), BrighTNess
  and I-SPY2 arms, GSE42127 and JBR.10 (adjuvant vs observation), CGGA (TMZ vs none).
- **(10) Randomised multi-arm platforms.** I-SPY2 (GSE194040; Wolf 2022, *Cancer Cell*,
  10.1016/j.ccell.2022.05.005) has 10 arms with a shared control and expression data, which
  makes it the best public resource for drug-specific HTE. Use it as the gold-standard Stage-2
  set for breast (pCR endpoint, randomised, so the interaction is causal).
- **(11) Negative controls and specificity.** Three checks for each claimed pathway–drug pair:
  - (a) the pathway × *unrelated-class* drug interaction should be about 0;
  - (b) the pathway's association in untreated or comparator patients should be about 0
    (otherwise it is prognostic);
  - (c) a negative-control outcome (for example non-cancer death, where recorded) should show
    no interaction.

  Report the empirical null distribution of interaction z across 1,000 random gene sets of
  matched size. ATLAS already does this for response.
- **(12) Data expansion with treatment lines.**

  | Source | Contents | Access |
  |---|---|---|
  | MMRF CoMMpass | about 1,100 myeloma patients with RNA-seq, line-of-therapy, PFS; bortezomib vs lenalidomide contrasts | GDC; clinical-treatment tier to be verified |
  | POG570 (Pleasance 2020) | metastatic WGTA plus per-drug time to discontinuation | expression downloadable from BC GSC; verify |
  | Hartwig (Priestley 2019) | metastatic WGS, RNA for a subset, treatment and response | data request |
  | CTR-DB (Liu 2022, *NAR*, 10.1093/nar/gkab860; DOI not re-checked) | aggregated clinical transcriptomic response cohorts | public |

- **(13) Paired pre/on-treatment patient biopsies.** Window-of-opportunity and neoadjuvant
  studies with biopsies before and during treatment show drug engagement *in human tumours with
  TME*. Use them as hypothesis source C, a stronger mechanistic anchor than cell lines.
  Catalogue them through Seqout/GEO (query: drug name + "pre-treatment" + "on-treatment" +
  "paired").

---

## 2. Coverage map: which ideas fix which problems

| Problem | Primary fix | Secondary |
|---|---|---|
| P1 noisy, scarce labels | 1, 2, 13 (label-free hypotheses) | 4 (time endpoints) |
| P2 prognostic ≠ predictive | 3 + 9 (HTE with active comparator) | 11 (negative controls), 10 |
| P3 pre-clinical non-transfer | 1, 2 (transfer mechanism, not weights) | 6 (malignant-only scores) |
| P4 untestable / underpowered | 5 (meta-analysis), 4, 12 | 8 (power targets, stop early) |
| P5 selection inflation | 8 (pre-registration), 5 (cross-fitting) | 1/2 (small hypothesis space) |
| P6 composition | 6 | 11 |

---

## 3. Power: what is detectable

For a standardised continuous pathway score × binary treatment, an approximately
Schoenfeld-type calculation (Schoenfeld 1983) gives the number of **events** D needed for 80%
power to detect an interaction HR per SD:

D ≈ (z₁₋α/₂ + z₀.₈)² / [p(1−p) · (log HR_int)²]

where p is the fraction treated with the drug of interest.

| two-sided α | p | HR_int 0.5 | 0.6 | 0.7 | 0.8 |
|---|---|---|---|---|---|
| 0.05 | 0.5 | 65 | 120 | 247 | 631 |
| 0.05 | 0.3 | 78 | 143 | 294 | 751 |
| 0.05 | 0.15 | 128 | 236 | 484 | 1,236 |
| 0.005 | 0.5 | 111 | 204 | 419 | 1,069 |
| 0.005 | 0.3 | 132 | 243 | 498 | 1,273 |

**What this means.**
- No single TCGA cancer × drug stratum in `FINAL_BIOMARKER_SURVIVAL.md` has more than about 70
  events. Single-stratum tests can only detect HR_int ≤ 0.5. That is implausibly large for
  chemotherapy, and it explains the "not testable / not supported" verdicts.
- Pooling 6–10 cancers, using PFI, and adding CoMMpass/POG570/Hartwig can reach 150–400 events
  per drug class. That is enough for HR_int ≈ 0.6–0.7 with a small, pre-specified hypothesis
  set.
- **Power target used in the gates:** at least 150 events pooled per tested hypothesis family
  at Stage 1, and at least 100 events (or n ≥ 150 for a binary endpoint) in each Stage-2
  cohort. A test below target is labelled *descriptive* and cannot pass or fail a gate.

---

## 4. Recommended integrated strategy: **ANCHOR**

**ANCHOR** = **A**nchored mechanism → **N**egative-control-calibrated → **C**ross-fitted
**H**eterogeneous-effect estimation → **O**ut-of-sample meta-analysis → **R**egistered
replication.

The logic is to separate where hypotheses come from (label-free mechanism) from where they are
tested (causal contrasts in patients), and never to let the same data do both.

```
 Stage 0  HYPOTHESES (no patient outcomes)
   L1000/Tahoe footprints (E) ─┐
   DepMap/PRISM co-dependency (C)  ├─► triangulate ─► ≤5 pathways/drug, signed ─► FREEZE + OSF
   SL/SR (SELECT-style), KG, paired biopsies ─┘

 Stage 1  SCREEN IN TCGA + open RWD (all patients, cross-fitted)
   target-trial protocols (active comparator) ─► DR pseudo-outcomes (OS/PFI/TTNT)
   ─► per-stratum interaction slopes ─► random-effects meta per pathway×drug(class)
   ─► negative controls / purity / malignant-only sensitivity

 Stage 2  CONFIRM (pre-registered, one-sided, frozen scores)
   randomised: I-SPY2 arms, BrighTNess, GSE41998, GSE72970, JBR.10/GSE42127, GSE103479
   comparator RWD: CoMMpass, POG570, Hartwig (DAR)

 Stage 3  PROSPECTIVE-RETROSPECTIVE (controlled access trial specimens; EGA plan in docs/data_access)
```

### Phased roadmap (about 12 weeks, 4 CPUs)

| Phase | Weeks | Work | Deliverables |
|---|---|---|---|
| **0. Infrastructure** | 1–2 | `obd/endpoints.py` (TTD/TTNT/PFI audit), `obd/target_trial.py` (about 15 protocols), purity/ESTIMATE covariates; simulation harness with known HTE for type-I error and power | endpoint audit table; simulation report |
| **1. Mechanism hypotheses** | 2–5 | `obd/mechanism.py` (L1000 + Tahoe footprints, death-axis removal, decoupler), `obd/dependency.py` (co-dependency, dependency markers, SELECT-style scores), KG tie-break; pre-clinical within-lineage check | `results/revamp/hypotheses_v1.tsv` frozen and registered |
| **2. Patient screen** | 5–8 | `obd/hte.py` (cross-fitted DR, BLP, RATE, overlap), `obd/meta.py`; InstaPrism for 6 priority cancers; negative controls | per-hypothesis forest plots; Stage-1 verdict table |
| **3. Confirmation** | 8–11 | frozen scores applied to Stage-2 randomised and comparator cohorts; CoMMpass/POG570 ingestion | Stage-2 verdict table (REMARK/TRIPOD+AI checklist) |
| **4. Iterate / extend** | 11–12 | exploratory causal-forest hypotheses feed `hypotheses_v2` (next cycle only); EGA requests for prospective–retrospective sets | v2 registration; decision memo |

### Decision gates (go / no-go)

**Gate G0, end of Phase 0: does the machinery work?**
- GO if both hold on simulated data that mimics TCGA:
  - type-I error of the interaction test is ≤ 0.06 at nominal 0.05;
  - power is ≥ 0.7 for HR_int = 0.6 at 150 events.
- GO also needs at least 6 target-trial protocols with ≥ 50 events each, on any endpoint.
- NO-GO: fix the estimator, or drop TTD and use PFI only.

**Gate G1, end of Phase 1: do the mechanistic hypotheses have pre-clinical support?**
- GO if at least 8 drugs have at least 1 hypothesis that:
  - (a) comes from at least 2 independent label-free sources (E, C, SL/SR, paired biopsy, KG);
  - (b) predicts within-lineage GDSC/PRISM sensitivity after deconfounding, at permutation
    p < 0.05;
  - (c) is drug-specific: the effect for the drug exceeds its mechanism-class siblings.
- NO-GO for a drug: it gets no confirmatory hypothesis and is reported as having no
  mechanistic anchor. This is not a failure.
- Global NO-GO (fewer than 4 drugs pass): mechanism-first is not informative for
  chemotherapy. Pivot the effort to targeted and immunotherapy cohorts (the ENLIGHT set,
  `data/external/trials`), where SELECT has shown signal.

**Gate G2, end of Phase 2: is there a predictive signal in patients?** Per hypothesis family,
with Holm correction over registered hypotheses at family α = 0.05:
- PROMOTE to Stage 2 if all of these hold:
  - pooled interaction is in the pre-specified direction at Holm-adjusted p < 0.05, or nominal
    p < 0.01 with ≥ 150 events;
  - the prediction interval does not cross into a harmful interaction of more than 1.2;
  - the effect survives purity adjustment, with the estimate attenuated by less than 50%;
  - all negative controls are null (\|z\| < 2);
  - the untreated/comparator-arm association is weaker than in the treated arm.
- HOLD: right direction but below the power target. Collect data; no claim.
- DROP: wrong direction with ≥ 150 events, or any negative control fails.
- **Programme-level NO-GO:** if after Stage 1 *none* of about 40 hypotheses is promoted, and the
  RATE/AUTOC for every drug is consistent with 0, conclude that bulk transcriptomic pathway
  state carries no detectable chemotherapy effect modification at public-data scale. Publish
  that as a negative result, and redirect to regimen-level, immunotherapy and targeted-therapy
  questions and to the EGA organoid–patient paired data.

**Gate G3, end of Phase 3: does it replicate?**
- VALIDATED requires:
  - one-sided p < 0.025 for the interaction (randomised) or the comparator-adjusted effect
    (real-world data) in at least one independent cohort that meets its power target;
  - the same direction in all other Stage-2 cohorts, with pooled Stage-2 p < 0.05.
- Any cohort below its power target is descriptive only.
- Validated biomarkers go to the Stage-3 prospective–retrospective request list
  (`docs/data_access/EGA_application_plan.md`).

### What success looks like (expected, not promised)

- Roughly 3–6 replicated pathway × drug(-class) effect modifiers, versus 1 today. The most
  likely are those with prior support:
  - DNA repair / Fanconi → anthracycline and topoisomerase II;
  - IFN/RIG-I → platinum resistance;
  - B-cell / immune → taxane benefit (or its TME counterpart);
  - innate immune / IRAK–TAK1 → oxaliplatin (already replicated);
  - SLFN11 / replication stress → platinum (a mechanistic candidate from DepMap).
- Every claim is labelled with its estimand: interaction HR or ΔRMST, cohort, endpoint and
  events. Prognostic-only pathways are explicitly marked as such.

### What to stop doing

- Training pathway weights on TCGA RECIST labels as the discovery engine. Keep ATLAS for
  description only.
- Transferring learned organoid/PDX weights. Keep organoids for **mechanism**: pre/post
  treatment footprints, and engagement checks for each hypothesis.
- Single-stratum log-rank tests on median splits.

---

## 5. References

DOIs marked † were checked against PubMed during this review. The others are from the
authors' records and should be spot-checked before external use.

**Mechanism and perturbation**
- Lamb J et al. 2006. *Science* 313:1929. The Connectivity Map. doi:10.1126/science.1132939
- Subramanian A et al. 2017. *Cell* 171:1437. LINCS L1000. doi:10.1016/j.cell.2017.10.049
- Szalai B et al. 2019. *Nucleic Acids Res* 47:10010. Cell death/proliferation signatures in perturbation data. doi:10.1093/nar/gkz805 †
- Schubert M et al. 2018. *Nat Commun* 9:20. PROGENy. doi:10.1038/s41467-017-02391-6
- Garcia-Alonso L et al. 2019. *Genome Res* 29:1363. DoRothEA. doi:10.1101/gr.240663.118
- Badia-i-Mompel P et al. 2022. *Bioinform Adv* 2:vbac016. decoupleR. doi:10.1093/bioadv/vbac016
- Dry JR et al. 2010. *Cancer Res* 70:2264. MEK pathway signature predicts selumetinib response. doi:10.1158/0008-5472.CAN-09-1577 †
- Pabon NA et al. 2018. *PLoS Comput Biol* 14:e1006651. Drug–knockdown signature correlation reveals targets. doi:10.1371/journal.pcbi.1006651 †
- Zhang J et al. 2026. *Cell* 189:5945. Tahoe-100M. doi:10.1016/j.cell.2026.08.035 †
- Srivatsan SR et al. 2020. *Science* 367:45. sci-Plex. doi:10.1126/science.aax6234
- McFarland JM et al. 2020. *Nat Commun* 11:4296. MIX-Seq. doi:10.1038/s41467-020-17440-w
- Weichselbaum RR et al. 2008. *PNAS* 105:18490. IFN-related DNA damage resistance signature. doi:10.1073/pnas.0806187105

**Dependency and synthetic lethality**
- Tsherniak A et al. 2017. *Cell* 170:564. DepMap. doi:10.1016/j.cell.2017.06.010
- Corsello SM et al. 2020. *Nat Cancer* 1:235. PRISM repurposing. doi:10.1038/s43018-019-0018-6
- Dempster JM et al. 2021. *Genome Biol* 22:343. Chronos. doi:10.1186/s13059-021-02540-7
- Gonçalves E et al. 2020. *Mol Syst Biol* 16:e9405. Drug MoA from pharmacological + CRISPR screens. doi:10.15252/msb.20199405 †
- Pacini C et al. 2024. *Cancer Cell* 42:301. Clinically informed dependency map. doi:10.1016/j.ccell.2023.12.016 †
- Lee JS et al. 2018. *Nat Commun* 9:2546. ISLE. doi:10.1038/s41467-018-04647-1
- Lee JS et al. 2021. *Cell* 184:2487. SELECT. doi:10.1016/j.cell.2021.03.030
- Dinstag G et al. 2023. *Med* 4:15. ENLIGHT. doi:10.1016/j.medj.2022.11.001
- Warren A et al. 2021. *Nat Commun* 12:22. Celligner. doi:10.1038/s41467-020-20294-x
- Geeleher P et al. 2017. *Genome Res* 27:1743. Imputing drug response in TCGA. doi:10.1101/gr.221077.117 †

**Causal inference and HTE**
- Tian L et al. 2014. *J Am Stat Assoc* 109:1517. Modified covariate method. doi:10.1080/01621459.2014.951443
- Wager S, Athey S. 2018. *J Am Stat Assoc* 113:1228. Causal forests. doi:10.1080/01621459.2017.1319839
- Athey S, Tibshirani J, Wager S. 2019. *Ann Stat* 47:1148. Generalized random forests. doi:10.1214/18-AOS1709
- Nie X, Wager S. 2021. *Biometrika* 108:299. R-learner. doi:10.1093/biomet/asaa076
- Kennedy EH. 2023. *Electron J Stat* 17:3008. DR-learner. doi:10.1214/23-EJS2157
- Cui Y et al. 2023. *J R Stat Soc B* 85:179. Causal survival forests. doi:10.1093/jrsssb/qkac001
- Yadlowsky S et al. 2025. *J Am Stat Assoc*. RATE for treatment prioritisation rules. doi:10.1080/01621459.2024.2393466 (verify)
- Chernozhukov V et al. 2018. *Econom J* 21:C1. Double/debiased ML. doi:10.1111/ectj.12097
- Chen S et al. 2017. *Biometrics* 73:1199. Subgroup identification framework. doi:10.1111/biom.12676
- Ballarini NM et al. 2018. *PLoS One* 13:e0205971. Predicted individual treatment effect. doi:10.1371/journal.pone.0205971 †
- Hernán MA, Robins JM. 2016. *Am J Epidemiol* 183:758. Target trial emulation. doi:10.1093/aje/kwv254
- Suissa S. 2008. *Am J Epidemiol* 167:492. Immortal time bias. doi:10.1093/aje/kwm324

**Endpoints and real-world data**
- Liu J et al. 2018. *Cell* 173:400. TCGA Clinical Data Resource. doi:10.1016/j.cell.2018.02.052
- Blumenthal GM et al. 2019. *Ann Oncol* 30:830. Time to treatment discontinuation in NSCLC trials. doi:10.1093/annonc/mdz060 †
- Singal G et al. 2019. *JAMA* 321:1391. Clinico-genomic real-world database. doi:10.1001/jama.2019.3241
- Pleasance E et al. 2020. *Nat Cancer* 1:452. POG570, therapy–genome interactions. doi:10.1038/s43018-020-0050-6 †
- Jee J et al. 2024. *Nature* 636:728. MSK-CHORD. doi:10.1038/s41586-024-08167-5 †
- Lavery JA et al. 2022. *JCO Clin Cancer Inform* 6:e2100105. GENIE BPC curation. doi:10.1200/CCI.21.00105 †
- Priestley P et al. 2019. *Nature* 575:210. Hartwig metastatic WGS. doi:10.1038/s41586-019-1689-y †
- Putter H et al. 2007. *Stat Med* 26:2389. Multi-state models tutorial. doi:10.1002/sim.2712

**Validation, meta-analysis and power**
- DerSimonian R, Laird N. 1986. *Control Clin Trials* 7:177. doi:10.1016/0197-2456(86)90046-2
- Riley RD et al. 2011. *BMJ* 342:d549. Random-effects interpretation. doi:10.1136/bmj.d549
- Bates S, Hastie T, Tibshirani R. 2024. *J Am Stat Assoc* 119:1434. What CV estimates. doi:10.1080/01621459.2023.2197686
- Tibshirani R, Efron B. 2002. *Stat Appl Genet Mol Biol* 1:1. Pre-validation. doi:10.2202/1544-6115.1000
- Simon RM et al. 2011. *Brief Bioinform* 12:203. Cross-validated survival risk classifiers. doi:10.1093/bib/bbr001
- Freidlin B, Jiang W, Simon R. 2010. *Clin Cancer Res* 16:691. Cross-validated adaptive signature design. doi:10.1158/1078-0432.CCR-09-1357 †
- Brookes ST et al. 2004. *J Clin Epidemiol* 57:229. Power of interaction tests. doi:10.1016/j.jclinepi.2003.08.009 †
- Schoenfeld DA. 1983. *Biometrics* 39:499. Sample size for proportional hazards. doi:10.2307/2531021
- Simon RM, Paik S, Hayes DF. 2009. *J Natl Cancer Inst* 101:1446. Prospective–retrospective designs. doi:10.1093/jnci/djp335
- McShane LM et al. 2005. *J Natl Cancer Inst* 97:1180. REMARK. doi:10.1093/jnci/dji237
- Collins GS et al. 2024. *BMJ* 385:e078378. TRIPOD+AI. doi:10.1136/bmj-2023-078378
- Ballman KV. 2015. *J Clin Oncol* 33:3968. Predictive vs prognostic. doi:10.1200/JCO.2015.63.3651

**Deconvolution and cell state**
- Chu T et al. 2022. *Nat Cancer* 3:505. BayesPrism. doi:10.1038/s43018-022-00356-3
- Hu M, Chikina M. 2024. *Bioinformatics* 40:btae440. InstaPrism. doi:10.1093/bioinformatics/btae440 †
- Newman AM et al. 2019. *Nat Biotechnol* 37:773. CIBERSORTx. doi:10.1038/s41587-019-0114-2
- Yoshihara K et al. 2013. *Nat Commun* 4:2612. ESTIMATE. doi:10.1038/ncomms3612
- Aran D et al. 2015. *Nat Commun* 6:8971. Tumour purity. doi:10.1038/ncomms9971
- Luca BA et al. 2021. *Cell* 184:5482. EcoTyper. doi:10.1016/j.cell.2021.09.014
- Bagaev A et al. 2021. *Cancer Cell* 39:845. TME subtypes. doi:10.1016/j.ccell.2021.04.014
- Gavish A et al. 2023. *Nature* 618:598. Malignant meta-programs (3CA). doi:10.1038/s41586-023-06130-4 †

**Foundation models and knowledge graphs**
- Cui H et al. 2024. *Nat Methods* 21:1470. scGPT. doi:10.1038/s41592-024-02201-0
- Theodoris CV et al. 2023. *Nature* 618:616. Geneformer. doi:10.1038/s41586-023-06139-9
- Kedzierska KZ et al. 2025. *Genome Biol* 26:101. Zero-shot limits of single-cell foundation models. doi:10.1186/s13059-025-03574-x †
- Ahlmann-Eltze C et al. 2025. *Nat Methods* 22:1657. Perturbation DL vs linear baselines. doi:10.1038/s41592-025-02772-6 †
- Himmelstein DS et al. 2017. *eLife* 6:e26726. Hetionet. doi:10.7554/eLife.26726
- Chandak P et al. 2023. *Sci Data* 10:67. PrimeKG. doi:10.1038/s41597-023-01960-3
- Huang K et al. 2024. *Nat Med* 30:3601. TxGNN. doi:10.1038/s41591-024-03233-x

**Repo anchors**
- Kong J et al. 2020. *Nat Commun* 11:5485. doi:10.1038/s41467-020-19313-8
- Wolf DM et al. 2022. *Cancer Cell* 40:609. I-SPY2 predictive biomarkers across 10 therapies. doi:10.1016/j.ccell.2022.05.005
