# Revamped strategy: ANCHOR

This plan combines three reports:
- `DATA_LANDSCAPE.md`: data resources;
- `METHODS_STRATEGY.md`: methods;
- `AUDIT_AND_POWER.md`: root-cause audit and power analysis. Scripts are in `robust/revamp_audit/`.

## Why the current approach cannot reach the goal (audit)

| Problem | Evidence |
|---|---|
| Response labels mostly measure recurrence, not tumour shrinkage | 64% of TCGA labels come from adjuvant or first-line treatment; 90% of "responders" there are CR, i.e. no recurrence. Agreement with 12-month PFI is weak (κ 0.29–0.49). A true AUROC of 0.65 would appear as about 0.575. |
| Drug-specific effects are not identifiable | 66% of labelled records are combination regimens; partner drugs share the label 99% of the time. |
| Too few events | The median held-out row has 9 deaths (smallest detectable HR 0.39). Testing a treatment × biomarker interaction with HR 0.7 needs about 250 events (about 370 with FDR control). No single TCGA cancer × drug pair reaches that. |
| Multiple testing | 34 of 282 p-values are below 0.05, against 14 expected by chance. The best row (colon oxaliplatin) has q = 0.16. Significant HRs from small n are probably inflated 1.6–3.8-fold. |
| Composition confounding | Immune hits track tumour purity (BCR ρ = −0.47, IFN-γ ρ = −0.53). |
| Discovery and validation share patients | Labelled patients are used up in discovery, which leaves untestable rows. |

Realistic targets from the literature: a predictive effect of about HR 0.8–0.9 per SD, and a response AUROC of about 0.60–0.65.

## The revamp: ANCHOR

**Anchored mechanism hypotheses → cross-fitted interaction screen → registered replication.**

Hypotheses come from data without patient outcomes, so every patient remains available for testing.

### Stage 0: mechanism hypotheses (no patient outcomes)
- **Drug-induced transcriptional footprints:** LINCS L1000 and Tahoe-100M, with the generic proliferation/death axis removed. This gives the pathways each drug must engage, and its resistance programmes.
- **Genetic dependency evidence:** DepMap CRISPR × PRISM/GDSC co-dependency, SELECT-style synthetic-lethal/rescue scores, and organoid pre/post-treatment data (existing NIT/BRIDGE sets used as evidence, not as labels).
- **Output:** at most 5 signed pathways per drug, frozen and pre-registered. This is about 100× fewer tests than 673 pathways × drugs. A hypothesis needs support from 2 independent sources and must be specific to its drug compared with related drugs.

### Stage 1: patient screen with all patients and a drug-specific design
- **Units:** drug or regimen pooled across cancers (stratified), not single cancer × drug pairs. Pooled TCGA groups have 150–506 PFI events, enough for interaction tests.
- **Comparisons:**
  - between treatments where possible: FOLFOX vs 5-FU (this design produced the one external success), cisplatin vs carboplatin, taxane vs no taxane;
  - otherwise treated vs untreated, with target-trial emulation, propensity/IPW and cross-fitted doubly-robust estimators.
- **Endpoints:** PFI timed from drug start, OS, and time to next treatment. The RECIST free-text label is dropped as the primary endpoint.
- **Required checks:** tumour purity (ABSOLUTE) and deconvolution adjustment; nulls from random signatures, proliferation and purity; negative-control outcomes; random-effects meta-analysis across cancers.
- **Gate G2:** the pooled interaction points the pre-specified way (Holm-adjusted p < 0.05); the result survives purity adjustment; the negative controls are null; the effect is weaker in comparator patients.

### Stage 2: registered confirmation in randomised or controlled cohorts already available
- **Breast neoadjuvant:** I-SPY2 (10 arms with a shared control), BrighTNess, GSE20271, GSE41998.
- **Lung adjuvant, chemo vs observation:** JBR.10, GSE42127, GSE68465, GSE37745. Pooled, these give 90% power for an interaction.
- **Colorectal:** GSE72970 (FOLFOX vs FOLFIRI), GSE103479 (5-FU vs surgery), FOCUS (GSE156915).
- **Other:** CheckMate-025 (nivolumab vs everolimus), CGGA (temozolomide), ICON7 (GSE140082, bevacizumab), APEX (GSE9782, bortezomib).
- **Gate G3:** one-sided p < 0.025 in at least one adequately powered cohort, and the same direction in the others.

### Stage 3: new data (start applications now)
- Hartwig (metastatic RNA + treatment + response).
- POG570 (open download of RNA + treatment durations; fetch locally).
- Roche EGA trials (IMmotion150/151, POPLAR/OAK).
- Paired organoid–patient sets (Vlachogiannis, Ooft, de Witte; EGA plan in `docs/data_access`).
- CTR-DB 2.0 (10,856 patients, 346 regimens) for breadth.

## What changes compared with now

| Now | ANCHOR |
|---|---|
| Discover on TCGA RECIST labels, then validate on what is left | Discover without patient outcomes; test on all patients |
| 673 pathways × many pairs | At most 5 pre-registered pathways per drug |
| Single cancer × drug pairs (≤ 70 events) | Drug or regimen pooled across cancers (150–500 events) + meta-analysis |
| Response or median-split survival | Treatment × biomarker interaction (predictive) as the primary test |
| Purity ignored | Purity and deconvolution adjustment, plus negative controls |
| Claims for every row | Power-gated: an underpowered test is reported as descriptive, not pass or fail |
| Drug-specific claims for combination partners | Regimen-level claims unless a comparator arm separates the drugs |

## Roadmap (about 12 weeks) and decision gates

| Weeks | Work | Gate |
|---|---|---|
| 1–2 | Infrastructure: PFI from drug start, regimen units, purity, cross-fitting, simulation of type-I error and power | **G0:** type-I error ≤ 0.06; ≥ 6 comparisons with ≥ 50 events |
| 2–5 | Stage 0 hypotheses (L1000, Tahoe, DepMap/PRISM, organoids) | **G1:** ≥ 8 drugs with a two-source hypothesis (if < 4, pivot to targeted/IO) |
| 5–8 | Stage 1 pooled interaction screen | **G2** |
| 8–11 | Stage 2 registered confirmation | **G3** |
| 11–12 | Write-up; start the next cycle with Stage 3 data | – |

If no hypothesis passes G2, the honest outcome is a negative result. That result is still publishable: public transcriptomes cannot identify drug-specific pathway biomarkers at current sample sizes.
