# Audit and power analysis of the organoid → patient biomarker programme

Scope: the robust / NIT / TRIAD / ATLAS / BRIDGE layers and the final per-cancer survival table
(`robust/results/survival_validation/FINAL_BIOMARKER_SURVIVAL.*`). Existing code and results were not
changed. Every number below comes from the scripts in `robust/revamp_audit/`, which write their outputs to
`robust/revamp_audit/out/`. The exception is the literature anchors in §3.

| Script | Question | Main outputs |
|---|---|---|
| `a01_labels_regimens.py` | How noisy are the TCGA response labels, and how often do several drugs share one regimen label? | `label_vs_pfi_*.tsv`, `regimen_confounding.tsv`, `regimen_identical_pairs.tsv`, `label_summary.json` |
| `a02_events_power.py` | How many events exist per cancer × drug, how much power do they give, and which rows are testable? | `power_targets.tsv`, `tcga_events_by_pair.tsv`, `tcga_events_by_drug_pooled.tsv`, `final_rows_testability.tsv`, `followup_by_cancer.tsv`, `power_simulation.tsv` |
| `a03_confounding.py` | Purity, prognostic confounding (LGG), multiple testing, list instability, LOCO noise, external cohorts dominated by one axis | `purity_*.tsv`, `lgg_tmz_confounding.tsv`, `multiple_testing.tsv`, `loco_noise.tsv`, `external_dominant_axis.tsv`, `confounding_summary.json` |
| `a04_external_power.py` | How much power do the 77 curated external cohorts and arms give? | `external_cohort_power.tsv` |
| `a05_calibration.py` | Literature effect sizes, and how much the significant TCGA hits are likely exaggerated (Gelman–Carlin) | `literature_effect_sizes.tsv`, `retrodesign_final_rows.tsv` |
| `a06_class_pooling.py` | How much does pooling by mechanism class within a cancer help? | `tcga_events_by_cancer_class.tsv` |

Tumour purity comes from the TCGA PanCanAtlas ABSOLUTE calls (GDC file `4f277128-f793-4354-a13d-30cc7fe9f6b5`,
saved as `robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt`, sha256 `f430a975…`).
Purity is available for 8,866 of the 10,203 TCGA patients with expression.

Re-run: `cd robust/revamp_audit && for s in a0*.py; do python3 $s; done` (about 1 minute).

---

## 0. Bottom line

1. **The project failed mainly because the data cannot answer its questions. The modelling is a secondary cause.**
   - Of the 43 final survival rows, 31 have less than 50% power, 7 have 50–80% power and 1 (GBM temozolomide)
     has at least 80% power to detect a held-out main effect of HR 0.7 per SD.
   - The median held-out row has **9 OS events (11 PFI events)**. 62 are needed at α = 0.05 and 92 at an FDR-level α.
   - **No row has 80% power for the treatment × biomarker interaction.** The median power for that interaction is 22%.
     Detecting an interaction HR of 0.7 needs about **247 events**, or about 368 with FDR control.
2. **The TCGA "RECIST" label mostly measures recurrence (a prognostic outcome), and it is attached to the regimen rather than to a single drug.**
   - 64% of labels come from the early (adjuvant / first-line) setting. There, 90% of "responders" are "CR", which
     usually means "no recurrence".
   - 66% of labelled records belong to multi-drug regimens. A drug's label is identical to its partner drug's label
     in 99% of shared regimens (median).
   - So the label cannot carry drug-specific information. A model trained on it will learn a *shared, prognostic*
     programme, which is exactly what ATLAS reported.
3. **The weak AUROC (0.56) is real, but the per-pair spread is noise.** Across the 64 LOCO pairs, the observed SD of
   AUROC is 0.164. Sampling alone predicts 0.144, and the estimated between-pair heterogeneity is τ ≈ 0.04
   (Cochran p = 0.26). Pair-specific "wins" cannot be told apart from noise.
4. **Every apparently positive survival finding is explained by multiple testing, effect-size inflation
   (type M error) or confounding.**
   - The final table contains 282 p-values. 34 are nominally significant, against about 14 expected if every null were true.
   - After BH correction, only the LGG temozolomide signal survives with q < 0.01, and it did not replicate in CGGA.
   - Detected with 8–22 events, the significant HRs (0.25–0.53) are 1.6–3.8-fold exaggerations of any plausible
     true effect.
5. **Realistic targets:**
   - Expect a per-SD HR of about 0.75–0.85 for prognostic effects. For predictive effects, expect an interaction ratio
     of about 0.8–0.9 per SD.
   - Expect a response AUROC of about 0.60–0.70 at best (§3).
   - At those sizes, single cancer × drug cells in TCGA are untestable by construction. Only **drug-level pan-cancer
     pooling** (300–500 PFI events per drug) and a few **pooled randomised external cohorts** can reach adequate power.

---

## 1. Root-cause analysis

### 1.1 Failure → dominant cause (with evidence)

| Failure | Dominant cause | Contributing causes | Key evidence (this audit) |
|---|---|---|---|
| **Weak held-out AUROC (ATLAS 0.56, BRIDGE 0.50–0.54)** | **Label construction.** The TCGA label is a regimen-level recurrence surrogate, not drug-specific tumour shrinkage. | Small n per pair (median 27.5); class imbalance; noise in the pair-level AUROC | 90% of early-setting responders are CR; PR (69% PFI-free at 12 m) is indistinguishable from SD (71%), while CR is 91% and PD is 22% (`label_vs_pfi_by_category.tsv`). Label vs 12-month PFI from drug start: κ = 0.49 (early setting), κ = 0.29 (post-progression, using OS). 14 of 52 pairs with n ≥ 20 have a responder rate above 85% or below 15%. LOCO τ ≈ 0.04 (`loco_noise.tsv`). |
| **No drug specificity (ATLAS shared layer dominates; randomised add-on designs null)** | **Regimen confounding.** Combination drugs share labels and patients. | Trial designs also combine drugs (AC backbone in BrighTNess / I-SPY2) | 66% of labelled records come from multi-drug regimens (71% weighted over the top 21 drugs). Median identity between a drug's label and its partner's label is 0.99. 14 cancer × drug pairs have patient-set Jaccard ≥ 0.6, including BLCA methotrexate = vinblastine (J = 1.0; identical top-10 lists in `TOP10_RELAXED.md`), TGCT cisplatin/etoposide (0.96), SARC docetaxel/gemcitabine (0.81), UCEC/UCS carboplatin/paclitaxel (0.75), and COAD 5-FU/oxaliplatin (0.67). |
| **Organoid → patient non-transfer (NIT, BRIDGE)** | **Unmeasurable on the patient side.** The patient test has 15–50 deaths per drug; the patient label is prognostic; organoids lack stroma and immune cells. | Platform / culture shift; organoid n 14–50 | Patient-side MDHR (§2) is 0.37–0.65 for the NIT studies. Even a perfect organoid model with true HR 0.8 would have had <25% power in every NIT row. BRIDGE pathway concordance is \|ρ\| ≤ 0.17. The ATLAS "responder" hits are largely composition axes: BCR ρ(purity) = −0.47, RIG-I/IFN −0.32, prostanoid −0.32, semaphorin −0.41, CTLA4 −0.42, IFN-γ −0.53 (median within-cancer Spearman; `purity_pathway_rho.tsv`). These are signals organoids cannot express. |
| **Survival non-replication (FINAL table: 1 "validated", 2 "TCGA only")** | **Small n (event-starved), plus multiple testing and winner's curse.** | Prognostic confounding; non-randomised "untreated" comparator; short follow-up; platform shift (CGGA, microarrays) | Median 9 held-out OS events. 34 of 282 p-values are below 0.05 (14 expected under the null); BH q < 0.1 for 7 tests, 4 of them from one LGG row. COAD oxaliplatin, the only "validated externally" row, has external interaction p = 0.010, which becomes q = 0.16 across the 139 primary tests, and its FOLFOX arm has only 27 PFS events. Retrodesign: HNSC cisplatin (22 events, HR 0.53) has 18% power if the true HR is 0.8, and a significant estimate is expected to be inflated 2.4× in log scale. LUAD pemetrexed (8 events, HR 0.25): 10% power, 3.8× inflation. |
| **"Contradicted" rows (BRCA doxorubicin HR 1.69, LUAD cisplatin HR 2.91)** | **Noise with few events, plus multiplicity.** | Prognostic confounding (purity) | 18 and 16 events. About 2 wrong-sign "significant" rows are expected among roughly 80 held-out tests under the null. LUAD carboplatin PFI HR moves from 1.60 (p = 0.04) to 1.18 (p = 0.57) after purity adjustment; ρ(score, purity) = 0.45 (`purity_biomarker_rows.tsv`). |
| **Untestable rows (20 of 43)** | **Discovery consumed the labelled patients**, and the remaining held-out patients have fewer than 8 events. | Short follow-up in several cancers (median potential follow-up 22–35 months: STAD, ESCA, PAAD, READ, COAD) | 15 rows had to fall back on "all-treated", which includes the discovery patients. Held-out events for these rows are 1–10 (`final_rows_testability.tsv`). |
| **LGG temozolomide "TCGA-only" HR 0.50** | **A prognostic signal within TCGA, not a demonstrated predictive one.** | Platform shift to CGGA; untreated comparator with 23 events | Across all LGG the score is prognostic (OS HR 0.61, p = 1e-7). The held-out HR survives adjustment for grade and histology (0.56, p = 0.001), and the score is unrelated to grade (AUC 0.49). The untreated HR of 0.99 has 23 events (95% CI 0.63–1.55) and cannot rule out the same effect; interaction power is 29% (OS). TMZ-treated patients are 70% grade 3 vs 29% of untreated patients, so treatment is strongly confounded by indication. In CGGA, HR is 1.03 among TMZ-treated grade II/III patients. IDH / 1p19q status, the obvious candidate confounder, is not in the CDR. |
| **External "biomarkers" with extreme z (TACE GSE104580 z ≈ 5.6; STORM)** | **Cohort-level structure (batch or composition).** | — | In GSE104580, 29% of all 673 pathways individually reach AUROC > 0.65 for response, and PC1 and PC2 each give AUROC 0.68–0.69. With n = 147 the null SD of AUROC is about 0.05, so that fraction should be about 0.1%. The response signal is diffuse, so no specific pathway can be credited. GSE32646 behaves the same way (24%). STORM GSE109211: PC1 AUROC 0.96, already excluded by the repo's QC (`external_dominant_axis.tsv`). |

### 1.2 Quantified diagnostics

**Label agreement (TCGA, earliest label per patient × drug; n = 3,371 pairs, 1,615 patients).**

| Setting (share of labels) | Endpoint measured from drug start | n with known 12-m status | κ | PPV (responder → event-free) | NPV (non-responder → event) | C-index (label vs time) |
|---|---|---|---|---|---|---|
| Early (64%) | PFI-free at 12 m | 1,843 | **0.49** | 0.88 | 0.59 | 0.72 |
| Post-progression (26%) | Alive at 12 m | 732 | **0.29** | 0.92 | 0.55 | 0.63 |
| All | PFI-free at 12 m | 2,125 | 0.51 | 0.87 | 0.63 | 0.72 |

12-month PFI-free rate by category (all settings): CR 0.91, PR 0.63, SD 0.66, PD 0.20.
The CR/PR vs SD/PD split therefore cuts *between* CR and PR, not where the biology changes.

The agreement with the CDR first-course outcome (κ = 0.82) is high only because both fields come from the same
case-report source; the label is no more valid for it. The same applies to the Ding 2016 overlap (99.3% agreement)
quoted in `data/external/tcga_response/SUMMARY.md`: it shows the labels are reproducible, not that they are valid.

**Regimen confounding (top 21 drugs).** Multi-drug fraction by drug:

| Drug | Multi-drug fraction |
|---|---|
| Cyclophosphamide | 0.96 |
| Epirubicin | 0.98 |
| Etoposide | 0.93 |
| Oxaliplatin | 0.89 |
| Pemetrexed | 0.88 |
| Carboplatin | 0.82 |
| Doxorubicin | 0.80 |
| Paclitaxel | 0.68 |
| Cisplatin | 0.68 |
| Temozolomide | 0.05 |

Temozolomide is essentially the only monotherapy label with real n, and it is 90% LGG
(`regimen_confounding.tsv`).

**Purity / composition.**
- 20.5% of Reactome pathways have a median within-cancer \|ρ\| with ABSOLUTE purity above 0.3, and 5.2% above 0.5.
- In 18 of 39 scorable final rows, the biomarker score has \|ρ\| > 0.3 with purity.
  Examples: GBM TMZ −0.70; LUSC cisplatin/carboplatin 0.62; BLCA gemcitabine 0.55; DLBC 0.55.
- Adjusting held-out Cox models for purity changes log-HR by a median of only 0.05. Purity therefore confounds the
  *interpretation* (immune and stromal pathways are proxies for composition) more than the survival numbers.
  Exceptions:
  - GBM TMZ: OS HR moves from 0.91 to 0.76.
  - LUAD carboplatin: PFI HR moves from 1.60 to 1.18.
- Purity itself does not predict the TCGA label: median per-pair AUROC is 0.53, and the pooled within-cancer AUROC
  is 0.48.

**Discovery stability.**

| | Standard top-10 lists | Relaxed lists |
|---|---|---|
| Median bootstrap top-10 frequency | 0.03 | 0.08 |
| Entries with local q < 0.1 | 0 | 1.3% |

In 90% of standard entries, |ATLAS prior z| exceeds |local z|. The per-cancer lists are mostly the pan-cancer
prior re-labelled, not local discoveries.

**Discovery/validation overlap.**
- 15 of 43 rows report "all-treated" estimates that include the discovery patients.
- The GBM row was discovered in CGGA / GSE7696, and those cohorts appear again on the external side; the script flags
  them as `DISCOVERY_COHORTS`.
- ATLAS TMZ effects are 90% LGG-derived. The LGG "pan-cancer prior" is therefore the same cohort's RECIST patients.
  It is not a leak into the held-out patients, but it gives no independent support.

**Follow-up.** Median potential follow-up (reverse KM, OS) is 22–35 months for STAD, ESCA, PAAD, READ, COAD, LIHC,
LGG, LUAD and CESC. Long-tail endpoints such as LGG OS (24% event rate) and BRCA (14%) accrue events slowly.
PFI gives roughly 1.2–2× the events of OS in most cancers and should be the default endpoint.

---

## 2. Power analysis

Formulas:
- **Main effect.** Schoenfeld/Hsieh: D = (z₁₋α/₂ + z₁₋β)² / (log HR)², for a per-SD effect.
- **Interaction.** Var(β_int) ≈ 1/D_treated + 1/D_control.
- **Binary endpoints.** Binormal equal variance, δ = √2·Φ⁻¹(AUROC), Var ≈ 1/n₁ + 1/n₀.

**Check by simulation** (statsmodels PHReg, 300 replicates; `power_simulation.tsv`):

| Events | Simulated power | Formula power |
|---|---|---|
| 237 | 0.74 | 0.78 |
| 473 | 0.97 | 0.97 |

The formula is slightly optimistic, so treat the requirements below as lower bounds.

FDR column: α* = 0.0099 (Jung 2005; 40 rows, 4 true, q = 0.10, power 0.8). Bonferroni (0.05/40) requires about 1.4× more.

| Target (80% power) | α = 0.05 | FDR α* = 0.0099 | Bonferroni 0.00125 |
|---|---|---|---|
| Main effect HR 0.7 / SD (events) | **62** | **92** | 130 |
| Main effect HR 0.7 / SD, marker R² = 0.3 with covariates | 88 | 131 | 186 |
| Main effect HR 0.8 / SD | 158 | 235 | 333 |
| **Interaction HR 0.7 / SD, events 50:50** | **247** | **368** | 521 |
| Interaction HR 0.7, events 20:80 | 385 | 575 | 813 |
| Interaction HR 0.8 / SD, 50:50 | 630 | 941 | 1,331 |
| Binary response AUROC 0.65, prevalence 0.6 (patients) | 110 | 164 | — |
| Binary response AUROC 0.60, prevalence 0.6 | 255 | — | — |
| Binary response AUROC 0.56 (ATLAS level), prevalence 0.6 | **717** | — | — |

### 2.1 TCGA: what is theoretically testable?

There are 131 cancer × drug pairs with at least 10 treated patients that have expression (`tcga_events_by_pair.tsv`).

| Analysis unit | Main effect HR 0.7 testable (α 0.05 / FDR) | Interaction HR 0.7 testable | RECIST AUROC 0.65 testable |
|---|---|---|---|
| Single cancer × drug (131 pairs) | **11 / 7** | **0** | **0** |
| Cancer × mechanism class (59 cells) | 12 | 0 | — |
| Drug, pooled across cancers (stratified) | 9 drugs at power ≥ 0.85 | **8 drugs** at power ≥ 0.95 (carboplatin, paclitaxel, cisplatin, gemcitabine, doxorubicin, 5-FU, docetaxel, bevacizumab) | cisplatin (391 labels), 5-FU (265), carboplatin (232), paclitaxel (226) reach ≥ 0.8 for AUROC 0.65 pooled |

The 11 cancer × drug pairs where a main effect is testable:

| Cancer | Drug | Treated | OS events | PFI events | PFI events in untreated | Interaction power |
|---|---|---|---|---|---|---|
| OV | carboplatin | 347 | 217 | 266 | 19 | 0.32 |
| OV | paclitaxel | 347 | 214 | 265 | 19 | 0.32 |
| LGG | temozolomide | 256 | 90 | 131 | 47 | 0.55 |
| OV | doxorubicin | 115 | 96 | 112 | 19 | 0.30 |
| OV | cisplatin | 123 | 76 | 101 | 19 | 0.30 |
| OV | gemcitabine | 100 | 82 | 99 | 19 | 0.30 |
| OV | topotecan | 97 | 93 | 97 | 19 | 0.30 |
| GBM | temozolomide | 100 | 72 | 79 | 33 | 0.41 |
| UCEC | paclitaxel | 163 | 43 | 71 | 45 | 0.46 |
| UCEC | carboplatin | 163 | 43 | 70 | 45 | 0.46 |
| PAAD | gemcitabine | 104 | 51 | 69 | 30 | 0.37 |

Most of these are ovarian recurrent-disease drugs (doxorubicin, gemcitabine, topotecan). Those drugs are given after
progression, so they need a landmark at drug start, and they have no untreated comparator (OV has 19 untreated PFI
events).

Drug-level pooling (`tcga_events_by_drug_pooled.tsv`), with an MDHR of 0.82–0.88 per SD:

| Drug | PFI events, treated | PFI events, untreated |
|---|---|---|
| Carboplatin | 506 | 612 |
| Paclitaxel | 475 | 589 |
| Cisplatin | 368 | 670 |
| Gemcitabine | 296 | 441 |
| Doxorubicin | 233 | 311 |
| Temozolomide | 210 | 80 |

**This is the only TCGA unit with adequate power for interactions.** The "untreated" arm is not randomised, so a
TCGA interaction remains confounded by indication. Treat it as a screening statistic.

**FINAL table rows** (`final_rows_testability.tsv`):
- Held-out OS events have a median of 9 and a maximum of 72 (GBM).
- The median minimum detectable HR (MDHR) is **0.39 per SD**, smaller than any credible transcriptomic effect (§3).
- The minimum detectable *interaction* HR has a median of about 0.4.

| Row | Held-out events | Best main-effect power at HR 0.7 |
|---|---|---|
| GBM TMZ | 72 OS | 0.86 (the only adequate row) |
| OV carboplatin | 42 OS | 0.64 |
| OV paclitaxel | 39 OS | 0.61 |
| LGG TMZ | 60 PFI | 0.79 |
| UCEC carboplatin / paclitaxel | 47–49 PFI | 0.69–0.70 |
| BRCA doxorubicin | 31 PFI | 0.51 |
| BRCA cyclophosphamide | 40 PFI | 0.62 |

Every other row has less than 50% power.

### 2.2 External cohorts (`external_cohort_power.tsv`, 77 cohorts / arm sets)

| Capability (80% power) | Count | Cohorts |
|---|---|---|
| Response AUROC 0.65 | 14 | I-SPY2 (987), BrighTNess (482), GSE25066 (488), GSE41998, GSE20194, GSE9891, GSE169455, GSE72970, IMvigor210, IMmotion150, BRAUN/CheckMate, Liu2019, GSE104580\*, GSE109211\* |
| Response AUROC 0.60 | 3 | the three largest breast cohorts |
| Main effect HR 0.7 (time-to-event) | 12 | GSE72970, GSE9891, GSE63885, GSE42127, CGGA693, BRAUN, GSE25066 (DRFS), IMmotion150, IMvigor210, Liu2019, GSE68465, GSE37745 |
| **Single-cohort interaction HR 0.7** | **2** | CGGA693 (TMZ vs none, 296/83 deaths, non-randomised, 0.82); CheckMate-025 (nivolumab vs everolimus, 159/117, 0.83) |
| Binary interaction (AUROC 0.65 in one arm vs 0.5 in the other) | 0 | best are BrighTNess 0.65 and I-SPY2 0.55 |

\* QC-flagged or diffuse-signal cohorts; see §1.1.

**Pooled meta-analyses that reach adequate interaction power** (inverse-variance pooled information of the per-cohort
interaction terms):

| Pool | Interaction power, α 0.05 | Interaction power, FDR α* | Note |
|---|---|---|---|
| **Lung adjuvant platinum vs observation** (JBR.10 30/30, GSE42127 55/9, GSE68465 161/64, GSE37745 51/22) | **0.90** | 0.75 | JBR.10 is randomised; the others are not |
| Ovarian GSE9891 (platinum + taxane vs platinum, 141/41) + GSE63885 (36/30) | 0.70 | 0.46 | |
| GBM CGGA WHO IV + GSE7696 (RT + TMZ vs RT, 35/24) | 0.27 | — | already used for discovery |

---

## 3. Realistic effect sizes (calibration)

Sources were retrieved from PubMed (`literature_effect_sizes.tsv`). TAILORx is cited from memory and has not been
re-checked here.

| Marker | Setting | Reported effect | What it implies here |
|---|---|---|---|
| 21-gene RS (Oncotype DX) – NSABP B-20 | ER+ N0 breast, tamoxifen ± chemo, RCT, n = 651 | Interaction p = **0.038**. Chemo RR 0.26 if RS ≥ 31 vs 1.31 if RS < 18 ([doi](https://doi.org/10.1200/JCO.2005.04.7985)) | The best-known transcriptomic *predictive* marker, with an extreme contrast between categories, only just reached significance with 651 patients. |
| 21-gene RS – TAILORx | RS 11–25, about 6,700 randomised | Chemo adds nothing overall; benefit is limited to age ≤ 50, RS 16–25 ([doi](https://doi.org/10.1056/NEJMoa1804710)) | In the intermediate range the predictive effect is near zero. |
| 70-gene MammaPrint – MINDACT | Clinical-high / genomic-low, n = 1,550 | 5-yr DMFS difference of 1.5 percentage points with vs without chemo ([doi](https://doi.org/10.1056/NEJMoa1602253)) | The marker is mainly prognostic; its predictive component is small. |
| Stromal TILs (continuous) | Early TNBC, pooled IPD n = 2,148, 736 iDFS events | iDFS HR **0.87 per 10% TILs** (SD 20%), i.e. about **0.76 per SD** ([doi](https://doi.org/10.1200/JCO.18.01010)) | A strong, well-replicated *prognostic* immune effect, measured with more than 700 events. |
| MGMT methylation | GBM, RT vs RT + TMZ, n = 206 | Prognostic HR **0.45**; TMZ benefit mainly in methylated tumours ([doi](https://doi.org/10.1056/NEJMoa043331)) | A large binary prognostic effect; the interaction is modest even for the canonical predictive marker. |
| T-cell-inflamed GEP + TMB | Pembrolizumab, 22 tumour types, > 300 patients | Independently predictive; per-tumour AUROCs about 0.6–0.75 ([doi](https://doi.org/10.1126/science.aar3593)) | The ceiling for transcriptomic response prediction is about AUROC 0.65–0.75, and that is for IO, where biology and expression align. |
| ERCC1 for adjuvant cisplatin | NSCLC, discovery IALT, validation JBR.10 / CALGB 9633 (n = 494) | Predictive interaction **failed to validate** ([doi](https://doi.org/10.1056/NEJMoa1214271)) | Plausible DNA-repair predictive markers often do not replicate. |
| Random signatures | Breast outcome | More than 90% of random signatures with > 100 genes are prognostic; the effect is driven by proliferation ([doi](https://doi.org/10.1371/journal.pcbi.1002240)) | Any broad pathway score will look prognostic unless benchmarked against random and proliferation nulls. |

**Calibrated targets for the revamp:**

| Quantity | Credible range | Implausible in this setting |
|---|---|---|
| Prognostic per-SD HR | 0.75–0.85 | — |
| Predictive interaction ratio (continuous score) | 0.80–0.90 per SD | HR ≤ 0.6 per SD, unless it is a single-gene, target-defined marker |
| Response AUROC (chemotherapy) | 0.58–0.65 | ≥ 0.7, unless driven by subtype (ER) or cohort structure |

At an interaction HR of 0.85, about **1,190 events** (50:50 split) are needed at α 0.05. No public dataset in the
repo offers that for a single drug in a single cancer.

---

## 4. Recommendations

### 4.1 Answerable now (with the available data)

| Question | Design | Power (target) |
|---|---|---|
| Is there a shared, chemotherapy-modulated prognostic programme? | Pan-cancer TCGA, stratified Cox on PFI with a landmark at drug start, cancer-stratified, purity- and proliferation-adjusted. Score × any-chemo interaction | Already done by ATLAS (6,902 patients, 2,460 events): interaction HR 0.93, p = 0.064. That design can detect an interaction of about 0.90; the observed 0.93 is at the edge. Report it as is. |
| Drug-level biomarker main effect, pooled across cancers | Per drug: carboplatin, paclitaxel, cisplatin, gemcitabine, doxorubicin, 5-FU, docetaxel, TMZ. Stratified by cancer, PFI from drug start | 210–506 PFI events, MDHR 0.82–0.88 → adequate for HR 0.8 |
| Breast neoadjuvant taxane/anthracycline response | Frozen score; meta-analysis of I-SPY2 control arm, BrighTNess, GSE25066, GSE41998, GSE20194, GSE20271, GSE32646 (about 2,000 patients, about 550 pCR). Adjust for ER/HER2 and proliferation | Adequate for AUROC 0.58–0.60 pooled. This is the one place the existing signal (0.60–0.66) can be confirmed or refuted. |
| Platinum predictive interaction in NSCLC | Frozen score × arm, meta-analysis of JBR.10 + GSE42127 + GSE68465 + GSE37745 | 0.90 power for interaction HR 0.7; about 0.53 for HR 0.8. Use JBR.10 as the randomised anchor. |
| Nivolumab vs everolimus predictive biomarker (RCC) | CheckMate-025 arm interaction | 0.83 for HR 0.7 |
| TMZ predictive (glioma) | CGGA693, all grades, score × TMZ, adjusted for IDH / 1p19q / MGMT (available in CGGA) | 0.82 for HR 0.7. Not randomised, so confounding by indication must be modelled. |

### 4.2 Not answerable with available data: retire or relabel

- **All 43 per-cancer × drug rows as predictive claims.** Drop the "validated externally" label for COAD oxaliplatin:
  q = 0.16, 27 events in the FOLFOX arm, and the same FOLFOX/FOLFIRI cohort is also used as a 5-FU validation.
  Relabel HNSC cisplatin and LUAD pemetrexed as "hypothesis (expected inflation 2–4×)".
  Relabel LGG TMZ as "prognostic in TCGA; predictive status unknown".
- **Per-cancer top-10 lists.** Their stability is about 0.03–0.08, so present them as prior-driven shortlists, not
  discoveries.
- **Drug-specific effects for any drug that is more than 80% given in combination**: cyclophosphamide, epirubicin,
  etoposide, oxaliplatin, pemetrexed, carboplatin, doxorubicin. Model them at the **regimen** level (AC, FOLFOX,
  carboplatin–paclitaxel, MVAC/GC, BEP) instead.
- **Organoid → patient transfer tested against TCGA overall survival.** The patient side has an MDHR of about 0.4–0.6,
  so this test cannot succeed even with a perfect model.

### 4.3 Needs new data

| Data | Why |
|---|---|
| Paired PDO response ↔ same patient's clinical response (EGA: Ooft 2019, Vlachogiannis 2018, de Witte 2020; `docs/data_access/`) | The only direct test of transfer. For a correlation of 0.4 at 80% power, about 47 pairs are needed. |
| Randomised single-agent or add-on trials with expression: GeparTrio / GeparSixto, CALGB 40603 (carboplatin add-on), STORM with clean processing, CALGB / NSABP colon adjuvant | Interaction-focused designs need 250+ events at HR 0.7. Only multi-trial IPD pools reach that. |
| Centrally reviewed RECIST or pCR labels, rather than TCGA site-reported best response | Label κ is 0.3–0.5. Under non-differential misclassification, the observed AUROC is 0.5 + (AUC_true − 0.5)(PPV + NPV − 1). With PPV 0.87 and NPV 0.63 (factor 0.50), a true AUROC of 0.65 appears as about **0.575**, close to the reported 0.56. |
| IDH / 1p19q / MGMT and HPV status joined to TCGA (PanCanAtlas subtype tables) | Needed to separate prognostic from predictive effects in LGG, GBM and HNSC. |

### 4.4 Design rules that maximise power

1. **Pool, then stratify.**
   - Use one model per drug (or regimen) across cancers, with cancer-stratified baselines and cancer-centred features.
     Partial pooling (random cancer slopes) replaces per-cancer fits.
   - Pre-specify about 10 hypotheses, not 43 × 6 tests. With m = 10 the FDR α* is about 0.02–0.06
     (2–4 true effects), and 70–80 events suffice for a main effect of HR 0.7.
2. **Make the interaction the primary estimand, and only where a comparator exists.**
   - Valid comparators are randomised arms or a reasonable non-randomised control (landmark plus IPTW on
     stage / age / grade / subtype).
   - Run TCGA "untreated vs treated" only at drug level, and only as a screening step.
3. **Use the right endpoint and time origin.** PFI or DFI measured from drug start (landmark) gives about 1.5× more
   events than OS and removes immortal-time bias. Model labels as ordinal (CR > PR ≈ SD > PD), not CR/PR vs SD/PD.
4. **Cross-fitting / sample splitting for every selection step.** Select pathways in fold k, score fold −k, and pool
   the out-of-fold scores, so that 100% of patients contribute to validation instead of the 30–60% held out now.
   Use permutation of the whole pipeline (already present) as the null.
5. **Benchmark against nulls that kill prognostic artefacts**:
   - random signatures of the same size (Venet);
   - proliferation and purity covariates (ABSOLUTE is now available in `robust/revamp_audit/data/`);
   - an immune-composition covariate, such as an IFN-γ / immunoregulatory pathway score.
6. **Freeze, then meta-analyse externally.** Use random-effects pooling over cohorts with the same regimen, rather
   than reporting the "best independent cohort". Report the pooled estimate, I² and the prediction interval.
7. **Report the MDHR beside every estimate.** Any row with MDHR < 0.6 should be shown as "not testable", not as
   "not supported".

---

## 5. Caveats of this audit

- The power formulas assume a per-SD continuous marker, proportional hazards, and independent tests. Pathway scores
  are correlated, so the effective number of tests is smaller than the nominal count. Simulation shows the
  interaction formula overstates power by about 4 points at around 240 events.
- `label_vs_pfi` takes PFI as the reference. PFI is itself imperfect, so κ is a lower bound on agreement with "true"
  response.
- The external interaction power uses the two largest arms of each cohort (for I-SPY2, control vs AMG386). Choosing
  other arms changes the numbers.
- The LGG analysis lacks IDH / 1p19q status. The prognostic interpretation is the most parsimonious one, not a proven one.
- The TAILORx effect description is from memory. All other literature numbers were checked against PubMed abstracts.
