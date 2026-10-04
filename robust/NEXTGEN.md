# NIT: Network-Informed, deconfounded Transfer framework

NIT is the next-generation model built on the Kong et al. (2020) framework. It keeps that
framework's foundations:
- drug-target network proximity;
- Reactome pathway features;
- an organoid-trained linear model;
- a patient survival readout.

It adds one fix for each failure mode found while making the pipeline robust. Code:
`robust/obd/nextgen.py`, driver `robust/run_nextgen.py`.

## Why the original framework failed, and the fix for each failure

| Failure (evidence) | Fix in NIT |
|---|---|
| **Proliferation and general-sensitivity confounding.** In viability assays, fast-growing organoids look sensitive to most drugs. LICOB sorafenib AUC correlates with proliferation (rho = -0.35, p = 0.012) and with PC1 of all 76 drugs (rho = 0.40, p = 0.004). In patients proliferation is prognostic, so a raw-AUC model learns "slow tumours do better" and the HR flips. | `deconfound`: learn from the residual of the response after regressing out a Hallmark E2F/G2M proliferation score and leave-drug-out general sensitivity (PC1 of the other drugs). |
| **Immortal-time bias.** 93% of TCGA-LIHC sorafenib patients started more than 90 days after diagnosis (median 346 days, at recurrence), but survival was counted from diagnosis. | `landmark`: survival re-based at treatment start (`tcga_biotab_drug_start`). |
| **Prognostic, not predictive, signal.** | Patient proliferation added as a Cox covariate. The score x treatment interaction is reported, and a signature must beat random signatures (Venet 2011). |
| **Hard proximity cut-off is brittle.** STRING version alone changes 30-50% of proximal pathways. | `soft_prior`: logistic weight in the proximity z for every pathway, equivalent to a per-pathway ridge penalty of 1/w^2. |
| **Tiny n (15-50 organoids).** | `transfer`: a GDSC2 pan-cancer cell-line model of the same drug (about 930 lines) gives prior coefficients. Organoid fits shrink toward it instead of toward zero. |
| **Mechanisms organoids cannot model.** Sorafenib is largely anti-angiogenic (Liu 2006; Wilhelm 2008). Published patient predictors are immune-related (Pinyol 2019, STORM). | Documented as a limit of the organoid paradigm. Such drugs need tissue-level features or co-culture models. |

## Validation hierarchy

Every model passes through the same robust runner:
- whole-pipeline permutation of the organoid response (500x);
- organoid jackknife;
- patient bootstrap;
- cut-off sensitivity;
- a 5-check pathway scorecard;
- external validation of frozen signatures with a 1,000-signature random null.

## Ablation: adding one fix at a time

Adjusted HR per SD of predicted resistance in drug-treated TCGA patients (> 1 = expected).
The steps are cumulative from left to right.

| Study (model) | v1 | +landmark | +prolif_cov | +deconfound | +soft_prior | +transfer |
|---|---|---|---|---|---|---|
| LIHC sorafenib (LICOB PDO) | **0.39** (reversed, p = 0.035) | 0.59 | 0.63 | **1.28** | **1.55** | 1.42 |
| COAD 5-FU (vdW PDO) | 1.22 | 1.05 | 1.06 | 0.92 | 0.88 | 1.01 |
| PAAD 5-FU (Tiriac PDO) | 1.25 | 1.30 | 1.38 | 1.37 | 1.10 | 1.10 |
| PAAD gemcitabine (Tiriac PDO) | 0.75 | 0.72 | 0.81 | 0.80 | 0.80 | 0.82 (interaction p = 0.016) |
| BRCA doxorubicin (GDSC proxy) | 0.84 | 1.10 | 2.32 | 1.18 | 1.19 | - |
| LUAD cisplatin (GDSC proxy) | 0.95 | 0.93 | 1.25 | 1.05 | 0.71 | - |
| BLCA gemcitabine (GDSC proxy) | 0.75 | 0.74 | 0.72 | 0.60 | 0.69 | - |

Full tables: `results/nextgen/ablation_*.tsv`. Reports: `results/nextgen/<STUDY>/REPORT.md`.
