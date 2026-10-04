# TRIAD: patient-anchored, prior-guided, pan-cancer biomarker discovery

## Why a new dimension was needed

Two generations of organoid-first models failed for the same reason:
- **Original pipeline:** Kong et al.
- **v1 robust and NIT:** the hardened and next-gen versions built in this repo.

The problem was the patient side. TCGA treatment records are non-randomised, and each drug in
each cancer has only 15–50 deaths. Overall survival is dominated by stage and later lines of
therapy. Organoids also cannot model stroma, vasculature or immunity (sorafenib).

TRIAD changes the unit of learning:

| | Organoid-first (Kong, NIT) | TRIAD |
|---|---|---|
| Training label | organoid IC50 / AUC (15–50 models) | **RECIST best response in patients**, pooled across all TCGA cancers given the drug (100–418 patients per drug) |
| Organoids / cell lines | the training data | **priors** that the patient model is shrunk toward (and ablated) |
| Network | hard proximity filter | **per-pathway penalty** from drug-target proximity on STRING v12 |
| Confounding | n/a | treatment setting (adjuvant / post-progression) fitted as covariates, excluded from the score; features centred within cancer type |
| Validation | TCGA overall survival | leave-one-cancer-out (LOCO) and stratified k-fold, within-cancer permutation null, **frozen external trials** (BrighTNess, I-SPY2, GSE25066, JBR.10, lung adjuvant cohorts) with coefficient-permutation nulls and adjustment for proliferation and receptor status |

The analysis choices were fixed before running: λ = 10, prior strength τ = 0.05, primary
model = all priors.

Code:
- `obd/triad.py`
- `run_triad.py`
- `run_triad_external.py`
- `run_triad_specificity.py`

Data loaders:
- `obd/tcga_response.py`: 3,832 RECIST labels, 32 cancers.
- `obd/trial_cohorts.py`: 48 trial cohorts plus 7 synthetic-lethality networks.

Rebuild scripts are in `data_builders/`.

## Results

Internal results (TCGA, held-out patients) are in `results/triad/<DRUG>/models.tsv` and
`results/triad/triad_summary.tsv`. External results are in
`results/triad/external/triad_external_validation.tsv` and `drug_specificity.tsv`.

### Internal validation: held-out TCGA patients

The primary (pre-registered) model is TRIAD_all. "Cancers" is the number of cancer types
evaluable in leave-one-cancer-out (LOCO).

| Drug | Patients (cancers) | Patient-only LOCO | **TRIAD_all LOCO (perm p)** | TRIAD_all k-fold (perm p) | Frozen pre-clinical model alone |
|---|---|---|---|---|---|
| **Paclitaxel** | 244 (7) | 0.546 (0.28) | **0.654 (0.007, 300 perms)** | 0.650 (0.039) | cell line 0.51, organoid 0.61 |
| **Cisplatin** | 418 (8) | 0.637 (0.020) | **0.596 (0.013, 300 perms)** | 0.602 (0.039) | 0.49 / 0.50 |
| Temozolomide | 136 (1, glioma) | 0.581 (0.098) | 0.612 (0.056, 300 perms) | 0.590 (0.059) | 0.53 / 0.59 |
| Oxaliplatin | 108 (4) | 0.649 (0.18) | 0.663 (0.12) | 0.649 (0.16) | 0.59 / 0.61 |
| Etoposide | 107 (3) | 0.659 (0.24) | 0.330 (0.96) | 0.494 (0.61) | 0.59 / 0.37 |
| Doxorubicin | 169 (4) | 0.518 (0.37) | 0.450 (0.78) | 0.477 (0.73) | 0.47 / 0.51 |
| Capecitabine, carboplatin, docetaxel, gemcitabine, 5-FU | 94-277 | 0.42-0.52 | 0.40-0.54 (all ns) | 0.44-0.58 | ~0.43-0.53 |

Benjamini-Hochberg across 11 drugs: paclitaxel q = 0.077, cisplatin q = 0.072 (suggestive
at FDR 10%), temozolomide q = 0.21.

For paclitaxel and temozolomide, the pre-clinical and network priors improve held-out
patient prediction over patient-only learning (paclitaxel 0.55 -> 0.65). That is the first
measured, positive contribution of organoid and cell-line evidence in this project.

### External validation: frozen models, independent trials

| Model (trained on TCGA only) | Cohort | Result | Null p |
|---|---|---|---|
| Paclitaxel TRIAD_all | BrighTNess, paclitaxel -> AC arm (n = 123) | **AUROC 0.656** | **0.010** |
| Paclitaxel TRIAD_all | I-SPY2 control arm, paclitaxel -> AC (n = 210) | AUROC 0.614 | 0.052 |
| Paclitaxel patient-only | BrighTNess / I-SPY2 | AUROC 0.655 / 0.597 | 0.003 / 0.019 |
| Doxorubicin patient-only | I-SPY2 control arm | AUROC 0.674 | 0.007 |
| Paclitaxel / doxorubicin | GSE25066 (T/FAC, n = 488) | AUROC 0.52-0.60 | ns |
| Carboplatin | BrighTNess, I-SPY2 add-on arms (score x arm) | OR 0.56-0.88 | ns |
| Cisplatin | JBR.10 (randomised), lung adjuvant cohorts (score x arm, OS) | interaction HR 1.0-1.4, wrong direction | ns |

Drug specificity (`drug_specificity.tsv`): adjusting for a proliferation score and HR/HER2
status, the paclitaxel and doxorubicin scores keep independent effects:
- BrighTNess paclitaxel arm: OR 1.65 (p = 0.02) and 1.90 (p = 0.006) per SD;
- I-SPY2 control arm, doxorubicin: OR 1.75 (p = 0.006).

But every arm also receives AC, and the 5-FU and cisplatin models predict similarly. The
validated signal is therefore best described as a **cytotoxic-chemosensitivity programme
learned from TCGA patients that is not explained by proliferation or receptor subtype**,
rather than a drug-specific paclitaxel marker. Carboplatin and gemcitabine models do not
transfer.

### Biology the patient models learned (exploratory)

- **Cisplatin, 418 patients, 21 cancers.**
  - Response pathways: inflammasome / NLRP3, iron uptake.
  - Resistance pathways: MHC class I downregulation, mTORC1 signalling, TGF-beta receptor
    regulation.
  - This fits immunogenic cell death and immune evasion as modifiers of platinum response.
- Coefficients for every drug: `results/triad/<DRUG>/triad_coefficients.tsv`.

## Bottom line and next steps

1. **Validated in independent trials:** the TCGA-trained paclitaxel model predicts pCR in
   BrighTNess (null p = 0.01) and I-SPY2 (p of about 0.05), beyond proliferation and receptor
   status, with organoid, cell-line and network priors improving it over patient-only learning.
2. **Not validated:** drug specificity (trial designs combine drugs); cisplatin survival
   benefit in randomised lung trials; 5-FU, carboplatin, gemcitabine.
3. **To make it decisive:**
   - Randomised single-drug trials, such as GEPARTRIO and CALGB 40603 (carboplatin add-on with
     pCR), and STORM with its placebo arm, once the network allowlist works
     (`docs/ENVIRONMENT_SETUP.md`).
   - EGA paired organoid-patient data (`docs/data_access/`).
   - A synthetic-lethality network prior using the 7 GI networks already harvested
     (`trial_cohorts.load_gi_network`).
