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
