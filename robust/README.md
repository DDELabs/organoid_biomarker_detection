# Robust organoid-to-patient biomarker pipeline (`robust/`)

A Python 3 extension of the Kong et al. (2020) network-based pipeline. The original
code in `python/` and `utilities/` is left in place (Python 2.7) and is reproduced
exactly by this package. Every candidate biomarker additionally has to pass a
**robustness scorecard** before it is reported.

## Why

The original readout is fragile. Swapping the 2019 TCGA files for the current
GDC release moves pathway scores only slightly (median per-pathway correlation
0.91), yet the 5-FU log-rank p for k=5 pathways moves from 0.13 to 0.72. Picking
the best of 3 models x 9 values of k also inflates significance: the paper's best
p = 0.076 has a whole-pipeline permutation p of about 0.5. See
`results/COAD_FLUOROURACIL/REPORT.md`.

## What it does

| Step | Original | Added here |
|---|---|---|
| Patient expression | HTSeq FPKM-UQ (GENCODE v22, no longer distributed) | Current GDC STAR FPKM-UQ (GENCODE v36) from the Xena GDC mirror, primary tumours only, one aliquot per patient (`cohorts.tcga_star_fpkm_uq`) |
| Pathway scores | ssGSEA, z-scored within each cohort | ssGSEA raw ES, plus a single-sample rank score (singscore-style) that is invariant to units (FPKM/FPKM-UQ/TPM) and to cohort composition (`scoring.py`) |
| Network step | closest-distance proximity z on STRING > 700, one cut-off | Python 3 port, plus random-walk-with-restart propagation, consensus of both, cut-off sensitivity (`network.py`) |
| Organoid model | Ridge / SVR / OLS fitted once, top-k by \|coef\| | Stability selection across subsamples x 4 learners (Ridge, linear SVR, elastic net, Spearman), sign agreement, LOOCV with selection inside folds, optional PRECISE organoid/tumour alignment (`models.py`) |
| Patient test | median split, log-rank | + continuous score, Cox adjusted for stage/age/sex, C-index, score x treatment interaction (predictive vs prognostic), patient bootstrap (`study.py`, `survival.py`) |
| Significance | nominal log-rank p | Whole-pipeline permutation of organoid IC50, including for the original best-of-k procedure (`Study.permutation_null`) |
| Fragility | none | Data-version comparison (paper / current ssGSEA / current rank), leave-one-organoid-out, cut-off sensitivity (`runner.py`) |

### Scorecard

Each pathway gets five checks: selection frequency >= 0.5, sign agreement >= 0.9,
leave-one-organoid-out retention >= 0.7, consistent in every data version, and a
patient hazard ratio in the direction implied by the organoid coefficient.
5/5 = **robust**, 3-4/5 = **supported**, otherwise **fragile**.

The signature as a whole is **validated** only if the permutation p < 0.05, the
direction holds in every data version and the treatment-interaction p < 0.1. The
report also flags cohorts with fewer than 30 events in treated patients as underpowered.

## Run

```bash
pip install pandas numpy scipy scikit-learn statsmodels gseapy networkx matplotlib joblib pytest
python -m pytest robust/tests -q          # exact reproduction of the paper + invariance tests
python robust/run_coad_5fu.py --perm 500  # writes robust/results/COAD_FLUOROURACIL/
```

The first run downloads `TCGA-COAD.star_fpkm-uq.tsv.gz`, `.clinical.tsv.gz` and
`.survival.tsv.gz` (~130 MB) from
`https://gdc-hub.s3.us-east-1.amazonaws.com/download/` into `data/external/`
(git-ignored).

## Adding a cancer type

You need:
- organoid expression (genes x samples) and drug response (long table: sample, drug, response), loaded with `cohorts.generic_organoids`;
- TCGA expression via `cohorts.tcga_star_fpkm_uq("<PROJECT>")` and survival via `cohorts.tcga_xena_clinical`;
- the set of treated patients (`{drug: set(barcodes)}`), parsed from the BCR biotab `clinical_drug` file with `cohorts.tcga_biotab_drugs`;
- drug targets (`reference.drug_targets` / `drugbank_targets`) and a PPI edge list for `network.proximity_z` / `rwr_z`, or the precomputed table.

Build a `Study` per data version and call `runner.run`. See `run_coad_5fu.py`.
