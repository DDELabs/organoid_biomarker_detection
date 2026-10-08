# BRAUN2020_CHECKMATE

- Accession: Braun2020 NatMed suppl (dbGaP phs001493 / EGA for raw)
- Paper: Braun DA et al. 2020 Nat Med, PMID 32472114
- Cancer: kidney renal clear cell carcinoma (advanced); setting: metastatic
- Platform: RNA-seq (Illumina)
- Samples: 311; genes: 19487
- Responders (responder==1): 44; non-responders: 237; unlabelled: 30
- Randomised: True; control/comparator arm: True

## URLs

- https://www.nature.com/articles/s41591-020-0839-y
- https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/MediaObjects/41591_2020_839_MOESM2_ESM.xlsx

## Downloaded files (sha256)

- `41591_2020_839_MOESM2_ESM.xlsx`: a1f8683676d416b255291dcaa007aa0fa48c736c23e6a7a566664ba528bf205b

## Processing

- Braun DA et al. 2020 Nat Med Supplementary Table S1 (clinical) and S4A (normalised RNA expression matrix, 311 samples) downloaded from the Springer static-content server.
- Expression kept as deposited (log2-scale normalised, batch-corrected per paper; values ~0-68, median ~24 so absolute levels are shifted vs log2(TPM+1)); gene symbols restricted to UniProt gene names -> 19487 genes.
- responder: ORR CR/PR/CRPR = 1, SD/PD = 0, NE = NaN (RECIST 1.1 per trial).
- clinical_benefit per paper: CB (CR/PR or SD with tumour shrinkage and PFS >= 6 mo) = 1, NCB (PD with PFS < 3 mo) = 0, ICB (intermediate) = NaN.
- PFS/OS in months; *_CNSR = 1 is an event (871/1006 PFS events) as in the paper's code.
- CM-025 is the randomised phase III (nivolumab vs everolimus); CM-009 and CM-010 are nivolumab-only phase I/II cohorts. Everolimus arm = comparator for treatment-interaction analyses.

## Arms

- NIVOLUMAB: n=181 (responders 39/172 labelled)
- EVEROLIMUS: n=130 (responders 5/109 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.090, PC2 0.047, PC3 0.036

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.513 | 0.446 | 0.489 | ok |
| arm | 0.527 | 0.507 | 0.528 | ok |

Built by `data/curated/trials/build_extra_trials.py BRAUN2020_CHECKMATE`.
