# GSE78220

- Accession: GSE78220
- Paper: Hugo W et al. 2016 Cell, PMID 26997480
- Cancer: skin cutaneous melanoma (metastatic); setting: metastatic
- Platform: RNA-seq (GPL11154)
- Samples: 28; genes: 19412
- Responders (responder==1): 15; non-responders: 13; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/

## Downloaded files (sha256)

- `GSE78220_series_matrix.txt.gz`: b2d7829dac364e7773390255252adb5c42f07104136a1f93d2a767e9522c086f
- `GSE78220_PatientFPKM.xlsx`: ae3b044f23a0a4cd2859da36726a220856c35216a169c17f93dfc3bd20b04de6

## Processing

- Series matrix GSE78220_series_matrix.txt.gz (no expression) + supplementary GSE78220_PatientFPKM.xlsx (FPKM, 28 patients); log2(FPKM + 1), symbols restricted to UniProt gene names -> 19412 genes.
- responder: irRECIST CR/PR = 1, PD = 0 (no SD in this cohort).
- Treatment from the GEO 'treatment' field (pembrolizumab for all 28); one biopsy is annotated on-treatment (biopsy_time; Pt16), the rest pre-treatment. Pt27A/Pt27B are two baseline lesions of one patient (same `patient`).
- OS days / 30.4375; vital status Dead = event.

## Arms

- ANTI-PD-1: n=28 (responders 15/28 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.156, PC2 0.126, PC3 0.088

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.364 | 0.272 | 0.405 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_extra_trials.py GSE78220`.
