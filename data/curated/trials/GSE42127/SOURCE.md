# GSE42127

- GEO accession: GSE42127
- Paper: Tang H et al. 2013 Clin Cancer Res, PMID 23357979
- Cancer: NSCLC; setting: adjuvant
- Platform: GPL6884
- Samples: 176; genes: 19405
- Responders (responder==1): 0; non-responders: 0; unlabelled: 176
- Randomised: False; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42127
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42127/

## Downloaded files (sha256)

- `GSE42127_series_matrix.txt.gz`: 0bcc8f1b641db90b76488a2e173444e1a9cdf6cd74eb183b774cfb7ad8fcef43
- `GPL6884.annot.gz`: 8bef75ebe5b7e28bf61cf398a31e28d63d3a9884405a3510f60ac534a0b59abb

## Processing

- Series matrix GSE42127_series_matrix.txt.gz (48803 probes); already log scale (range 2.41..15.26).
- Probes mapped with GPL6884.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 19405 genes.
- Resected NSCLC (MDACC); adjuvant chemotherapy yes/no, non-randomised. Regimen not annotated (platinum doublets per Tang et al.) -> 'PLATINUM'. No response labels; OS only.

## Arms

- observation: n=127
- adjuvant chemo: n=49

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.178, PC2 0.102, PC3 0.053

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | NA | NA | NA | not computable |
| arm | 0.453 | 0.441 | 0.510 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE42127`.
