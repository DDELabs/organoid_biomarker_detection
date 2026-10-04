# GSE103479

- GEO accession: GSE103479
- Paper: Allen WL et al. 2018 JCO Precis Oncol, PMID 30088816
- Cancer: colorectal (stage II/III); setting: adjuvant
- Platform: GPL23985
- Samples: 156; genes: 24457
- Responders (responder==1): 0; non-responders: 0; unlabelled: 156
- Randomised: False; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103479
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103479/

## Downloaded files (sha256)

- `GSE103479_series_matrix.txt.gz`: 2d045df7f49352d32f4cd32ae51fcaf7791df6873ad387c21af0c0c309b4c4dd
- `GPL23985.soft.txt`: 399ee434dcda1d22906824cc2ba29122560fe48533c0a40eaea58cabf6b4b68b

## Processing

- Series matrix GSE103479_series_matrix.txt.gz (110961 probes); already log scale (range 0.74..15.01).
- Probes mapped with GPL23985.soft.txt; probes with no / multiple symbols dropped; median over probes per symbol -> 24457 genes.
- Stage II/III colorectal cancer, adjuvant chemotherapy yes/no (non-randomised); regimen not annotated (fluoropyrimidine +/- oxaliplatin per Allen et al.) -> 'FLUOROURACIL'.
- pfs_event from 'recurrence'; os_event from vital status; times as deposited (months).
- Almac Xcel array (GPL23985) symbols from the GEO platform table.

## Arms

- surgery only: n=87
- adjuvant chemo: n=66

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.096, PC2 0.086, PC3 0.066

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | NA | NA | NA | not computable |
| arm | 0.496 | 0.516 | 0.543 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE103479`.
