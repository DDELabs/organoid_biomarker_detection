# GSE67501

- Accession: GSE67501
- Paper: Ascierto ML et al. 2016 Cancer Immunol Res, PMID 27491898
- Cancer: renal cell carcinoma (metastatic); setting: metastatic
- Platform: GPL14951
- Samples: 11; genes: 20794
- Responders (responder==1): 4; non-responders: 7; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67501/

## Downloaded files (sha256)

- `GSE67501_series_matrix.txt.gz`: e48b9c3b01a160ff59233d185a17df1c1210ba244b637a0a85ea08aaf4b5bd2a
- `GPL14951.soft.txt`: c04d2950c404c8124ba971114697f8c2c1b9f658c7e71e45ca237559f76a5906

## Processing

- Series matrix GSE67501_series_matrix.txt.gz (29377 probes); already log scale (range 5.71..14.61).
- Probes mapped with GPL14951.soft.txt; probes with no / multiple symbols dropped; probe with the highest mean kept per symbol -> 20794 genes.
- responder: CR/PR = 1, SD/NR = 0 (matches the series' response/no_response split).
- Archival tumour (collected 2-81 months before nivolumab). n=11 only; no survival data.

## Arms

- NIVOLUMAB: n=11 (responders 4/11 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.195, PC2 0.176, PC3 0.137

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.786 | 0.714 | 0.607 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_extra_trials.py GSE67501`.
