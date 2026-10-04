# GSE52219

- GEO accession: GSE52219
- Paper: Choi W et al. 2014 Cancer Cell, PMID 24525232
- Cancer: bladder (MIBC); setting: neoadjuvant
- Platform: GPL14951
- Samples: 23; genes: 20794
- Responders (responder==1): 6; non-responders: 17; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52219
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52219/

## Downloaded files (sha256)

- `GSE52219_series_matrix.txt.gz`: daabc6a1f8b8d79e7b7d3fb351ca8e67d6f5f8fb4ea2f9005648f13bdb98477f
- `GPL14951.soft.txt`: c04d2950c404c8124ba971114697f8c2c1b9f658c7e71e45ca237559f76a5906

## Processing

- Series matrix GSE52219_series_matrix.txt.gz (29377 probes); already log scale (range 3.51..16.13).
- Probes mapped with GPL14951.soft.txt; probes with no / multiple symbols dropped; median over probes per symbol -> 20794 genes.
- Pre-treatment TURBT of MIBC; responder = downstaged to pT0 or pT1 after neoadjuvant MVAC.

## Arms

- MVAC: n=23 (responders 6/23 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.156, PC2 0.122, PC3 0.074

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.471 | 0.471 | 0.480 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE52219`.
