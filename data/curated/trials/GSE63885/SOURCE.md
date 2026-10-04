# GSE63885

- GEO accession: GSE63885
- Paper: Lisowska KM et al. 2014 Front Oncol, PMID 24478986 (also PMID 27028324)
- Cancer: ovarian; setting: adjuvant
- Platform: GPL570
- Samples: 101; genes: 20848
- Responders (responder==1): 65; non-responders: 10; unlabelled: 26
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63885
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63885/

## Downloaded files (sha256)

- `GSE63885_series_matrix.txt.gz`: 9815d81e390dd2e429c5857cb35da57a072118ce68bdea21282c2db130b15f94
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE63885_series_matrix.txt.gz (54675 probes); already log scale (range 3.01..14.44).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- Surgical tumour samples before first-line chemotherapy (taxane/platinum or platinum/cyclophosphamide); platinum agent (cisplatin vs carboplatin) not specified -> 'PLATINUM'.
- response = clinical status after 1st line chemo (CR/PR/SD/P->PD); responder = 1 for CR/PR.
- OS/DFS days converted to months (/30.4375); os_event = 1 for DOD (dead of disease), 0 for AWD/NED. DFS deposited as 0 for non-CR patients, so dfs_months is only meaningful for CR patients; no DFS event indicator deposited.
- platinum_sensitivity (resistant/moderately/highly sensitive) kept as extra column.

## Arms

- taxane/platinum: n=41 (responders 38/41 labelled)
- platinum/cyclophosphamide: n=34 (responders 27/34 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.119, PC2 0.089, PC3 0.059

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.615 | 0.548 | 0.600 | ok |
| arm | 0.491 | 0.355 | 0.575 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE63885`.
