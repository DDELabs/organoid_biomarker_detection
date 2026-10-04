# GSE62080

- GEO accession: GSE62080
- Paper: Del Rio M et al. 2007 J Clin Oncol, PMID 17327601 (also PMID 30863148)
- Cancer: colorectal (metastatic); setting: metastatic
- Platform: GPL570
- Samples: 21; genes: 20848
- Responders (responder==1): 9; non-responders: 12; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62080
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62080/

## Downloaded files (sha256)

- `GSE62080_series_matrix.txt.gz`: 8e830c1cb0e5c2930a24257f78b03223c79ca657c4e15b591e40d550d58c1424
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE62080_series_matrix.txt.gz (54675 probes); already log scale (range 2.59..14.82).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- First-line FOLFIRI in metastatic CRC; sensitive (objective response) -> responder 1, resistant -> 0.

## Arms

- FOLFIRI: n=21 (responders 9/21 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.188, PC2 0.122, PC3 0.097

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.620 | 0.537 | 0.667 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE62080`.
