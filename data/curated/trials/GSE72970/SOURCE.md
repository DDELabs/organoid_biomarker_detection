# GSE72970

- GEO accession: GSE72970
- Paper: Del Rio M et al. 2017 J Clin Oncol, PMID 28284171 (also PMID 28659146, 30863148)
- Cancer: colorectal (metastatic); setting: metastatic
- Platform: GPL570
- Samples: 124; genes: 20848
- Responders (responder==1): 63; non-responders: 61; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72970
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72970/

## Downloaded files (sha256)

- `GSE72970_series_matrix.txt.gz`: c1b1e66ab63f67f9ec121e7b2250cfd883d15de94abcd89ee8d1ca3405996d9e
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE72970_series_matrix.txt.gz (54675 probes); already log scale (range 2.32..14.93).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- responder from 'response status' (R = CR/PR, NR = SD/PD; RECIST best response 'response category').
- 'pfs censored'/'os censored' = 1 taken as event (1 for 114/124 PFS, i.e. progression observed).
- PFS/OS units: months as deposited.
- Regimens expanded to generic names (FOLFIRI = FLUOROURACIL;LEUCOVORIN;IRINOTECAN, ERBITUX = CETUXIMAB, XELIRI = CAPECITABINE;IRINOTECAN).

## Arms

- FOLFIRI: n=60 (responders 27/60 labelled)
- FOLFOX: n=32 (responders 20/32 labelled)
- FOLFIRI+BEVACIZUMAB: n=22 (responders 9/22 labelled)
- FOLFOX+BEVACIZUMAB: n=4 (responders 1/4 labelled)
- FOLFIRINOX: n=3 (responders 3/3 labelled)
- FOLFIRI+ERBITUX: n=1 (responders 1/1 labelled)
- XELIRI+BEVACIZUMAB: n=1 (responders 1/1 labelled)
- FOLFIRINOX+BEVACIZUMAB: n=1 (responders 1/1 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.148, PC2 0.099, PC3 0.089

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.464 | 0.508 | 0.468 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE72970`.
