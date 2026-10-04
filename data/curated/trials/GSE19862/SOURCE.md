# GSE19862

- GEO accession: GSE19862
- Paper: Watanabe T et al. (Teikyo Univ.); no PMID in GEO
- Cancer: colorectal (advanced); setting: metastatic
- Platform: GPL570
- Samples: 14; genes: 20848
- Responders (responder==1): 7; non-responders: 7; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19862
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19862/

## Downloaded files (sha256)

- `GSE19862_series_matrix.txt.gz`: 334be6cd123d7b5db435a6c8393a5171e9a7d549e4e10c8ec2c53d4822e15831
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE19862_series_matrix.txt.gz (54675 probes); already log scale (range -4.61..7.38).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- Advanced CRC, response to bevacizumab-based therapy (BV_Responder / BV_Non-responder). The chemotherapy backbone is not annotated in GEO; drugs lists BEVACIZUMAB only.
- Deposited values are on an unusual scale (about -4.6..7.4, unlike GSE19860 from the same group); kept as deposited - standardise per gene before use. n=14 only.

## Arms

- bevacizumab-containing: n=14 (responders 7/14 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.138, PC2 0.132, PC3 0.105

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.061 | 0.510 | 0.347 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE19862`.
