# GSE19860

- GEO accession: GSE19860
- Paper: Watanabe T et al. 2011 Int J Cancer (FOLFOX response signature); no PMID in GEO
- Cancer: colorectal (advanced); setting: metastatic
- Platform: GPL570
- Samples: 40; genes: 20848
- Responders (responder==1): 15; non-responders: 25; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19860
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19860/

## Downloaded files (sha256)

- `GSE19860_series_matrix.txt.gz`: 6a8b15d5844a5cde3ebe166436623e389470a54cce4f2d49e5f31eff94fa79be
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE19860_series_matrix.txt.gz (54675 probes); linear values (max 16801) -> log2 (floor 1).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- Advanced CRC treated with mFOLFOX6; responder from 'FL_Responder' / 'FL_Non_responder'.
- Some patients later received bevacizumab; that response is kept in bevacizumab_response.

## Arms

- mFOLFOX6: n=40 (responders 15/40 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.994, PC2 0.001, PC3 0.001

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.405 | 0.595 | 0.477 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE19860`.
