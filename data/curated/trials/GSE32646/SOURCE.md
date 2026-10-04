# GSE32646

- GEO accession: GSE32646
- Paper: Miyake T et al. 2012 Cancer Sci, PMID 22320227
- Cancer: breast; setting: neoadjuvant
- Platform: GPL570
- Samples: 115; genes: 20848
- Responders (responder==1): 27; non-responders: 88; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32646
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32646/

## Downloaded files (sha256)

- `GSE32646_series_matrix.txt.gz`: 1a0f6c79a73982c51386da46cd9623a8daad0c124f89ccdd9d8a773d8212cc40
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE32646_series_matrix.txt.gz (54613 probes); already log scale (range -5.25..16.53).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- Single-arm: paclitaxel (80 mg/m2 weekly x12) followed by FEC x4 (Miyake et al.).
- responder = 1 for pCR, 0 for nCR (non-pCR).

## Arms

- P->FEC: n=115 (responders 27/115 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.080, PC2 0.065, PC3 0.042

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.222 | 0.429 | 0.670 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE32646`.
