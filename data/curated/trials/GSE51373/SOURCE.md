# GSE51373

- GEO accession: GSE51373
- Paper: Koti M et al. 2013 BMC Cancer, PMID 24237932
- Cancer: ovarian (high-grade serous); setting: adjuvant
- Platform: GPL570
- Samples: 28; genes: 20848
- Responders (responder==1): 16; non-responders: 12; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51373
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51373/

## Downloaded files (sha256)

- `GSE51373_series_matrix.txt.gz`: 5b1391902cb392c3423bd7d98a7b7aec7083f5912ca9cd707220baffb7ecb91d
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE51373_series_matrix.txt.gz (54675 probes); linear values (max 89446) -> log2 (floor 1).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- High-grade serous ovarian cancer; chemo sensitive vs resistant (platinum-free interval based, Koti et al.).
- Regimen per paper: carboplatin/paclitaxel after primary debulking (not per-sample annotated).

## Arms

- platinum-taxane: n=28 (responders 16/28 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.131, PC2 0.076, PC3 0.062

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.620 | 0.349 | 0.583 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE51373`.
