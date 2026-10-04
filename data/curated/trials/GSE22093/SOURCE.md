# GSE22093

- GEO accession: GSE22093
- Paper: Iwamoto T et al. 2011 J Natl Cancer Inst, PMID 21191116
- Cancer: breast; setting: neoadjuvant
- Platform: GPL96
- Samples: 103; genes: 12502
- Responders (responder==1): 28; non-responders: 69; unlabelled: 6
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22093
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22093/

## Downloaded files (sha256)

- `GSE22093_series_matrix.txt.gz`: 5c71ee1057b37094db1af5b39b002c81e7b1167813a9d55a9e19818cae2d7389
- `GPL96.annot.gz`: 88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4

## Processing

- Series matrix GSE22093_series_matrix.txt.gz (22283 probes); already log scale (range -2.99..18.62).
- Probes mapped with GPL96.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 12502 genes.
- Anthracycline-based FAC/FEC (no taxane); per-sample anthracycline (doxorubicin vs epirubicin) not annotated, drugs lists DOXORUBICIN.
- responder = 1 for pCR, 0 for RD; 'NA' -> NA.
- Overlaps with other MDACC series possible.

## Arms

- FAC/FEC: n=103 (responders 28/97 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.164, PC2 0.087, PC3 0.060

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.286 | 0.594 | 0.595 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE22093`.
