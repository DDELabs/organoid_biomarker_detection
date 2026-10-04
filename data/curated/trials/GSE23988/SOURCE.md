# GSE23988

- GEO accession: GSE23988
- Paper: Iwamoto T et al. 2011 J Natl Cancer Inst, PMID 21191116
- Cancer: breast; setting: neoadjuvant
- Platform: GPL96
- Samples: 61; genes: 12502
- Responders (responder==1): 20; non-responders: 41; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23988
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23988/

## Downloaded files (sha256)

- `GSE23988_series_matrix.txt.gz`: 7b13b1601c4e5eea54d1b8333983ab58c0e6c3ffd69c110bc685bc39b3d3cbc5
- `GPL96.annot.gz`: 88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4

## Processing

- Series matrix GSE23988_series_matrix.txt.gz (22283 probes); already log scale (range -2.43..18.70).
- Probes mapped with GPL96.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 12502 genes.
- Single regimen (US Oncology trial): FAC x4 followed by docetaxel + capecitabine x4.
- responder = 1 for pCR, 0 for RD.

## Arms

- FAC->TX: n=61 (responders 20/61 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.128, PC2 0.083, PC3 0.056

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.379 | 0.417 | 0.674 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE23988`.
