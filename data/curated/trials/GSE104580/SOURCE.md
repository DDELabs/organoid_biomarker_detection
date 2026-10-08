# GSE104580

- Accession: GSE104580
- Paper: unpublished (GEO submitter: Kam Hui, National Cancer Centre Singapore)
- Cancer: hepatocellular carcinoma; setting: locoregional
- Platform: GPL570
- Samples: 147; genes: 20848
- Responders (responder==1): 81; non-responders: 66; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104580
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104580/

## Downloaded files (sha256)

- `GSE104580_series_matrix.txt.gz`: d28ad11495a5861c3ccd8c92069975b79b13b167c461e4c3c595a11bc172d09a
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE104580_series_matrix.txt.gz (54675 probes); already log scale (range -2.40..18.64).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; probe with the highest mean kept per symbol -> 20848 genes.
- Deposited data: GCRMA, log2, batch-corrected on scan date in Partek (per GEO).
- responder from 'subject subgroup' (TACE responders = 1 / non-responders = 0). GEO and the submitter (National Cancer Centre Singapore) give no response criterion (presumably mRECIST/EASL radiological response) and no publication is linked; no survival data.
- Chemotherapeutic agent(s) used in TACE are not stated; drugs = 'TACE' (placeholder).

## Arms

- TACE: n=147 (responders 81/147 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.150, PC2 0.135, PC3 0.056

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.777 | 0.591 | 0.500 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_extra_trials.py GSE104580`.
