# GSE7696

- Accession: GSE7696
- Paper: Murat A et al. 2008 J Clin Oncol, PMID 18565887
- Cancer: glioblastoma; setting: adjuvant
- Platform: GPL570
- Samples: 70; genes: 20848
- Responders (responder==1): 0; non-responders: 0; unlabelled: 70
- Randomised: False; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7696
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7696/

## Downloaded files (sha256)

- `GSE7696_series_matrix.txt.gz`: 6b4cb5d333e12f934e88e926026ca6a4aa161f1f37466ae802f3111b9c0caa25
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE7696_series_matrix.txt.gz (54675 probes); already log scale (range 3.00..14.38).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; probe with the highest mean kept per symbol -> 20848 genes.
- Primary GBM only (disease status 'GBM', n=70); recurrent / re-recurrent and non-tumoral samples dropped.
- arm: RT+TMZ (concomitant/adjuvant temozolomide) vs RT alone. Patients came from the EORTC 26981/NCIC CE.3 trial and associated studies; the expression subset is not a randomised sample, so arm is treated as non-randomised.
- No response endpoint (responder = NaN); OS in months, survival status 1 = dead.

## Arms

- RT+TMZ: n=43 (responders 0/0 labelled)
- RT: n=27 (responders 0/0 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.154, PC2 0.095, PC3 0.088

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | NA | NA | NA | not computable |
| arm | 0.481 | 0.506 | 0.476 | ok |

Built by `data/curated/trials/build_extra_trials.py GSE7696`.
