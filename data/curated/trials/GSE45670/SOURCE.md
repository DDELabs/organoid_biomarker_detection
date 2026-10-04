# GSE45670

- GEO accession: GSE45670
- Paper: Wen J et al. 2014 Dis Esophagus / Oncotarget, PMID 24907633
- Cancer: oesophageal squamous cell carcinoma; setting: neoadjuvant
- Platform: GPL570
- Samples: 28; genes: 20848
- Responders (responder==1): 11; non-responders: 17; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45670
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45670/

## Downloaded files (sha256)

- `GSE45670_series_matrix.txt.gz`: 2c18d12865956901ec52aa77adcd4efe456b8df0798db2ba3df37cc723acbbab
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394

## Processing

- Series matrix GSE45670_series_matrix.txt.gz (54675 probes); linear values (max 104780) -> log2 (floor 1).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20848 genes.
- Pre-treatment endoscopic ESCC biopsies (n=28); 10 normal-epithelium samples dropped.
- Preoperative chemoradiotherapy (vinorelbine + cisplatin + RT); responder = 1 for pCR.

## Arms

- CRT: n=28 (responders 11/28 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.115, PC2 0.074, PC3 0.061

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.663 | 0.278 | 0.674 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE45670`.
