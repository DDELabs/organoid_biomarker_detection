# GSE20271

- GEO accession: GSE20271
- Paper: Tabchy A et al. 2010 Clin Cancer Res, PMID 20829329 (also PMID 23185353)
- Cancer: breast; setting: neoadjuvant
- Platform: GPL96
- Samples: 178; genes: 12502
- Responders (responder==1): 26; non-responders: 152; unlabelled: 0
- Randomised: True; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20271
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20271/

## Downloaded files (sha256)

- `GSE20271_series_matrix.txt.gz`: 711fe86a702893a68290b3b6d46f31baa2257d07a3bf69197f079ef662becca3
- `GPL96.annot.gz`: 88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4

## Processing

- Series matrix GSE20271_series_matrix.txt.gz (22283 probes); linear values (max 31169) -> log2 (floor 1).
- Probes mapped with GPL96.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 12502 genes.
- arm = randomised arm ('randomized (1=fac, 2=t/fac)'); treatment_received and the free-text 'preoperative treatment' (drug) are also kept. drugs derived from the regimen actually given (FEC -> EPIRUBICIN instead of DOXORUBICIN; T = weekly PACLITAXEL x12). drugs include preoperative switches after non-response (14 FAC-arm patients also got paclitaxel before surgery); use 'arm' for ITT.
- responder = 1 for pCR, 0 for RD.
- Possible sample overlap with other MDACC series (GSE20194/GSE22093) - deduplicate before pooling.

## Arms

- FAC: n=94 (responders 8/94 labelled)
- T/FAC: n=84 (responders 18/84 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.188, PC2 0.130, PC3 0.073

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.573 | 0.649 | 0.668 | ok |
| arm | 0.534 | 0.438 | 0.453 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE20271`.
