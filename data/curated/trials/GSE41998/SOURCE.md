# GSE41998

- GEO accession: GSE41998
- Paper: Horak CE et al. 2013 Clin Cancer Res, PMID 23340299
- Cancer: breast; setting: neoadjuvant
- Platform: GPL571
- Samples: 279; genes: 12502
- Responders (responder==1): 69; non-responders: 184; unlabelled: 26
- Randomised: True; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41998
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41998/

## Downloaded files (sha256)

- `GSE41998_series_matrix.txt.gz`: 2199e27f41415f4089f97c77c6d38d0c71cbcc6fe464c7bf9d771eb25b2ea4c2
- `GPL571.annot.gz`: 9f5a482f9f3d478b47fc80124600319e569f1e9ab447f3f256d61b1fbacce4fa

## Processing

- Series matrix GSE41998_series_matrix.txt.gz (22277 probes); already log scale (range 2.50..14.77).
- Probes mapped with GPL571.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 12502 genes.
- All patients received 4x AC, then were randomised to ixabepilone or paclitaxel; 'none' = not randomised / AC only (kept, randomised=False).
- responder from 'pcr' (Yes=1, No=0); value '0' (n=20) and blanks are ambiguous and set to NA. pcr_rcb1 (pCR or RCB-I) and clinical response to AC (ac_response) kept as extra columns.
- control_arm = paclitaxel (standard) arm.

## Arms

- AC->ixabepilone: n=138 (responders 35/132 labelled)
- AC->paclitaxel: n=127 (responders 34/121 labelled)
- AC only: n=14 (responders 0/0 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.203, PC2 0.104, PC3 0.059

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.346 | 0.601 | 0.444 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE41998`.
