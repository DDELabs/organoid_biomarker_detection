# GSE15622

- GEO accession: GSE15622
- Paper: Ahmed AA et al. 2007 Cancer Cell, PMID 18068629 (also PMID 25560085)
- Cancer: ovarian (high-grade serous); setting: neoadjuvant
- Platform: GPL8414
- Samples: 35; genes: 11196
- Responders (responder==1): 22; non-responders: 13; unlabelled: 0
- Randomised: False; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15622
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15622/

## Downloaded files (sha256)

- `GSE15622_series_matrix.txt.gz`: 4f234c7f571ef376a62dc3705da14cf10660770187d9ff8ce57a93facf52c390
- `GPL8414.soft.txt`: cfae5be98742af02cb48b8ce67801fab379b4e1dbcff4849f57daa15b0162d97

## Processing

- Series matrix GSE15622_series_matrix.txt.gz (11896 probes); already log scale (range 5.19..13.54).
- Probes mapped with GPL8414.soft.txt; probes with no / multiple symbols dropped; median over probes per symbol -> 11196 genes.
- CTCR-OV01: 3 cycles of single-agent paclitaxel or carboplatin before debulking; only pre-treatment biopsies kept (post-treatment samples dropped). Allocation is not described as randomised in GEO, so randomised=False; the two single-agent arms are comparators.
- response sensitive/resistant (CA-125 based, per GEO); responder = 1 for sensitive.
- Almac Ovarian Cancer DSA (GPL8414) probesets annotated with ENSG ids; mapped to HGNC symbols with data/ENSG_GENESYMBOL.txt.

## Arms

- paclitaxel: n=20 (responders 13/20 labelled)
- carboplatin: n=15 (responders 9/15 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.200, PC2 0.165, PC3 0.091

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.573 | 0.524 | 0.500 | ok |
| arm | 0.473 | 0.530 | 0.600 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE15622`.
