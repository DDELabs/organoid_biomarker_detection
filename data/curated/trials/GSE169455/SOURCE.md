# GSE169455

- GEO accession: GSE169455
- Paper: Sjodahl G et al. 2022 Eur Urol, PMID 34782206
- Cancer: bladder (MIBC); setting: neoadjuvant
- Platform: GPL6244
- Samples: 149; genes: 13624
- Responders (responder==1): 48; non-responders: 101; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169455
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169455/

## Downloaded files (sha256)

- `GSE169455_series_matrix.txt.gz`: 7fcf4a34fc1543d2000e285d769fbd59aae7ba7696222996a543cf6e8cb92590
- `GSE169455_normalized_by_gene.txt.gz`: e228f59d595dab6b3b5fe46bbe89ae1600ce7158be46dac2eb930e30ed1ab1a7

## Processing

- Expression from supplementary GSE169455_normalized_by_gene.txt.gz (gene-level RMA, Affymetrix Human Gene 1.0 ST, GPL6244; rows already gene symbols; duplicated symbols median-collapsed); already log scale (range 0.66..13.54). Columns renamed from sample titles to GSM ids.
- arm = neoadjuvant vs induction (cN+) cisplatin-based chemotherapy (Sjodahl et al.: mostly gemcitabine/cisplatin or (dd)MVAC; per-sample regimen not deposited, drugs lists CISPLATIN only).
- responder = 1 for pT0N0 at cystectomy (pCR), 0 otherwise; 'downstaged' = <=pT1N0 (incl. pTa/pTis).
- Labeling kit / batch kept (two amplification kits) - check batch before pooling.

## Arms

- neoadjuvant: n=125 (responders 42/125 labelled)
- induction: n=24 (responders 6/24 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.175, PC2 0.096, PC3 0.084

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.483 | 0.555 | 0.415 | ok |
| arm | 0.462 | 0.527 | 0.496 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE169455`.
