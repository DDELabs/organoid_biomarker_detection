# GSE140901

- Accession: GSE140901
- Paper: Hsu CL et al. 2021 Liver Cancer, PMID 34414122
- Cancer: hepatocellular carcinoma (advanced); setting: metastatic
- Platform: GPL19965 (NanoString IO360, ~770 genes)
- Samples: 24; genes: 784
- Responders (responder==1): 6; non-responders: 18; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140901
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140901/

## Downloaded files (sha256)

- `GSE140901_series_matrix.txt.gz`: 8aed1c712c8a33c2154b5ad241d3a5b04e9fec7dd8d9805ab6c50b17b8a9ee9b
- `GSE140901_processed_data.txt.gz`: acb4888efbd4c15339445f9130243f28075de2dcf971fe1a2c2a0a70ed4e72bf

## Processing

- Supplementary GSE140901_processed_data.txt.gz: NanoString PanCancer IO 360 panel (GPL19965), TMM-normalised log2 CPM as deposited; 784 genes only (targeted immune panel, not genome-wide).
- Archival pre-treatment tumour; anti-PD-1/PD-L1-based ICI (agent/combination per patient not deposited) -> drugs = 'ICI_ANTI_PD1_PDL1' placeholder.
- responder: best response PR = 1, SD/PD = 0 (no CR); clinical_benefit as deposited.
- pfs_time / os_time are in weeks (shortest PFS 5.1 = first restaging) -> months = weeks * 7 / 30.4375.

## Arms

- ICI: n=24 (responders 6/24 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.381, PC2 0.090, PC3 0.075

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.796 | 0.713 | 0.435 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_extra_trials.py GSE140901`.
