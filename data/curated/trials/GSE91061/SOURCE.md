# GSE91061

- Accession: GSE91061
- Paper: Riaz N et al. 2017 Cell, PMID 29033130
- Cancer: skin cutaneous melanoma (metastatic); setting: metastatic
- Platform: RNA-seq (GPL9052)
- Samples: 51; genes: 18130
- Responders (responder==1): 10; non-responders: 39; unlabelled: 2
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/
- https://raw.githubusercontent.com/riazn/bms038_analysis/master/data/bms038_clinical_data.csv

## Downloaded files (sha256)

- `GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz`: 987c7637e72ff5ad9c9630b0ee959064713cd8f57348d05003fb85a2f4cc43d5
- `Homo_sapiens.gene_info.gz`: 6e2f79000e79ecb430f2613b878d514302df7f81263b5db979010d032c3235a2
- `bms038_clinical_data.csv`: a70ae5a0b1156b6d85ec89f0d58f3f67538bfb969d53ec5d98232411d19869ab
- `GSE91061_series_matrix.txt.gz`: f3ebd04d019afab9f12af0b92cfcafa516bb2b170f5b45d22d5e330e2eebdca0

## Processing

- Supplementary GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz (FPKM, Entrez gene IDs) mapped to symbols with NCBI Homo_sapiens.gene_info; pre-treatment biopsies only; log2(FPKM + 1); symbols restricted to UniProt gene names -> 18130 genes.
- Clinical data (BOR, OS, PFS in days) from the authors' repository riazn/bms038_analysis (data/bms038_clinical_data.csv); *_SOR = 1 is censored, so event = 1 - SOR.
- responder: RECIST BOR CR/PR = 1, SD/PD = 0, NE = NaN. arm: NIV3-NAIVE (ipilimumab-naive) vs NIV3-PROG (progressed on ipilimumab); all received nivolumab (CheckMate 038).

## Arms

- NIV3-PROG: n=26 (responders 4/26 labelled)
- NIV3-NAIVE: n=25 (responders 6/23 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.196, PC2 0.110, PC3 0.079

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.351 | 0.531 | 0.626 | ok |
| arm | 0.475 | 0.537 | 0.529 | ok |

Built by `data/curated/trials/build_extra_trials.py GSE91061`.
