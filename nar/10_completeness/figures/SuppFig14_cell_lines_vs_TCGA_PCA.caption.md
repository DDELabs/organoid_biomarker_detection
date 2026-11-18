# Supp. Fig. 14 (with the explained-variance numbers; extended to 8 cancers)

Each row is one cancer: Cell Model Passports cancer cell lines of that type (red; RNA-seq, log2 TPM) and the TCGA
tumours of the matching cohort (blue). Columns: gene level raw, gene level after scaling each gene within each data
set (z-score), pathway level raw (mean expression of the Reactome set's genes), pathway level scaled. Axis labels give
the % of variance explained by PC1 and PC2 (the numbers the paper reports). Raw: cell lines and tumours form two
separate clouds (a technical and biological batch). Scaled: the clouds overlap, which is why the method scales each
data set before training. Numbers: tables/SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv. The paper used ssGSEA NES for
the pathway level; here a simpler mean-expression score is used for the cell lines because no ssGSEA run of all
cell lines exists in the project (stated plainly; the gene-level panels are identical in kind to the paper's).
