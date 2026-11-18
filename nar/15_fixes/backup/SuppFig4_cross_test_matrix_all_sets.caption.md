# Supp. Fig. 4 (extended to every set): is each biomarker specific to its own cancer and drug?

Each row is one biomarker: the rank-1 pathway of one lab-model set and drug, with the sign of its organoid Ridge
weight. Each column is one group of drug-treated TCGA patients. Every cell applies the row's biomarker to the column's
patients (median split, log-rank p printed). Red frames = the biomarker in its own patients (the main result).
Darker blue = smaller p. The paper's cross-tests are: colon 5-FU patients scored with the bladder biomarker
(row "Bladder, paper's amino acid pathway", column "Colon organoids, paper data | 5-FU") and bladder cisplatin patients
scored with the colon biomarker (row "Colon organoids, paper data", column "Bladder organoids, paper data").
Paper: p = 0.16 and 0.86. Blank = pathway not scored in that cohort. Numbers incl. HR and 95% CI:
tables/SuppFig4_cross_test_matrix.tsv.
