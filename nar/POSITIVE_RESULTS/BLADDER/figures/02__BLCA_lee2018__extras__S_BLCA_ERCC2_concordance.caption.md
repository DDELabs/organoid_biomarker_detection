# ERCC2 mutation and predicted cisplatin resistance (paper Fig. 5b and Supp. Fig. 13)

ERCC2 mutations make bladder tumours more sensitive to cisplatin in clinical studies. Each panel scores every TCGA-BLCA
patient with mutation data (open GDC somatic mutations, Xena GDC hub; non-silent ERCC2 variants = mutant) using a model trained
on the bladder organoids' cisplatin response, then asks whether mutant tumours get LOWER predicted resistance
(one-sided Mann-Whitney U test, as in the paper).

- whole transcriptome (17717 features): mutant n = 36, WT n = 366, p = 0.0647.
- whole pathway (paper data) (488 features): mutant n = 36, WT n = 367, p = 0.0103.
- whole pathway (current data) (668 features): mutant n = 36, WT n = 366, p = 0.0222.
- neighbours up to 1 (317 features): mutant n = 36, WT n = 366, p = 0.0814.
- neighbours up to 2 (5569 features): mutant n = 36, WT n = 366, p = 0.0345.
- neighbours up to 3 (11784 features): mutant n = 36, WT n = 366, p = 0.0517.
- network rank 1 (paper data) (1 features): mutant n = 36, WT n = 367, p = 0.878.
- paper biomarker AA (paper data) (1 features): mutant n = 36, WT n = 367, p = 0.00302.
- network rank 1 (current data) (1 features): mutant n = 36, WT n = 366, p = 0.33.
- paper biomarker AA (current data) (1 features): mutant n = 36, WT n = 366, p = 0.000546.

Paper: ERCC2 mut n = 37, WT n = 367; own model p = 0.002; whole transcriptome 0.04, whole pathway 0.037, neighbours 0.017 / 0.021 / 0.028.
