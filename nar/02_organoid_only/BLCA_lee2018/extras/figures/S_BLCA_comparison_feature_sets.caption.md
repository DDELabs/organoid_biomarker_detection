# Other ways of choosing features, cisplatin-treated TCGA-BLCA (paper Fig. 3b-d, Supp. Fig. 4b, 8b, 9)

Each panel trains Ridge on the 9 bladder organoids with a different feature set (or uses one gene or the colorectal
biomarker), scores the cisplatin-treated patients and splits them at the median.

- whole transcriptome: n = 79, deaths = 28, p = 0.705, HR = 1.16 (0.55-2.44).
- whole pathway (paper data): n = 79, deaths = 20, p = 0.592, HR = 1.27 (0.53-3.07).
- whole pathway (current data): n = 79, deaths = 28, p = 0.863, HR = 0.94 (0.44-1.98).
- neighbours up to 1: n = 79, deaths = 28, p = 0.696, HR = 1.16 (0.55-2.46).
- neighbours up to 2: n = 79, deaths = 28, p = 0.942, HR = 0.97 (0.46-2.06).
- neighbours up to 3: n = 79, deaths = 28, p = 0.986, HR = 1.01 (0.47-2.13).
- network rank 1 (paper data): n = 79, deaths = 20, p = 0.0435, HR = 2.63 (0.99-6.96).
- network rank 1 (current data): n = 79, deaths = 28, p = 0.0672, HR = 2.03 (0.94-4.38).
- cross-test: colorectal 5-FU biomarker BH3 (paper data): n = 79, deaths = 20, p = 0.982, HR = 0.99 (0.41-2.41).
- cross-test: colorectal 5-FU biomarker BH3 (current data): n = 79, deaths = 28, p = 0.445, HR = 0.75 (0.35-1.58).
- target gene A2M (high = resistant): n = 79, deaths = 28, p = 0.74, HR = 0.88 (0.42-1.85).
- target gene ATOX1 (high = resistant): n = 79, deaths = 28, p = 0.365, HR = 0.71 (0.33-1.50).
- target gene MPG (high = resistant): n = 79, deaths = 28, p = 0.681, HR = 0.86 (0.41-1.80).
- target gene TF (high = resistant): n = 79, deaths = 28, p = 0.653, HR = 0.84 (0.40-1.78).

Paper (77 patients): whole transcriptome p = 0.79, whole pathway 0.82, first neighbours 0.37, 2nd/3rd neighbours 0.79/0.79, BH3 cross-test 0.86; target genes A2M 0.92, ATOX1 0.18, MPG 0.83, TF 0.48.
