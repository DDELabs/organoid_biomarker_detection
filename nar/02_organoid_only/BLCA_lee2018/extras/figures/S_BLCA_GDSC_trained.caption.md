# Models trained on cancer cell lines instead of organoids (paper Supp. Fig. 15c,d)

Same procedure, but the training data are GDSC2 cisplatin ln IC50 values of cell lines that also have CCLE (DepMap) RNA-seq.
Top: all such lines; bottom: bladder lines only. Each panel splits cisplatin-treated TCGA-BLCA patients (current GDC data) by one of the top-3 pathways.

- all GDSC lines, rank 1: lipoprotein metabolism; n = 79, p = 0.819, HR = 0.92 (0.43-1.93).
- all GDSC lines, rank 2: base free sugar phosphate removal via the single nucleotide replaceme…; n = 79, p = 0.0118, HR = 0.36 (0.16-0.83).
- all GDSC lines, rank 3: resolution of ap sites via the single nucleotide replacement pathway; n = 79, p = 0.0114, HR = 2.76 (1.22-6.27).
- bladder GDSC lines, rank 1: bmal1 clock npas2 activates circadian expression; n = 79, p = 0.213, HR = 1.61 (0.75-3.45).
- bladder GDSC lines, rank 2: digestion of dietary carbohydrate; n = 79, p = 0.694, HR = 1.16 (0.55-2.48).
- bladder GDSC lines, rank 3: gluconeogenesis; n = 79, p = 0.235, HR = 1.58 (0.74-3.38).

Paper: 831 lines -> p = 0.102, 0.0092, 0.84; 14 bladder lines -> p = 0.13, 0.23, 0.0092 (77 patients).
