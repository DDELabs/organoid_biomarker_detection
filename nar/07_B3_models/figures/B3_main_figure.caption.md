**Main figure: the paper's network method next to its comparison methods and modern models, on identical splits.**
(a) Summary heatmap (rows = models, columns = tasks; numbers as in the full benchmark heatmap). Lab-model task = held-out Spearman;
patient tasks = C-index minus 0.5 and log HR per SD in drug-treated TCGA patients; cell-line tasks = cell-blind and drug-blind CV.
Purple = expected direction. (b) Mean C-index with 95% CI over the 41 set x drug pairs. (c) Mean held-out Spearman
in lab models with 95% CI. (d-g) Kaplan-Meier curves for the paper's two headline pairs (colorectal organoids -> 5-FU-treated TCGA-COAD
patients; bladder organoids -> cisplatin-treated TCGA-BLCA), paper's single network pathway versus the modern model with the
highest mean C-index across pairs (GCN on STRING (genes)). Row colours: blue = paper, grey = paper's comparison methods,
red = modern, green = with alignment.
