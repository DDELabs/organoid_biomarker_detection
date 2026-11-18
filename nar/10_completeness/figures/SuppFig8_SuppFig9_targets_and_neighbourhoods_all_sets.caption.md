# Supp. Figs. 8 and 9 (extended to every set): simpler network features instead of the proximal pathway

For each lab-model set and drug (x axis) the treated TCGA patients are split at the median of different scores and
their survival compared (y = -log10 log-rank p; higher = stronger split; dashed line = p 0.05).
Red star: the paper's method (rank-1 network-proximal pathway).
Crosses (Supp. Fig. 8): each drug-target gene alone. Its direction comes from its correlation with drug response in the
lab models.
Circles, squares, triangles (Supp. Fig. 9): ridge regression trained on lab models using all genes within 1, 2 or 3
steps of the drug targets in the STRING network (first neighbours = paper Fig. 3d; 2 and 3 steps = Supp. Fig. 9).
Numbers: tables/SuppFig8_target_genes_alone.tsv and tables/SuppFig9_network_neighbourhood_1_2_3_steps.tsv.
Paper (bladder): neighbours 0.70 / 0.94 / 0.99 type values; colon: target TYMS alone did not split patients.
