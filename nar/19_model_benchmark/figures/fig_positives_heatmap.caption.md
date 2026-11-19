# All models on every positive result (heatmap)

Each column is one positive result and each row is one model. The number in a cell is the main score on held-out data: the C-index
(for survival: 0.5 = chance, 1 = perfect ordering of who lives longer) or the AUC (for mutation and response tasks: 0.5 = chance, 1 = perfect).
Purple-brown means above chance, orange means below chance. The first row is the original method of that result (its name is in the
bar-plot of the task and in the report). A plus or minus sign next to a number means the model differs from the original method with a
paired bootstrap p below 0.05 (plus = better, minus = worse); no sign means no clear difference. 'RandomSets median' is the median of
200 random pathway sets of the same size as the network shortlist. The two olaparib columns marked context are extra cohorts where the original signature was not
positive. Source table: tables/positives_all_models.tsv.
