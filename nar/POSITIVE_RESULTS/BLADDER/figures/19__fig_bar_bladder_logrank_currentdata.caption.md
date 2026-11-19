# Bladder cisplatin survival (current data): every model against the original method

Bladder cisplatin (current data): organoids -> TCGA-BLCA treated, survival. Each bar is one model. The bar shows the main score on held-out data (C-index; 0.5 = chance), the black line is its 95% confidence interval
(bootstrap over patients). The red bar is the original method (paper amino-acid pathway); the orange bars are other published-style methods of the same study;
blue bars are the models added here (trained with settings chosen only inside the training data); grey bars are simple baselines and the median of 200 random pathway sets (with its 2.5-97.5% range as the line).
The text at the right gives the paired difference to the original method and its bootstrap p-value (same patients for both).
Original method: 0.572. Highest bar: Paper rank1 at 0.590. Source: tables/positives_all_models.tsv.
