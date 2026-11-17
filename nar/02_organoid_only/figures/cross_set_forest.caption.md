# Cross-set forest plot

One row per lab-model set x patient cohort x drug (23 rows). Each row uses the paper's single-pathway biomarker
(rank-1 pathway of a Ridge model trained on the lab models) and splits TCGA patients at the median prediction.
Filled blue circles: drug-treated patients; open grey squares: untreated patients of the same cancer with the same
biomarker and rule. Lines: 95% confidence intervals of the hazard ratio (Cox model). A hazard ratio above 1 means
predicted non-responders died sooner, which is the direction the method expects. "Shuffle" is the share of 1,000
analyses with the lab-model response shuffled that gave a p-value as small or smaller. * = 12-14 treated patients
(below the ~15 target, kept and flagged). Numbers: `tables/cross_set_summary.tsv`.
