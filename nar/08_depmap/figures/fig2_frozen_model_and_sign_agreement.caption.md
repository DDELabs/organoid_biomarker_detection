**Figure 2. The organoid model applied unchanged to colorectal cell lines.**
(c, top row) Each dot is a colorectal cell line. x = predicted 5-FU resistance from the frozen organoid signature
"robust_v1" (BH3-only activation with a negative weight, hyaluronan uptake and CO2 hydration with positive weights);
y = observed 5-FU response (higher = more resistant). The title gives Spearman rho and the 2.5-97.5% range of rho for
1,000 random signatures with the same number of Reactome pathways and random signs ("random-pathway null band").
(c, bottom left) Same test for all five frozen signatures from the earlier external-validation work. Grey bar =
random 95% band; coloured dot = the real signature. A dot above the band would mean better than random.
(d) Sign agreement for the 37 pathways the network step kept for 5-FU. Left: organoid weight (from the robust
stability model) against the cell-line Spearman rho in GDSC2 colorectal lines. Right: fraction of pathways whose
sign agrees, per screen, with a one-sided binomial test against 50%. Numbers: `tables/A_c_*.tsv`, `tables/A_d_*.tsv`.
