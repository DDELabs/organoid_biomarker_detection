# BH3-only pathway score and drug response in the 19 colorectal organoids (van de Wetering 2015)

**What is shown.** Each dot is one organoid. The x-axis is the score of the Reactome pathway "Activation of BH3-only
proteins" (ssGSEA, recomputed here from the published expression). The y-axis is the drug response as area under
the dose-response curve (AUC): lower AUC means the organoid is more sensitive. In cell lines, a higher BH3-only score
went with more sensitivity to navitoclax (a negative rho).

- **a** navitoclax (blocks BCL-2, BCL-XL and BCL-W): rho = 0.05 (95% CI -0.44 to 0.51), p = 0.850.
- **b** obatoclax (pan-BCL-2-family blocker with other off-target effects): rho = -0.38, p = 0.113.
- **c** 5-FU, for contrast: rho = -0.41, p = 0.083.
- **d** All drug x score x metric combinations with 95% confidence intervals (Fisher z). Filled dots use our
  recomputed score, open dots the authors' own organoid scores (NES). "ln IC50" is the log of the half-maximal
  inhibitory concentration (lower = more sensitive).
- **e** Random-pathway check: 2000 random gene sets of the same size (8 genes); grey =
  their correlations with navitoclax, blue line = the BH3-only set. The share of random sets with a correlation at
  least as negative is 0.583; among all 533 Reactome pathways it is 0.612.
- **f** The BH3-only score correlated with every drug in the screen; bars are sorted from most negative to most
  positive. Numbers in brackets are ranks.

Note on range: navitoclax was tested up to 2 uM only; organoid AUCs lie between 0.94 and 1.00, so almost no
organoid was killed in the tested range (see tables/G3_response_range_organoids.tsv). This leaves little
variation for any marker to explain.

With 19 organoids, a 95% interval spans roughly +/-0.45 around the estimate, so only large correlations can be told
apart from zero.
