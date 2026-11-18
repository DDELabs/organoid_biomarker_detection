# Does the BH3-only pathway tell us more than BCL2-family genes and proliferation?

Cell-line screens: GDSC1, GDSC2 (response = ln IC50) and CTRPv2 (response = AUC); lower = more sensitive. Navitoclax
and ABT-737 block BCL-2, BCL-XL and BCL-W. All values were rank-normalised; "partial r" is the correlation that is
left after removing the listed covariates (always including cancer type).

- **a** BH3-only score vs response, adjusted for cancer type only (grey), and then also for BCL2L1 (the gene for
  BCL-XL), for the three main anti-apoptotic genes, for a proliferation score (Hallmark E2F targets), and for all of
  them together (red). Bars = 95% CI.
- **b** Single markers against navitoclax (cancer type removed): each gene alone, the proliferation score, the six
  core BH3-only genes alone, and the BH3-only set without E2F1 and TFDP1.
- **c-e** How special is the BH3-only set? Grey = 1,000 random gene sets of the same size, blue = all 672 Reactome
  pathways, red line = BH3-only. Panel titles give the rank among pathways, the share of random sets that do as
  well, and the extra variance explained (R2) by adding the BH3-only score to a model with BCL2, BCL2L1, MCL1, E2F
  and cancer type.

- GDSC1 navitoclax (n = 926): added R2 = 0.010, F-test p = 0.00029; BH3-only beta -0.153 (-0.236 to -0.071); rank 1 of 672 pathways; random-set share 0.000
- GDSC2 navitoclax (n = 961): added R2 = 0.010, F-test p = 4.5e-05; BH3-only beta -0.158 (-0.234 to -0.083); rank 1 of 672 pathways; random-set share 0.000
- GDSC2 abt-737 (n = 951): added R2 = 0.022, F-test p = 7.6e-09; BH3-only beta -0.226 (-0.302 to -0.150); rank 1 of 672 pathways; random-set share 0.000
- CTRP navitoclax (n = 772): added R2 = 0.018, F-test p = 3.4e-05; BH3-only beta -0.189 (-0.279 to -0.100); rank 1 of 672 pathways; random-set share 0.000
- CTRP abt-737 (n = 706): added R2 = 0.019, F-test p = 3.2e-05; BH3-only beta -0.194 (-0.286 to -0.103); rank 1 of 672 pathways; random-set share 0.000
