# Cell lines (DepMap/GDSC/CTRPv2): the BH3-only pathway score tracks sensitivity to navitoclax and ABT-737

## What was tested
The colon 5-FU model's most stable pathway, "Activation of BH3-only proteins", was scored (single-sample rank score) in cancer cell lines and correlated with drug response, with cancer type removed (lineage-adjusted partial Spearman). Pre-stated reading: a higher score means the cell is closer to the death threshold, so more sensitive; negative rho for response metrics where higher = more resistant.

## Data
GDSC1 (378 drugs), GDSC2 (286), CTRPv2 (460), PRISM (1,442); about 500-960 lines per drug; CCLE/Cell Model Passports RNA-seq.

## Results (BH q over all drugs of the screen)
| screen | drug | n | partial rho | p | BH q | drug rank in screen | pathway rank of 672 (after BCL2, BCL2L1, MCL1, E2F removed) |
|---|---|---|---|---|---|---|---|
| GDSC2 | ABT-737 | 951 | -0.246 | 4.9e-14 | 1.4e-11 | 1 | 1 (added R2 0.0216, F-test p 7.6e-09) |
| GDSC2 | navitoclax | 961 | -0.202 | 6.2e-10 | 8.8e-08 | 2 | 1 (R2 0.0104, p 4.5e-05) |
| GDSC1 | navitoclax | 926 | -0.206 | 6.6e-10 | 2.5e-07 | 2 | 1 (R2 0.0098, p 2.9e-04) |
| CTRPv2 | ABT-737 | 706 | -0.238 | 3.3e-10 | 1.5e-07 | 2 | 1 (R2 0.0189, p 3.2e-05) |
| CTRPv2 | navitoclax | 772 | -0.211 | 5.8e-09 | 1.3e-06 | 5 | 1 (R2 0.0178, p 3.4e-05) |
Share of 1,000 random same-size gene sets as strong: 0.000 in all five. Other BH3 mimetics: venetoclax GDSC2 -0.104 (p 0.0017), CTRPv2 -0.226 (p 1.2e-06); WEHI-539 GDSC2 beta -0.212 (p 1.1e-06); obatoclax does not follow (GDSC1 +0.046, GDSC2 +0.036), which fits its off-target action. Sources: `tables/08__A_e_bh3_vs_all_drugs.tsv`, `tables/08__A_g_bh3_vs_bh3_mimetics.tsv`, `tables/13__G7_regression_added_value.tsv`, `tables/13__G7_partial_correlations.tsv`, `tables/08__S4_permutation_random_controls.tsv`.

## How robust is it
Seen in three screens (GDSC1, GDSC2, CTRPv2) and across five drug-screen pairs; survives removal of BCL2, BCL2L1, MCL1 and proliferation (partial r -0.202 -> -0.134 in GDSC2 navitoclax); ranks first of 672 pathways every time. PRISM is weaker (navitoclax -0.112, p 0.017, q 0.42; ABT-737 p 0.19).

## Caveats
- The effect size is modest (partial r about -0.2, added R2 about 1-2%).
- GDSC1, GDSC2 and CTRPv2 share many cell lines, so the three screens are not fully independent. The mimetic link was found by scanning all drugs (nar/08 part e) and then examined in more detail (nar/13 G7); q values correct for the drug scan within each screen.
- It is not a 5-FU result: 5-FU sits in the middle of the ranking (pan-cancer rho 0.02), and within colorectal lines the score does not follow 5-FU (rho 0.05, 0.02, 0.26, 0.41). Within colorectal lines alone navitoclax shows no link, except PRISM at the opposite sign (rho +0.405, n 27, p 0.036; `../BORDERLINE.md`). The 19 colorectal organoids cannot test it (navitoclax tested only to 2 uM; rho 0.05).

## What I think
A real, replicated biological readout: the score behaves like an apoptotic-priming marker for BCL-2/BCL-XL inhibitors. It supports the biology of the pathway, not the organoid-to-patient claim.

## Where things are
`figures/`, `tables/`, `source_reports/08__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
