# Colon: the paper's 5-FU result, isogenic cells, cetuximab vs BRAF V600E, and the response cohort GSE62080

## What was tested
The Kong et al. colon results re-run from the paper's own files and from current data: 5-FU model (19 van de Wetering organoids; rank-1 pathway "Activation of BH3-only proteins"), isogenic HCT-8 cells, cetuximab model vs BRAF V600E tumours; plus external colorectal cohorts.

## Results
| id | test | n | effect | p | tier |
|---|---|---|---|---|---|
| CR1 | 5-FU model, 5-FU-treated TCGA-COAD, authors' 2019 files, median split | 114 patients, 8 deaths (0 vs 8) | HR not estimable; continuous HR 1.52 (0.94-2.46) | **0.014** | nominal (BH q 0.273 over 57; shuffle p 0.034; bootstrap share 0.0072) |
| CR2 | isogenic HCT-8 sensitive vs 5-FU-resistant, BH3-only higher in sensitive | 3 vs 9 samples | t 3.27 | **0.0085** | strong direction, tiny n (recomputed from raw GSE81005: p 0.039) |
| CR3 | cetuximab model vs BRAF V600E, authors' data | 49 mutant / 369 wild type | one-sided Mann-Whitney | 0.031 | nominal |
| CR3 | same, current data (colon) | 50 / 377 | | **0.00072** | strong |
| CR4 | same, colon + rectum, current data | 51 / 522 | | **0.00048** | strong |
| CR5 | benchmark: L1 on shortlist vs original rank-1, AUC | 573 / 51 V600E | 0.704 vs 0.640, diff +0.064 (0.011-0.113) | paired 0.018, q 0.036 | q < 0.10 |
| CR6 | response, GSE62080 (metastatic FOLFIRI), trained on other colorectal cohorts | 21 / 9 responders | AUC 0.806 (0.605-1.000) | 0.021, q 0.056 | nominal, tiny n |
| CR7 | Gao CRC PDX 5-FU top-10 signature in GSE39582, OS, x 5-FU interaction | 139 / 40 | HR 1.46 (1.01-2.10) | 0.042 | nominal (q 0.650 over 196) |
| CR8 | organoid irinotecan top-10 signature, GSE72970 FOLFIRI vs FOLFOX PFS interaction | 84 / 80 | HR 1.65 (1.02-2.67) | 0.042 | nominal (q 0.650) |
Sources: `tables/02__cross_set_summary.tsv`, `tables/02__Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv`, `tables/02__Fig5a_cetuximab_BRAF_V600E.tsv`, `tables/10__Fig5_mutation_concordance_other_cancers.tsv`, `tables/13__G8_BRAF_null_and_adjustment.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/18__part1_v2.tsv`, `tables/05__B1_external_results.tsv`, `tables/05__A3_master_results.tsv`.

## How robust is it
- **5-FU (CR1) is fragile.** It reproduces exactly on the authors' 2019 inputs, but p moves to 0.866 on current GDC data (n 129, 25 deaths): the main drivers are the new organoid scores, the new patient expression pipeline and longer follow-up (nar/13 G2). With the 2019 setup but current survival the frozen BH3-only pathway gives p 0.002 (114 patients, 20 deaths), so the direction is not always lost. With the rank score instead of ssGSEA p becomes 0.32 (nar/14 G18). 94 of the 114 "5-FU" patients also received oxaliplatin; the split is almost entirely a FOLFOX group. In independent 5-FU cohorts the BH3-only pathway does not carry over (GSE39582 RFS HR 0.95, p 0.70). No cell-line support: BH3-only vs 5-FU in colorectal lines rho 0.05 (GDSC1, n 44), 0.02 (GDSC2), 0.26 (CTRPv2), 0.41 (PRISM).
- **Cetuximab/BRAF (CR3-CR5) is strong but not specific.** After adjusting for MSI, side, stage and TMB p is 0.026 (current) / 0.093 (authors); within MSS tumours p 0.007 (10 mutants). 17.7% of random pathway sets of the same size do as well (current data), so it is not special to the network. Refits leaving k organoids out still give p < 0.05 in 100% (k=1), 93%, 94%, 87% (k=4). KRAS mutations are not picked up (p 0.88).
- **Isogenic HCT-8 (CR2)**: same direction in both versions, 3 vs 9 samples.
- The paper's colon 5-FU result ranks q 0.27 across the 57 organoid-only analyses; none of the 57 passes BH.

## What I think
Two things hold up: the cetuximab/BRAF concordance (strong p, but partly a tumour-subtype signal) and the isogenic cell direction. The headline colon 5-FU survival split reproduces on the paper's data but disappears on current data, so I would present it as a reproduction with that limit stated. GSE62080 has 21 patients; treat it as a lead.

## Where things are
`figures/`, `tables/`. Checksums in `../COPY_MANIFEST.tsv`. Rectal cancer is in `../RECTAL/`.
