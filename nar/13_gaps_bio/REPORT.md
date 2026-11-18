# Biology gap analyses (G2, G3, G5, G7, G8, G13, G14)

Every number in this report was copied by `code/make_report.py` from a file in `tables/`. Figures are in `figures/` (PNG + PDF + a plain-language `.caption.md`). "What I think:" lines are my reading; the decision is yours.

Words used: **rho** = Spearman rank correlation (-1 to +1). **HR** = hazard ratio (above 1 = shorter survival). **CI** = 95% confidence interval. **Frozen** = the pathway and direction fixed from the 2019 paper, with no re-selection. **AUC** = area under the dose-response curve (lower = more sensitive).

## G2. Colon 5-FU: why p moved from 0.014 (2019) to 0.866 (now)

Files: `tables/G2_decomposition_colon_5FU.tsv`, `tables/G2_data_versions.tsv`, `tables/G2_gene_universe_sizes.tsv`. Figure: `figures/G2_decomposition_colon_5FU.png/.pdf`.

**a. One change at a time, adding current inputs step by step (orange = model re-picks its top pathway; blue = the paper's BH3-only pathway kept fixed).**

| step | what was changed | n | deaths | re-picked top pathway | log-rank p | frozen BH3-only log-rank p | frozen HR per SD (95% CI) |
|---|---|---|---|---|---|---|---|
| 0 | authors' setup (2019) | 114 | 8 | activation of bh3 only proteins | 0.014 | 0.014 | 1.52 (0.94 to 2.46) |
| 1 | + only patients with current STAR data | 114 | 8 | activation of bh3 only proteins | 0.014 | 0.014 | 1.52 (0.94 to 2.46) |
| 2 | + organoid scores recomputed | 114 | 8 | endosomal sorting complex required for transport escrt | 0.396 | 0.014 | 1.52 (0.94 to 2.46) |
| 3 | + patient scores from current STAR FPKM-UQ | 114 | 8 | endosomal sorting complex required for transport escrt | 0.839 | 0.319 | 1.63 (0.84 to 3.17) |
| 4 | + shared gene universe (02 pipeline) | 117 | 9 | endosomal sorting complex required for transport escrt | 0.228 | 0.596 | 1.21 (0.60 to 2.43) |
| 5 | + organoid response from curated table | 117 | 9 | purine ribonucleoside monophosphate biosynthesis | 0.414 | 0.596 | 1.21 (0.60 to 2.43) |
| 6 | + current drug-record parser | 129 | 11 | purine ribonucleoside monophosphate biosynthesis | 0.413 | 0.774 | 1.17 (0.63 to 2.17) |
| 7 | + TCGA-CDR survival = current setup | 129 | 25 | purine ribonucleoside monophosphate biosynthesis | 0.866 | 0.620 | 1.20 (0.82 to 1.75) |

The re-picked pathway starts as BH3-only (p = 0.014) and ends as purine ribonucleoside monophosphate biosynthesis (p = 0.866). The frozen BH3-only pathway goes from p = 0.014 to p = 0.620, HR per SD 1.20 (0.82 to 1.75).

**b. Only one of the authors' inputs swapped in (everything else as in 2019).**

| swapped input | n | deaths | re-picked top pathway | log-rank p (re-picked) | frozen BH3-only p |
|---|---|---|---|---|---|
| organoid scores recomputed | 114 | 8 | endosomal sorting complex required for transport escrt | 0.396 | 0.014 |
| patient scores from current STAR | 117 | 9 | activation of bh3 only proteins | 0.512 | 0.512 |
| shared gene universe, both sides | 117 | 9 | endosomal sorting complex required for transport escrt | 0.228 | 0.596 |
| organoid response from curated table | 114 | 8 | activation of bh3 only proteins | 0.014 | 0.014 |
| current drug-record parser | 125 | 10 | activation of bh3 only proteins | 0.091 | 0.091 |
| TCGA-CDR survival | 114 | 20 | activation of bh3 only proteins | 0.002 | 0.002 |
| TCGA-CDR survival capped at 2019 follow-up | 114 | 8 | activation of bh3 only proteins | 0.014 | 0.014 |
| only patients with current STAR data | 114 | 8 | activation of bh3 only proteins | 0.014 | 0.014 |

**c. Starting from the current setup, put back one of the authors' inputs.**

| authors' input put back | n | deaths | re-picked top pathway | log-rank p (re-picked) | frozen BH3-only p |
|---|---|---|---|---|---|
| authors' organoid NES | 129 | 25 | activation of bh3 only proteins | 0.620 | 0.620 |
| authors' patient NES | 125 | 24 | purine ribonucleoside monophosphate biosynthesis | 0.597 | 0.005 |
| authors' organoid + patient NES | 125 | 24 | activation of bh3 only proteins | 0.005 | 0.005 |
| authors' organoid response file | 129 | 25 | endosomal sorting complex required for transport escrt | 0.161 | 0.620 |
| authors' drug-record parser | 117 | 21 | purine ribonucleoside monophosphate biosynthesis | 0.805 | 0.634 |
| authors' 2019 survival | 129 | 11 | purine ribonucleoside monophosphate biosynthesis | 0.413 | 0.774 |
| TCGA-CDR capped at 2019 follow-up | 129 | 11 | purine ribonucleoside monophosphate biosynthesis | 0.413 | 0.774 |
| current Reactome release (both sides) | 129 | 25 | unwinding of dna | 0.340 | 0.089 |

**d. Frozen BH3-only pathway alone, by regimen (current data).** `5-FU + oxaliplatin` = FOLFOX-type.

| patient group | n | deaths | log-rank p | HR (95% CI), predicted non-responder vs responder | HR per SD (95% CI) |
|---|---|---|---|---|---|
| all 5-FU patients | 129 | 25 | 0.620 | 1.22 (0.56 to 2.68) | 1.20 (0.82 to 1.75) |
| 5-FU + oxaliplatin (FOLFOX-type) | 106 | 18 | 0.764 | 0.87 (0.34 to 2.20) | 1.22 (0.77 to 1.92) |
| 5-FU without oxaliplatin | 23 | 7 | 0.135 | 4.59 (0.51 to 41.20) | 1.50 (0.65 to 3.48) |
| 5-FU alone | 12 | 0 | 1.000 | NA (0 deaths) | NA |

Same, with current follow-up cut at the 2019 date: all 5-FU patients n = 129, deaths = 11, frozen log-rank p = 0.774, HR per SD 1.17 (0.63 to 2.17).

**What I think:** The patient-side scores are the main driver. With the paper's pathway fixed, putting the authors' 2019 patient scores back into the current setup gives p = 0.005 (n = 125, 24 deaths), while putting back only the authors' organoid scores gives p = 0.620. Going the other way, swapping only the current STAR patient scores into the 2019 setup moves the fixed pathway to p = 0.512. The organoid side matters for which pathway the model picks: recomputing organoid scores alone moves the re-picked pathway from 0.014 to 0.396 (step 2). Survival data and follow-up are not the cause: with the authors' scores, the current TCGA-CDR survival gives p = 0.002 (20 deaths), and the authors' scores plus current survival give p = 0.005 (HR 3.5 (1.4 to 8.8), n = 125, 24 deaths). The 2019 result rested on 8 deaths. In short, the 2019 patient expression scores and today's STAR scores rank the same 5-FU patients differently for this pathway.

**Not done (G2, bladder).** The request also asked for the same frozen-pathway report on bladder. The script `g2_decompose.py` is colon-only and I did not write a bladder version. The bladder ERCC2 side is covered in G8.

## G3. BH3-only score vs navitoclax in the 19 colorectal organoids

Files: `tables/G3_bh3_vs_drugs_organoids.tsv`, `G3_single_genes_vs_bh3_mimetics_organoids.tsv`, `G3_bh3_vs_all_drugs_organoids_AUC.tsv`, `G3_response_range_organoids.tsv`, `G3_data_versions.tsv`. Figure: `figures/G3_bh3_vs_bh3_mimetics_organoids.png/.pdf`.

Spearman rho between the BH3-only score and drug response (AUC). A negative rho means a higher score goes with more sensitivity. Null = share of random gene sets (same size, 2,000 sets) and of all Reactome pathways with rho at least as negative.

| drug | score used | n | rho (95% CI) | p | share of random sets as negative | share of all pathways as negative |
|---|---|---|---|---|---|---|
| navitoclax | recomputed | 19 | 0.05 (-0.44 to 0.51) | 0.850 | 0.583 | 0.612 (of 533) |
| navitoclax | authors' NES | 19 | -0.04 (-0.51 to 0.45) | 0.872 | NA | 0.503 (of 612) |
| obatoclax | recomputed | 19 | -0.38 (-0.72 to 0.12) | 0.113 | 0.059 | 0.098 (of 533) |
| obatoclax | authors' NES | 19 | -0.20 (-0.62 to 0.30) | 0.403 | NA | 0.327 (of 612) |
| fluorouracil | recomputed | 19 | -0.41 (-0.74 to 0.09) | 0.083 | 0.042 | 0.069 (of 533) |
| fluorouracil | authors' NES | 19 | -0.57 (-0.82 to -0.13) | 0.011 | NA | 0.015 (of 612) |
| oxaliplatin | recomputed | 19 | -0.33 (-0.70 to 0.17) | 0.162 | 0.088 | 0.116 (of 533) |
| oxaliplatin | authors' NES | 19 | -0.34 (-0.70 to 0.17) | 0.159 | NA | 0.111 (of 612) |

Tested range: navitoclax was tested up to 2 uM; organoid AUC runs from 0.944 to 0.996, and 19 of 19 IC50 values lie above the top dose. For 5-FU, AUC runs from 0.623 to 0.995 (SD 0.109) vs navitoclax SD 0.012.

Single genes vs navitoclax (rho, 95% CI, p): BCL2 0.09 (-0.41 to 0.54), p = 0.729; BCL2L1 -0.18 (-0.61 to 0.32), p = 0.452; MCL1 -0.17 (-0.60 to 0.33), p = 0.486; BAK1 0.19 (-0.31 to 0.61), p = 0.424; PMAIP1 0.01 (-0.47 to 0.49), p = 0.954; BID -0.04 (-0.51 to 0.44), p = 0.864; BMF -0.16 (-0.59 to 0.34), p = 0.513; E2F1 -0.25 (-0.65 to 0.26), p = 0.298; TP53 -0.32 (-0.69 to 0.19), p = 0.183.

Across the whole screen (83 drugs, AUC), the BH3-only score ranks navitoclax 70, 5-FU 16 and obatoclax 19 from the most negative correlation.

**What I think:** In these 19 organoids the navitoclax link seen in cell lines is not visible (rho = 0.05, CI -0.44 to 0.51, p = 0.850). The data cannot settle it either way: navitoclax was only tested up to 2 uM, nearly every organoid stayed at AUC near 1, and the interval is wide. The BH3-only score points the 5-FU way in the same organoids (recomputed rho = -0.41, p = 0.083; authors' NES rho = -0.57, CI -0.82 to -0.13, p = 0.011). Obatoclax, the planned control, shows rho = -0.38, so it did not behave as a flat control here. A drug panel with a higher top dose is the clean way to test the apoptotic-priming claim in organoids.

## G5. FOLFOX question: is the 5-FU marker a regimen (oxaliplatin) marker?

Files: `tables/G5_cell_lines_bh3_vs_oxaliplatin.tsv`, `G5_folfox_patients.tsv`, `G5_data_versions.tsv`, `G3_bh3_vs_drugs_organoids.tsv`. Figure: `figures/G5_folfox_question.png/.pdf`.

**a. 19 organoids (AUC, recomputed score).** 
- oxaliplatin: rho = -0.33 (-0.70 to 0.17), p = 0.162
- fluorouracil: rho = -0.41 (-0.74 to 0.09), p = 0.083

**b. Cancer cell lines.** rho of BH3-only score vs response (colorectal lines only, and all lines with cancer type removed).

| screen | drug | lines | n | rho (95% CI) | p |
|---|---|---|---|---|---|
| GDSC1 | oxaliplatin | colorectal | 0 | no data in table | NA |
| GDSC1 | fluorouracil | colorectal | 44 | 0.05 (-0.26 to 0.36) | 0.730 |
| GDSC1 | fluorouracil | pan-cancer, lineage-adjusted | 903 | 0.02 (-0.05 to 0.08) | 0.605 |
| GDSC2 | oxaliplatin | colorectal | 48 | -0.02 (-0.32 to 0.28) | 0.903 |
| GDSC2 | oxaliplatin | pan-cancer, lineage-adjusted | 961 | -0.05 (-0.12 to 0.01) | 0.114 |
| GDSC2 | fluorouracil | colorectal | 48 | 0.02 (-0.29 to 0.31) | 0.916 |
| GDSC2 | fluorouracil | pan-cancer, lineage-adjusted | 962 | 0.02 (-0.04 to 0.08) | 0.515 |
| CTRP | oxaliplatin | colorectal | 44 | 0.08 (-0.24 to 0.38) | 0.621 |
| CTRP | oxaliplatin | pan-cancer, lineage-adjusted | 768 | -0.12 (-0.19 to -0.05) | 0.001 |
| CTRP | fluorouracil | colorectal | 43 | 0.26 (-0.07 to 0.53) | 0.096 |
| CTRP | fluorouracil | pan-cancer, lineage-adjusted | 775 | 0.01 (-0.06 to 0.08) | 0.763 |
| PRISM | oxaliplatin | colorectal | 26 | 0.37 (-0.04 to 0.68) | 0.060 |
| PRISM | oxaliplatin | pan-cancer, lineage-adjusted | 448 | 0.03 (-0.06 to 0.12) | 0.504 |
| PRISM | fluorouracil | colorectal | 22 | 0.40 (-0.05 to 0.72) | 0.062 |
| PRISM | fluorouracil | pan-cancer, lineage-adjusted | 429 | 0.07 (-0.02 to 0.17) | 0.139 |

**c. TCGA colon 5-FU patients, frozen BH3-only score (HR per SD; above 1 = the paper's direction).**

| data | group | n | events | HR per SD (95% CI) | p |
|---|---|---|---|---|---|
| authors' 2019 files | all 5-FU patients | 114 | 8 | 1.52 (0.94 to 2.46) | 0.091 |
| authors' 2019 files | 5-FU + oxaliplatin (FOLFOX-type) | 94 | 6 | 1.67 (0.94 to 2.98) | 0.082 |
| authors' 2019 files | 5-FU without oxaliplatin | 20 | 2 | 2.87 (0.25 to 33.19) | 0.399 |
| current data | all 5-FU patients | 129 | 25 | 1.20 (0.82 to 1.75) | 0.357 |
| current data | 5-FU + oxaliplatin (FOLFOX-type) | 106 | 18 | 1.22 (0.77 to 1.92) | 0.398 |
| current data | 5-FU without oxaliplatin | 23 | 7 | 1.50 (0.65 to 3.48) | 0.341 |

Score x oxaliplatin interaction (current data): HR 0.88 (0.33 to 2.31), p = 0.789.

**d. GSE39582 (stage II-III), frozen BH3-only score by adjuvant regimen.**

| endpoint | group | n | events | HR per SD (95% CI) | p | interaction HR (95% CI), p |
|---|---|---|---|---|---|---|
| RFS | 5FU alone | 76 | 28 | 0.87 (0.59 to 1.28) | 0.481 | |
| RFS | FUFOL (5-FU + leucovorin) | 50 | 13 | 1.12 (0.66 to 1.88) | 0.677 | |
| RFS | 5FU or FUFOL (no oxaliplatin) | 126 | 41 | 0.96 (0.70 to 1.31) | 0.788 | |
| RFS | FOLFOX | 12 | 7 | 0.69 (0.30 to 1.63) | 0.400 | |
| RFS | no chemotherapy | 258 | 62 | 1.20 (0.93 to 1.54) | 0.153 | |
| RFS | FOLFOX vs 5FU/FUFOL (interaction) | 138 | 48 | | | 0.59 (0.15 to 2.29), p = 0.443 |
| RFS | 5FU/FUFOL vs no chemotherapy (interaction) | 384 | 103 | | | 0.79 (0.53 to 1.19), p = 0.262 |
| OS | 5FU alone | 76 | 19 | 1.00 (0.62 to 1.59) | 0.987 | |
| OS | FUFOL (5-FU + leucovorin) | 51 | 16 | 1.18 (0.73 to 1.89) | 0.495 | |
| OS | 5FU or FUFOL (no oxaliplatin) | 127 | 35 | 1.08 (0.77 to 1.51) | 0.666 | |
| OS | FOLFOX | 12 | 5 | 1.15 (0.49 to 2.66) | 0.749 | |
| OS | no chemotherapy | 262 | 89 | 0.98 (0.79 to 1.20) | 0.814 | |
| OS | FOLFOX vs 5FU/FUFOL (interaction) | 139 | 40 | | | 1.15 (0.31 to 4.27), p = 0.835 |
| OS | 5FU/FUFOL vs no chemotherapy (interaction) | 389 | 124 | | | 1.11 (0.75 to 1.64), p = 0.618 |

**What I think:** The score points the same way for oxaliplatin as for 5-FU in the 19 organoids (oxaliplatin rho = -0.33, p = 0.162), so the organoid data cannot tell the two drugs apart. In colorectal cell lines the numbers are close to zero or positive (higher score with more resistance) (for example CTRPv2 oxaliplatin rho = 0.08, n = 44); across all CTRPv2 lines with cancer type removed rho = -0.12 (p = 0.001). In TCGA, the HR per SD is 1.22 with oxaliplatin and 1.50 without, and the interaction p = 0.789: the groups are too small (23 patients without oxaliplatin) to separate a 5-FU effect from a FOLFOX effect. In GSE39582 (stage II-III) the HR per SD is 0.96 for relapse-free survival in 5-FU/FUFOL patients (n = 126) and 0.69 in FOLFOX (n = 12); the FOLFOX vs 5-FU interaction p = 0.443. Most of the 5-FU patients in TCGA got FOLFOX, so the paper's "5-FU" claim is in practice a FOLFOX claim.

## G7. Does the BH3-only pathway add beyond BCL2, BCL2L1, MCL1 and proliferation?

Files: `tables/G7_regression_added_value.tsv`, `G7_partial_correlations.tsv`. Figure: `figures/G7_bh3_beyond_simple_markers.png/.pdf`. Cell-line screens GDSC1, GDSC2 (ln IC50), CTRPv2 (AUC); all variables rank-normalised; cancer type always removed.

| screen | drug | n | partial r, cancer type only (95% CI) | partial r, + BCL2, BCL2L1, MCL1 + E2F (95% CI) | extra R2 from BH3-only | F-test p | rank among pathways | share of random sets as strong |
|---|---|---|---|---|---|---|---|---|
| GDSC1 | navitoclax | 926 | -0.206 (-0.268 to -0.142) | -0.122 (-0.186 to -0.056) | 0.0098 | 0.00029 | 1 of 672 | 0.000 |
| GDSC2 | navitoclax | 961 | -0.202 (-0.263 to -0.139) | -0.134 (-0.197 to -0.070) | 0.0104 | 4.5e-05 | 1 of 672 | 0.000 |
| GDSC2 | abt-737 | 951 | -0.246 (-0.306 to -0.184) | -0.190 (-0.252 to -0.127) | 0.0216 | 7.6e-09 | 1 of 672 | 0.000 |
| CTRP | navitoclax | 772 | -0.211 (-0.278 to -0.141) | -0.151 (-0.221 to -0.080) | 0.0178 | 3.4e-05 | 1 of 672 | 0.000 |
| CTRP | abt-737 | 706 | -0.238 (-0.307 to -0.165) | -0.159 (-0.232 to -0.085) | 0.0189 | 3.2e-05 | 1 of 672 | 0.000 |

Other BH3-mimetics (BH3-only beta in the full model, 95% CI, extra R2):
- GDSC1 venetoclax (n = 901): beta -0.105 (-0.185 to -0.024), extra R2 0.0047, p = 0.011
- GDSC1 obatoclax (n = 858): beta 0.068 (-0.029 to 0.165), extra R2 0.0019, p = 0.171
- GDSC2 venetoclax (n = 952): beta -0.125 (-0.192 to -0.058), extra R2 0.0066, p = 0.00027
- GDSC2 wehi-539 (n = 956): beta -0.212 (-0.296 to -0.127), extra R2 0.0188, p = 1.1e-06
- GDSC2 azd5991 (n = 711): beta -0.079 (-0.159 to 0.002), extra R2 0.0026, p = 0.055
- GDSC2 obatoclax (n = 782): beta -0.023 (-0.118 to 0.073), extra R2 0.0002, p = 0.640
- CTRP venetoclax (n = 473): beta -0.198 (-0.310 to -0.086), extra R2 0.0197, p = 0.00054
- CTRP obatoclax (n = 784): beta -0.091 (-0.179 to -0.004), extra R2 0.0041, p = 0.040

Single markers vs navitoclax (partial r with cancer type removed, GDSC2): BH3-only pathway -0.202 (-0.263 to -0.139); BCL2 -0.179 (-0.241 to -0.116); BCL2L1 -0.035 (-0.100 to 0.030); MCL1 0.248 (0.187 to 0.308); E2F -0.072 (-0.136 to -0.007); BH3core -0.163 (-0.225 to -0.099); BH3_noE2F -0.189 (-0.250 to -0.126).

**What I think:** The BH3-only score keeps most of its link after removing BCL2L1, BCL2, MCL1 and proliferation, so it is not just BCL-XL or E2F. For example in GDSC2 navitoclax the partial r goes from -0.202 to -0.134. It adds a small but consistent slice of variance (extra R2 0.010 in GDSC2 navitoclax) and ranks 1 of 672 pathways in all three screens. Removing E2F1 and TFDP1 from the set leaves the link intact (see `BH3_noE2F` in the partial-correlation table). The effect is modest in size, so I would describe it as an added signal, not a replacement for BCL2-family genes.

## G8. Null and confounder control for the mutation checks

Files: `tables/G8_ERCC2_null_and_adjustment.tsv`, `G8_BRAF_null_and_adjustment.tsv`, and their `_data_versions.tsv`. Figure: `figures/G8_mutation_concordance_null_and_adjustment.png/.pdf`. One-sided p values (expected direction).

**ERCC2 (bladder, cisplatin).**

| data | model | mutant / wild type | p | share of all pathways (both directions) as strong | share of random proximal-sized sets as strong | p after TMB adjustment (OLS) | OR per SD after TMB (95% CI) |
|---|---|---|---|---|---|---|---|
| paper | paper biomarker (amino acid synthesis) | 36 / 367 | 0.003 | 0.017 | 0.063 | 0.024 | 0.63 (0.42 to 0.93) |
| paper | rank-1 pathway | 36 / 367 | 0.878 | 0.820 | 0.987 | 0.922 | 1.41 (0.96 to 2.06) |
| current | paper biomarker (amino acid synthesis) | 36 / 366 | 0.00055 | 0.002 | 0.000 | 0.008 | 0.54 (0.36 to 0.82) |
| current | rank-1 pathway | 36 / 366 | 0.330 | 0.360 | 0.602 | 0.276 | 0.90 (0.63 to 1.28) |

**BRAF V600E (colon, cetuximab).**

| data | mutant / wild type | p | share of pathways as strong | share of random sets as strong | p adj. MSI | p adj. side | p adj. TMB | p adj. all (MSI, side, stage, TMB) | MSS tumours only (mutants, p) |
|---|---|---|---|---|---|---|---|---|---|
| authors | 49 / 369 | 0.031 | 0.272 | 0.352 | 0.087 | 0.041 | 0.011 | 0.093 | 9, 0.034 |
| current | 50 / 377 | 0.00072 | 0.163 | 0.177 | 0.015 | 0.00082 | 0.001 | 0.026 | 10, 0.007 |

**What I think:** ERCC2 looks specific: only 0.002 of all pathways and 0.000 of random sets separate mutants as well on current data (paper data: 0.017 and 0.063), and it holds after TMB adjustment (p = 0.008 current, 0.024 paper). BRAF V600E is less specific: 0.163 of pathways and 0.177 of random sets do as well on current data (0.272 and 0.352 on the authors' data). After adjusting for MSI, side, stage and TMB the current-data p is 0.026 (authors' data 0.093), and within MSS tumours p = 0.007 with 10 mutants. So the BRAF check is best described as a consistent direction whose size partly tracks tumour subtype, while the ERCC2 check stands on its own.

## G13. The 'HIV' frozen pathways: what genes are in them?

Files: `tables/G13_frozen_pathway_composition.tsv`, `G13_most_overlapping_pathways.tsv`, `G13_purity_immune_correlation.tsv`, `G13_genes_*.txt` (gene lists), `G13_data_versions.tsv`. Figure: `figures/G13_odd_pathway_composition.png/.pdf`.

| pathway | genes | in proteasome set | in nuclear export/pore set | in immune-system set | in adaptive-immune set | in cell-cycle set | host-factor genes (PPIA, BANF1, XRCC*...) |
|---|---|---|---|---|---|---|---|
| apobec3g mediated resistance to hiv1 infection | 12 | 0 (0.00) | 0 (0.00) | 0 (0.00) | 0 (0.00) | 0 (0.00) | 5 |
| hiv infection | 207 | 40 (0.19) | 27 (0.13) | 104 (0.50) | 76 (0.37) | 64 (0.31) | 9 |

Closest other Reactome pathways (Jaccard overlap of gene lists):
- apobec3g mediated resistance to hiv1 infection: integration of provirus (Jaccard 0.65, 11 shared genes); early phase of hiv life cycle (Jaccard 0.50, 11 shared genes); binding and entry of hiv virion (Jaccard 0.47, 7 shared genes)
- hiv infection: host interactions of hiv factors (Jaccard 0.64, 132 shared genes); hiv life cycle (Jaccard 0.60, 125 shared genes); late phase of hiv life cycle (Jaccard 0.50, 104 shared genes)

Correlation with tumour purity and an adaptive-immune score, median over 32 TCGA cancer types (10th to 90th percentile):

| pathway | vs purity | vs adaptive immune score |
|---|---|---|
| apobec3g mediated resistance to hiv1 infection | -0.23 (-0.41 to -0.01) | 0.48 (0.27 to 0.67) |
| hiv infection | -0.06 (-0.32 to 0.25) | 0.48 (0.22 to 0.63) |

For all 34 distinct frozen pathways the median (over pathways) of the per-pathway median rho with the immune score is 0.07.

**What I think:** These two sets are human genes that viruses use, not virus genes. The 12-gene APOBEC3G set is built from HIV life-cycle pathways (closest: integration of provirus, 11 shared genes). The 207-gene HIV infection set is 19% proteasome, 13% nuclear pore/export, and 50% immune system genes. Both scores rise with an adaptive-immune score (median rho 0.48 and 0.48), and the APOBEC3G set also falls with purity (median rho -0.23). A reader can therefore take these as proteasome, nuclear-transport and immune-infiltration biology under a virus-themed name. Note: the frozen list has 35 rows but 34 distinct pathways; the table covers the distinct ones.

## G14. Supp. Fig. 15 for pancreas, liver, ovary and sarcoma (cell-line-trained models)

Files: `tables/G14_SuppFig15_cell_line_models.tsv`, `G14_data_versions.tsv`. Figure: `figures/G14_SuppFig15_cell_line_models.png/.pdf`. Models trained on GDSC2 cell lines (all lines, or lines of the matching cancer type), applied to treated TCGA patients; each of the top 3 pathways splits patients at the median. HR above 1 = predicted non-responders die sooner (expected direction).

| setting | training | lines | patients | deaths | rank | pathway | log-rank p | HR (95% CI) |
|---|---|---|---|---|---|---|---|---|
| liver | sorafenib | all GDSC2 lines | 952 | 30 | 15 | 1 | fgfr1 ligand binding and activation | 0.045 | 3.49 (0.95 to 12.84) |
| liver | sorafenib | all GDSC2 lines | 952 | 30 | 15 | 2 | fgfr2c ligand binding and activation | 0.138 | 0.43 (0.14 to 1.35) |
| liver | sorafenib | all GDSC2 lines | 952 | 30 | 15 | 3 | signaling by activated point mutants of fgfr1 | 0.255 | 0.52 (0.16 to 1.63) |
| liver | sorafenib | lineage-matched lines | 15 | 30 | 15 | 1 | negative regulation of the pi3k akt network | 0.084 | 0.35 (0.10 to 1.21) |
| liver | sorafenib | lineage-matched lines | 15 | 30 | 15 | 2 | platelet adhesion to exposed collagen | 0.017 | 0.22 (0.06 to 0.85) |
| liver | sorafenib | lineage-matched lines | 15 | 30 | 15 | 3 | insulin receptor recycling | 0.161 | 2.19 (0.71 to 6.75) |
| ovary | paclitaxel | all GDSC2 lines | 952 | 348 | 215 | 1 | g0 and early g1 | 0.369 | 0.88 (0.68 to 1.16) |
| ovary | paclitaxel | all GDSC2 lines | 952 | 348 | 215 | 2 | prolonged erk activation events | 0.433 | 1.11 (0.85 to 1.45) |
| ovary | paclitaxel | all GDSC2 lines | 952 | 348 | 215 | 3 | arms mediated activation | 0.978 | 1.00 (0.76 to 1.30) |
| ovary | paclitaxel | lineage-matched lines | 40 | 348 | 215 | 1 | activation of chaperone genes by xbp1s | 0.222 | 0.84 (0.65 to 1.11) |
| ovary | paclitaxel | lineage-matched lines | 40 | 348 | 215 | 2 | signaling by fgfr1 mutants | 0.186 | 1.20 (0.92 to 1.57) |
| ovary | paclitaxel | lineage-matched lines | 40 | 348 | 215 | 3 | dcc mediated attractive signaling | 0.00049 | 0.62 (0.47 to 0.81) |
| ovary | platinum (carboplatin or cisplatin) | all GDSC2 lines | 757 | 376 | 228 | 1 | base free sugar phosphate removal via the single nucleotide  | 0.754 | 0.96 (0.74 to 1.24) |
| ovary | platinum (carboplatin or cisplatin) | all GDSC2 lines | 757 | 376 | 228 | 2 | dcc mediated attractive signaling | 0.00041 | 1.60 (1.23 to 2.08) |
| ovary | platinum (carboplatin or cisplatin) | all GDSC2 lines | 757 | 376 | 228 | 3 | resolution of ap sites via the single nucleotide replacement | 0.572 | 1.08 (0.83 to 1.40) |
| ovary | platinum (carboplatin or cisplatin) | lineage-matched lines | 30 | 376 | 228 | 1 | metabolism of polyamines | 0.339 | 1.14 (0.87 to 1.47) |
| ovary | platinum (carboplatin or cisplatin) | lineage-matched lines | 30 | 376 | 228 | 2 | initial triggering of complement | 0.650 | 0.94 (0.73 to 1.22) |
| ovary | platinum (carboplatin or cisplatin) | lineage-matched lines | 30 | 376 | 228 | 3 | gluconeogenesis | 0.846 | 1.03 (0.79 to 1.33) |
| pancreas | gemcitabine | all GDSC2 lines | 954 | 104 | 51 | 1 | cdc6 association with the orc origin complex | 0.056 | 1.73 (0.98 to 3.06) |
| pancreas | gemcitabine | all GDSC2 lines | 954 | 104 | 51 | 2 | lagging strand synthesis | 0.076 | 1.66 (0.94 to 2.94) |
| pancreas | gemcitabine | all GDSC2 lines | 954 | 104 | 51 | 3 | e2f enabled inhibition of pre replication complex formation | 0.124 | 0.65 (0.37 to 1.13) |
| pancreas | gemcitabine | lineage-matched lines | 30 | 104 | 51 | 1 | purine salvage | 0.088 | 0.61 (0.35 to 1.08) |
| pancreas | gemcitabine | lineage-matched lines | 30 | 104 | 51 | 2 | purine catabolism | 0.810 | 1.07 (0.61 to 1.87) |
| pancreas | gemcitabine | lineage-matched lines | 30 | 104 | 51 | 3 | synthesis and interconversion of nucleotide di and triphosph | 0.033 | 1.84 (1.04 to 3.25) |
| sarcoma | gemcitabine | all GDSC2 lines | 954 | 42 | 20 | 1 | cdc6 association with the orc origin complex | 0.114 | 2.03 (0.83 to 5.00) |
| sarcoma | gemcitabine | all GDSC2 lines | 954 | 42 | 20 | 2 | lagging strand synthesis | 0.126 | 1.99 (0.81 to 4.90) |
| sarcoma | gemcitabine | all GDSC2 lines | 954 | 42 | 20 | 3 | e2f enabled inhibition of pre replication complex formation | 0.111 | 0.48 (0.19 to 1.21) |
| sarcoma | gemcitabine | lineage-matched lines | 53 | 42 | 20 | 1 | cell cycle mitotic | 0.114 | 2.03 (0.83 to 5.00) |
| sarcoma | gemcitabine | lineage-matched lines | 53 | 42 | 20 | 2 | g1 s specific transcription | 0.191 | 0.55 (0.22 to 1.36) |
| sarcoma | gemcitabine | lineage-matched lines | 53 | 42 | 20 | 3 | cell cycle | 0.114 | 2.03 (0.83 to 5.00) |

**What I think:** All 5 settings (x2 training sets) are now filled in, so the paper's comparison exists for every cancer. Across the rank-1 pathways, the smallest p is 0.045 (liver | sorafenib, all GDSC2 lines, n = 30, HR 3.49 (0.95 to 12.84)); the ovarian platinum and paclitaxel sets (n = 376 and 348) have the most patients and give HRs close to 1 (see table). Pancreas and liver have 104 and 30 patients, so their intervals are wide. The sign flips between the all-lines and matched-lines models for pancreas and liver, which suggests the pathway choice is driven by which lines are in the training set. The ovarian platinum model is trained on cisplatin only (no carboplatin screen in GDSC), as stated in the figure caption.

## Files and data versions

Data versions (file name, bytes, sha256) for each gap are in `tables/G*_data_versions.tsv`. Scripts are in `code/` (`python3 <script>` from inside `code/`; `PROGRESS.md` gives the order). Random seeds are set in `gcommon.py` (`SEED`).

## Fixes made while running

- `g7_bh3_beyond_markers.py`: skip a drug that is absent from a screen (WEHI-539 is not in GDSC1); the script had stopped there.
- `g5_folfox.py`: GSE39582 stage is stored as text, so the stage II-III filter returned no patients; it now compares text values.
