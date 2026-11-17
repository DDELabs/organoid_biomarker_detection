# 08 DepMap / cell-line mechanism: report

Agent folder: `nar/08_depmap/`. Code: `code/` (run in order with `python -I` from the repo root:
`01_prepare.py`, `02_partA.py`, `03_supp.py`, `04_proximity.py`, `05_stage0_assoc.py`, `06_freeze.py`).
Figures: `figures/` (PNG 300 dpi + PDF + `.caption.md`). Tables: `tables/`. Every number below comes from a table
written by this code. Seed 20261117; random-signature nulls use fixed seeds in the code.

## What was asked, in one paragraph

The paper's colon 5-FU model was trained on 19 organoids. Its most stable pathway is **"Activation of BH3-only
proteins"**. (BH3-only proteins such as BIM, PUMA and NOXA are the cell's pro-death switches.) In organoids, a higher
score went with a lower IC50, which means more 5-FU sensitivity. Part A checks this in cancer cell lines from four
drug screens (GDSC1, GDSC2, CTRPv2, PRISM), using RNA-seq, CRISPR and mutation data. Part B builds the ANCHOR
"Stage 0" list for 10 drugs: at most 5 signed pathways per drug, from cell-line and network data only. The list is
frozen before anyone looks at patient outcomes.

## Data used (versions)

All downloads are in `nar/_cache/depmap/`. Each file is in its own new folder and was read with `python -I`.
File names, sizes and sha256 hashes are in `tables/data_versions.tsv`. Main files:
- GDSC release 8.5: `GDSC1_fitted_dose_response_27Oct23.xlsx` (29,353,210 bytes),
  `GDSC2_fitted_dose_response_27Oct23.xlsx` (21,330,376 bytes) and `screened_compounds_rel_8.5.csv`.
- Cell Model Passports:
  - `model_list_latest.csv.gz` (cancer type, and MSI status);
  - `rnaseq_all_20220624.zip` (940,707,972 bytes; the TPM matrix was used as log2(TPM+1));
  - `mutations_summary_latest.csv.gz` (driver mutations).
- Project Score, combined Sanger v1 + Broad 20Q2 (20210311): the batch-corrected fold-change matrix. Only QC-pass
  screens were kept: 961 screens, which became 811 models after duplicates were averaged, with 17,303 genes.
- CTRPv2 and PRISM (with CCLE 21Q4 RNA-seq) were reused from `data/curated/preclinical/`, not downloaded again.
- Bladder organoids (Lee 2018) from `data/curated/preclinical/bladder_lee2018`.
- Other files:
  - STRING v12 human links, score >= 700 (15,187 proteins, 217,561 edges);
  - the current Reactome GMT (for S8);
  - GEO GSE81005, GSE196900 and GSE56322, with the platform tables GPL15207 and GPL14550 (for S7).
- Pathway scores: Reactome gene sets from the repo's MSigDB v6.1 file, scored with the single-sample rank score
  (`robust/obd/scoring.rank_score`). The earlier external validation used the same scoring.
- Mutations: the full CMP mutation file (`mutations_all_latest`, 309 MB) was reachable, but the 0.6 MB driver
  summary already covers TP53, KRAS, BRAF and APC, so it was used instead.

**S1, matching flow** (`tables/S1_data_matching_flow.tsv`, `figures/S1_data_matching_flow.*`):

| Screen | Colorectal lines with 5-FU and RNA-seq | ...also with CRISPR |
|---|---|---|
| GDSC1 | 44 | 35 |
| GDSC2 | 48 | 36 |
| CTRPv2 | 43 | 31 |
| PRISM | 22 | 15 |

All lines with 5-FU and RNA-seq: 903, 962, 775 and 429.

---

## Part A. Does the organoid 5-FU biomarker hold in cell lines?

Direction convention:
- Response is ln IC50 (GDSC) or AUC (CTRPv2, PRISM). **Higher = more resistant.**
- The organoid model predicts a **negative** correlation between the BH3-only score and 5-FU response.

### (a) Workflow (Figure 1a)
The steps are:
1. Start from the organoids.
2. Keep the 37 network-proximal pathways.
3. Train the organoid model and freeze it.
4. Score it, unchanged, in cell lines.
5. Run tests (b) to (h).

### (b) BH3-only score vs 5-FU in colorectal lines (Figure 1b, `tables/A_b_bh3_vs_5fu_colon.tsv`)

| Screen | n | Spearman rho (95% CI) | p |
|---|---|---|---|
| GDSC1 (ln IC50) | 44 | 0.054 (−0.247 to 0.345) | 0.73 |
| GDSC2 (ln IC50) | 48 | 0.016 (−0.270 to 0.298) | 0.92 |
| CTRPv2 (AUC) | 43 | 0.257 (−0.047 to 0.518) | 0.096 |
| PRISM (AUC) | 22 | 0.405 (−0.020 to 0.706) | 0.062 |

What I think: in colorectal cell lines, the BH3-only score does not follow 5-FU sensitivity the way it did in
organoids. GDSC shows no relationship. CTRPv2 and PRISM lean the other way (lines with a high BH3-only score are
slightly more resistant). With only 22–43 lines the confidence intervals are wide.

### (c) The organoid model applied unchanged (Figure 2c, `tables/A_c_frozen_model_colon_lines.tsv`)
The main signature is `robust_v1`. Its weights are: BH3-only −0.347, hyaluronan uptake +0.334, CO2 hydration
+0.263. We expect a positive rho between predicted resistance and observed response. The random band is the
2.5–97.5% range of rho for 1,000 random signatures, each made of 3 pathways with random signs.

| Screen | n | rho | random 95% band | random-null p (one-sided) |
|---|---|---|---|---|
| GDSC1 | 44 | −0.111 | −0.414 to 0.386 | 0.67 |
| GDSC2 | 48 | −0.076 | −0.431 to 0.377 | 0.60 |
| CTRPv2 | 43 | −0.202 | −0.396 to 0.382 | 0.81 |
| PRISM | 22 | −0.298 | −0.424 to 0.421 | 0.89 |

We tested all five frozen signatures (paper_top7, robust_v1, v1_rank, nit_soft, nit_deconf) in all four screens.
Every result fell inside the random band: rho ranged from −0.360 to 0.105, and the random-null p from 0.34 to 0.94.

What I think: the frozen organoid model ranks colorectal cell lines by 5-FU response no better than random pathway
sets. This matches the earlier patient result (random-null p of about 0.2 in GSE39582).

### (d) Sign agreement of the 37 network-proximal pathways (Figure 2d, `tables/A_d_sign_agreement_37.tsv`)

| Screen | same sign | fraction | binomial p (one-sided, > 50%) |
|---|---|---|---|
| GDSC1 | 24/37 | 0.649 | 0.049 |
| GDSC2 | 14/37 | 0.378 | 0.95 |
| CTRPv2 | 14/37 | 0.378 | 0.95 |
| PRISM | 12/37 | 0.324 | 0.99 |
| all four (not independent) | 64/148 | 0.432 | 0.96 |

What I think: only GDSC1 agrees more often than chance would give (p = 0.049). The other three screens lean the
other way. As a set, the organoid pathway weights do not carry over to colorectal lines.

### (e) Drug specificity: BH3-only score vs every drug (Figure 3e, `tables/A_e_*.tsv`)
We used a lineage-adjusted partial Spearman correlation across all lines, for every drug tested on at least 30 lines.
"Lineage-adjusted" means that differences between cancer types are removed first. In the rank column, rank 1 is
the drug where lines with a high BH3-only score are most sensitive.

| Screen | drugs | 5-FU pan-cancer rho | 5-FU rank | 5-FU colorectal rho (rank) |
|---|---|---|---|---|
| GDSC1 | 378 | 0.018 | 203/378 | 0.054 (322/378) |
| GDSC2 | 286 | 0.021 | 204/286 | 0.016 (184/273) |
| CTRPv2 | 460 | 0.011 | 409/460 | 0.257 (386/452) |
| PRISM | 1,442 | 0.073 | 1179/1442 | 0.405 (1348/1435) |

The drugs most strongly linked to the BH3-only score are **BH3 mimetics** (drugs that block the BCL-2 family of
survival proteins):
- GDSC2: ABT-737 rho −0.246 and navitoclax −0.202 (both q < 0.001), plus the BCL-XL inhibitor WEHI-539 −0.159
  (q = 0.0001).
- GDSC1: navitoclax −0.206 (rank 2).
- CTRPv2: ABT-737 −0.238, venetoclax −0.226 and navitoclax −0.211 (all q < 0.001).
- PRISM: no drug reaches q < 0.05.

What I think: the BH3-only score is a real biological readout. It measures how close cells are to the death
threshold in general ("apoptotic priming"), not something specific to 5-FU. It picks out exactly the drugs that act
on the BCL-2 family, while 5-FU sits in the middle of the ranking.

### (f) Bladder: organoid-trained cisplatin model (Figure 3f, `tables/A_f_*.tsv`)
- Training data: 9 bladder organoids from Lee 2018 (cisplatin AUC), using the 47 cisplatin-proximal pathways and the
  robust stability learner.
- Resulting signature: glucose transporters (+0.293), gene expression in beta cells (−0.240), organic anion
  transport (−0.355).

| Screen | drug | bladder lines | rho | random 95% band | random-null p |
|---|---|---|---|---|---|
| GDSC1 | cisplatin | 18 | 0.191 | −0.684 to 0.664 | 0.35 |
| GDSC2 | cisplatin | 17 | 0.162 | −0.596 to 0.527 | 0.33 |
| CTRPv2 | carboplatin (CTRPv2 has no cisplatin) | 24 | −0.033 | −0.365 to 0.387 | 0.59 |
| PRISM | cisplatin | 22 | 0.316 | −0.457 to 0.461 | 0.11 |

What I think: three of the four screens point the expected way (positive rho), but every value lies inside the
random band. With 9 training organoids and 17–24 bladder lines, this test cannot separate signal from noise.

### (g) BH3-only score vs BH3-mimetic sensitivity (Figure 4g, `tables/A_g_bh3_vs_bh3_mimetics.tsv`)
Pan-cancer, lineage-adjusted results (negative rho = lines with a high BH3-only score are more sensitive):

| Drug | GDSC1 | GDSC2 | CTRPv2 | PRISM |
|---|---|---|---|---|
| Navitoclax (ABT-263) | −0.206 (n = 926) | −0.202 (n = 961) | −0.211 (n = 772) | −0.112 (n = 471) |
| ABT-737 | not screened | −0.246 | −0.238 | −0.064 (p = 0.19) |
| Venetoclax | −0.087 | −0.104 | −0.226 | +0.028 (p = 0.55) |
| Obatoclax | +0.046 | +0.036 | −0.092 | +0.004 |

All four navitoclax correlations have p < 0.02. Within colorectal lines alone (n = 25–48), every confidence interval
crosses 0 (rho from −0.164 to +0.405).

What I think: this is the clearest positive result in Part A. Across about 500–960 lines, the score behaves like a
priming marker for the BCL-2 and BCL-XL inhibitors navitoclax and ABT-737. Obatoclax does not follow, which fits its
known weak and off-target action. There are too few colorectal lines to see the effect within colon alone.

### (h) CRISPR dependency vs 5-FU (Figure 4h, `tables/A_h_crispr_vs_5fu.tsv`)
We compared the knockout effect of MCL1, BCL2L1, BCL2 and the BH3-only genes with 5-FU response.
- Pan-cancer, lineage-adjusted: rho ranged from −0.092 to +0.110. Only two had p < 0.05, both in CTRPv2 only:
  BBC3 (PUMA) +0.110 (p = 0.030) and BMF +0.105 (p = 0.038).
- Colorectal lines (n = 15–36): rho ranged from −0.311 to +0.309, and none had p < 0.05. The smallest p was HRK in
  GDSC2 (rho 0.307, p = 0.068).

What I think: how much a line depends on BCL-2 family genes does not predict its 5-FU response. In 2-D cell lines,
5-FU killing does not seem to be controlled by the BCL-2 family.

---

## Supplementary analyses (S1–S8)

**S1 data matching:** see the Data section above.

**S2 GDSC1 vs GDSC2** (`tables/S2_gdsc1_vs_gdsc2_agreement.tsv`):

| Comparison of 5-FU response | all lines: rho (n) | colorectal: rho (n) |
|---|---|---|
| GDSC1 vs GDSC2 (ln IC50) | 0.523 (902) | 0.677 (44) |
| GDSC2 vs CTRPv2 | 0.570 (548) | 0.447 (36) |

What I think: the screens share only about a quarter to a half of their ranking variance for 5-FU. That limits how
well any marker can replicate from one screen to another.

**S3 colorectal vs pan-cancer** (`tables/S3_colon_vs_pancancer.tsv`):
- Without lineage adjustment, the pan-cancer correlation runs in the organoid direction: GDSC1 −0.241, GDSC2 −0.233,
  CTRPv2 −0.234 (all p < 1e-12).
- After cancer type is removed, it is close to 0: GDSC1 0.018, GDSC2 0.021, CTRPv2 0.011, PRISM 0.073.
- Within single cancer types, rho scatters around 0. Examples from CTRPv2: blood −0.389, upper aerodigestive
  −0.381, oesophagus +0.470.

What I think: the raw pan-cancer signal comes from differences between cancer types, not from differences within
one tumour type. For example, blood-cancer lines are both BH3-high and 5-FU-sensitive.

**S4 permutation and random controls** (`tables/S4_permutation_random_controls.tsv`):

| Control | GDSC1 | GDSC2 | CTRPv2 | PRISM |
|---|---|---|---|---|
| Label-permutation p (two-sided) | 0.73 | 0.91 | 0.096 | 0.059 |
| BH3-only percentile among 672 Reactome pathways | 0.60 | 0.43 | 0.85 | 0.96 |
| Random 16-gene sets, one-sided p | 0.63 | 0.53 | 0.91 | 0.96 |

A low percentile means the pathway is among the most linked to 5-FU sensitivity.

What I think: in colorectal lines the BH3-only score is an ordinary pathway for 5-FU, not an outlier.

**S5 adjustment for TP53, MSI, KRAS and BRAF** (`tables/S5_mutation_adjustment.tsv`):

| BH3-only effect per SD (negative = organoid direction) | GDSC1 | GDSC2 | CTRPv2 | PRISM |
|---|---|---|---|---|
| Unadjusted | 0.062 | 0.042 | 0.270 | 0.319 |
| Adjusted | 0.238 (p 0.12) | 0.231 (p 0.10) | 0.374 (95% CI 0.020–0.727, p 0.039) | 0.434 (p 0.33) |

In GDSC, MSI lines were more sensitive to 5-FU: β −1.09 (p 0.002) in GDSC1 and −0.95 (p 0.004) in GDSC2.

What I think: adjusting for the driver mutations moves the BH3 estimate further from the organoid direction. MSI
is the strongest 5-FU covariate in colorectal lines, which fits the biology.

**S6 each BH3-only gene** (`tables/S6_single_genes.tsv`):
- No core BH3-only gene has a consistent negative rho.
- Some go the other way: BMF in GDSC1 (0.398, p 0.007) and GDSC2 (0.285, p 0.050), and BAD in CTRPv2 (0.335, p 0.028).
- Among the other pathway genes, E2F1 (GDSC2 0.377, p 0.008; CTRPv2 0.438, p 0.003) and YWHAB (GDSC2 0.319;
  CTRPv2 0.389) are positive.

What I think: in cell lines the score leans towards "BH3-high = resistant". That lean comes from the proliferation
and E2F-type members of the gene set, not from the death genes themselves.

**S7 isogenic 5-FU-resistant vs parental lines** (`tables/S7_isogenic_pairs.tsv`, `figures/S7_isogenic_pairs.*`):
The table shows the resistant-minus-parental BH3-only score (MSigDB set). The organoid model expects a negative value.

| Dataset | Comparison | Difference | p |
|---|---|---|---|
| GSE81005 HCT-8 | 3 vs 3, before 5-FU (0 h) | −0.008 | 0.18 |
| GSE81005 HCT-8 | 3 vs 3, after 24 h of 5-FU | −0.020 | 0.007 |
| GSE81005 HCT-8 | 3 vs 3, after 48 h of 5-FU | −0.016 | 0.019 |
| GSE196900 HCT116 | 2 vs 2 | +0.009 | 0.12 |
| GSE196900 SW480 | 2 vs 2 | −0.011 | 0.16 |
| GSE56322 | 1 vs 1 | −0.003 | no test |

All four gene-set versions give a negative difference in HCT-8 at 24 h and at 48 h (p 0.0002–0.026).

What I think: this is the one cell-line setting that supports the organoid direction. When exposed to 5-FU, the
resistant HCT-8 line switches on the BH3-only programme less than its parent. It is one isogenic pair, so it is a
mechanistic hint, not a validation. Before treatment (0 h) the difference is small.

**S8 IC50 vs AUC, and gene-set version** (`tables/S8_metric_and_geneset_robustness.tsv`):
- We tried both response metrics and all four versions of the BH3-only gene set.
- Colorectal rho with 5-FU ranged from −0.086 to 0.135 in GDSC, and from 0.072 to 0.405 in CTRPv2 and PRISM.
- Pan-cancer adjusted rho ranged from −0.025 to 0.113.

What I think: the answer does not depend on the response metric or on which version of the gene set is used.

### Part A overall
What I think: the cell-line data show that the BH3-only score measures how primed cells are to die through the
BCL-2 family. Its link to BH3 mimetics is strong and is seen in three screens. But this priming does not set 5-FU
sensitivity in 2-D colorectal lines, so the organoid association is not reproduced there. One isogenic pair (HCT-8)
does show the expected change after 5-FU exposure. Organoids and cell lines differ in several ways: 3-D growth,
longer assays, and the organoid IC50 range. Any of these could explain the gap, and the cell-line data cannot tell
them apart.

---

## Part B. ANCHOR Stage 0: frozen, outcome-free hypotheses

**No patient data were opened by this agent.** Every hypothesis comes from cell lines and the protein network only.

### How the list was built (rules fixed in `code/06_freeze.py` before the result was produced)
Three independent label-free sources were used.
1. **Expression association** (`code/05_stage0_assoc.py`). For each drug and screen, every Reactome pathway score is
   correlated with drug response after removing cancer type and a "general drug sensitivity" factor. That factor is
   the line's average z-scored response to all other drugs, which takes out "this line dies from everything".
   - GDSC1 and GDSC2 share lines, so their z values are averaged and counted as one source.
   - The Sanger, CTRPv2 and PRISM results are then combined with Stouffer's method.
   - A pathway passes if the combined BH q < 0.05 and at least 2 screens have the same sign.
2. **CRISPR co-dependency.** The same adjusted drug response is correlated with the Project Score knockout effect
   of every gene. The gene results are then summarised per pathway against 2,000 random gene sets of the same size.
   A pathway passes if q < 0.05.
3. **Network proximity.** STRING v12 (score >= 700) closest-distance z-score between the pathway and the drug's
   DrugBank targets, using `robust/obd/network.proximity_z` with 1,000 degree-matched randomisations. A pathway passes
   if z <= −1.2816.
   - Carboplatin and oxaliplatin have no protein target in DrugBank, so they use the cisplatin target set, as in
     `robust/run_survlearn.py`.
   - The paper's own proximity table is also listed for comparison (`tables/B_network_proximity_z.tsv`).

A pathway is **frozen** for a drug only if all of these hold:
- at least 2 of the 3 sources pass;
- the sign is clear (expression p < 0.05). The sign comes from the expression association:
  "+" means a higher pathway score goes with resistance, "−" with sensitivity;
- it is **specific**: the drug's directional z is in the top 25% of all drugs in the same screens, and is larger
  than the mean for related drugs of the same mechanism class (for example, the other platinum drugs for a platinum
  drug).

Among eligible pathways, the ones with more sources and a larger |z| come first. Near-duplicate pathways (more than
50% shared genes) are skipped, and at most 5 are kept per drug.

**Drug coverage.**
- Irinotecan uses SN-38 (its active metabolite) in GDSC1 and CTRPv2, where irinotecan itself was not screened.
- Pemetrexed is only in GDSC1. Doxorubicin is only in GDSC1 and PRISM (CTRPv2 has no doxorubicin under that name).
- Carboplatin is only in CTRPv2 and PRISM (GDSC has no carboplatin).

### Result (`tables/frozen_hypotheses.tsv`, `tables/B_stage0_summary.tsv`, Figure 5)
Frozen at **2026-11-17T13:03:40Z**. The file's sha256 is
`a0bdd620293e2559667cf7f4c218894ae3e70cf1eb404558a084d28aee878ab7` (recorded in `tables/frozen_hypotheses.sha256`).
The script refuses to overwrite the frozen file: a rerun writes `frozen_hypotheses.rerun.tsv` instead.

| Drug | pathways passing expression / CRISPR / network | 2+ sources | 2+ sources, signed and specific | frozen |
|---|---|---|---|---|
| 5-FU | 99 / 81 / 30 | 16 | 9 | 5 |
| Oxaliplatin | 312 / 95 / 70 | 78 | 72 | 5 |
| Irinotecan | 296 / 96 / 52 | 91 | 14 | 5 |
| Cisplatin | 0 / 82 / 70 | 4 | 0 | 0 |
| Carboplatin | 9 / 29 / 70 | 2 | 0 | 0 |
| Gemcitabine | 257 / 69 / 36 | 79 | 54 | 5 |
| Paclitaxel | 275 / 71 / 85 | 71 | 46 | 5 |
| Doxorubicin | 54 / 10 / 8 | 0 | 0 | 0 |
| Pemetrexed | 80 / 53 / 33 | 15 | 9 | 5 |
| Sorafenib | 122 / 40 / 154 | 51 | 31 | 5 |

In total, **35 hypotheses for 7 of the 10 drugs**. The ones supported by all 3 sources are:
- **oxaliplatin**: TGF-β receptor signalling (+) and down-regulation of SMAD2/3:SMAD4 transcription (+);
- **gemcitabine**: metabolism of nucleotides (−) and APOBEC3G-mediated resistance to HIV-1 (−);
- **sorafenib**: SEMA4D semaphorin signalling (+) and RHO GTPase signalling (+).

The 5-FU list is:
1. generic transcription pathway (+);
2. sulfur amino-acid metabolism (−);
3. respiratory electron transport (−);
4. metabolism of vitamins and cofactors (−);
5. CDC6 association with the ORC origin complex (+).

The BH3-only pathway is **not** among them, which agrees with Part A.

What I think:
- Seven drugs have at least one two-source, drug-specific, signed hypothesis. The plan's gate G1 asks for 8 drugs
  to "GO" and treats fewer than 4 as a global stop, so 7 sits between the two thresholds.
- Cisplatin, carboplatin and doxorubicin get no frozen hypothesis.
  - Cisplatin: no pathway passed the expression source after general sensitivity was removed (0 of 672). This fits
    the idea that, in cell lines, cisplatin response is mostly a "general sensitivity" signal.
  - Carboplatin: 2 pathways had 2+ sources, but none was specific against the other platinum drugs.
  - Doxorubicin: no pathway had 2 or more sources.
- By the plan's own rule, these three drugs have no mechanistic anchor. That is a result in itself, not a failure.
- Several frozen pathways are biologically sensible:
  - nucleotide metabolism for gemcitabine;
  - mitotic and cell-cycle pathways for paclitaxel;
  - TGF-β/SMAD for oxaliplatin;
  - integrin and extracellular-matrix (ECM) signalling for pemetrexed resistance.
- Others are hard to interpret, such as "APOBEC3G-mediated resistance to HIV-1" and "HIV infection". These Reactome
  sets are mostly made of proteasome, nuclear-pore or immune genes.
- Stage 1 should test the list exactly as frozen, and check those odd sets for composition effects (purity, immune
  cells).

Also written:
- `tables/B_candidates_two_or_more_sources.tsv`: all 407 drug–pathway pairs with 2 or more sources, with every flag,
  for transparency;
- `tables/B_source1_expression_association.tsv.gz`, `tables/B_source2_crispr_codependency.tsv.gz`: the full source
  tables;
- `tables/B_drug_targets_used.tsv`: the drug targets used.

---

## Figures and tables

| Figure | Content |
|---|---|
| `fig1_workflow_bh3_vs_5fu` | (a) workflow, (b) BH3-only score vs 5-FU in colorectal lines, 4 screens |
| `fig2_frozen_model_and_sign_agreement` | (c) frozen organoid model vs observed, with random null band; (d) sign agreement of the 37 pathways |
| `fig3_drug_specificity_and_bladder` | (e) BH3-only score vs every drug; (f) bladder organoid cisplatin model in bladder lines |
| `fig4_bh3_mimetics_and_crispr` | (g) BH3 mimetics; (h) CRISPR dependency vs 5-FU |
| `fig5_anchor_stage0_frozen_hypotheses` | Part B heatmap of the frozen hypotheses with evidence counts, plus the evidence funnel |
| `S1`–`S8` | matching flow; GDSC1 vs GDSC2; colon vs pan-cancer; permutation and random controls; mutation adjustment; single genes; isogenic pairs; metric and gene-set robustness |

Each figure has a PNG, a PDF and a `.caption.md` file in `figures/`. Each table is in `tables/`, named `A_*`, `S*_*`
or `B_*`.

## Limits and substitutions (one line each)
- In CTRPv2 there is no cisplatin, so carboplatin was used for the bladder test (f).
- In GDSC there is no carboplatin. Pemetrexed is only in GDSC1. ABT-737 is only in GDSC2 among the GDSC screens.
- PRISM has only 22 colorectal lines with expression, and 15 with CRISPR.
- The bladder organoid model rests on 9 organoids, because only these have both cisplatin AUC and RNA-seq.
- For S7, GEO had three 5-FU resistant/parental colorectal sets with usable expression: GSE81005, GSE196900 and
  GSE56322.
