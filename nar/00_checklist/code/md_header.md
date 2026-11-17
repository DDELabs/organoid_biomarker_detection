# Paper checklist: Kong et al. 2020, Nature Communications 11:5485

"Network-based machine learning in colorectal and bladder organoid models predicts anti-cancer drug efficacy in patients".
DOI [10.1038/s41467-020-19313-8](https://doi.org/10.1038/s41467-020-19313-8), PMID 33127883, PMC7599252 (open access, CC BY).
Full text was retrieved via PubMed (PMC) and the NCBI PMC open-access copy.

This file lists **every** main figure panel, supplementary figure panel, supplementary table, source-data sheet and
text-only number in the paper. It was built by `code/build_checklist.py`. The same rows are in `checklist.tsv`
(machine-readable, tab-separated).

## 0. How to use this checklist (for the agents that come later)

- **Name every file you make with the checklist ID at the start**, e.g. `Fig2d-left_KM_COAD_5FU.png`,
  `SuppFig8b_targets_BLCA.pdf`, `SuppTable3_isogenic.tsv`. For other cancers add the cancer after the ID,
  e.g. `Fig2d-left_LUAD_cisplatin_KM.png`. The completeness checker searches `nar/` for the `search_token`
  column of `checklist.tsv` (the ID without the `-left`/`-COAD` suffix; ranges like `SuppFig1a-d` become `SuppFig1`). Match the token followed by a non-digit, so `SuppFig1` does not match `SuppFig10`.
- `status` column starts as `todo`. The checker should set `done` / `substitute` / `blocked` with a one-line reason.
- Numbers tagged **[derived]** are not printed in the paper; they were computed from the paper's own source data
  by `code/derive_source_data_facts.py` (output: `tables/source_data_facts.tsv`). Numbers tagged
  **[read from plot, approx]** were read by eye from a plotted axis. Every other number is printed in the paper.

## 1. Sources used (all cached under `nar/_cache/`)

| file | what | sha256 (first 16) |
|---|---|---|
| `41467_2020_Article_19313.pdf` | main article PDF (PMC open-access S3 copy, `PMC7599252.1.pdf`) | 1d4d34cad1b88adc |
| `PMC7599252_jats.xml` | JATS full text from NCBI efetch (figure legends, text with figure refs) | f031c8ab4074db96 |
| `main_text_PMC7599252.txt` | PubMed MCP full text (figure numbers stripped; used only as a cross-check) | - |
| `41467_2020_19313_MOESM1_ESM.pdf` | Supplementary Information: 18 Supp Figs + 4 Supp Tables (28 pages) | baf1c7f4d603eb2c |
| `41467_2020_19313_MOESM2_ESM.pdf` | Reporting summary (no results) | bfc626963bae8f30 |
| `41467_2020_19313_MOESM3_ESM.xlsx` | Source data, 13 sheets | c1e065aa9490694a |

Note: nature.com and media.springernature.com are blocked here. Main-figure numbers were read from the article PDF
pages (rendered to `nar/_cache/main_pages/`). Supplementary figures were extracted as images into `nar/_cache/supp_imgs/`
(img-000 = Supp Fig 1 ... img-017 = Supp Fig 18, img-018..021 = Supp Tables 1-4).

## 2. The paper's key claims in plain language (with the statistics it reports)

1. **Network filter.** Of 674 Reactome pathways, 37 are "close" to the 5-fluorouracil (5FU) target and 30 are close
   to cisplatin targets in the STRING protein network (13,824 proteins, 323,774 links; close = z-score <= -1.2816,
   from 1,000 random gene sets).
2. **Organoid biomarker.** Ridge regression on 19 colorectal organoids picks "activation of BH3-only proteins" as the
   5FU biomarker (rank 1). On bladder organoids (text: 9) it picks "amino acid synthesis and interconversion"
   for cisplatin (it is rank 2 in the bar chart, chosen because it was consistent across 3 algorithms).
   Coefficients agree across ridge, linear and SVR (R2 0.993-0.998).
3. **Patient survival.** In TCGA, predicted 5FU responders live longer than predicted non-responders
   (colon, n = 114, log-rank P = 0.014; responders had no deaths). Predicted cisplatin responders live longer
   (bladder, n = 77, P = 0.01; 5-year survival about 75% vs 50%).
4. **Not just prognosis.** In untreated patients the same split gives P = 0.16 (colon, n = 298) and P = 0.066
   (bladder, n = 294).
5. **Specific to cancer and drug.** Swapping the two biomarkers gives P = 0.16 (colon) and P = 0.86 (bladder).
   Predicted groups do not differ in age, sex, T stage or N stage (all chi-squared P >= 0.093).
6. **Other methods do worse.** Whole transcriptome (P = 0.69 / 0.79), whole pathway (0.81 / 0.82), direct target
   neighbours (0.96 / 0.37), 2nd neighbours (0.24 / 0.79), 3rd neighbours (0.36 / 0.79), Bolis gene filter (0.15 / 0.96),
   deep learning (0.91 / 0.82), network-centrality pathways (all 18 P between 0.091 and 0.997), drug-target
   genes alone (0.18-0.92), and GDSC cell-line-trained models (top pathway P = 0.17, 0.86, 0.102, 0.13)
   (colon / bladder).
7. **Multi-pathway.** Using the top 2 pathways keeps the bladder result (P = 0.022) but not colon (P = 0.88).
8. **Isogenic cells agree.** The BH3 pathway is higher in 5FU-sensitive than resistant HCT-8 cells
   (n = 3 vs 9, t-test P = 0.0085). The amino-acid pathway is higher in cisplatin-sensitive than resistant T24 cells
   (n = 3 vs 3, P = 0.00022).
9. **Random pathways rarely do as well.** 10,000 random pathway sets: empirical P = 0.0012 (colon) and 0.014 (bladder).
   Pathway size does not explain survival P (R2 = 0.0026 and 0.00034).
10. **Agreement with known mutation biomarkers.** A cetuximab model from colon organoids gives higher predicted
    resistance in BRAF V600E tumours (n = 47 vs 328 wild-type, one-sided Mann-Whitney P = 0.037; other feature sets
    P = 0.63-0.99). The cisplatin model gives lower predicted resistance in ERCC2-mutant bladder tumours
    (n = 37 vs 367, P = 0.002; other feature sets P = 0.017-0.04). With 1-3 organoids removed, >70% of subsets still
    recover the BRAF result.
11. **Within-organoid fit.** A 60/10/30 split of the organoids gives observed-vs-predicted R2 = 0.98 (colon) and
    0.89 (bladder).

## 3. Generalisation targets (used in the last column)

"GEN" in the table means the set of new cancer-drug pairs below. Patient counts are from the project brief
(TCGA drug tables) unless marked "to count".

| cancer | TCGA project, treated n | lab model with drug response (repo: `data/curated/preclinical/`) | isogenic pair idea (check GEO) | known mutation biomarker for a Fig5-style test |
|---|---|---|---|---|
| colorectal (external) | GSE39582: 5FU 82, FUFOL 54, FOLFOX 23, untreated 326 | van de Wetering 2015 organoids; Isella 2017 PDX | HCT-8/5FU (GSE81008, as paper) | BRAF V600E / KRAS vs cetuximab |
| bladder (other drugs) | BLCA cisplatin 77, methotrexate 12, gemcitabine to count | Lee 2018 organoids | cisplatin-resistant bladder lines | ERCC2 mut vs cisplatin |
| lung adenocarcinoma | LUAD cisplatin 87, carboplatin 77, pemetrexed 60, gemcitabine 25 | GDSC/CTRP lung lines; GSE165611 patient-derived cells (drug table not in repo) | A549 vs A549/DDP | EGFR mut vs EGFR TKI; KEAP1/STK11 vs platinum |
| ovarian | OV carboplatin 464, paclitaxel 448, cisplatin 156, olaparib 0 | Vias 2023 organoids; GDSC ovary lines | PEO1/PEO4, A2780/A2780cis | BRCA1/2 mut vs platinum |
| pancreas | PAAD gemcitabine (to count) | Tiriac 2018, Shi 2022 organoids | gemcitabine-resistant PANC-1 / MIA PaCa-2 | BRCA/ATM (HRD) vs platinum |
| liver | LIHC sorafenib (to count, likely small) | Broutier 2017, Ji 2023 organoids | sorafenib-resistant Huh7 / HepG2 | none established; use external sorafenib-response cohort if one exists in GEO |
| pan-cancer | all of the above | CTRPv2/CCLE, PRISM, Gao 2015 PDX | - | - |

## 4. Checklist table
