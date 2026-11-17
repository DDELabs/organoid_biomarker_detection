## 5. Things in the paper that do not match each other (worth knowing before reproducing)

1. **Supp Fig 5 labels are swapped in the main text.** Text: "bladder P = 0.022, Supp Fig 5a; colon P = 0.88, Supp Fig 5b".
   Legend and figure: a = colon (0.88), b = bladder (0.022). Source data agrees with the figure.
2. **Supp Fig 3 legend** says bladder panels are "(d)-(e)" but the figure has d, e and f.
3. **The bladder biomarker is rank 2, not rank 1**, in ridge (Fig 2f), linear and SVR (Supp Fig 2c,d). The bootstrap
   (Fig 4f) also places the observed result at rank 2 [derived: 140 of 10,000 random sets reach rank <= 2 -> P = 0.014].
4. **The colon biomarker only just passes the proximity cut-off**: z = -1.2857 vs cut-off -1.2816 [derived].
5. **Proximal-pathway counts depend on a hidden filter.** The source z-score matrix gives 38 (5FU) and 47 (cisplatin)
   pathways with z <= -1.2816; only 37 and 30 of these have organoid ssGSEA scores [derived]. So "37" and "30" are
   after removing pathways missing from the organoid data (the bladder organoid sheet has only 488 of 674 pathways).
6. **Bladder organoid count.** Text says 9 bladder organoids; source sheet 'Figure 2e' has 11 sample columns
   (some patients have two lines, e.g. SCBO-3 and SCBO-3.2); Supp Fig 17b shows 9 with cisplatin IC50.
   Colon sheet 'Figure 2b' has 22 columns for 19 organoids.
7. **Source-data sheet names do not match panels.** 'Figure 2g' and 'Figure 2h' hold TCGA COAD/BLCA ssGSEA scores
   (inputs to Fig 2d/2h KM), not the network drawings. Sheet 'FigureS5' is titled "Figure 17" inside.
8. **Fig 5a** model was trained on 16 organoids (3 hypermutated removed: P19a, P19b, P24a), while Supp Fig 11 shows 19.
9. **Supp Fig 13**: all five alternative feature sets also give P < 0.05 for ERCC2 (0.017-0.04); the paper
   describes them as "reduced statistical power" compared with P = 0.002.
10. **Supp Fig 15c,d**: cell-line-trained models give P = 0.0092 for their rank-2 (all lines) or rank-3 (bladder lines)
    pathway; the paper reports these models as not strongly predictive because the top-ranked pathway was not significant.
11. **TCGA COAD "5FU-treated" n = 114** in the paper, but the project brief says only ~5 TCGA-COAD patients got 5FU alone
    and most got FOLFOX. The paper's 114 almost certainly includes 5FU-containing combinations.
12. **Isogenic BLCA data (T24S/T24R, Yeon et al.)** were obtained by personal request; only ssGSEA scores are public
    (MOESM3 'Figure 4c').
13. Text says pathway mean proximity count 4.9; the source matrix gives 4.86 [derived] (rounding, consistent).
14. STRING sheet lists each link twice (647,548 rows = 2 x 323,774) and has 13,825 distinct UniProt IDs vs 13,824 in the text [derived].

## 6. What the original code in `python/` already produces

| script | produces | checklist items it covers |
|---|---|---|
| `python/run_ssGSEA.py` (+ `utilities/`) | Reactome ssGSEA NES for organoids and TCGA (`python/results/COAD/...`) | inputs for Fig2b, Fig2d, Fig2f, Fig2h, Fig3c, SuppFig1 |
| `python/single_pathway_prediction.py` | per-ML (Ridge/SVR/Linear) coefficients + KM plot of the rank-k pathway; `single_pathway_predictions.txt` (COAD 5FU only configured) | Fig2b, Fig2d-left, Fig2f/Fig2h-left (if BLCA enabled), SuppFig2, SuppFig3 (coefficients) |
| `python/multiple_pathway_prediction.py` | top-k (2-10) multi-pathway KM plots + `multi_pathway_predictions.txt` (present in `python/results/multi_pathway_predictions/`) | SuppFig5, SourceData-multipathway |
| `modified_codes_py3/proximity_code.py` | network proximity z-scores (user's Python 3 port) | Fig2a, Fig2e, SuppFig16 |

Not in the original code (must be written new): untreated-control KM (Fig2d/h right), all Fig3 baselines,
isogenic t-tests (Fig4b,c, SuppTable3), bootstrap (Fig4e,f, SuppTable4), mutation concordance (Fig5, SuppFig12, 13, 18),
cross-biomarker test (SuppFig4), centrality (SuppFig6, 7), target-gene KM (SuppFig8), k-hop neighbours (SuppFig9),
pathway-size check (SuppFig10), cetuximab heatmap (SuppFig11), GDSC models + PCA (SuppFig14, 15), split-organoid test
(SuppFig17), clinical tables (SuppTable1, 2), network drawings (Fig2c, g), PCA (SuppFig1).

## 7. Counts

Main-figure items: 31. Supplementary-figure items: 37. Supplementary tables: 4. Source-data-only items: 3. Text-only items: 2.
Total: 77 checklist rows (see `checklist.tsv`).
