# Fig. 4b/c equivalents for other cancers: biomarker pathways in drug-resistant cell lines

Paper Fig. 4b/c compared the biomarker pathway between a drug-sensitive parental cell line and a drug-resistant copy of
the same line. We found two public pairs that match drugs of our other sets:
- GSE15709: ovarian cancer line A2780, parental (cisplatin-sensitive, 5 replicates) vs Round5 cisplatin-resistant
  (5 replicates), Affymetrix arrays. Tested: the platinum biomarkers of bladder, ovarian and lung sets
  (carboplatin uses the cisplatin network, as in the main analyses).
- GSE322742: liver cancer lines HepG2 and Huh7, parental vs sorafenib-resistant (2 vs 2 each, RNA-seq). Tested: the
  sorafenib biomarkers of the Ji and Broutier liver organoids. Note: the resistant cells were kept in sorafenib, so drug
  exposure and resistance are mixed.
Pathways were scored with ssGSEA (as in the paper). Bars: -log10 Welch t-test p; green = the difference goes the way
the organoid Ridge weight predicts (positive weight -> higher in resistant cells), red = opposite. Dashed line = p 0.05.
With 2 vs 2 samples a t-test has little power. No public isogenic pair with expression data was found for gemcitabine
in pancreatic organoid drugs within this work, nor for sarcoma, melanoma (dacarbazine) or the PDX drugs.
Numbers: tables/Fig4bc_isogenic_other_cancers.tsv.
