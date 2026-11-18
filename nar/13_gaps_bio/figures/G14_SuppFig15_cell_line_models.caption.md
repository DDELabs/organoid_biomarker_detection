# Supp. Fig. 15 for the other cancers: models trained on cancer cell lines instead of organoids

The paper compared its organoid-trained model with the same pipeline trained on cancer cell lines (GDSC). We repeat
this for pancreas (gemcitabine), liver (sorafenib), ovary (platinum and paclitaxel) and sarcoma (gemcitabine).
Training data: GDSC2 cell lines (response = log IC50, lower = more sensitive), either all lines with RNA-seq or only
lines from the matching cancer type. Ovarian platinum: the model is trained on cisplatin (no carboplatin screen in
GDSC) and applied to patients who received carboplatin or cisplatin. Each of the top three pathways splits treated
TCGA patients at the median prediction.

- **a** Hazard ratio (predicted non-responders vs predicted responders, 95% CI) and log-rank p for every setting and
  each of the top 3 pathways. HR > 1 = the expected direction.
- **b onward** Kaplan-Meier curves for the rank-1 pathway of each setting.

- pancreas | gemcitabine, all GDSC2 lines (954 lines): rank-1 cdc6_association_with_the_orc_origin_complex; treated n = 104, deaths = 51, p = 0.056, HR 1.73 (0.98-3.06); inside the lines leave-one-out rho = nan
- pancreas | gemcitabine, lineage-matched lines (30 lines): rank-1 purine_salvage; treated n = 104, deaths = 51, p = 0.088, HR 0.61 (0.35-1.08); inside the lines leave-one-out rho = 0.37
- liver | sorafenib, all GDSC2 lines (952 lines): rank-1 fgfr1_ligand_binding_and_activation; treated n = 30, deaths = 15, p = 0.045, HR 3.49 (0.95-12.84); inside the lines leave-one-out rho = nan
- liver | sorafenib, lineage-matched lines (15 lines): rank-1 negative_regulation_of_the_pi3k_akt_network; treated n = 30, deaths = 15, p = 0.084, HR 0.35 (0.10-1.21); inside the lines leave-one-out rho = 0.23
- ovary | platinum (carboplatin or cisplatin), all GDSC2 lines (757 lines): rank-1 base_free_sugar_phosphate_removal_via_the_single_nucleotide_replacement_pathway; treated n = 376, deaths = 228, p = 0.754, HR 0.96 (0.74-1.24); inside the lines leave-one-out rho = nan
- ovary | platinum (carboplatin or cisplatin), lineage-matched lines (30 lines): rank-1 metabolism_of_polyamines; treated n = 376, deaths = 228, p = 0.339, HR 1.14 (0.87-1.47); inside the lines leave-one-out rho = 0.72
- ovary | paclitaxel, all GDSC2 lines (952 lines): rank-1 g0_and_early_g1; treated n = 348, deaths = 215, p = 0.369, HR 0.88 (0.68-1.16); inside the lines leave-one-out rho = nan
- ovary | paclitaxel, lineage-matched lines (40 lines): rank-1 activation_of_chaperone_genes_by_xbp1s; treated n = 348, deaths = 215, p = 0.222, HR 0.84 (0.65-1.11); inside the lines leave-one-out rho = 0.12
- sarcoma | gemcitabine, all GDSC2 lines (954 lines): rank-1 cdc6_association_with_the_orc_origin_complex; treated n = 42, deaths = 20, p = 0.114, HR 2.03 (0.83-5.00); inside the lines leave-one-out rho = nan
- sarcoma | gemcitabine, lineage-matched lines (53 lines): rank-1 cell_cycle_mitotic; treated n = 42, deaths = 20, p = 0.114, HR 2.03 (0.83-5.00); inside the lines leave-one-out rho = 0.21
