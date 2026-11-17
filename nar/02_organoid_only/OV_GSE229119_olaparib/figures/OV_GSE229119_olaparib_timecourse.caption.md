# Olaparib time course in ovarian organoids (GSE229119) and what it says about TCGA-OV

What the data are: 3 high-grade serous ovarian cancer organoids (2414, 2445, 3008) and 4 cell lines (Hey, MCF7 [breast],
OVCA432, SKOV3), each treated with olaparib (a PARP inhibitor) and sampled at day 3, day 7 and after the drug was removed.
There is one DMSO control per model. Every model got the drug, and no model has a 'how strongly did it respond' number,
so the paper's regression step (which needs a response value per model) cannot be run on these data.

**a.** Change in Reactome pathway score (ssGSEA) after olaparib, minus the same model's DMSO control. Rows are the
15 pathways that went up and 15 that went down at day 7 in all 3 organoids (the 'olaparib-response signature').
**b.** The pathways closest to olaparib's targets PARP1/2/3 in the paper's network (115 pass z <= -1.2816) and how much
they changed at day 7.
**c.** Mean time course in the 3 organoids for the 3 strongest up and down pathways.
**d.** TCGA-OV patients given carboplatin or cisplatin, split at the median of the signature score (score = mean z of the up
pathways minus mean z of the down pathways). n = 375, deaths = 229, log-rank p = 0.539, HR (high vs
low score) = 1.08 (0.84-1.41).
**e.** Same split, progression-free interval (TCGA-CDR): n = 374, events = 283, p = 0.273,
HR = 1.14 (0.90-1.44).
**f.** Score in patients who progressed within 6 months (platinum-resistant proxy, n = 23) vs later (n = 326):
medians 0.257 vs 0.001, Mann-Whitney p = 0.35.
**g.** Score by BRCA1/2 somatic mutation (open data; germline BRCA status is controlled access): mutant n = 17, wild-type
n = 271, medians -0.008 vs 0.042, p = 0.683.
**h.** Similarity of the day-7 change between models.

What this can answer: which pathways olaparib moves in ovarian organoids, how fast, and whether they move back after the
drug is removed; and whether tumours that already 'look like' an olaparib-treated organoid differ in platinum outcome.
What it cannot answer: which organoid or patient responds to olaparib (no response values, and no TCGA-OV patient got
olaparib), so no organoid-trained biomarker and no olaparib-treated survival test is possible.
