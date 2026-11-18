# Fig. 5 equivalents for other cancers, Supp. Fig. 12 and Supp. Fig. 18

The idea (paper Fig. 5): a drug model should call tumours with a known resistance (or sensitivity) mutation more
resistant (or more sensitive), even without survival data. Score = rank-1 pathway x its Ridge weight (higher =
predicted more resistant), in all TCGA tumours with open somatic mutation calls (Xena GDC hub). One-sided
Mann-Whitney test in the expected direction. Dashed line = p 0.05.
a. Mutation pairs with an established clinical meaning: cetuximab vs KRAS (codons 12, 13, 61, 146) and BRAF V600E in
colon + rectum (mutants expected more resistant); erlotinib vs activating EGFR mutations (L858R, L861Q, G719X, exon 19
deletions) in lung adenocarcinoma (mutants expected more sensitive); platinum vs somatic BRCA1/2 mutations in ovarian
cancer (expected more sensitive; germline BRCA status is controlled access, so only somatic mutations are used).
No widely accepted mutation marker exists for gemcitabine (pancreas, sarcoma), sorafenib (liver) or dacarbazine
(melanoma), so no such test is made there. Bladder ERCC2: nar/02_organoid_only/BLCA_lee2018/extras.
b. Supp. Fig. 12: colon cetuximab model trained on the 19 van de Wetering organoids (current data) with other feature
sets, tested on BRAF V600E. The network rank-1 pathway refit here is Gastrin creb signalling pathway via pkc and mapk (same as the organoid-only agent).
c. Supp. Fig. 18: the network model refit after leaving 1-4 organoids out (100 random draws each); boxes show the
spread of the BRAF V600E p. Share with p < 0.05 and share keeping the same rank-1 pathway:
tables/SuppFig18_leave_k_out_summary.tsv. Numbers: tables/Fig5_mutation_concordance_other_cancers.tsv,
tables/SuppFig12_cetuximab_BRAF_feature_sets.tsv.
