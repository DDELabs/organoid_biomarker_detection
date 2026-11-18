**Paper Fig. 3 (b-f) replica for CRC PDX (Gao) | fluorouracil.** Each panel trains one model on the 41 lab models and applies it to the
129 TCGA patients given fluorouracil (25 deaths). Patients are split at the median predicted resistance;
blue = predicted responders, red = predicted non-responders; log-rank p and hazard ratio (non-responders vs responders, 95% CI).
(a) the paper's network method: the network-proximal pathway with the largest ridge weight (hyaluronan uptake and degradation).
(b) ridge on every gene. (c) ridge on every Reactome pathway. (d) ridge on the drug targets and their direct STRING neighbours.
(e) ridge on genes picked by the Bolis et al. correlation filter (Spearman p < 0.05 with response in all 10 random halves of
the lab models; 10 genes). (f) MOLI-style deep network (encoder + classifier, cross-entropy + triplet loss),
lab models split at the median response; best of 32 hyper-parameter settings by 3-fold CV AUC (0.61; Supp Table 4).
