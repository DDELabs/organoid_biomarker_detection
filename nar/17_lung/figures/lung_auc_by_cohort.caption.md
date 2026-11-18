**Can pathway scores predict drug response in lung patients?** One row per lung cohort with a response label. Each dot is the area under the ROC curve
(AUC: the chance that a random responder gets a higher score than a random non-responder; 0.5 is a coin flip, 1 is perfect) with its 95% confidence interval.
Four ridge-regression models were each tested on the cohort without ever seeing it: one trained on the other cohorts of the same treatment group
(immunotherapy cohorts together, TCGA-LUAD and TCGA-LUSC together), one trained on all breast response cohorts, one trained on every other cohort of any cancer,
and one trained on every cancer except lung. GSE33072 has no same-group partner. GSE126044 (16 patients) and GSE207422 (24 patients, pathological response after neoadjuvant treatment) are very small.
