# survlearn: patient-trained drug-specific pathway model (survival + treatment + response notes)

Code: `robust/obd/survlearn.py`, `robust/run_survlearn.py`, `robust/run_survlearn_external.py`.

## Design
- Patients: 10,032 TCGA patients (33 cancers) with pathway scores and progression-free interval (3,419 events).
- Features: 673 Reactome rank scores, standardised within cancer. Per drug, only the network-proximal shortlist
  (STRING v12, z <= -1.2816; carboplatin/oxaliplatin use cisplatin's targets; drugs without a protein target use all pathways).
- Survival model per drug: Cox stratified by cancer; age, stage, purity, other systemic therapy, a cross-fitted general
  risk score, treated, pathway main effects and treated x pathway terms (ridge; penalty chosen inside the training folds).
  Treated = drug started before the progression-free interval ended.
- Notes model per drug: ridge logistic of recorded response (CR/PR vs SD/PD) on the same shortlist.
- Combined score: mean of the two standardised scores.
- All scores are out-of-fold (5 folds by patient). Control: 10 random shortlists of the same size per drug.
- External: full-fit coefficients frozen and applied to the curated trial cohorts; AUROC within cohort, pooled.

## Results
- 37 drugs have >= 30 treated patients. 18 have enough events for the survival model, 20 enough labels for the notes model.
- General risk score (out-of-fold): HR 1.42 per SD, p = 1e-88.
- Drug-specific survival score, treated x score interaction: right direction in 9 of 18 drugs; nominal p < 0.05 in the
  right direction for 2 (carboplatin HR 0.86, p = 0.010; vincristine HR 0.53, p = 0.039); smallest BH q = 0.135.
- Median within-cancer AUROC for recorded response: 0.51 (survival score), 0.48 (notes score), 0.51 (combined).
- Network shortlist vs random shortlists: better than the random sets in 57% of comparisons (no clear advantage).
  For several drugs the random sets are not centred on zero (e.g. paclitaxel mean z -2.2, capecitabine +3.1), so single-drug
  p-values should be read against the random-set distribution, not at face value.
- External cohorts: notes-trained scores predict response for breast chemotherapy drugs (cyclophosphamide AUROC 0.67,
  doxorubicin 0.65, epirubicin 0.64, paclitaxel 0.62; n = 271-1,007), but equally well in patients who did not receive the
  drug (paclitaxel 0.64, epirubicin 0.66, doxorubicin 0.61). The fluorouracil survival score reaches 0.62 pooled, driven by
  breast cohorts (0.58-0.76) and not by two of three colorectal cohorts (0.45, 0.36).

## Conclusion
No drug-specific pathway biomarker passes multiple-testing correction in cross-validation, and the external signal is a
general chemotherapy-response signal rather than a drug-specific one.

## Known limits
- Treatment was not randomised; drugs given in combination cannot be separated.
- Time zero is diagnosis, not drug start (mild immortal-time bias for treated patients).
- The non-zero centre of the random-shortlist control for some drugs is unexplained.
