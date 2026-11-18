**Figure G16. Does adding organoid data to a patient model improve prediction in unseen patients?**
Each row is one cancer x drug unit and outcome type (survival: C-index; response: AUROC; both 0.5 = chance). Dot = the
difference between two models in held-out patients; line = 95% paired bootstrap interval (patients resampled within
each cohort, the same resample for both models, 2,000 times). (a) joint organoid+patient model minus patient-only
model; (b) joint model minus organoid-only model; (c) as (a) after letting the inner cross-validation choose the lab
weight from a wider grid (w = 0, 0.05, 0.1, 0.25, 0.5, 0.75; w = weight of the lab-model loss, 0 = patients only).
A line that crosses the dashed zero line means the data cannot tell the two models apart.
