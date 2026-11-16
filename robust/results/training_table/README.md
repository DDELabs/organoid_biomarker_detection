# Combined training table

Built by `robust/fetch_geo_matrices.py`, `robust/build_training_cohorts.py` and `robust/build_training_table.py`.

| File | Content |
|---|---|
| `patients.tsv.gz` | one row per tumour sample (20,766): cohort, cancer, drugs, treated, response, survival, age, sex, stage |
| `pathway_scores_tcga.tsv.gz`, `pathway_scores_geo.tsv.gz` | activity of 673 Reactome pathways per sample (within-sample rank score), same row ids |
| `cohorts.tsv` | usable numbers per cohort |
| `drug_counts.tsv` | treated patients with an outcome, per cancer group and drug |
| `skipped.tsv` | downloaded cohorts that could not be used, with the reason |

Conventions
- `drugs`: names joined by `;`. `NONE` = known to have had no drug. Empty = not recorded.
  `CHEMOTHERAPY`, `ENDOCRINE`, `TAXANE`, `ANTHRACYCLINE` = the record gives only the class, not the drug.
- `responder`: 1/0 as defined by each study (complete response at surgery, or tumour shrinkage on scans).
- Survival columns are in months; `*_event` = 1 when the event happened.
- `possible_duplicate_of`: the sample's pathway profile correlates above 0.95 with a sample of an earlier
  cohort from the same hospital. Drop these rows before pooling cohorts.
- Pathway scores are comparable within a cohort. Standardise within cohort before pooling, because cohorts
  were measured on different platforms.

Known limits
- Untreated groups outside randomised trials differ from treated ones (stage, age), so they are a rough comparison.
- Response coding of GSE123845 and event coding of GSE14333 were inferred from the data (see the builder).
