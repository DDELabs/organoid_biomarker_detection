# Import report: ovary_gse229119

Folder read: `_data/inputs/ovary_gse229119`

**Result: imported.** 7 organoid samples, 2 drugs, 671 pathways scored.

## Problems (must be fixed or checked)

- None.

## Warnings (imported, but please check)

- Values include negatives, so they look already centred (z-scores or log-ratios). Kept as they are; pathway scores use only ranks within each sample.
- Only 671 of the 673 Reactome pathways the repo uses could be scored (some pathways have fewer than 5 genes in your file). Scores for the missing pathways are not available.
- Only 7 organoids have drug response. Models trained on so few organoids are very unstable; treat results as exploratory.

## What I found

- Samples file: samples.tsv (0.00 MB).
- Response file: drug_response.tsv (0.00 MB).
- Expression file: expression.tsv.gz (0.54 MB).
- Expression table is laid out as: genes in rows, samples in columns (31374 genes x 7 samples).
- Values look like log-scale expression (range -9.00 to 16.05).
- Genes are gene symbols.
- 31374 genes kept; 17478 of them are in the repo's gene-name system and are used for scoring.
- 7 samples in the expression file: Hey, MCF7, Organoid 2414, Organoid 2445, Organoid 3008, OVCA432, SKOV3.
- Drug-response table is long format (sample column 'sample', drug column 'drug', value column 'value', metric column 'metric').
- Metric: ln_IC50; lower value = more sensitive.
- Samples table: 7 rows with columns sample, patient, cancer_type, line_type.
- Cancer type given by you: ovarian.
- Scored 671 Reactome pathways in 7 samples (rank score).
- 2 drugs. 0 of them have a network-proximity score in the paper's table and can go through the paper method.

## What I changed in your data

- Kept values as they are (already on a log scale).
- Took the natural log of the IC50 values (ln scale, as in the paper).

## Files written

- `expression.tsv.gz`: genes x samples, gene symbols, log scale.
- `response.tsv`: one row per sample x drug (lower or higher = more sensitive is stated in the `direction` column).
- `samples.tsv`, `gsm_map.tsv`: sample information.
- `pathway_scores.tsv.gz`: samples x 671 Reactome pathway scores.
- `organoid_set.json`: all settings and the sha256 of your input files.
