# Import report: luad_gse165611

Folder read: `_data/inputs/luad_gse165611`

**Result: imported.** 102 organoid samples, 30 drugs, 673 pathways scored.

## Problems (must be fixed or checked)

- None.

## Warnings (imported, but please check)

- None.

## What I found

- Samples file: samples.tsv (0.00 MB).
- Response file: drug_response.tsv (0.10 MB).
- Expression file: expression.tsv.gz (3.61 MB).
- Expression table is laid out as: genes in rows, samples in columns (19861 genes x 102 samples).
- Values look like TPM/FPKM (non-negative, up to 360,832, not on a log scale).
- Genes are gene symbols.
- 19861 genes kept; 19505 of them are in the repo's gene-name system and are used for scoring.
- 102 samples in the expression file: NCCLU-009, NCCLU-012, NCCLU-015, NCCLU-017, NCCLU-019, NCCLU-025, NCCLU-027, NCCLU-041 ....
- Drug-response table is long format (sample column 'sample', drug column 'drug', value column 'value', metric column 'metric').
- Metric: LOG2FC; lower value = more sensitive (as you said).
- Samples table: 102 rows with columns sample, patient, cancer_type, cluster.
- Cancer type given by you: lung.
- Scored 673 Reactome pathways in 102 samples (rank score).
- 30 drugs. 12 of them have a network-proximity score in the paper's table and can go through the paper method: CAMPTOTHECIN, GEFITINIB, GEMCITABINE, NAVITOCLAX, PACLITAXEL, PALBOCICLIB, PF-562271, RUXOLITINIB, SELUMETINIB, TANESPIMYCIN, TRAMETINIB, XAV939.

## What I changed in your data

- Converted TPM/FPKM to log2(value + 1).
- Matched 102 drug-response sample name(s) to expression names ignoring upper/lower case, spaces and dashes (e.g. 'NCCLu-041' -> 'NCCLU-041').
- Matched 102 samples-table sample name(s) to expression names ignoring upper/lower case, spaces and dashes (e.g. 'NCCLu-041' -> 'NCCLU-041').

## Files written

- `expression.tsv.gz`: genes x samples, gene symbols, log scale.
- `response.tsv`: one row per sample x drug (lower or higher = more sensitive is stated in the `direction` column).
- `samples.tsv`, `gsm_map.tsv`: sample information.
- `pathway_scores.tsv.gz`: samples x 673 Reactome pathway scores.
- `organoid_set.json`: all settings and the sha256 of your input files.
