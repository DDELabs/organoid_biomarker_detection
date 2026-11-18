# Use your own organoids

This folder lets you drop in your own organoid expression and drug-response files. The code checks them, scores the
Reactome pathways the rest of this repo uses, writes them in the same layout as `data/curated/preclinical/<set>/`, and
then runs the paper's two methods on them.

## One command

```
python3 nar/21_user_organoids/code/import_organoids.py <input_folder> <set_name> --cancer colorectal
```

- `<input_folder>`: a folder that holds your files (see below).
- `<set_name>`: a short name without spaces, for example `mylab_crc2025`. It must not match a set in `data/curated/preclinical/`.
- `--cancer`: colorectal, bladder, pancreas, liver, ovarian, lung, breast, sarcoma, stomach, kidney, melanoma ... or a TCGA code such as `COAD`.
  If you leave it out, the `cancer_type` column of `samples.tsv` is used.

Useful options: `--import-only` (stop after the check), `--drugs FLUOROURACIL,OXALIPLATIN` (only these drugs),
`--fast` (1,000 bootstraps and 100 permutations instead of 10,000 and 1,000), `--metric AUC`, `--direction lower|higher`,
`--value-type counts|tpm|log`, `--features shortlist|all|both`.
It uses one process on one CPU. The slow part is the paper method (about 2.5 minutes per drug and cohort with `--fast`; much longer without it).

## What to put in the folder

| File | Required | What it contains |
|---|---|---|
| `expression.*` | yes | One expression table. |
| `drug_response.*` | yes | The drug results. |
| `samples.*` | no | One row per sample: `sample`, `patient`, `cancer_type`, `organoid_line`, `replicate`. |

Accepted file types for all three: `.tsv`, `.csv`, `.txt`, `.xlsx`, and any of them gzipped (`.gz`). The first sheet of an Excel file is read.
The file name only needs to contain the word (`expression`, `drug_response`, `samples`); if it does not, I guess from size and say so in the report.

### Expression table
- Genes in rows and samples in columns, **or** samples in rows and genes in columns. I work out which one.
- The first column holds the labels. Gene labels may be gene symbols (`TP53`) or Ensembl gene ids (`ENSG00000141510`, with or without `.7`).
  Entrez numbers and probe ids are not accepted (no offline map): convert them first.
- Values may be raw counts (converted to log2 counts per million), TPM or FPKM (converted to log2(x+1)), or log values (kept).
  Detection: whole numbers with a maximum above 100 = counts; non-negative with a maximum above 50 = TPM/FPKM; otherwise log. Override with `--value-type`.
- Repeated gene labels are combined by the median. Excel date-mangled genes (`1-Sep`) are dropped and counted.

### Drug-response table
- **Long**: one row per sample and drug. Columns `sample`, `drug`, `value` (or `auc` / `ic50` / `response` ...), optionally `metric`.
- **Wide**: first column = sample, other columns = drugs (or the reverse; I flip it using the sample names). Put the metric in the file
  name (`drug_response_AUC.csv`) or give `--metric`.
- Drug names are matched to the repo's drug-name map (synonyms and brand names, for example `5-FU` becomes `FLUOROURACIL`). The original name is kept in `response.tsv`.
- IC50, GI50 and EC50 are converted to natural log, as in the paper. AUC and viability are used as they are.
- Direction: AUC, IC50, viability = lower means more sensitive. Names such as inhibition or kill = higher means more sensitive.
  If the name is unknown I assume lower and tell you in the report; set `--direction` to change it.
- Several metrics in one table: one is kept (AUC first), the rest go to `response_other_metrics.tsv`. Metrics are never mixed.

### Samples
Sample names are matched between files ignoring case, spaces and dashes (`P-19 TA` = `P19TA`). Names that do not match are listed in the report.
If `samples.*` has an `organoid_line` column and several samples share a line, their expression is combined (median) into one column per line.

## What comes out

Everything is written inside this folder; nothing else in the repo is changed.

`imports/<set_name>/` (same layout as `data/curated/preclinical/<set>/`):
- `IMPORT_REPORT.md`: what I found, what I changed, warnings, problems. **Read this first.**
- `expression.tsv.gz` (genes x samples, symbols, log scale), `response.tsv`, `samples.tsv`, `gsm_map.tsv`, `SOURCE.md`
- `pathway_scores.tsv.gz`: samples x Reactome pathways, scored with the same rank score as the patient table. 673 pathways exist; you get as many as have at
  least 5 genes in your file (for a 10,912-gene table, 533).
- `organoid_set.json`: settings, detected formats, sha256 of your input files.

`runs/<set_name>/`:
- `ANALYSIS_SUMMARY.md`: plain-language tables.
- `organoid_only/`: paper method (reuses `nar/02_organoid_only/code`). `organoid_only_summary.tsv` has hazard ratio with 95% CI, log-rank p, number of
  treated patients and events, bootstrap and permutation p, leave-one-organoid-out correlation. Only drugs that have an entry in the paper's
  network-proximity table can be run (oxaliplatin uses the cisplatin entry, SN-38 uses irinotecan). Others are listed in `drugs_not_analysed.tsv`.
- `organoid_plus_patient/`: organoid + patient method (reuses `nar/03_organoid_plus_patient/code`) for each cancer-drug unit defined there
  (colorectal, bladder, pancreas, liver, ovarian, lung, sarcoma), with your organoids as the lab models. `joint_summary.tsv` and `tables/`.
  Units that need a drug your organoids do not have, or that have fewer than 4 organoid responses, are listed in `units_not_run.tsv`.

The TCGA tables are downloaded once into `_cache/` when a cancer has not been used before (needs internet). This folder is safe to delete.

## Tests

```
python3 nar/21_user_organoids/tests/test_loader.py
```
It writes `data/curated/preclinical/colorectal_vandewetering2015` in three user formats and checks that the imported pathway scores equal the original
ones (see `TEST_RESULTS.md`), plus unit tests for orientation, Ensembl ids, counts / TPM / log, wide / long and name mismatches.

## Limits to keep in mind
- With fewer than about 15 organoids per drug, models are unstable. The report warns below 8.
- Rank scoring needs about 10,000 or more genes to be comparable with patient data. A small gene panel gives a warning.
- Patient cohorts come from TCGA and the repo's curated GEO table; they depend on the cancer you name.
