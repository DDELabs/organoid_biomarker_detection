# Test results (all numbers from `tests/test_output.txt` and the runs described here)

Source set: `data/curated/preclinical/colorectal_vandewetering2015` (22 organoids, 10,912 genes, 83 drugs).

## 1. Pathway scores from three user formats equal the original (13 tests pass)
| Format | Genes | Expression layout | Values | Drug table | Max abs. difference in pathway scores |
|---|---|---|---|---|---|
| A | symbols, 10,912 | genes in rows, csv | log | long, with metric column | 9.9e-17 |
| B | Ensembl ids with `.7` suffix, 10,478 (the ones that have an id in the repo's table) | samples in rows, tsv.gz | log | wide, metric AUC in file name, plus samples.tsv | 1.0e-16 |
| C | symbols, 10,912 | genes in rows, xlsx | TPM (not log) | wide, drugs in rows, csv | 9.9e-17 |

Each gives 22 samples x 533 pathways (673 exist; 533 have >= 5 genes in this 10,912-gene table). For B the reference was computed on the same 10,478 genes.
The imported AUC per sample and drug also equals the original (difference < 1e-9).

## 2. Paper method on the imported set equals nar/02
Format A-style import with ln IC50 as the metric, drug FLUOROURACIL (`--fast`): hazard ratio, log-rank p, number of close pathways and treated patients are identical to
`nar/02_organoid_only/_data/results/COAD_vandewetering__*__FLUOROURACIL.pkl`:
COAD HR 1.069985, p 0.865938 (129 treated); READ HR 2.123732, p 0.408821 (59); COADREAD HR 1.209941, p 0.602825 (188).

## 3. Organoid + patient method equals nar/03
Format A import (AUC), unit CRC_5FU, shortlist features: all 60 rows of the evaluation table match `nar/03_organoid_plus_patient/tables/CRC_5FU__shortlist__eval.tsv`
(max difference in C = 0, same n).
(With ln IC50 instead of AUC as the response, C differs by up to 0.019, as expected: different response metric.)
