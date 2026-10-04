# PRISM Repurposing secondary screen + CCLE RNA-seq

**Paper**: Corsello SM et al. *Discovering the anticancer potential of non-oncology drugs by systematic viability
profiling.* Nat Cancer 2020;1:235-248. doi:10.1038/s43018-019-0018-6. Expression: Ghandi M et al. Nature 2019
(CCLE) / DepMap Public 21Q4.

**Accessions / URLs**: PRISM Repurposing 19Q4 secondary screen (figshare article 9393293, file 20237739);
DepMap Public 21Q4 (figshare article 16924132: `CCLE_expression.csv` file 31315882, `sample_info.csv` file 31316011).
The DepMap portal itself sits behind a Cloudflare challenge here, so the 21Q4 figshare release
was used instead of the newer `OmicsExpressionProteinCodingGenesTPMLogp1.csv`.

| file | URL | sha256 |
|---|---|---|
| `secondary-screen-dose-response-curve-parameters.csv` | https://ndownloader.figshare.com/files/20237739 | `88d1013506e0cd6f191a51c5f3fdd3fb2be54f8afb4e19a5d1f8538e81fbfec8` |
| `CCLE_expression_21Q4.csv` | https://ndownloader.figshare.com/files/31315882 | `94a1ce7231ce6f839f13a4140ffa5bd6820755bb9189b7783c1eef7026994413` |
| `sample_info_21Q4.csv` | https://ndownloader.figshare.com/files/31316011 | `6c4b4d9c5fec5fe21a6d721beaa85ba33f71d2e983b1fdef65fce616a78b59e7` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Derivation
* Response: `secondary-screen-dose-response-curve-parameters.csv` (701004 curves) -> kept curves with
  `passed_str_profiling == TRUE` and a fitted `auc`. For a cell line x compound (broad_id) measured in several
  screens one curve is kept, preferring MTS010 > MTS006 > MTS005 > HTS002 (as in DepMap's secondary AUC matrix).
  Compound `name` upper-cased = `drug`; if several broad_ids share a name, the median AUC is reported
  (`n_compounds`, `broad_id`, `screen_id` columns record this).
* `metric` = AUC of the fitted 4-parameter log-logistic viability curve over the screened dose range
  (8 doses, 3-fold, max 10 uM; the curve is normalised so AUC = 1 is no effect; values can exceed 1).
  **Lower = more sensitive** (no conversion needed). `lineage` column copied from DepMap sample_info.
* Expression: `CCLE_expression.csv` (RSEM log2(TPM+1), protein-coding, columns `SYMBOL (ENTREZ)`) -> Entrez IDs
  mapped to current HGNC symbols (NCBI gene_info, protein-coding only), restricted to cell lines with PRISM
  secondary data. Units: **log2(TPM + 1)**.
* `samples.tsv`: DepMap annotations (CCLE name, lineage, subtype, COSMIC/Sanger IDs) for every PRISM line.

## Sample-ID matching
Both tables use DepMap IDs (ACH-xxxxxx). PRISM lines (STR-passed): 480;
with CCLE RNA-seq: **476** (expression columns = 476). Genes: 19115. Drugs: 1448.
