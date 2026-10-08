# Colorectal cancer PDOs - living organoid biobank (van de Wetering et al. 2015)

**Paper**: van de Wetering M, Francies HE, Francis JM, et al. *Prospective derivation of a living organoid biobank
of colorectal cancer patients.* Cell 2015;161(4):933-945. doi:10.1016/j.cell.2015.03.053 (PMID 25957691).

**Accessions**: drug screen = Table S2 (`mmc3.xlsx`, Elsevier supplementary-content CDN), sheet S2b (GDSC-style
screen at the Sanger Institute, 83 compounds, technical/biological replicates per organoid); expression = GEO
**GSE64392** (Affymetrix HuGene 2.0 ST, RMA-sketch, tumour and matched normal organoids).

| file | URL | sha256 |
|---|---|---|
| `vdWetering2015_mmc3.xlsx` | https://ars.els-cdn.com/content/image/1-s2.0-S0092867415003736-mmc3.xlsx | `9605067b702790c5e937238b40326df58d359cd1c51f2aa47e87e665d03bc4f9` |
| `GSE64392_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64392/matrix/GSE64392_series_matrix.txt.gz | `ee44313ca90eb818cf7308b366c506fa2bfb40104ca17842b37bd4be288a9b6d` |
| `GPL16686_full.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16686&targ=self&view=data&form=text | `09ad24b3a7768bc8c6ca495636687768d37e10fd0d46727a0111b5b391559582` |
| `hgnc_complete_set.txt` | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | `8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `6e2f79000e79ecb430f2613b878d514302df7f81263b5db979010d032c3235a2` |

## Derivation
* Response: **AUC** of the fitted dose-response curve as published (fraction of the tested range, 0-1; screen
  used a 5-dose, 256-fold range up to `max_conc_uM`). Median over replicate screens (`n_screens`). **Lower = more
  sensitive**. The published ln(IC50 / uM) median is kept as `ln_ic50` (not used by the loader). Drug names ->
  generic via DrugBank vocabulary (e.g. `Nutlin-3a (-)` -> NUTLIN-3A); research codes are kept.
  This replaces the degenerate CoderData 2.2 colorectal refit noted in REPORT_PRECLINICAL.md.
* Expression: GEO series-matrix RMA values (already log2), **tumour organoids only** (`p<n>t`, `p19ta/b`,
  `p24ta/b`; normal organoids `p<n>n` dropped). Transcript clusters -> RefSeq (`GB_ACC` of GPL16686) -> HGNC
  symbol (HGNC complete set, RefSeq/MANE columns), highest-mean cluster per gene, protein-coding only. Units:
  **log2 RMA**. Only ~13.5k genes map (GPL16686 lacks a gene-symbol column).

## Sample-ID matching
Drug-screen IDs `P10`, `P19a`, `P24a` = GEO titles `p10t`, `p19ta`, `p24ta` (written `P10`, `P19TA`, `P24TA`).
Expression organoids: 22; screened: 19; **overlap n = 19**.
Genes: 10912; drugs: 83. No matched patient clinical response is public.
