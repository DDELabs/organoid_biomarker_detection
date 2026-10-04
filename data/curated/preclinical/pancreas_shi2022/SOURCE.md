# Pancreatic cancer PDOs (Shi et al. 2022)

**Paper**: Shi X, Li Y, Yuan Q, et al. *Integrated profiling of human pancreatic cancer organoids reveals chromatin
accessibility features associated with drug sensitivity.* Nat Commun 2022;13:2169.
doi:10.1038/s41467-022-29857-6 (PMC9023604, open access).

**Accessions**: RNA-seq GEO **GSE194249** (`GSE194249_PDPCOs_FPKM.txt.gz`); drug screen = Supplementary Data 8
(normalised AUC, 39 exocrine PDPCOs x 59 compounds + 5 chemotherapies) with compound list in Supplementary Data 7,
fetched as the Europe PMC supplementary bundle.

| file | URL | sha256 |
|---|---|---|
| `GSE194249_PDPCOs_FPKM.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194249/suppl/GSE194249_PDPCOs_FPKM.txt.gz | `479517fbed1d8f5bb9d141b6b932b6ab56ba2ba72d4394d13aea2a3da0a50f27` |
| `PMC9023604_supplementary.zip` | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9023604/supplementaryFiles | `890fc937488f5bbbb017c9c79b54c7a18c1565a355f905e6dd6a7227c95f636d` |

## Derivation
* Response: Supplementary Data 8 "Normalized AUC" (area under the viability curve, normalised to the dose range;
  0-1). **Lower = more sensitive** (no conversion). Rows are Selleck catalogue numbers mapped to names with
  Supplementary Data 7 (`catalog` column kept). Chemotherapy abbreviations: GEM gemcitabine, 5-FU fluorouracil,
  PTX paclitaxel, OXA oxaliplatin, IRI irinotecan.
* Expression: FPKM -> log2(FPKM + 1); GEO gene symbols -> current HGNC symbols (NCBI gene_info symbol or unique
  synonym), protein-coding only; duplicates -> highest-mean row. Units: **log2(FPKM + 1)**.

## Sample-ID matching
Organoid IDs (`CAS-DAC-1`, `CAS-IPMN-1`, `CAS-NEN-*` ...) identical in both tables. Expression models:
87 (includes neuroendocrine/normal lines without drug data); screened: 39;
**overlap n = 39**. Genes: 14661; drugs: 63.
