# Colorectal cancer liver-metastasis PDX, cetuximab (Isella et al. 2017 / Bertotti et al.)

**Papers**: Isella C, Brundu F, Bellomo SE, et al. *Selective analysis of cancer-cell intrinsic transcriptional traits
defines novel clinically relevant subtypes of colorectal cancer.* Nat Commun 2017;8:15107. doi:10.1038/ncomms15107
(PMID 28561063). PDX cetuximab trials: Bertotti A, et al. Cancer Discov 2011;1:508 and Nature 2015;526:263.

**Accessions**: expression GEO **GSE76402** (Illumina HumanHT-12 v4, GPL10558, 529 arrays of 244 PDX models,
lumi/loess-normalised, human-specific probes); response: Supplementary Data 4 (`ncomms15107-s5.xlsx`, Europe PMC
open-access supplement of PMC5499209), keyed by array barcode (= GEO `Sample_description`).

| file | URL | sha256 |
|---|---|---|
| `GSE76402_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76402/matrix/GSE76402_series_matrix.txt.gz | `1f2d8e09cad23d149a8fccf3cf9ce4c965414c03133febe132d5b214c2f05940` |
| `GPL10558.annot.gz` | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz | `c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45` |
| `PMC5499209_supplementary.zip` | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5499209/supplementaryFiles | `4cc169298231067d8d556d50d7ef0c98d9320da7d9aad4aecb8eaa8dfe0f28bb` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `6e2f79000e79ecb430f2613b878d514302df7f81263b5db979010d032c3235a2` |

## Derivation
* Response: **tumour volume change after 3 weeks of cetuximab** (20 mg/kg twice weekly), as a fraction of the
  volume at treatment start (-0.5 = 50% shrinkage; the paper's PR is < -0.5, PD > +0.35). **Lower = more
  sensitive**. The 6-week value (`vol_change_6w`, when available) and the paper's response class
  (`response_class`: PD 72, SD 36, PR 28, SD-PD 25, PD-mild 19, SD-PR 12) are kept. PR and SD were the paper's
  "cetuximab-sensitive". KRAS/NRAS/BRAF status from the same table is kept.
* Expression: series-matrix values (linear) -> log2; arrays of the same PDX model (regions A/B, replicate
  hybridisations) averaged on the log scale (`gsm_map.tsv`); probes -> HGNC symbols via the GPL10558 GEO
  annotation, highest-mean probe per symbol, protein-coding only. Units: **log2 (lumi-normalised intensity)**.

## Caveats
* Single agent (cetuximab, anti-EGFR antibody) only. Mostly KRAS-wild-type selected for EGFR biology; response is
  strongly driven by RAS/BRAF status.
* No matched patient clinical response is public for these models (no `patient_response.tsv`).

## Sample-ID matching
PDX model IDs (`CRC0014LM` ...) from Supplementary Data 4; GEO arrays matched by barcode. Expression models:
248; models with a cetuximab response: 192; **overlap n = 192**.
Genes: 18375; drugs: 1.
