# Pancreatic cancer PDOs (Tiriac et al. 2018; HCMI)

**Paper**: Tiriac H, Belleau P, Engle DD, et al. *Organoid Profiling Identifies Common Responders to Chemotherapy
in Pancreatic Cancer.* Cancer Discov 2018;8(9):1112-1129. doi:10.1158/2159-8290.CD-18-0349.

**Accessions**: dose-response from Tiriac 2018 (AACR figshare 39996295, refit by CoderData); RNA-seq of the same
organoids from the NCI Human Cancer Models Initiative (GDC, HCMI). Both harmonised in **CoderData 2.1.0**
(`pancpdo_*`, figshare files 53779676, 53779637, 53779610, 53779682).

| file | URL | sha256 |
|---|---|---|
| `coderdata_pancpdo_experiments.tsv.gz` | https://ndownloader.figshare.com/files/53779676 | `203fc9572f975a061ee4972ed1326db6ac3fcd61d79eb58f451da95baba3f947` |
| `coderdata_pancpdo_samples.csv` | https://ndownloader.figshare.com/files/53779637 | `fc51cce27d27fa02ae6741b4f427d114f1eda345a61bbd028bf2c07a02b7229e` |
| `coderdata_pancpdo_drugs.tsv.gz` | https://ndownloader.figshare.com/files/53779610 | `858de031c6e354373e956f77b6e0a9d16370bd48508401a6e86c6a6e9f05d911` |
| `coderdata_pancpdo_transcriptomics.csv.gz` | https://ndownloader.figshare.com/files/53779682 | `3e71f737906648be476f7241a6f474b4a73ac4c1a160b7bfd3d00a1d836c5339` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Derivation
* Response: CoderData `fit_auc` (AUC of a fitted Hill curve on fraction viability, 0-1, 120 h assay).
  **Lower = more sensitive**. Drug names from DrugBank via InChIKey (SN-38 kept as SN-38, the active irinotecan
  metabolite; `1-ohp` = oxaliplatin).
* Expression: CoderData transcriptomics (GDC STAR-count TPM of the HCMI organoid RNA-seq) -> log2(TPM + 1);
  Entrez -> HGNC symbols, protein-coding only. Units: **log2(TPM + 1)**.

## Sample-ID matching
Models are named by the Tiriac organoid IDs (`hF32`, `hM1A`, ...), taken from CoderData's `experimentId` other_id.
Drug-screened organoids: 58; with RNA-seq: 46; **overlap n = 36**.
Genes: 18975; drugs: 5.
