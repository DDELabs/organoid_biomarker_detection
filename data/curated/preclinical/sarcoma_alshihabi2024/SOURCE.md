# Sarcoma PDOs (Al Shihabi et al. 2024)

**Paper**: Al Shihabi A, Tebon PJ, Nguyen HTL, et al. *The landscape of drug sensitivity and resistance in
sarcoma.* Cell Stem Cell 2024;31(10):1524-1542.e4. doi:10.1016/j.stem.2024.08.010.

**Accessions**: Synapse syn61892224 (drug table), syn61894695/syn61894699 (omics); harmonised in
**CoderData 2.1.0** (`sarcpdo_*`, figshare files 53956238, 53956244, 53956235, 53956247).

| file | URL | sha256 |
|---|---|---|
| `coderdata_sarcpdo_experiments.tsv.gz` | https://ndownloader.figshare.com/files/53956238 | `d381d9f1195e226e1f808aec504b100bba75856e2223dae140f54b3fbe2dffe1` |
| `coderdata_sarcpdo_samples.csv` | https://ndownloader.figshare.com/files/53956244 | `e750c47f1eeb2cb5d97181dee7116a56b461d207673974a50554e9c5285d03ff` |
| `coderdata_sarcpdo_drugs.tsv.gz` | https://ndownloader.figshare.com/files/53956235 | `a2e8a26dbcc832b933632253f8d70d25c0bb36a593444ef14923c1791212891e` |
| `coderdata_sarcpdo_transcriptomics.csv.gz` | https://ndownloader.figshare.com/files/53956247 | `a2582e095960f357eb49a20ca4fb1892a4df03279342fad599e8ca81ba4463d8` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Derivation and caveats
* Response: the published per-drug **viability score** (`Viability_Score` in syn61892224; CoderData calls it
  `published_auc` but notes that full dose-response data were not released). It is % viability relative to
  vehicle, so `metric = viability_pct`, *not* an AUC. Higher = more viable, so **lower = more sensitive** without
  conversion. Values > 100 occur (growth above control).
  CoderData attaches these scores to the `_Tumor` sample IDs, but the screens were run on the PDOs of the
  same patient. They are keyed here by patient/model ID (`SARC0065`, `SARC0139_1`, ...).
* Expression: CoderData transcriptomics (TPM) of the **organoid** samples only -> log2(TPM + 1); Entrez -> HGNC
  symbols, protein-coding. Units: **log2(TPM + 1)**. Histology per model in `samples.tsv`.
* Drug names: DrugBank via InChIKey, then synonyms.

## Sample-ID matching
Organoid RNA-seq models: 16; screened models: 17; **overlap n = 15**.
Genes: 17851; drugs: 33.
