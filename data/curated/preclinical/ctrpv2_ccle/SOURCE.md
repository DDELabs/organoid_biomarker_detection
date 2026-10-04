# CTRPv2 (CoderData refit) + CCLE RNA-seq

**Papers**: Seashore-Ludlow B et al. Cancer Discov 2015;5:1210 and Rees MG et al. Nat Chem Biol 2016;12:109 (CTRPv2);
curves refit with PharmacoGx and harmonised by CoderData 2.1.0. Expression: DepMap Public 21Q4 CCLE RNA-seq.

**Overlap with the repo**: `robust/obd/preclinical_sets.py` already pairs CTRPv2 AUC with **GDSC RMA microarray**
expression keyed by COSMIC ID (data/external/gdsc). This set is complementary: CoderData's uniform refit
(`fit_auc` on a 0-1 scale rather than CTRP's unnormalised area) paired with **CCLE RNA-seq** keyed by DepMap ID.
Use one or the other, not both, in a pooled analysis. The original CTD2 portal
(ctd2-data.nci.nih.gov) is blocked by the build host's egress policy.

| file | URL | sha256 |
|---|---|---|
| `coderdata_ctrpv2_experiments.tsv.gz` | https://ndownloader.figshare.com/files/53779400 | `5eaaa53a7dcfea7b866a2d810f216497241124ae933046af936ededa9910cfea` |
| `coderdata_ctrpv2_samples.csv` | https://ndownloader.figshare.com/files/53779412 | `649fc0b7825d7c6e795e07246f51056fd981b2d4372ce4d4435055d5cbda04cc` |
| `coderdata_ctrpv2_drugs.tsv.gz` | https://ndownloader.figshare.com/files/53779397 | `6742a1f6dd7351a25d47ce0f8d25dd137f70d57f671ad987e8f84640ef2dd24b` |
| `CCLE_expression_21Q4.csv` | https://ndownloader.figshare.com/files/31315882 | `94a1ce7231ce6f839f13a4140ffa5bd6820755bb9189b7783c1eef7026994413` |
| `sample_info_21Q4.csv` | https://ndownloader.figshare.com/files/31316011 | `6c4b4d9c5fec5fe21a6d721beaa85ba33f71d2e983b1fdef65fce616a78b59e7` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Derivation
* Response: CoderData `fit_auc` (Hill-curve AUC on fraction viability over the tested range, 0-1, 72 h).
  **Lower = more sensitive**. Drug names: DrugBank via InChIKey; else the PRISM Repurposing name of the same
  Broad compound (BRD-K ID among CoderData's synonyms); else CoderData's source name / a DrugBank synonym / the
  shortest readable CoderData synonym; a few tool compounds keep IUPAC-like or BRD-ID names (`drug_map.tsv` lists improve_drug_id -> name). If two improve_drug_ids share a name, the
  median is used (`n_compounds`).
* Expression: as for `prism_repurposing` (CCLE 21Q4 log2(TPM + 1), protein-coding, Entrez -> HGNC), restricted to
  CTRPv2 lines.

## Sample-ID matching
CoderData sample -> DepMap ID (other_id_source == DepMap). Screened lines with a DepMap ID: 844;
with RNA-seq: **819**. Genes: 19115; drugs: 460.
