# Bladder cancer PDOs (Lee et al. 2018)

**Paper**: Lee SH, Hu W, Matulay JT, et al. *Tumor Evolution and Drug Response in Patient-Derived Organoid Models
of Bladder Cancer.* Cell 2018;173(2):515-528.e17. doi:10.1016/j.cell.2018.03.017 (PMID 29625057).

**Accessions**: RNA-seq GEO **GSE103990** (SRA SRP118077). Drug screen raw dose-response files: Synapse
folder **syn64765430** ("Lee Bladder PDO Datasets", deposited by CoderData), curve-refit and released in
**CoderData 2.1.0** (figshare files 53956220 experiments, 53956226 samples, 53956214 drugs).

| file | URL | sha256 |
|---|---|---|
| `GSE103990_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` | https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file&acc=GSE103990&file=GSE103990_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz | `c48a0a3bc55757d1c8af26779741b465a116e9773054d60f94754a4eaa8b6d41` |
| `GSE103990_gsm_brief.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103990&targ=gsm&form=text&view=brief | `164ca193308f984b04bff819b7abeef48a9ba12fe36e69a46ecfd266465a5647` |
| `coderdata_bladderpdo_experiments.tsv.gz` | https://ndownloader.figshare.com/files/53956220 | `df6ef1f69105cde515aae44d2354a6fed48662e7e46538800297c4bac38511cc` |
| `coderdata_bladderpdo_samples.csv` | https://ndownloader.figshare.com/files/53956226 | `d4c951dad463ec6b13b6bb43d963bf2aefd10934311a8418b37b70c07352f587` |
| `coderdata_bladderpdo_drugs.tsv.gz` | https://ndownloader.figshare.com/files/53956214 | `90ad73cf60998492b365aa2996f96257c766b07d106b8e014a6e02c83432df44` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Access notes
* The author-supplied `GSE103990_Normalized_counts.txt.gz` (DESeq2 VST) answers HTTP 403 at GEO. NCBI's
  uniform re-quantification of the SRA runs (GRCh38.p13, `..._norm_counts_TPM_...`) was used instead.
* The Cell supplementary tables (PMC5890941) are not open access, and PMC is not reachable from the build host.
  Synapse file *contents* need a login, but folder metadata is public; the per-drug file names
  (`1) Gemcitabine.txt` ... `50) GSK126.txt`) were read anonymously to assign generic drug names to CoderData's
  `improve_drug_id`s (50/50 matched 1:1 by synonym).

## Derivation
* Response: CoderData `dose_response_metric == fit_auc`: the AUC of a fitted Hill curve on fraction viability over
  the tested dose range (6-day assay), 0-1 scale. **Lower = more sensitive**. CoderData labels screens as
  `<line>_Organoid_P<passage>`, `<line>_Parental`, `<line>_XenoOrganoid_P<n>` and `<line>_Xenograft`.
  Only organoid screens (`Organoid_P*`, and `Parental*`, which we take to be the parental organoid line before
  xenografting; CoderData does not define it) are kept (78 screens in total before
  filtering), and per line the **median** AUC over passages is reported (`n_screens`, `screens` columns).
  SCBO-3.2 / SCBO-11.2 etc. (organoids from recurrences) are separate lines written as `SCBO-3_2`, `SCBO-11_2`.
* Expression: NCBI TPM -> log2(TPM + 1); organoid samples only (GEO titles `SCBO-x_orgP<n>`, tumour tissue dropped);
  multiple passages of one line averaged on the log scale (`gsm_map.tsv` lists the GSMs per line). Entrez IDs ->
  HGNC symbols, protein-coding only. Units: **log2(TPM + 1)**.

## Sample-ID matching
Line IDs (`SCBO-5`, `SCBO-3_2`, ...) shared by both tables. Expression lines: 21;
drug-screened lines: 11; **overlap n = 11** (SCBO-1, SCBO-10, SCBO-11, SCBO-11_2, SCBO-2, SCBO-3, SCBO-3_2, SCBO-4, SCBO-5, SCBO-6, SCBO-8).
Genes: 19401; drugs: 50. Note that the screened passage is not always the sequenced passage.
