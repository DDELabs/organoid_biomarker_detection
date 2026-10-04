# Primary liver cancer PDOs (Ji et al. 2023)

**Paper**: Ji S, Feng L, Fu Z, et al. *Pharmaco-proteogenomic characterization of liver cancer organoids for
precision oncology.* Sci Transl Med 2023;15(706):eadg3358. doi:10.1126/scitranslmed.adg3358 (PMID 37494474).

**Accessions**: data deposited on Synapse by CoderData (syn66401300-syn66401303, syn66593307; login needed for
file contents), harmonised in **CoderData 2.2.x** and taken from the full release zip of figshare article
29923646 (version 4). Individual file IDs of that release are not listable here (api.figshare.com is blocked by the
egress policy), but the whole-article download works.

| file | URL | sha256 |
|---|---|---|
| `coderdata_v2.2_article29923646_v4.zip` | https://ndownloader.figshare.com/articles/29923646/versions/4 | `f5f81bf339ec07922907b19c1c8312e1fd5b0689a9ea542b9c943f9a9b230c22` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |
| `liver_experiments.tsv.gz` (zip member) | | `a45d8bd3b9e09f44bfcc25cc299c3184039612862d497e2479636a351f396d09` |
| `liver_samples.csv` (zip member) | | `c8a5470b08e674833d9eb808fa967244866a0b128c168d034b4644c4cb9fe30e` |
| `liver_drugs.tsv.gz` (zip member) | | `2c3e1221034e282d6a318e03618b7911bae7db54acce2b2ce097510ce9d525b6` |
| `liver_transcriptomics.csv.gz` (zip member) | | `2ccfbbac4dc704016e8a919aff7f9b7b8ff4e7cb1116a338b387e25fa104dac3` |

## Derivation
* Response: CoderData `fit_auc` (AUC of a fitted Hill curve on fraction viability over the tested range, 0-1,
  72 h). **Lower = more sensitive**. The fit R^2 is kept as `fit_r2` (median 0.68; filter on it if needed).
  Drug names: DrugBank via InChIKey, else CoderData's source name / synonyms.
* Expression: CoderData transcriptomics (Synapse RNA-seq table, TPM-scale: ~1.07e6 per sample over all genes) ->
  log2(x + 1); Entrez -> HGNC symbols, protein-coding only. Units: **log2(TPM + 1)** (TPM-scale input).
  Histology per model (HCC, ICC, combined HCC-CC, hepatoblastoma) in `samples.tsv`.

## Sample-ID matching
Organoid IDs (`HCCO4`, `ICCO1`, `CHCO1`, `HBO1` ...) identical in both tables. Expression models: 62;
screened: 61; **overlap n = 61**. Genes: 16303; drugs: 73.
