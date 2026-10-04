# ovarian_vias2023: high-grade serous ovarian cancer organoids (Cambridge, OV04)

Paper: Vias M, Morrill Gavarró L, et al. "High-grade serous ovarian cancer organoids as models of chromosomal instability." eLife 2023;12:e83867.
Code and data repo: https://github.com/lm687/Organoids_Compositional_Analysis (branch master, commit `6dd6895402b48ee6e9b58bfbad966d874c47433a`, cloned with `git clone --depth 1`).
Loader: `robust.obd.preclinical_sets.load_organoid_set_extra('ovarian_vias2023')`.

## Files

| file | content | sha256 |
|---|---|---|
| `expression.tsv.gz` | log2(TPM+1), 19,267 gene symbols x 14 organoids (columns are organoid IDs). TPM is summed over Ensembl genes sharing a symbol. Genes that are 0 in every sample are dropped. | `7d50233ca7cfb33f87751cb04e38c1fca37d5f9572f4f683a27304bf959274d0` |
| `response.tsv` | 196 rows (14 organoids x 14 drugs), columns `sample, drug, response, metric`; metric = `AUC_LL5_viability` | `755e61cc70efa08fb89078939478d263a2b5486bfd771904cfc4fde197c2e1a3` |

## Raw inputs (repo paths, raw.githubusercontent.com/lm687/Organoids_Compositional_Analysis/master/...)

| raw | sha256 |
|---|---|
| `RNASeq_DE_resistant_sensitive/files/20191218_ViasM_BJ_orgaBrs_tpm.csv` (24 RNA-seq libraries JBLAB-199xx, TPM) | `63d573bf8c997cd62894578e62dfdf07b973627d836b0c2f8f9a1de8a8b9a300` |
| `survival_analysis/data/20200419-AUC-organoids.csv` (14 organoids x 14 drugs, `auc_ll5.*`) | `f08077a76eec57e7fedd67f0e9ec435bc274018ba04889fb627933c56fc0c986` |
| `survival_analysis/data/20200419-AUC-organoids_PFI.csv` (paclitaxel/oxaliplatin AUC + patient platinum-free-interval class) | `0c2c699cac4f06bf42435c745a16f770fe3df124cd16e76ed0f6f2e7d53927f1` |
| `RNASeq_DE_resistant_sensitive/files/deObject_SampleGroup_sensitive_vs_resistant.RData` (sample sheet linking libraries to organoids) | `0d59fc0bf935db3515f8d87e789d471c2ae2e65d83cefe7cb91bbcb864cb79bb` |

Library-to-organoid map, read from the colData of the DESeq object above (columns SampleName / organoid):
19902→119148, 19904→54327, 19905→118976, 19906→119178, 19907→119127, 19917→23868, 19921→119058, 19925→54288, 19937→54059, 19938→151723, 19939→151761, 19940→54276, 19941→32077, 19942→151773.
The other 10 libraries (19903, 19916, 19920, 19936, 19950–19955) have no organoid label in the repo and are not used.

Drugs: AZD8186, VISTUSERTIB (AZD2014), CAPIVASERTIB (AZD5363), OLAPARIB (AZD2281), AZD0156, CERALASERTIB (AZD6738), ADAVOSERTIB (AZD1775), AZD8835, ELESCLOMOL (spelled "Elescamol" in the source), PACLITAXEL, OXALIPLATIN, DOXORUBICIN, GEMCITABINE, EPRENETAPOPT (APR-246).

## Direction check (lower = more sensitive)
The response is the area under a 5-parameter log-logistic viability curve, normalised to untreated (about 1 means no effect; values slightly above 1 occur). A lower value means more sensitive.
Check against the patients: the median paclitaxel AUC is 0.83 for organoids from platinum-resistant patients (PFI class) and 0.39 for platinum-sensitive ones, so the direction holds. Oxaliplatin AUC is flat (0.93 vs 0.94) and close to 1 for most organoids, so it gives little signal.

## Caveats
- n = 14 organoids. Three patients contribute two organoids each (OV04 466: 118976 and 119058; OV04 627: 119127 and 119148; OV04 366: 32077 and 54059). Use patient-aware resampling.
- No platinum compound other than oxaliplatin, and no carboplatin, was screened.
- Expression was quantified by Kallisto (Ensembl/GRCh37 per the repo scripts), as distributed. RNA-seq passage and drug-screen passage may differ.
- The drug screen is a single (or few) biological replicate. The repo also has replicate and passage AUC tables (`AUCbiolrepW.csv`, `AUCdifferentPassagesW.csv`); these are not used here.
