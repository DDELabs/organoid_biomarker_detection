# TCGA-READ (rectal adenocarcinoma; RNA-seq)

Retrieved 2026-10-04. Build: `scratchpad/work/build_geo.py::tcga_read` (same conventions as `obd.cohorts.tcga_star_fpkm_uq`
and `tcga_xena_clinical`).

| file | URL | sha256 |
|---|---|---|
| STAR FPKM-UQ log2(x+1) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.star_fpkm-uq.tsv.gz | e58ee370b8199f7398b2eecdaf65ea2989e548f6c7e0f63661e376748f9e56a8 |
| Xena clinical | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.clinical.tsv.gz | d0e4b4cd84979811f5cc8d9db84ce2e72e7e117725d2dde363618590c38836d0 |
| Xena survival (OS) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.survival.tsv.gz | 3dc70834cc28e876daeaadb197fc4c2e0f922c21e93a2a66825d8ca8e9ec00d2 |
| BCR biotab drug table | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_read.txt | 33e3be6f805a73df5b0d17f0082598d9e5e369911d99fd5a4b835b26991f5292 |

## expression.tsv.gz
Primary tumours (-01), one aliquot per patient (first sorted barcode), columns = 12-char patient barcodes. Ensembl IDs
(version stripped) mapped to symbols with the repo's `data/2017_07_31_biomart_protein_coding_genes.txt` (protein-coding
only), median per symbol -> 19,708 x 166. Values log2(FPKM-UQ+1).

## clinical.tsv (166 rows)
- `os_months` = OS.time/30.44, `os_event` = OS (Xena GDC survival); RFS not derived (NA).
- `stage` from `ajcc_pathologic_stage.diagnoses`; `age` = `age_at_index`; `sex` 1 = male.
- `chemo` = 1 if the patient has any Chemotherapy record in the biotab drug table, else NA (absence of a record is not
  evidence of no chemotherapy). `regimen` = '+'-joined drug names as entered; `fu_based` = 1 if any of 5-FU /
  fluorouracil / capecitabine (Xeloda) / FOLFOX / FOLFIRI etc. appears. `setting` = biotab `therapy_regimen`
  (ADJUVANT, NEOADJUVANT/chemoradiation, ...).
- `response` = best `treatment_best_response` over the patient's chemotherapy records mapped to CR/PR/SD/PD
  (mostly recorded for adjuvant/neoadjuvant treatment, so "CR" largely means no evidence of disease; interpret with care).
