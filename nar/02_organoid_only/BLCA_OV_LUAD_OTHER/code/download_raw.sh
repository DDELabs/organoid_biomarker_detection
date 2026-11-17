#!/bin/bash
# Raw downloads for the BLCA / OV / LUAD / other organoid-only analyses (re-runnable; skips existing files).
D=$(dirname "$0")/../_data/raw; mkdir -p "$D"; cd "$D"
X=https://gdc-hub.s3.us-east-1.amazonaws.com/download
K="https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA"
get(){ [ -s "$2" ] || { curl -sS -L --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"; }; }
get https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp Survival_SupplementalTable_S1_20171025_xena_sp.tsv
for P in BLCA OV LUAD SARC LUSC BRCA PAAD SKCM STAD COAD; do
  get $K/nationwidechildrens.org_clinical_drug_${P,,}.txt nationwidechildrens.org_clinical_drug_${P,,}.txt
  get $K/nationwidechildrens.org_clinical_patient_${P,,}.txt nationwidechildrens.org_clinical_patient_${P,,}.txt
  get $X/TCGA-$P.survival.tsv.gz TCGA-$P.survival.tsv.gz
  get $X/TCGA-$P.clinical.tsv.gz TCGA-$P.clinical.tsv.gz
done
get $X/TCGA-BLCA.somaticmutation_wxs.tsv.gz TCGA-BLCA.somaticmutation_wxs.tsv.gz
get $X/TCGA-OV.somaticmutation_wxs.tsv.gz TCGA-OV.somaticmutation_wxs.tsv.gz
for P in BLCA OV LUAD SARC LUSC PAAD SKCM STAD BRCA COAD; do get $X/TCGA-$P.star_fpkm-uq.tsv.gz TCGA-$P.star_fpkm-uq.tsv.gz; done
echo done > done.flag
# also used (fetched by hand in the first session, same URLs):
#  GEO: GSE229119_log2TPM.gene.expr.table.txt.gz, GSE229119_series_matrix.txt.gz, GSE165611_RAW.tar, GSE165611 gsm brief text
#  Springer: 13046_2023_2606_MOESM1-3 (lung paper supplement) -> geo/JECCR_MOESM*.{docx,xlsx}
#  raw.githubusercontent.com/lm687/Organoids_Compositional_Analysis/6dd6895... (Vias TPM + AUC csv) -> vias/
#  cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC{1,2}_fitted_dose_response_27Oct23.xlsx -> GDSC{1,2}_fitted.tsv.gz
