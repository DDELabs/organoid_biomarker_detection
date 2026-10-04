# GSE9891

- GEO accession: GSE9891
- Paper: Tothill RW et al. 2008 Clin Cancer Res, PMID 18698038
- Cancer: ovarian (serous/endometrioid); setting: adjuvant
- Platform: GPL570
- Samples: 267; genes: 20838
- Responders (responder==1): 146; non-responders: 83; unlabelled: 38
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9891
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9891/
- https://bioconductor.org/packages/release/data/experiment/src/contrib/curatedOvarianData_1.50.0.tar.gz

## Downloaded files (sha256)

- `GSE9891_series_matrix.txt.gz`: bc6d88539b6c77d3bab3b51410b096b05bf91115d3c0eca890c414f605ac2ab6
- `GPL570.annot.gz`: d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394
- `curatedOvarianData_1.50.0.tar.gz`: dba53fbc9a116770628fb422ad5b81620f349f12780ab198bd1048dfc2adebfb

## Processing

- Series matrix GSE9891_series_matrix.txt.gz (54621 probes); already log scale (range 2.22..14.81).
- Probes mapped with GPL570.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 20838 genes.
- Clinical annotation from Bioconductor curatedOvarianData 1.50.0 (GSE9891_eset phenoData; Tothill et al. clinical data curated by Ganzfried et al. 2013); borderline/LMP tumours (n=18) excluded.
- arm from pltx/tax flags (platinum agent not specified -> 'PLATINUM'; taxane assumed paclitaxel).
- rfs/os from days_to_tumor_recurrence / days_to_death (/30.4375) with recurrence_status / vital_status.
- responder is a DERIVED platinum-response proxy (no RECIST/CA-125 response deposited): 1 = recurrence-free >=12 months after surgery, 0 = recurrence <12 months, NA if censored <12 months or not platinum-treated.
- Array hybridisation batch (date) kept in 'batch'.

## Arms

- platinum+taxane: n=194 (responders 120/183 labelled)
- platinum: n=49 (responders 26/46 labelled)
- no platinum: n=21 (responders 0/0 labelled)
- unknown: n=3 (responders 0/0 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.140, PC2 0.078, PC3 0.054

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.361 | 0.393 | 0.453 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE9891`.
