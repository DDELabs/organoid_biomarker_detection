# EMTAB3267

- Accession: E-MTAB-3267
- Paper: Beuselinck B et al. 2015 Clin Cancer Res, PMID 25593300
- Cancer: kidney renal clear cell carcinoma (metastatic); setting: metastatic
- Platform: A-AFFY-141 (HuGene 1.0 ST; GPL6244 annotation)
- Samples: 53; genes: 19036
- Responders (responder==1): 19; non-responders: 24; unlabelled: 10
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3267
- https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/267/E-MTAB-3267/Files/

## Downloaded files (sha256)

- `E-MTAB-3267.sdrf.txt`: 851fbeeda295b3397e052cdb078c513329a51dbf59eb2341a1551c96370116e7
- `GPL6244.annot.gz`: 99ee70f87e02c75ff1c8dc57b6c40c0482c7ef2760e004301315ac34344ac6dc

## Processing

- 59 CEL files (HuGene-1_0-st, A-AFFY-141; 53 tumours + 6 normals) from BioStudies/ArrayExpress.
- Core-target PM intensities read with Bioconductor oligo (pd.hugene.1.0.st.v1); RMA done in numpy over all 59 arrays: affy-style convolution background (KDE-mode estimates), quantile normalisation, log2, median polish per transcript cluster -> 33297 transcript clusters (approximates oligo::rma(target='core')).
- Transcript clusters mapped with GPL6244.annot.gz; probe with highest mean kept per symbol -> 19036 genes. Normals dropped (53 tumours kept).
- responder: RECIST best response PR = 1, SD/PD = 0; 'CLINICAL BENEFIT' (n=10, response category not given in ArrayExpress) = NaN. clinical_benefit: PR/SD/CLINICAL BENEFIT = 1, PD = 0.
- PFS in months; pfs_event = 'progression' (1 = progressed). No OS deposited.
- Beuselinck et al. derived ccrcc1-4 subtypes on these data; the paper's primary endpoints were PFS/OS.

## Arms

- SUNITINIB: n=53 (responders 19/43 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.165, PC2 0.096, PC3 0.066

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | 0.713 | 0.476 | 0.471 | FLAG |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_extra_trials.py EMTAB3267`.
