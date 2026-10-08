# CGGA693

- Accession: CGGA mRNAseq_693
- Paper: Zhao Z et al. 2021 Genomics Proteomics Bioinformatics (CGGA), PMID 33662628
- Cancer: glioma (WHO II-IV; incl. GBM); setting: primary/recurrent, adjuvant TMZ vs none
- Platform: RNA-seq (Illumina HiSeq)
- Samples: 693; genes: 16982
- Responders (responder==1): 0; non-responders: 0; unlabelled: 693
- Randomised: False; control/comparator arm: True

## URLs

- http://www.cgga.org.cn/download.jsp
- https://zenodo.org/records/8193658
- https://github.com/JackWJW/LGG_Prognosis_Prediction

## Downloaded files (sha256)

- `CGGA.mRNAseq_693.csv`: 32cc37e3cfd4a8cbb37c55ecd31a3745d423ccca7156190eae8b8db4bfb17f28
- `CGGA.mRNAseq_693_clinical.20200506.txt`: deca8eb6529903ac314640714dcda9b3ca459507cb9737d5e7da97053ce55bd0

## Processing

- CGGA mRNAseq_693 (release 20200506). cgga.org.cn is not reachable from the build environment (egress allow-list), so mirrors were used: expression = Zenodo record 8193658 'CGGA.mRNAseq_693.csv' (samples x genes; values are the CGGA RSEM FPKM, 23987 genes, despite the record text mentioning z-normalisation), clinical = the original CGGA clinical file mirrored in GitHub JackWJW/LGG_Prognosis_Prediction.
- log2(FPKM + 1); symbols restricted to UniProt gene names -> 16982 genes.
- No response endpoint: responder = NaN. arm = 'TMZ' (Chemo_status 1) vs 'no TMZ' (0); NA if unknown. Treatment was not randomised (TMZ given by clinical indication; confounded by grade/era/IDH) - use as a treatment-interaction cohort with OS, adjusting for grade/IDH/MGMT/radio.
- os_months = OS days / 30.4375; os_event = Censor (1 = dead). cancer = 'glioblastoma' for GBM/rGBM/sGBM histology else 'lower-grade glioma'. Primary and recurrent tumours included (prs_type); patients may contribute >1 sample.

## Arms

- TMZ: n=486 (responders 0/0 labelled)
- no TMZ: n=161 (responders 0/0 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.438, PC2 0.136, PC3 0.087

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (AUC<0.3 or AUC>0.7) |
|---|---|---|---|---|
| responder | NA | NA | NA | not computable |
| arm | 0.492 | 0.599 | 0.498 | ok |

Built by `data/curated/trials/build_extra_trials.py CGGA693`.
