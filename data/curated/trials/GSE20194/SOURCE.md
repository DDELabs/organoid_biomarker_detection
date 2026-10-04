# GSE20194

- GEO accession: GSE20194
- Paper: Popovici V et al. 2010 Breast Cancer Res (MAQC-II), PMID 20064235
- Cancer: breast; setting: neoadjuvant
- Platform: GPL96
- Samples: 278; genes: 12502
- Responders (responder==1): 56; non-responders: 222; unlabelled: 0
- Randomised: False; control/comparator arm: False

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20194/

## Downloaded files (sha256)

- `GSE20194_series_matrix.txt.gz`: ad0467234c5b7440f63e4c6dd96b17b4e930ce155985d8c421641b4f74ebb7cb
- `GPL96.annot.gz`: 88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4
- `GSE20194_MDACC_Sample_Info.xls.gz`: d489e24a8c38f80e9c4f7e55ec82f824c771f28ac59978c9523d8728dac1ec13

## Processing

- Series matrix GSE20194_series_matrix.txt.gz (22283 probes); already log scale (range -3.16..20.07).
- Probes mapped with GPL96.annot.gz; probes with no / multiple symbols dropped; median over probes per symbol -> 12502 genes.
- MAQC-II breast set (MDACC). Characteristics parsed per sample by key (GEO rows are misaligned).
- Treatment codes: T = paclitaxel, X = capecitabine (TXFAC: assumed paclitaxel+capecitabine then FAC), FAC/FEC, H = trastuzumab, /HT = hormone therapy (not listed); drugs derived from code, FAC->FEC when the free-text comment says FEC. Most patients received TFAC (non-randomised).
- responder = 1 for pCR, 0 for RD.
- Overlaps with GSE20271/GSE22093 MDACC samples are possible - deduplicate before pooling.
- Supplementary GSE20194_MDACC_Sample_Info.xls.gz downloaded for provenance (not needed: same fields as series matrix).

## Arms

- TFAC: n=210 (responders 45/210 labelled)
- TFEC: n=33 (responders 1/33 labelled)
- TXFAC: n=9 (responders 1/9 labelled)
- TH/FAC: n=6 (responders 2/6 labelled)
- FAC: n=3 (responders 0/3 labelled)
- TH/FEC: n=2 (responders 2/2 labelled)
- FECT: n=2 (responders 0/2 labelled)
- Tonly: n=1 (responders 0/1 labelled)
- FEC: n=1 (responders 0/1 labelled)
- TFAC/HT: n=1 (responders 0/1 labelled)
- FACT: n=1 (responders 1/1 labelled)
- FACT+XRT/X: n=1 (responders 0/1 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.165, PC2 0.063, PC3 0.037

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.511 | 0.211 | 0.436 | ok |
| arm | NA | NA | NA | not computable |

Built by `data/curated/trials/build_geo_trials.py GSE20194`.
