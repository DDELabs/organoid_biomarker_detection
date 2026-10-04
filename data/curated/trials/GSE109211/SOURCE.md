# GSE109211

- GEO accession: GSE109211
- Paper: Pinyol R et al. 2019 J Hepatol, PMID 30108162
- Cancer: hepatocellular carcinoma; setting: adjuvant
- Platform: GPL13938
- Samples: 140; genes: 20702
- Responders (responder==1): 42; non-responders: 98; unlabelled: 0
- Randomised: True; control/comparator arm: True

## URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109211
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109211/

## Downloaded files (sha256)

- `GSE109211_series_matrix.txt.gz`: 589134c7c62500439c396de8b2da6a7f34456cb331375d40caa1e355769309a4
- `GPL13938.soft.txt`: a50dca86c71de09110cd36f5a21abc77eda29aafcc76f6d4c0cd93c67eac282b

## Processing

- Series matrix GSE109211_series_matrix.txt.gz (29285 probes); linear values (max 18897) -> log2 (floor 34.3532).
- Probes mapped with GPL13938.soft.txt; probes with no / multiple symbols dropped; median over probes per symbol -> 20702 genes.
- Phase 3 STORM trial (adjuvant sorafenib vs placebo after resection/ablation of HCC); FFPE tumour, Illumina WG-DASL HumanHT-12 v4 (GPL13938).
- response/responder from the deposited 'outcome' (responder vs non-responder as defined by Pinyol et al.; recurrence-based - see paper). Labels exist in both arms, so the placebo arm is a prognostic control.
- No survival times in GEO.

## Arms

- placebo: n=73 (responders 21/73 labelled)
- sorafenib: n=67 (responders 21/67 labelled)

## QC: PCA of top-2000-variance genes (complete genes, centred)

PC variance explained: PC1 0.478, PC2 0.122, PC3 0.033

| label | PC1 AUC | PC2 AUC | PC3 AUC | flag (abs(AUC-0.5)>0.35) |
|---|---|---|---|---|
| responder | 0.037 | 0.515 | 0.439 | FLAG |
| arm | 0.503 | 0.516 | 0.573 | ok |

Built by `data/curated/trials/build_geo_trials.py GSE109211`.
