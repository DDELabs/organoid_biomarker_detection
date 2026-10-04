# GEO chemotherapy / treated cohorts - curation report

Built with `build_geo_trials.py` (all cohorts: `python build_geo_trials.py`; raw downloads cached in `$GEO_CACHE`). Each cohort folder has `expression.tsv.gz` (HGNC symbols x samples, log2; arrays median-collapsed over uniquely annotated probes using the GEO GPL annotation), `clinical.tsv` and `SOURCE.md` (URLs, sha256, derivations, PCA QC). Catalog: `CATALOG_GEO.tsv`.

QC: PCA on the top-2000-variance genes; AUC of PC1-PC3 vs `responder` (and vs `arm` for two-arm cohorts). `qc_flag` = any |AUC-0.5| > 0.35; flagged cohorts are `usable=False`.

## Obtained (20 cohorts, 2316 samples)

| cohort | cancer | setting | n | responders / labelled | drugs | randomised | control arm | endpoint | QC PC1 AUC (resp) | usable |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE72970 | colorectal (metastatic) | metastatic | 124 | 63/124 | BEVACIZUMAB;CAPECITABINE;CETUXIMAB;FLUOROURACIL;IRINOTECAN;LEUCOVORIN;OXALIPLATIN | False | False | RECIST response; PFS; OS | 0.464 | True |
| GSE109211 | hepatocellular carcinoma | adjuvant | 140 | 42/140 | PLACEBO;SORAFENIB | True | True | responder (recurrence-based, per paper) | 0.037 FLAG | False |
| GSE20271 | breast | neoadjuvant | 178 | 26/178 | CYCLOPHOSPHAMIDE;DOXORUBICIN;EPIRUBICIN;FLUOROURACIL;PACLITAXEL | True | True | pCR | 0.573 | True |
| GSE41998 | breast | neoadjuvant | 279 | 69/253 | CYCLOPHOSPHAMIDE;DOXORUBICIN;IXABEPILONE;PACLITAXEL | True | True | pCR | 0.346 | True |
| GSE32646 | breast | neoadjuvant | 115 | 27/115 | CYCLOPHOSPHAMIDE;EPIRUBICIN;FLUOROURACIL;PACLITAXEL | False | False | pCR | 0.222 | True |
| GSE20194 | breast | neoadjuvant | 278 | 56/278 | CAPECITABINE;CYCLOPHOSPHAMIDE;DOXORUBICIN;EPIRUBICIN;FLUOROURACIL;PACLITAXEL;TRASTUZUMAB | False | False | pCR | 0.511 | True |
| GSE22093 | breast | neoadjuvant | 103 | 28/97 | CYCLOPHOSPHAMIDE;DOXORUBICIN;FLUOROURACIL | False | False | pCR | 0.286 | True |
| GSE23988 | breast | neoadjuvant | 61 | 20/61 | CAPECITABINE;CYCLOPHOSPHAMIDE;DOCETAXEL;DOXORUBICIN;FLUOROURACIL | False | False | pCR | 0.379 | True |
| GSE45670 | oesophageal squamous cell carcinoma | neoadjuvant | 28 | 11/28 | CISPLATIN;VINORELBINE | False | False | pCR | 0.663 | True |
| GSE15622 | ovarian (high-grade serous) | neoadjuvant | 35 | 22/35 | CARBOPLATIN;PACLITAXEL | False | True | CA-125 response | 0.573 | True |
| GSE9891 | ovarian (serous/endometrioid) | adjuvant | 267 | 146/229 | PACLITAXEL;PLATINUM | False | False | RFS; OS; derived RFS>=12m response | 0.361 | True |
| GSE51373 | ovarian (high-grade serous) | adjuvant | 28 | 16/28 | CARBOPLATIN;PACLITAXEL | False | False | platinum sensitivity | 0.62 | True |
| GSE63885 | ovarian | adjuvant | 101 | 65/75 | CYCLOPHOSPHAMIDE;PACLITAXEL;PLATINUM | False | False | clinical response; OS; platinum sensitivity | 0.615 | True |
| GSE52219 | bladder (MIBC) | neoadjuvant | 23 | 6/23 | CISPLATIN;DOXORUBICIN;METHOTREXATE;VINBLASTINE | False | False | downstaging (<=pT1) | 0.471 | True |
| GSE169455 | bladder (MIBC) | neoadjuvant | 149 | 48/149 | CISPLATIN | False | False | pCR (pT0N0) | 0.483 | True |
| GSE42127 | NSCLC | adjuvant | 176 | 0/0 | NONE;PLATINUM | False | True | OS |  | True |
| GSE103479 | colorectal (stage II/III) | adjuvant | 156 | 0/0 | FLUOROURACIL;NONE | False | True | OS; PFS |  | True |
| GSE19860 | colorectal (advanced) | metastatic | 40 | 15/40 | FLUOROURACIL;LEUCOVORIN;OXALIPLATIN | False | False | response (R/NR) | 0.405 | True |
| GSE19862 | colorectal (advanced) | metastatic | 14 | 7/14 | BEVACIZUMAB | False | False | response (R/NR) | 0.061 FLAG | False |
| GSE62080 | colorectal (metastatic) | metastatic | 21 | 9/21 | FLUOROURACIL;IRINOTECAN;LEUCOVORIN | False | False | response (sensitive/resistant) | 0.62 | True |

## Notes and caveats

- **Randomised with comparator**: GSE20271 (FAC vs T/FAC), GSE41998 (AC->ixabepilone vs AC->paclitaxel), GSE109211 (STORM sorafenib vs placebo). Non-randomised treated-vs-untreated: GSE42127 (NSCLC adjuvant chemo vs observation, OS only), GSE103479 (CRC adjuvant chemo vs surgery, OS/PFS). GSE15622 has two single-agent comparator arms (paclitaxel vs carboplatin), allocation not stated as randomised.
- **Flagged**: GSE109211 - PC1 (48% of variance) separates responders from non-responders (AUC 0.04) in both arms, consistent with a processing/RNA-quality effect aligned with the outcome label; GSE19862 - n=14, PC1 AUC 0.06 and an unusual value scale. Both set usable=False.
- **Derived labels**: GSE9891 has no deposited response; `responder` is a proxy (recurrence-free >=12 months after surgery, platinum-treated only), clinical data from Bioconductor curatedOvarianData 1.50.0. GSE169455 responder = pT0N0; `downstaged` (<=pT1N0) also given. GSE52219 responder = downstaged to pT0/pT1.
- **Unspecified agents**: 'PLATINUM' is used where the platinum compound is not annotated (GSE9891, GSE63885, GSE42127). GSE103479 adjuvant regimen not annotated (listed FLUOROURACIL); GSE169455 lists CISPLATIN only (cisplatin-based combinations); GSE19862 lists BEVACIZUMAB only.
- **Overlap**: the MDACC breast series (GSE20194, GSE20271, GSE22093, GSE23988) may share patients; deduplicate (e.g. by expression correlation) before pooling.
- GSE63885 DFS is deposited as 0 for non-CR patients and has no event column; use OS and response.
- GSE20271 drugs include preoperative switches after non-response; use `arm` for intention-to-treat.

## Skipped

| accession | reason |
|---|---|
| GSE154524 (CALGB 40603) | GEO holds only pretreatment RNA-seq; no per-sample arm (carboplatin/bevacizumab) or pCR labels. The open-access JCO 2022 supplement (PMC9015203) has no per-patient table; labels need an Alliance/NCTN data request. |
| GSE62254 (ACRG) | Series matrix has no clinical annotation; adjuvant-chemo/survival data are only in the Nat Med 2015 supplement (not open access; nature.com unreachable from this environment). |
| GSE14210 | Only training/validation/post-treatment tags in GEO; per-patient cisplatin/5-FU response or TTP not deposited (PLoS One supplement PMC3041770 holds gene lists only). |
| GSE87304 | GEO has subtype and cT/cN only; no pathological response or survival per sample (Seiler 2017 Eur Urol, not open access). |
| GSE50081 | Survival only; adjuvant chemotherapy not annotated in GEO, so no treated/untreated split. |
| GSE57495 | Prognostic pancreatic cohort; no treatment (gemcitabine/FOLFIRINOX) annotation in GEO. |
