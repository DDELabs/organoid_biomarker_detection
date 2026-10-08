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

## Extra cohorts: liver, kidney, glioma, melanoma (`build_extra_trials.py`)

These cohorts were added because TCGA has almost no drug-response labels for LIHC, GBM and KIRC/KIRP. Build them with `GEO_CACHE=<dir> python build_extra_trials.py [ID ...]`. The script reuses the helpers in `build_geo_trials.py`, with these differences:
- Microarray probes are collapsed to one row per gene by keeping the probe with the highest mean.
- RNA-seq data are log2(x+1) and keep only symbols that resolve to UniProt gene names, which are the genes `obd.to_canonical` keeps.
- `clinical.tsv` has a per-sample `cancer` column.
- QC uses the stricter requested rule: `qc_flag` is set when any PC1-3 AUC vs `responder` is < 0.3 or > 0.7. Earlier cohorts used < 0.15 or > 0.85. Flagged cohorts get `usable=False`.

| cohort_id | cancer | drugs / arms | n | responders / endpoint | PC1-3 AUC | usable |
|---|---|---|---|---|---|---|
| GSE104580 | HCC | TACE (agent not stated) | 147 | 81/147 TACE responders (criterion not stated); no survival | 0.78 / 0.59 / 0.50 FLAG | False |
| GSE140901 | HCC (advanced) | anti-PD-1/PD-L1 ICI (agent not stated) | 24 | 6/24 PR vs SD/PD; PFS, OS | 0.80 FLAG | False |
| EMTAB3267 | ccRCC (metastatic) | SUNITINIB (1st line) | 53 | 19/43 PR vs SD/PD (10 "clinical benefit" = NaN); PFS | 0.71 / 0.48 / 0.47 FLAG | False |
| BRAUN2020_CHECKMATE | ccRCC (advanced) | NIVOLUMAB 181 / EVEROLIMUS 130 | 311 | 44/281 CR/PR vs SD/PD; clinical benefit; PFS; OS | 0.51 / 0.45 / 0.49 | True |
| GSE67501 | RCC (metastatic) | NIVOLUMAB | 11 | 4/11 CR/PR vs SD/NR | 0.79 FLAG | False |
| CGGA693 | glioma WHO II-IV (GBM + LGG) | TMZ 486 / no TMZ 161 / NA 46 | 693 | none (responder NaN); OS | arm PC1 0.49 | True |
| GSE7696 | GBM (primary) | RT+TMZ 43 / RT 27 | 70 | none (responder NaN); OS | arm PC1 0.48 | True |
| GSE78220 | melanoma (metastatic) | PEMBROLIZUMAB | 28 | 15/28 CR/PR vs PD; OS | 0.36 / 0.27 / 0.41 FLAG | False |
| GSE91061 | melanoma (metastatic) | NIVOLUMAB (ipi-naive 25 / ipi-progressed 26) | 51 | 10/49 CR/PR vs SD/PD; PFS; OS | 0.35 / 0.53 / 0.63 | True |

### Per-cohort notes

- **GSE104580** (GEO, GPL570, GCRMA plus batch correction by the submitter): these are pre-TACE tumour biopsies from National Cancer Centre Singapore. No paper is linked, and the submitter states neither the response criterion (likely mRECIST) nor the embolisation agent, so `drugs` is the placeholder `TACE`.
  - The flag is not a technical artefact. PC1 is a hepatocyte-differentiation axis: CYP7A1/CYP8B1/SLC10A1/GYS2 on one side and S100P/IGF2BP3/POSTN/MMP12 on the other.
  - Responders sit on the differentiated side. This is biologically plausible, so the parent may override `usable`.
- **GSE140901** (Hsu 2021 Liver Cancer): NanoString IO360 panel with 784 genes (763 after canonical mapping), so it is not genome-wide. ICI agents are not given per patient, so `drugs` is the placeholder `ICI_ANTI_PD1_PDL1`.
  - Times were deposited in weeks and converted to months.
  - The flag is expected: PC1 of an immune panel is immune infiltration, which is the paper's own finding.
- **EMTAB3267** (Beuselinck 2015 CCR; BioStudies/ArrayExpress has CEL files only):
  - Processing: the CEL files were read with Bioconductor oligo. RMA was re-implemented in numpy (convolution background, quantile normalisation, median polish on core transcript clusters) because preprocessCore's threaded C code fails in this container. The 6 normal kidneys were dropped.
  - Responder definition: `responder` is PR=1 and SD/PD=0. The 10 patients deposited only as "CLINICAL BENEFIT" are NaN in `responder`, and `clinical_benefit` gives a broader 1/0 label.
  - Survival: PFS only, no OS.
  - Checks: sex-gene and CA9 levels look as expected.
- **BRAUN2020_CHECKMATE** (Braun 2020 Nat Med Supplementary Tables S1 and S4A; open supplement, raw data controlled):
  - Expression and arms: the expression matrix is used as deposited, on a log2-like normalised scale with a median of about 24, so it is not on a TPM scale. CM-025 is the randomised nivolumab vs everolimus trial. CM-009 and CM-010 are nivolumab only.
  - Responder definitions: `responder` is CR/PR (including CRPR) = 1, SD/PD = 0, NE = NaN. `clinical_benefit` is CB=1, NCB=0, ICB=NaN.
  - Survival coding: `*_CNSR` = 1 is treated as an event, matching the paper's code (871 of 1006 PFS records are events).
  - Use: this is the best kidney cohort for treatment-interaction work.
- **CGGA693** (Chinese Glioma Genome Atlas, mRNAseq_693, release 20200506):
  - Source and expression: cgga.org.cn is blocked by the egress allow-list, so the files came from mirrors. Expression is RSEM FPKM from Zenodo record 8193658, transformed to log2(FPKM+1). The clinical file is the original CGGA file, mirrored on GitHub (JackWJW/LGG_Prognosis_Prediction).
  - Arms and endpoint: `arm` is TMZ vs no TMZ (Chemo_status). Treatment was not randomised and is confounded with grade and IDH status, so it is an OS-only treatment-interaction cohort; adjust for grade, IDH, MGMT and `radio_status`.
  - Samples: `cancer` is glioblastoma (GBM/rGBM/sGBM) or lower-grade glioma. Primary and recurrent tumours are both included (`prs_type`).
- **GSE7696** (Murat 2008 JCO, GPL570): primary GBM only (70 samples); recurrent and non-tumour samples were dropped. Patients come from the EORTC 26981/NCIC CE.3 setting, but this expression subset is not a randomised sample. Endpoint is OS only, with MGMT status included.
- **GSE78220** (Hugo 2016 Cell): FPKM converted to log2(FPKM+1). Pt27A and Pt27B are two lesions from one patient. One biopsy (Pt16) was taken on treatment. The flag comes from PC2 (AUC 0.27) with n=28.
- **GSE91061** (Riaz 2017 Cell): pre-treatment biopsies only. FPKM rows are Entrez IDs, mapped with NCBI gene_info. BOR, OS and PFS come from the authors' GitHub repository (riazn/bms038_analysis), where `*_SOR` = 1 means censored. `arm` separates ipilimumab-naive from ipilimumab-progressed patients.
- **GSE67501** (Ascierto 2016 CIR, Illumina GPL14951): n=11, archival tissue, no survival. Too small to rely on.

### Tried and not added

| source | reason |
|---|---|
| CGGA mRNAseq_325 | cgga.org.cn is not reachable. The only mirror found (Zenodo 8190371/8193658) holds a 520-gene subset, which is not genome-wide. |
| IMbrave150 / GO30140 (atezolizumab + bevacizumab HCC), IMmotion150/151, JAVELIN Renal 101 | Controlled access (EGA). |
| GSE235863 (anti-PD-1 + lenvatinib HCC) | The bulk RNA-seq is from post-treatment resections only (15 samples), so it is biased by outcome. |
| GSE202069 (anti-PD-1 HCC) | GEO has no per-sample response labels. |
| GSE109211 (sorafenib) | Already present and flagged as unusable. |
| Other sorafenib/lenvatinib, sunitinib/pazopanib and GEO GBM-TMZ cohorts | Not searched exhaustively. Candidates still to check: E-MTAB-1980 (ccRCC, prognostic only); GSE43378/GSE74187 (GBM, TMZ annotation unclear). |
