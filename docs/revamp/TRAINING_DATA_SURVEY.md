# Training-data survey: tumour gene activity paired with drug and outcome

Checked 2026-11-14 from this machine. Sample counts are live GEO counts (arrays or libraries),
which can be higher than the number of patients with a usable outcome. The clinical fields of the
datasets marked "new" have not been opened yet; arms and endpoints for those come from the GEO
title and the earlier desk survey in `DATA_LANDSCAPE.md`.

## A. Patient trials with a comparison arm (best for "which drug for which patient")

| Dataset | Cancer | Comparison | GEO samples | Status |
|---|---|---|---|---|
| GSE194040 (I-SPY2) | breast | 9 drug arms vs shared paclitaxel control | 988 | new |
| GSE279390 (VESPER) | bladder | two cisplatin regimens before surgery | 540 | new |
| GSE164458 (BrighTNess) | breast | paclitaxel with/without carboplatin | 482 | new |
| GSE319641 (NeoTRIP) | breast | chemo with/without atezolizumab, before and during | 401 | new |
| GSE140082 (ICON7) | ovary | chemo with/without bevacizumab | 380 | new |
| GSE156915 (FOCUS) | colorectal | 5-FU with/without oxaliplatin or irinotecan | 361 | new |
| GSE9782 (APEX) | myeloma | bortezomib vs dexamethasone | 528 | new |
| GSE116335 (CALGB 40601) | breast | trastuzumab vs lapatinib vs both | 265 | new |
| GSE234519 (NSABP B-41) | breast | trastuzumab vs lapatinib vs both | 245 | new |
| GSE107850 (EORTC 22033) | low-grade glioma | temozolomide vs radiotherapy | 195 | new |
| GSE109211 (STORM) | liver | sorafenib vs placebo | 140 | new |
| GSE14814 (JBR.10) | lung | cisplatin/vinorelbine vs observation | 133 | new |
| GSE20271, GSE41998 | breast | with/without paclitaxel; ixabepilone vs paclitaxel | 178, 279 | already used |

## B. Large patient cohorts with a treated/untreated flag or one known regimen

| Dataset | Cancer | What is recorded | GEO samples | Status |
|---|---|---|---|---|
| GSE96058 (SCAN-B) | breast | chemo / endocrine / none, survival | 3,409 | new |
| GSE39582 | colon | adjuvant chemo yes/no, relapse | 585 | wired in code, not in the main runs |
| GSE25066 | breast | taxane-anthracycline, relapse | 508 | new |
| GSE186102 | breast | eribulin before surgery | 473 | new |
| GSE68465 | lung | adjuvant chemo flag, survival | 462 | new |
| GSE87304 | bladder | cisplatin-based chemo before surgery, response, survival | 305 | new |
| GSE14333, GSE103479, GSE104645, GSE87211 | colorectal | chemo flag or regimen | 290, 156, 193, 363 | partly wired |
| GSE163882, GSE123845, GSE240671 | breast | chemo before surgery, response | 222, 227, 122 | new |
| GSE232569 | lung | pemetrexed-platinum | 95 | new |
| GSE235067, GSE48276, GSE70691 | bladder | cisplatin-based chemo before surgery | 84, 116, 49 | new |
| TCGA (33 cancers) | all | drug records, survival | 3,603 treated | already used |
| 24 curated response cohorts | several | response to one regimen | about 2,900 | already used |

## C. Lab models with gene activity and drug testing

| Dataset | Models | Drugs | Status |
|---|---|---|---|
| CTRP, PRISM (cell lines) | about 800 / 480 | hundreds | already used |
| GDSC1 and GDSC2 (cell lines) | about 970 | about 300 | reachable, not used yet |
| Beat AML (patient leukaemia samples tested in the dish) | several hundred patients | over 100 | reachable, not used yet |
| Organoids: colon 19, bladder 9-11, liver 61, pancreas 39 and 36, sarcoma 14 | | | already used |
| New organoid sets: GSE294511 colorectal liver metastases (169 samples), GSE232170 pancreas (90), GSE65253 colon (108), GSE202263 ovary (22), GSE235912 gastric (30), GSE280749 bladder (16), GSE294532 prostate (63) | | drug data must be taken from each paper | new, drug tables unverified |
| PDX: Novartis (already used), GSE291966 triple-negative breast (92), GSE142767 breast (66) | | | new |
| Sanger Cell Model Passports organoids | 74 | no public drug or RNA data flagged | not usable |

## D. Not reachable from this machine or needing an application

- CTR-DB 2.0 (about 10,900 curated patient samples, 346 regimens): site did not open from here.
- METABRIC via cBioPortal, POG570 (bcgsc.ca), PDMR: blocked by the network allow-list.
- Hartwig, Roche trials at EGA, ORIEN, Tempus, Caris: application or contract, weeks to months.
