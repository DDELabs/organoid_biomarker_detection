# Controlled-access PDO + patient-outcome data: application plan

Prepared 2026-10-04 for DDELabs. EGA, NCBI and journal sites were blocked from the research environment, so accessions come from search-engine snippets of EGA/dbGaP pages. **Open each EGA/dbGaP page before submitting.**
Status key: **V** = accession and DAC shown in a snippet of the EGA/dbGaP page itself; **P** = accession found but DAC/contact not confirmed; **U** = not confirmed.

## 1. Dataset table

| # | Paper | Accession | DAC (id, contact) | Data types | n PDOs with patient outcome (approx.) | Why it matters for us | Status |
|---|---|---|---|---|---|---|---|
| 1 | Ooft 2019 Sci Transl Med (TUMOROID, mCRC; irinotecan, 5-FU/oxaliplatin) | **No EGA deposit.** The paper says per-patient PDO, outcome and DNA-seq data are obtained through the NKI Institutional Review Board | NKI IRB, via corresponding author E. Voest (NKI). No EGAC | DNA-seq, PDO viability (supplement), RECIST response | ~61 patients enrolled; about 10-20 evaluable per regimen | Main paired chemo-response CRC set | P (no EGA), contact U |
| 2 | Vlachogiannis 2018 Science (mCRC/GOC, phase I/II trials) | EGAS00001002784 (targeted + WGS). No ArrayExpress found | Probably EGAC00001000807 "GI Organoid Genomic Profiling ICR DAC". Contact not found | Targeted panel (80 samples), WGS (2); drug screens in supplement | ~20 (71 patients, 110 PDOs total) | Response labels for targeted agents and chemo. No transcriptome | Accession V; DAC U |
| 3 | de Witte 2020 Cell Rep (ovarian, carboplatin/paclitaxel) | EGAD00001005707 (WGS, 23 samples; EGA titles it "Ovarian cancer organoid biobank - followup") | EGAC00001000432, Div. Biomedical Genetics UMC Utrecht, **DACDBG@umcutrecht.nl** | WGS | 36 PDOs; subset from interval debulking matched to clinical response | Platinum response in ovarian cancer; same DAC as #4 and #6 | P (paper-to-EGAD link only from snippet) |
| 4 | Kopper 2019 Nat Med (ovarian biobank) | EGAS00001003073; dataset EGAD00001004387 | EGAC00001000432 (as #3); access also needs an MTA and UMCU METC sign-off | DNA + RNA BAMs | Few. Mostly in vitro platinum response, plus recurrence chemoresistance | Organoid RNA-seq for training. Outcome labels are weak | V |
| 5 | Driehuis 2019 Cancer Discov (HNSCC; RT, cisplatin) | EGAS00001003628; EGAD00001005063 (WES, 8 patients) | Probably EGAC00001000432 / DACDBG@umcutrecht.nl | WES of tumour, normal and organoids | ~7 RT patients with follow-up | Radiotherapy response labels. Small n | Accession V; DAC U |
| 6 | Sachs 2018 Cell (breast) | EGAS00001002158; EGAD00001003751 (102 samples) | "Department of Biomedical Genetics UMC Utrecht", probably EGAC00001000432 | WGS (tumour/blood/PDO), PDO RNA-seq | ~1 patient (tamoxifen) | Transcriptome + drug-screen pretraining only | Accession V; DAC P |
| 7 | Yan 2018 Cell Stem Cell (gastric) | EGAS00001003145; EGAD00001004301 (WES, 130 samples), EGAD00001004302 (RNA-seq) | EGAC00001000981, HKU Gastric Cancer Organoids Genomics DAC, S.Y. Leung, **suetyi@hku.hk** | WES, RNA-seq, drug screen (supplement) | None reported | High-quality RNA + drug-response pairs for pretraining | V |
| 8 | Tiriac 2018 Cancer Discov (PDAC) | dbGaP **phs001611.v1.p1** | dbGaP DAC (probably NCI). PI D. Tuveson, CSHL | WGS, WES, RNA-seq. Gene-expression quantification is **open** (AWS Open Data "organoid-pancreatic") | ~9 patients with retrospective clinical correlation | Pharmacotranscriptomic signatures. Expression usable now with no application | Accession V; DAC U |
| 9 | Yao 2020 Cell Stem Cell (LARC, phase III CinClare) | Irradiation curves are open on Mendeley Data (10.17632/74b5j3dwf6.1). No controlled-access seq found | Corresponding authors (Fudan): Z. Zhang, G. Hua | Dose-response, clinical response | ~80 (96 PDOs from 112 patients) | Largest paired CRT cohort | P |
| 9b | Possibly related: GSE171680 (tissue) / GSE171681 (organoid), open GEO | Open | n/a | RNA-seq, 87 matched tumour-PDO pairs, survival | 87 | Open, matched, used for network biomarkers (Transl Oncol 2025, PMID 39754813) | **U** (source study not confirmed) |
| 10 | Ganesh 2019 Nat Med (rectal, MSK) | No EGA/dbGaP. Data "on reasonable request"; tumoroid data on cBioPortal | Corresponding author K. Ganesh (MSKCC) | MSK-IMPACT, CRT response | ~7 matched | Small CRT set | P |
| 11 | Smabers 2025 Clin Cancer Res (OPTIC, mCRC) | **Hartwig Medical Foundation** databank (not EGA). Apply at hartwigmedicalfoundation.nl/data/aanvragen-data | Hartwig DAC. Organoids via HUB Organoid Biobank | WGS, PDO 5-FU/oxaliplatin/irinotecan screens, RECIST, PFS/OS | ~42 in interim; 232 patients enrolled | Best modern paired set (AUROC up to ~0.88) | V (from data statement snippet) |
| 12 | Wensink 2024 J Exp Clin Cancer Res (OPTIC pilot) | Not confirmed | Corresponding author (UMCU) | Drug screens vs response | 6 / 10 / 11 (5-FU / irinotecan / oxaliplatin) | Screening protocol for OPTIC | U |
| 13 | Ooft 2021 ESMO Open (SENSOR) | Not confirmed. Likely NKI IRB, as #1 | NKI | PDO screens, off-label therapy response | 6 treated | A negative-result trial; useful for calibration | U |
| 14 | Narasimhan 2020 Clin Cancer Res (CRC peritoneal) | Supplement only | Corresponding author (Wake Forest/Ohio State) | PDO screens, treatment course | ~few | Peritoneal-metastasis context | U |

## 2. EGA process (~30 min of your time per DAC once the documents exist)

1. **EGA account.** Register at ega-archive.org (free). Use your institutional email. The account name must match the applicant named on the DAA.
2. **Find the DAC.** On each study or dataset page, note the EGAC id and the "Contact" email. Some DACs (UMCU, HKU) want an email first. Others take requests through the portal.
3. **Request access.** Log in, open the dataset (EGAD...) and click **"Request access"**. Paste the research-use statement (see `request_letter_template.md`). One request per EGAD. Bundle the UMCU datasets (#3-#6) into a single email to DACDBG@umcutrecht.nl.
4. **DAA / MTA.** The DAC sends its Data Access Agreement. It must be signed by an **institutional signatory** (legal or research office), not by the PI alone. UMCU also requires an MTA and approval by the UMCU medical ethics committee. Attach your local IRB/ethics number and the security checklist.
5. **Timelines (typical, unverified).** First reply 1-4 weeks. DAA negotiation 2-8 weeks. Total about 1-3 months. HKU and ICR are often slower.
6. **Download.** Once approved, the datasets appear under your EGA account.
   ```bash
   pip install pyega3
   # credentials.json: {"username": "...", "password": "..."}
   pyega3 -cf credentials.json datasets
   pyega3 -cf credentials.json files EGAD00001004387
   pyega3 -cf credentials.json -c 4 fetch EGAD00001004387 --output-dir /secure/ega/
   ```
   Files arrive decrypted over TLS. Store them only on encrypted, access-controlled storage (see `data_security_checklist.md`). Some newer datasets use Crypt4GH; keep any private key offline.

## 3. dbGaP process (phs001611)

1. The PI needs an **eRA Commons** account with the "PI" role. Non-US institutions need an institutional eRA/NIH registration (allow 2-4 weeks).
2. The institution's **Signing Official (SO)** must be registered in eRA Commons. An **IT Director** must be designated.
3. Submit a Data Access Request in dbGaP Authorized Access: research-use statement, non-technical summary, collaborators, and selection of phs001611.
4. The SO approves, then the NIH DAC reviews (usually 2-6 weeks). Approval lasts 1 year and is renewable. Annual progress reports are required.
5. Download with SRA Toolkit / fasterq-dump using the `.ngc` repository key. Follow the NIH Genomic Data Sharing security best practices (cloud use is allowed if the DAR declares it).
6. **First, use the open expression matrix (AWS Open Data) now.** Apply only if you need raw reads or variants.

## 4. Priority order

1. **UMCU bundle (EGAC00001000432: Kopper + de Witte + Driehuis + Sachs).** One DAC and one email cover 4 datasets with organoid DNA/RNA and the ovarian/HNSCC outcome pairs. The contact is confirmed.
2. **Hartwig / OPTIC (Smabers 2025).** This is the largest modern mCRC cohort with prospective PDO-vs-patient response and WGS. It is the closest match to the irinotecan/5-FU/oxaliplatin goal. Ask the corresponding author for the per-patient PDO AUC + RECIST table alongside it.
3. **NKI TUMOROID (Ooft 2019) via the NKI IRB**, plus SENSOR in the same request. These are paired chemo-response labels from the same group that set the field standard.

Do in parallel at zero cost: Tiriac open expression (AWS), Yao Mendeley curves, and GSE171680/1 (verify the source first). Yan (HKU) and Vlachogiannis (ICR) are useful but second-wave: they offer pretraining data or small n.
