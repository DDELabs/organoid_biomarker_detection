# Data landscape for drug-response pathway biomarkers (cancer type x drug)

Compiled 2026-10-11. This is a desk survey. Web search and NCBI E-utilities were used to check
the items. GEO sample counts in this file come from a live `esummary` call; they count GEO
samples (arrays/libraries), which is not always the number of patients.
Anything marked **(unverified)** could not be confirmed against a primary page during this survey.
Check it before relying on it.

## 0. Bottlenecks and what each one needs

| ID | Bottleneck | What fixes it |
|---|---|---|
| **B1** | TCGA has few, noisy RECIST labels and short follow-up. Using the labelled patients for discovery leaves too few for survival validation | (a) Discovery without patient labels (perturbation/mechanism data). (b) Separate large cohorts with long follow-up for validation |
| **B2** | Few public treated cohorts per drug that have pre-treatment expression plus survival | Aggregators (CTR-DB, ROCplotter), trial deposits, real-world clinico-genomic cohorts (Hartwig, POG) |
| **B3** | Organoid/PDX signatures do not transfer to patients | Patient-anchored data for defining pathways. Perturbation signatures across many contexts. Single-cell data to locate signal in the tumour compartment versus the TME |
| **B4** | Prognostic cannot be separated from predictive without untreated/control arms | Randomised trials with a control arm (treatment x biomarker interaction). As a second-best, observational cohorts with a treated/untreated flag analysed with propensity weighting |

Access codes: **O** = open download; **R** = registered/click-through (days); **C** = controlled application with a DAC/DUA (weeks to months); **M** = member-only or commercial.

---

## 1. Real-world clinico-genomic cohorts (tumour molecular data plus treatment lines plus outcome)

| Resource | Molecular | Treatment / endpoint | Cancers, n | Access, timeline | Fixes |
|---|---|---|---|---|---|
| **Hartwig Medical Foundation (HMF)** | WGS (tumour/normal) for all samples. RNA-seq for a subset, delivered as FASTQ | Pre-biopsy treatment history, post-biopsy treatment, RECIST response (first/best) **(field names unverified)**. OS/death dates not confirmed **(unverified)** | Metastatic pan-cancer: about 7,000+ patients (news reports). One user reports RNA for 394 mCRC, 332 mBC and 127 lung samples. A VHIO/HMF study built 56 treatment cohorts. TRANSACT used 226 RNA+response patients | **C**. Email dataaccess@hartwigmedicalfoundation.nl; online form; scientific review, then Data Access Board; institutional licence; delivered through GCP buckets (deleted after 1 month unless extended). Allow 1-3 months **(unverified)** | B2 (metastatic, drug-specific cohorts for platinum, taxanes, gemcitabine, fluoropyrimidines, ICB, endocrine). Response is measured after the biopsy, so expression is truly pre-treatment. No control arm (B4 not fixed) |
| **POG570 (BC Cancer Personalized OncoGenomics)** | WGS+RNA-seq (TPM) of advanced tumours | Treatment history with **treatment durations** (time on treatment as a response proxy) in Table S2 *Treatment* | 570 patients, more than 20 tumour types (Pleasance et al., *Nat Cancer* 1:452, 2020) | **O** for processed data: expression TPMs, small-mutation catalogue and README listing `Table_S2_Treatment.xlsx` at bcgsc.ca/downloads/POG570/ (bcgsc.ca did not resolve from this sandbox; listing confirmed through a search snippet). Raw data **C** at EGA EGAS00001001159 | B2. Time-on-treatment per drug class in metastatic disease, usable within days. Per-drug n is small (tens) **(unverified)** |
| **MSK-CHORD** (Jee et al., *Nature* 2024) | MSK-IMPACT DNA panel only. **No RNA** | NLP-derived progression, structured medications/timelines, OS | 24,950 patients: NSCLC 7,809, breast 5,368, CRC 5,543, prostate 3,211, pancreas 3,109 | **O** on cBioPortal (`msk_chord_2024`; study ID from a third-party R package, check it) | B2/B4 only for **mutation-defined** pathway activity (e.g. pathway-alteration scores). Very large n, so treated-versus-other-line contrasts are possible. Not an expression biomarker source |
| **AACR GENIE BPC** | DNA panel (multi-institution). No RNA found | Curated regimen-level real-world PFS, OS, line of therapy, imaging/pathology-based progression | NSCLC 1,846 (v2.0-public); CRC v2.0-public; Breast 1.0-public (2026); pancreas, prostate (2026); bladder, RCC expected 2026 | **R** through Synapse with a registered account and data-use terms (days). `genieBPC` R package | As MSK-CHORD: regimen-specific rwPFS for mutation-pathway biomarkers. Second-line versus first-line regimens give within-patient contrasts |
| **ORIEN Avatar** (M2GEN/Aster Insights) | WES tumour/normal plus RNA-seq | Lifetime clinical data (about 325 elements), treatment lines, outcome | More than 37,000 tumours (latest abstract); 18 US member centres | **M**. Member institutions only (project request, DUA, IRB; several months). Otherwise commercial licence | Ideal for B1/B2 at scale, but **not feasible within weeks** unless the team is at a member centre |
| **ICGC ARGO** | WGS/WXS, RNA-seq (STAR/HISAT2 pipelines) | Clinical model includes treatment and follow-up. Response fields sparse **(unverified)** | Release 14 (Mar 2026): 7,056 donors | **C** through DACO (grants 2 years). Open-tier clinical data browsable | B2 (moderate). Treatment/follow-up completeness varies by programme. Note POG is an ARGO programme (BC-POGP) |
| **TCGA-CDR** (Liu et al., *Cell* 2018) + GDC drug table | RNA-seq, all TCGA | OS/DSS/DFI/**PFI** (PFI recommended); drug table with sparse RECIST | About 11,000 patients, 33 cancers | **O** | B1 partly: use PFI, not OS, and restrict RECIST labels. Best used as **unlabelled** expression for unsupervised pathway structure, not for discovery labels |
| **CPTAC** (PDC/GDC) | Proteome, phosphoproteome, RNA, WGS/WXS | Mostly treatment-naive primaries. Follow-up and treatment forms sparse **(unverified)** | About 1,000+ tumours, about 10 cancers | **O** (processed), **C** (germline) | B3 (are pathways active at the protein level in patients?). Weak on outcome. Not a response cohort |
| **CGGA** (glioma) | mRNA-seq (cohorts of about 693 and 325) plus arrays | OS, **radiotherapy and TMZ chemotherapy flags** **(fields from memory; unverified)** | About 1,000+ gliomas | **O** | B4 (non-randomised): TMZ-treated versus untreated interaction for glioma pathway biomarkers |
| **METABRIC** | Illumina arrays, CNA | Chemo/hormone/radio flags; OS/DSS with **about 20-year follow-up** | About 1,980 breast | **O** (cBioPortal); raw at EGA (**C**) | B1 (long follow-up) and B4 (treated versus untreated, observational) |
| **SCAN-B** (GSE96058) | RNA-seq | Treatment group (endocrine, chemo, none), OS | About 3,273 breast **(count unverified here)** | **O** | B1/B4 as METABRIC, with RNA-seq on a modern platform |
| **Commercial: Tempus, Caris (WTS plus claims-derived time-on-treatment), Flatiron-FMI CGDB (DNA only)** | RNA (Tempus, Caris) | rwOS, time on treatment, time to next treatment | Tens of thousands | **M** (contract, months, fee) | Best scale for B2, but not feasible within weeks |

### Observational cohorts with a treated/untreated flag (proxy for B4, all **O**)
- **GSE39582** (CRC, 585 samples): adjuvant chemo yes/no, RFS/OS.
- **GSE14333** (CRC, 290): adjuvant chemo flag, DFS.
- **GSE68465** (lung adenocarcinoma, 462): adjuvant chemo/radiation flags, OS.
- **GSE42127** (NSCLC, adjuvant chemo versus none; already in the repo plan).
- **GSE6532** (ER+ breast, 741 arrays): tamoxifen-treated versus systemically untreated.
- **GSE10846** (DLBCL, 420): CHOP versus R-CHOP (historical, not randomised), so a rituximab interaction can be tested.
- METABRIC, SCAN-B, CGGA (above).

Caveat: these have confounding by indication. Use inverse-probability weighting on stage, grade and age, and report the **interaction** term, not the main effect.

---

## 2. Aggregated response-transcriptome resources

| Resource | Content | Endpoints | n | Access | Fixes |
|---|---|---|---|---|---|
| **CTR-DB 2.0** (Liu et al., *NAR* 2025, doi:10.1093/nar/gkae993) | Curated pre-treatment patient transcriptomes with response labels per regimen, from GEO/ArrayExpress and other sources | Mostly response (RECIST/pCR); survival for a subset | 10,856 patient samples, 39 cancers, 346 regimens | **O** (ctrdb.ncpsb.org.cn) | **B2** (main source of extra per-drug cohorts). B4 only where a source dataset has a control arm |
| **ROC Plotter** (Fekete & Gyorffy, *Int J Cancer* 2019) | Expression versus therapy response | pCR/RFS-at-5-years response classes | Breast 3,104 (chemo 2,108, endocrine 971, anti-HER2 267); ovarian 2,369; later GBM/CRC modules and an ICB module **(sizes for these unverified)** | **O** web, gene/signature input; no bulk download | B2: fast external validation of a *signature* (mean expression) with AUC. Cannot model interactions |
| **KM-plotter** (kmplot.com) | Expression versus survival | OS/RFS/PFS; breast/ovarian/gastric/lung subsets filterable by treatment (e.g. systemically untreated, endocrine-only, chemo) **(treatment filters from memory; unverified)** | More than 35,000 samples, 21 tumour types; 2024 colon module of 17 GEO cohorts (2,137) | **O** web | B1/B4 (untreated versus treated subgroup survival check of a signature in minutes) |
| **"BEST" (Biomarker Exploration of Solid Tumors)** | Curated GEO/TCGA plus ICB cohorts | Survival, ICB response | Tens of thousands **(unverified; no Gyorffy "BEST" paper found; likely a separate Chinese database)** | **O** web **(unverified)** | B2 screening only |
| **CRI iAtlas** | 15 ICB trials uniformly processed, plus TCGA/PCAWG | ICB response, OS/PFS | 1,142 ICB samples (12,677 total) | **O** (portal plus Synapse) | B2 for ICB drugs |
| **IMvigor210CoreBiologies** (R package) | Atezolizumab mUC RNA-seq | OS, RECIST | about 348 | **O** | B2 (ICB, single arm) |
| **PharmacoDB / PharmacoGx** | GDSC, CTRP, CCLE, gCSI, PRISM, NCI60, PDX dose-response | Cell-line/PDX AUC/IC50 | About 1,700 lines, more than 700k experiments **(unverified)** | **O** | B3: model-side only. Use to test whether a pathway's drug association is consistent across datasets |
| **cBioPortal timelines** | Treatment/progression timelines for MSK-CHORD, GENIE BPC and some ICB studies | Line/regimen, progression | See above | **O/R** | Exposes treatment timelines in the same format, so harmonisation is cheap |

---

## 3. Randomised trials with public expression **and** a control or comparator arm (key for B4)

Endpoints: pCR = pathological complete response; RFS/PFS/OS as usual. "n" is GEO samples (verified via E-utilities) unless noted.

| Trial | Cancer | Arms | Accession, n | Endpoint | Access | Notes |
|---|---|---|---|---|---|---|
| **JBR.10** | Resected NSCLC | Cisplatin/vinorelbine **vs observation** | GSE14814, 133 | OS, DSS | O | Classic chemo-versus-control interaction test |
| **STORM** (adjuvant) | HCC | Sorafenib **vs placebo** | GSE109211, 140 | RFS | O | Already listed in ENVIRONMENT_SETUP |
| **ICON7** (translational subset) | Ovarian HGSC | Carboplatin/paclitaxel **± bevacizumab** | GSE140082, 380 (ICON7-A pooled 745 with individual data) | PFS, OS | O | Bevacizumab x pathway (angiogenesis/mesenchymal) interaction |
| **APEX (039) + 025/040** | Relapsed myeloma | Bortezomib **vs dexamethasone** | GSE9782, 528 arrays (about 264 patients **(unverified)**) | Response, TTP, OS | O | Proteasome/UPR pathway interaction |
| **FOCUS** (MRC) | mCRC | FU **± oxaliplatin**/irinotecan strategy arms | GSE156915, 361 | PFS, OS | O | Oxaliplatin x DDIR pathway |
| **BrighTNess** | TNBC neoadjuvant | Paclitaxel **± carboplatin ± veliparib** | GSE164458, 482 | pCR (EFS reported separately **(unverified in GEO)**) | O | Platinum x DNA-repair pathway |
| **I-SPY2 (990 resource)** | Breast neoadjuvant | 10 arms with a shared paclitaxel **control** (210): veliparib/carbo 71, neratinib 114, MK2206 93, ganitumab 106, ganetespib 93, trebananib 134, T-DM1/pertuzumab 52, H/P 44, pembrolizumab 69 | GSE194040 (reported subseries; E-utilities lookup failed, check it), 987 patients | pCR (EFS published for 950, not deposited per arm **(unverified)**) | O | Best multi-drug control-arm resource; Agilent 44K plus RPPA (736) |
| **CALGB 40601** | HER2+ breast | Paclitaxel plus trastuzumab vs lapatinib vs both | GSE116335, 265 | pCR (RFS in follow-up papers) | O | Active-comparator, not untreated |
| **CHER-LOB** | HER2+ breast | Chemo plus trastuzumab vs lapatinib vs both | GSE66399, 156 | pCR | O | Active comparator |
| **GSE20271** | Breast neoadjuvant | FAC vs **TFAC** (± paclitaxel) | GSE20271, 178 | pCR | O | Paclitaxel x pathway |
| **GSE41998** | Breast neoadjuvant | AC then ixabepilone vs paclitaxel | GSE41998, 279 | pCR | O | Active comparator |
| **CheckMate 009/010/025** (Braun et al., *Nat Med* 2020) | ccRCC | Nivolumab **vs everolimus** (CM025) | Supplementary normalised RNA matrix (Table S4); about 170-311 patients **(unverified)** | PFS, OS, ORR | O | ICB versus active control |
| **IMmotion151** | ccRCC | Atezolizumab+bevacizumab **vs sunitinib** | EGAS00001004353 (about 823 RNA overall; 407 atezo+bev in one reuse) | PFS, OS | C (Roche DAC EGAC00001001748) | Randomised ICB+VEGF versus VEGF-TKI |
| **IMmotion150 / POPLAR / IMvigor210 / PCD4989g** | RCC / NSCLC / mUC | IMmotion150: atezo vs atezo+bev vs sunitinib; POPLAR: atezo **vs docetaxel** | EGAS00001004343 (FASTQ EGAD00001006631, TPM EGAD00001006741), 817 tumours | PFS/OS | C (same DAC) | OAK (atezo vs docetaxel) deposit not found **(unverified)** |
| **PETACC-3** | Stage II/III colon | FU/LV **± irinotecan** | ArrayExpress E-MTAB-990 **(accession and access unverified)**, about 688 | RFS, OS | ? | Worth checking |
| **S:CORT (FOCUS, COIN, New EPOC etc.)** | CRC | Cetuximab and oxaliplatin randomisations | About 1,300-2,000 patients, multi-omic | C/M (collaboration with Prof. T. Maughan, or commercial licence through Cancer Research Horizons) | Largest CRC randomised multi-omic set; slow |
| **CLASSIC** (gastric) | Gastric adjuvant | XELOX **vs observation** | NanoString classifier data; public deposit **(unverified)** | DFS | ? | Ask the authors |

Non-randomised but useful drug-specific cohorts (O): GSE25066 (taxane-anthracycline breast, 508, DRFS), GSE72970 (mCRC FOLFOX/FOLFIRI ± targeted, 124), GSE104645 (CRC CMS x chemo, 193), GSE5851 (cetuximab mCRC, 80), GSE91061/GSE78220 (ICB melanoma), GSE176307 (ICB mUC, 90), GSE68871 (VTD myeloma, 118).

**How to use for B4:** fit `outcome ~ arm * pathway_score + covariates` (Cox or logistic). A biomarker is *predictive* only if the interaction is significant and in the expected direction. A main effect seen equally in both arms means it is *prognostic*. Pre-register the score before unblinding these cohorts.

---

## 4. Perturbation and mechanism resources (no outcome labels needed, so they address B1 and B3)

| Resource | Content | Size | Access | Use |
|---|---|---|---|---|
| **LINCS L1000 / CMap (CLUE)** | 978-landmark transcriptional responses to compounds, shRNA/CRISPR and overexpression | About 3M signatures, about 30k perturbagens **(2020 release sizes unverified)**; GSE92742 and GSE70138 | O | Define a **drug-induced pathway response** per drug in many lines. Score patients by reversal or concordance. Discovery without patient labels, which saves labelled patients for validation |
| **LINCS P100 / GCP** | Targeted phospho- and chromatin proteomics after drug treatment | Hundreds of drugs x few lines **(unverified)** | O | Protein-level mechanism for pathway definition |
| **DepMap CRISPR (Chronos) + PRISM Repurposing** | Gene dependency; multiplexed drug viability | About 1,100+ lines (CRISPR); PRISM about 4,500 compounds x about 578 lines (primary) **(unverified exact)** | O (CC-BY) | Link the drug-sensitivity pathway to target dependency (orthogonal check of mechanism). Biomarker = pathway that predicts both PRISM AUC and CRISPR dependency on the drug target |
| **GDSC / CTRP / gCSI** (via PharmacoDB) | Dose-response plus baseline omics | About 1,000 lines x hundreds of drugs | O | Cross-dataset consistency filter for B3 |
| **Tahoe-100M** (Arc Virtual Cell Atlas) | scRNA-seq after drug treatment | More than 100M cells, 50 cancer lines, about 1,100 small-molecule perturbations (3 doses) | O (CC0; Hugging Face/Arc) | Cell-intrinsic pathway response per drug at single-cell resolution; separates responder and non-responder subpopulations |
| **sci-Plex** (Srivatsan 2020) | scRNA-seq chemical screens | 188 compounds x 3 lines | O | Smaller precursor to Tahoe |
| **DrugComb** | Combination screens | About 700k+ combinations **(unverified)** | O | Only if combinations are in scope |
| **decryptM** (Zecha et al., *Science* 2023) and drug-proteome atlases (e.g. Mitchell et al., *Nat Biotech* 2023) **(details unverified)** | Dose-resolved (phospho)proteome responses to drugs | About 30-900 drugs | O (PRIDE) | Protein-level pathway engagement. This is the realistic reading of "CPTAC drug-proteome". CPTAC itself profiles untreated tumours |

Recommended use: build **drug-response pathway signatures from perturbation data** (L1000/Tahoe/PRISM), freeze them, then test them in *all* labelled patient cohorts as pure validation. This removes the discovery/validation split that causes B1.

---

## 5. Single-cell / spatial atlases with treatment and response (pre/post)

| Dataset | Cancer, drug | Design | Accession | n | Access |
|---|---|---|---|---|---|
| Sade-Feldman 2018 | Melanoma, anti-PD-1/CTLA-4 | Pre/post, responder versus non-responder | GSE120575 | 48 samples, 32 patients | O |
| Jerby-Arnon 2018 | Melanoma, ICB resistance | Untreated versus post-ICB | GSE115978 | 7,186 cells (GEO samples); 33 tumours **(unverified)** | O |
| Yost 2019 | BCC, anti-PD-1 | Pre/post | GSE123813 | 86 GEO samples | O |
| Zhang 2021 (*Cancer Cell*) | TNBC, paclitaxel ± atezolizumab | Pre/on/post | GSE169246 | 83 GEO samples, 22 patients **(unverified)** | O |
| Bassez 2021 (*Nat Med*) | Breast, anti-PD-1 window | Pre/on | EGA plus lab website **(unverified)** | about 40 patients | O/C |
| Maynard 2020 (*Cell*) | NSCLC, targeted TKIs | Pre/residual/progression | PRJNA591860 **(unverified)** | 49 biopsies | O |
| NSCLC neoadjuvant PD-1 ± chemo | NSCLC | Pre/post | GSA-Human HRA001033 (15 patients), HRA006493 (more than 1M cells) | | C |
| **DRMref** (*NAR* 2024) | 13 cancers, 42 drugs, 35 regimens | Resistant versus sensitive; 14 paired pre/post | ccsm.uth.edu/DRMref | 42 datasets, 382 to 666 samples | O |
| **scDrugAtlas** (*Database* 2026) | 31 cancers, 77 drugs | Drug-response sc compendium | PMC12923164 | 1,023 samples | O |
| **scICB** | BRCA, HNSCC, CRC, SKCM, NSCLC; ICB | 247 patients with matched pre/post | (journal resource) | 338 patients | O **(unverified)** |
| TISCH2 | Pan-cancer, includes ICB datasets | Curated scRNA | tisch.comp-genomics.org | about 190 datasets **(unverified)** | O |

Value: B3. Tells whether an organoid pathway signal sits in tumour cells, which organoids model, or in stroma/immune cells, which they do not. Deconvolution of bulk trial data with these references makes compartment-specific pathway scores possible. These datasets have no survival power (n is too small).

---

## 6. Top-10 ranking (value x feasibility within weeks)

Scores from 1 to 5. Feasibility: O is about 5, R about 4, C about 2, M about 1.

| Rank | Resource | Value | Feas. | Score | Bottlenecks fixed | First action |
|---|---|---|---|---|---|---|
| 1 | **Randomised control-arm GEO trials** (GSE14814, GSE109211, GSE140082, GSE9782, GSE156915, GSE164458, GSE194040, GSE20271) | 5 | 5 | 25 | **B4**, B2 | Download, harmonise arm/endpoint, test arm x pathway interaction |
| 2 | **Perturbation signatures: LINCS L1000 + PRISM/DepMap + Tahoe-100M** | 5 | 5 | 25 | **B1**, B3 | Build frozen drug-response pathway scores without patient labels |
| 3 | **CTR-DB 2.0** | 4 | 5 | 20 | **B2** | Pull every cohort for target drugs; use it only for validation |
| 4 | **Observational treated/untreated cohorts** (SCAN-B, METABRIC, GSE39582, GSE68465, GSE6532, CGGA) | 4 | 5 | 20 | B1 (long follow-up), B4 (approximate) | IPW-adjusted interaction models |
| 5 | **POG570 open TPM + treatment-duration table** | 4 | 4 | 16 | B2 (metastatic, time on treatment) | Download from bcgsc.ca (blocked from this sandbox, so fetch locally) |
| 6 | **ICB harmonised sets** (CRI iAtlas, IMvigor210CoreBiologies, Braun CM025 with everolimus control) | 3 | 5 | 15 | B2, B4 (CM025) | Use only if ICB is in the drug list |
| 7 | **ROC Plotter / KM-plotter** | 3 | 5 | 15 | B2, B4 screen | Sanity screen of each final signature (hours) |
| 8 | **GENIE BPC / MSK-CHORD** | 3 | 4 | 12 | B2/B4 for mutation-pathway biomarkers only (DNA panel, no RNA) | Use if the biomarker can be expressed as pathway-alteration status |
| 9 | **Hartwig Medical Foundation** | 5 | 2 | 10 | B2 (largest metastatic RNA+treatment+response) | **Submit the DAR now** (1-3 months), in parallel with ranks 1-7 |
| 10 | **Roche EGA trials (IMmotion150/151, POPLAR)** | 4 | 2 | 8 | B4 (randomised, sunitinib/docetaxel control) | Apply to EGAC00001001748 now |

Not in the top 10: ORIEN and the commercial datasets (best scale, but member-only/contract; months). ICGC ARGO (DACO; modest treatment data). CPTAC (no outcomes). Single-cell atlases (mechanism and B3 only; too small for survival).

## 7. Suggested sequencing (weeks 1-6)
1. Week 1: submit the Hartwig DAR and the Roche EGA DAC request (long lead times). Download the ranks 1, 3, 4 and 5 open data.
2. Weeks 1-2: build drug-response pathway signatures from L1000/PRISM/Tahoe (no patient labels) and freeze them.
3. Weeks 2-4: validate in CTR-DB, POG570 and observational cohorts (survival/time on treatment).
4. Weeks 3-5: test predictiveness (arm x score interaction) in the randomised control-arm trials. Treat this as the primary claim.
5. Weeks 5+: add Hartwig/EGA data on arrival, as a held-out confirmation set.

## Sources consulted
- POG570 README (bcgsc.ca/downloads/POG570/POG570_README.txt) and EGAS00001001159; Pleasance et al., *Nat Cancer* 2020.
- Hartwig data access pages (hartwigmedicalfoundation.nl/en/data/data-access-request/; hartwigmedical.github.io documentation); CIBRA (arXiv 2403.03829); VHIO/HMF non-response study (ESMO RWD&DO; bioRxiv 10.1101/2025.10.06.680777).
- Jee et al., *Nature* 2024 (doi:10.1038/s41586-024-08167-5); chordR GitHub.
- GENIE BPC: AACR analytic data guides; mskcc.org GENIE BPC page; Choudhury et al., *Clin Cancer Res* 2023 (PMC10472103).
- ORIEN: UVA Cancer Center page; ORIEN/M2GEN abstracts.
- ICGC ARGO release notes (docs.icgc-argo.org), DACO (daco.icgc-argo.org).
- ROC Plotter, *Int J Cancer* 2019 (PMID 31020993); KM-plotter (Database Commons 9164).
- CRI iAtlas (ISB news; JITC 2023 abstract).
- I-SPY2 990 resource (PMC9426306); ICON7-A (medRxiv 10.1101/2023.02.14.23285833).
- Roche EGA: EGAC00001001748, EGAS00001004343, EGAS00001004353.
- S:CORT (cancer.ox.ac.uk; Cancer Research Horizons).
- Tahoe-100M (bioRxiv 10.1101/2025.02.20.639398; Arc Institute).
- DRMref (*NAR* 2024, doi:10.1093/nar/gkad1087); scDrugAtlas (PMC12923164).
- NCBI GEO E-utilities (esummary) for all GSE sample counts listed.
