# BRIDGE: research survey and proposed next strategy

**BRIDGE** = **B**ayesian **R**esponse **I**ntegration across **D**omains with **G**ene-set **E**xplanations.

Scope: methods published 2020-2026 (plus a few earlier anchors) that either
(i) train drug-response predictors on patient-derived models (cell lines, PDX, PDO) and transfer
them to patients, or (ii) learn directly from patient tumour transcriptomes with response labels.
Emphasis is on methods that work with small n and can be read at the pathway level.

Status: literature survey plus a design proposal. Nothing here has been run. Performance numbers
are as reported by the original authors. Where a number or accession could not be checked against
the primary source during this survey, it is marked **(unverified)**.

Repo context (from `robust/ATLAS.md`, `robust/NEXTGEN.md`, `robust/TRIAD.md`):
- Organoid-first (Kong 2020 / NIT) failed on the patient side: proliferation confounding,
  immortal-time bias, and signals that were prognostic rather than predictive.
- Patient-first ATLAS: about 2,769 TCGA patient x drug RECIST labels, 673 Reactome rank scores,
  held-out (leave-one-cancer-out) AUROC 0.558-0.561. External validation is reproducible only in
  breast neoadjuvant taxane/anthracycline. No per-drug within-cancer pathway passes FDR.
- Already curated: GDSC, CTRP v2, PRISM; PDOs for bladder (Lee 2018), liver (Broutier 2017,
  Ji 2023), pancreas (Tiriac 2018, Shi 2022) and sarcoma (Al Shihabi 2024); LICOB; 21 GEO trial
  cohorts plus the ENLIGHT/trial cohorts in `data/external/trials/`.

---

## 1. Survey of methods

### 1.1 Transfer from pre-clinical models to patients (domain adaptation)

| Method | Data needed | Reported patient-level performance | Code | Pathway biomarkers? |
|---|---|---|---|---|
| **PRECISE**, Mourragui et al. 2019, *Bioinformatics* 35:i510, doi:10.1093/bioinformatics/btz372 | Cell-line and/or PDX expression + response (source); unlabelled tumour expression (target) | Recovered known biomarkers in TCGA (e.g. ERBB2/lapatinib-type associations). Only modest separation of RECIST responders. | github.com/NKI-CCB/PRECISE | **Yes, indirectly.** Linear consensus factors (principal vectors between PCA subspaces) have gene loadings, so factor coefficients map back to genes or pathways. |
| **TRANSACT**, Mourragui et al. 2021, *PNAS* 118:e2106682118, doi:10.1073/pnas.2106682118 | As PRECISE; kernel PCA consensus space. The non-linearity was tuned on PDXE. | 23 drug challenges in TCGA plus 226 metastatic tumours (Hartwig). Significant gains over 4 competitors for 5 TCGA drugs and 3 HMF drugs (platinum, gemcitabine, paclitaxel). No patient labels used in training. Per-drug AUCs mostly 0.6-0.7 **(unverified)**. | github.com/NKI-CCB/TRANSACT | **Partial.** Consensus factors can be decomposed into linear and non-linear gene contributions; the paper shows this. |
| **CODE-AE**, He et al. 2022, *Nat Mach Intell* 4:879, doi:10.1038/s42256-022-00541-0 | Unlabelled CCLE + TCGA expression for pre-training; GDSC labels for fine-tuning | Beat AD-AE, TCRP, VAEN, DSN and Deep CORAL on TCGA chemotherapy response; screened 50 drugs over 9,808 patients. Per-drug AUROC: see the paper's supplement **(unverified)**. | github.com/XieResearchGroup/CODE-AE | **No.** Post-hoc attribution only. |
| **Velodrome**, Sharifi-Noghabi et al. 2021, *Nat Mach Intell* 3:962, doi:10.1038/s42256-021-00408-w | Labelled cell lines from several sources + unlabelled samples; target not seen during training (domain generalisation) | Evaluated on PDX and patient cohorts for a few drugs; modest gains **(unverified detail)**. | Ester lab GitHub (hosseinshn) **(unverified URL)** | No |
| **AITL**, Sharifi-Noghabi et al. 2020, *Bioinformatics* 36:i380, doi:10.1093/bioinformatics/btaa442 | GDSC + small labelled PDX/TCGA target | Adversarial input + output discrepancy; target AUROCs about 0.55-0.75 for docetaxel, paclitaxel, cisplatin and bortezomib **(unverified)** | github.com/hosseinshn/AITL | No |
| **TCRP** (few-shot meta-learning), Ma et al. 2021, *Nat Cancer* 2:233, doi:10.1038/s43018-020-00169-2 | Cell lines (meta-training); 5-10 labelled target samples | Rapid gains with ≤10 samples when moving to new tissues, PDTCs and PDX. Patient-level evaluation limited. | github.com/idekerlab/TCRP | No |
| **TUGDA**, Peres da Silva et al. 2021, *Bioinformatics* 37:i76, doi:10.1093/bioinformatics/btab299 | GDSC (source); unlabelled PDX/TCGA; multi-task over drugs with uncertainty weighting | Less negative transfer than single-task models on PDXE and TCGA | github.com/CSB5/TUGDA | No (multi-task design is relevant) |
| **WISER**, Shubham et al. 2024, *ICML*, PMLR 235 (arXiv:2405.04078) | Cell lines + unlabelled patients; weak supervision on patients | Better than CODE-AE and other baselines on TCGA patient response | github.com/kyrs/WISER | No |
| **TransDRP**, Liu et al. 2025, *AAAI* 39 | Cell lines + patients; drug-knowledge GNN; within- and cross-cancer alignment | Better than unsupervised domain adaptation baselines on TCGA | (check paper) | No |
| **PDXGEM**, Kim et al. 2020, *BMC Bioinformatics* 21:288, doi:10.1186/s12859-020-03633-z | PDXE expression + response; patient expression for a **co-expression concordance filter** | Separated responders in independent patient cohorts (paclitaxel breast, trastuzumab HER2+ pCR) and stratified OS for gemcitabine in pancreas | pdxgem.moffitt.org (R) | **Gene-level.** The concordance filter is directly portable to pathways. |
| **TG-LASSO**, Huang et al. 2020, *PLoS Comput Biol* 16:e1007607, doi:10.1371/journal.pcbi.1007607 | GDSC + TCGA tissue labels | Separated TCGA responders from non-responders for 7 of 13 drugs. Adding a network did **not** help, while tissue guidance did. | github (KnowEnG) | **Gene-level** sparse; maps to pathways |
| **PharmaFormer**, Zhou et al. 2025, *npj Precis Oncol*, doi:10.1038/s41698-025-01082-6 | GDSC pre-training, then PDO fine-tuning (colon, bladder); SMILES | Pan-cancer cell-line pre-training plus tumour-type organoid fine-tuning improved separation of responders in TCGA colon (5-FU/oxaliplatin), bladder (cisplatin/gemcitabine) and liver patients | github.com/zhouyuru1205/PharmaFormer; also in `drevalpy` | No |
| **Kong et al. 2020** (this repo's origin), *Nat Commun* 11:5485, doi:10.1038/s41467-020-19313-8 | PDO AUC (CRC, bladder), drug-target network, ssGSEA | 114 5-FU CRC and 77 cisplatin BLCA patients (survival split) | this repo | Yes (Reactome) |
| **Zhang et al. 2025**, *Transl Oncol* 52:102238, doi:10.1016/j.tranon.2024.102238 | Matched CRC tumour-organoid pairs + an independent organoid IC50 set; consensus WGCNA | Better than standard approaches on independent 5-FU datasets | (paper) | **Module level.** Consensus tumour/organoid modules are an alignment-by-design idea. |

**Alignment tools without a response model**

- **Celligner**: Warren et al. 2021, *Nat Commun* 12:22, doi:10.1038/s41467-020-20294-x.
  Contrastive PCA removes axes specific to one dataset (tumour stroma/immune, culture artefacts),
  then mutual-nearest-neighbour alignment. Most cell lines align to their own cancer type; several
  hundred mesenchymal or undifferentiated lines align poorly.
- **MOBER**: Dimitrieva et al. 2025, *Sci Adv*, doi:10.1126/sciadv.adn5596. Adversarial
  conditional VAE. Integrates 932 cell lines, 434 PDX and 11,159 tumours; can transform PDX/cell-line
  profiles to look like tumours. Code on Zenodo (14209839).

These tools are most useful as **fidelity filters**: drop pre-clinical models that do not
resemble any patient tumour of their labelled cancer.

**Optimal transport (OT)**

- Generic OT domain adaptation: Courty et al. 2017, *IEEE TPAMI* 39:1853,
  doi:10.1109/TPAMI.2016.2615921.
- CellOT, for single-cell perturbation maps applied to held-out patients: Bunne et al. 2023,
  *Nat Methods* 20:1759, doi:10.1038/s41592-023-01969-x.
- No peer-reviewed method was found that uses OT specifically for bulk model-to-patient
  drug-response transfer. This is a gap.
- Entropic OT (POT library: Flamary et al. 2021, *JMLR* 22(78)) is cheap and can be tested as an
  alignment variant inside BRIDGE.

### 1.2 Methods that learn on patients or encode prior knowledge (patient-side)

| Method | Data needed | Reported patient-level performance | Code | Pathway biomarkers? |
|---|---|---|---|---|
| **SELECT**, Lee et al. 2021, *Cell* 184:2487, doi:10.1016/j.cell.2021.03.030 | Tumour transcriptome + drug targets. Synthetic-lethal/rescue partners mined from TCGA/DepMap. **No response training.** | Predictive in about 80% of 35 targeted- and immunotherapy trials (10 cancer types) and in the WINTHER trial | Signatures and code for academic use (Ruppin lab) | **Gene-partner level**, interpretable; could be lifted to pathways |
| **ENLIGHT**, Dinstag et al. 2023, *Med* 4:15, doi:10.1016/j.medj.2022.11.001 | As SELECT, extended GI networks; unsupervised | 21 blinded cohorts, 697 patients, 15 drugs: OR = 2.59 for high match score (preprint figures). The repo already holds these cohorts as `enlight_*`. | Commercial (Pangea); cohorts public | Gene-partner level |
| **PERCEPTION**, Sinha et al. 2024, *Nat Cancer* 5:938, doi:10.1038/s43018-024-00756-7 | Cell-line bulk + scRNA (to learn); patient **scRNA** (to apply) | Validated in a multiple myeloma trial and a breast cancer trial, and for lung TKI resistance | github.com/ruppinlab/PERCEPTION | No (gene-level elastic net) |
| **CTR-DB 2.0**, Liu et al. 2025, *Nucleic Acids Res*, doi:10.1093/nar/gkae993 (v1: doi:10.1093/nar/gkab860) | Resource: 10,856 patient samples, 39 cancers, 346 regimens, pre-treatment transcriptome + response | Not a model; the best source of additional external cohorts | ctrdb.ncpsb.org.cn (free) | Has built-in GSEA/biomarker modules |
| **Ding et al. 2016**, *Bioinformatics* 32:2891, doi:10.1093/bioinformatics/btw344 | TCGA RECIST curation | Source of most "TCGA response" benchmarks (TRANSACT, CODE-AE, WISER). Labels are noisy. | supplement | n/a |
| **oncoPredict / pRRophetic**: Maeser et al. 2021, *Brief Bioinform* 22:bbab260, doi:10.1093/bib/bbab260; Geeleher et al. 2014, *Genome Biol* 15:R47 | Cell-line ridge, batch-corrected to patients | Classic baseline; weak but sometimes significant (e.g. docetaxel, bortezomib) | CRAN/GitHub | Gene-level |

### 1.3 Pathway- or hierarchy-structured predictors

- **DCell / DrugCell**: Ma et al. 2018, *Nat Methods* 15:290, doi:10.1038/nmeth.4627;
  Kuenzi et al. 2020, *Cancer Cell* 38:672, doi:10.1016/j.ccell.2020.09.014.
  - Visible neural network over the GO hierarchy, trained on cell lines (mutations + drug
    fingerprints).
  - Subsystem activations are interpretable.
  - Patient validation is limited; the combination predictions were validated in vitro and in PDX.
  - Code: github.com/idekerlab/DrugCell.
- **NeST-VNN**: Park et al. 2024, *Nat Cancer* 5:996, doi:10.1038/s43018-024-00740-1.
  - Visible neural network over NeST multiprotein assemblies, trained for palbociclib.
  - Eight assemblies predicted response in PDX and patients where single-gene markers failed.
- **PathDSP**: Tang & Gottlieb 2021, *Sci Rep* 11:3815, doi:10.1038/s41598-021-82612-7.
  Pathway-enrichment features (expression, mutation, CNV, drug-target pathway) in a deep network.
  Cell lines only. Code: github.com/TangYiChing/PathDSP.
- **Kernelised Bayesian matrix factorisation with pathway-response associations**:
  Ammad-ud-din et al. 2016, *Bioinformatics* 32:i455, doi:10.1093/bioinformatics/btw433.
  Pathway-by-drug association matrix learned jointly. A close ancestor of what BRIDGE needs.
- **Bayesian multitask multiple kernel learning (BMTMKL)**: Costello et al. 2014, *Nat Biotechnol*
  32:1202, doi:10.1038/nbt.2877. Won the DREAM challenge; pathway kernels shared across drugs.
- **scDEAL**: Chen et al. 2022, *Nat Commun* 13:6494, doi:10.1038/s41467-022-34277-7.
  Bulk-to-single-cell transfer at the cell level. Not a patient-response method, and not suitable
  for the bulk data here.

### 1.4 Does the patient-derived model predict the patient? (meta-analyses and key studies)

| Study | Design | Result |
|---|---|---|
| Wensink et al. 2021, *npj Precis Oncol* 5:30, doi:10.1038/s41698-021-00168-1 | Systematic review of 17 PDO studies | Pooled sensitivity 0.81 (0.69-0.89), specificity 0.74 (0.64-0.82) |
| Sakshaug et al. 2023, *Sci Rep* 13, doi:10.1038/s41598-023-45297-8 | CRC-only systematic review | PPV 68%, NPV 78% |
| Romero et al. 2025, *medRxiv*, doi:10.1101/2025.08.10.25333051 (PROSPERO CRD42023387229) | 411 patient-model pairs (267 PDX, 144 PDO) | 70% concordance; PDX = PDO; PDO response associated with patient PFS |
| Vlachogiannis et al. 2018, *Science* 359:920, doi:10.1126/science.aao2774 | GI metastases PDO | PPV 88%, NPV 100% (n = 21 patients) |
| Ooft et al. 2019, *Sci Transl Med* 11:eaay2574, doi:10.1126/scitranslmed.aay2574 | mCRC (TUMOROID) | Predicted irinotecan response; failed for 5-FU/oxaliplatin |
| Yao et al. 2020, *Cell Stem Cell* 26:17, doi:10.1016/j.stem.2019.10.010 | 80 LARC patients, chemoradiation | Accuracy about 84% **(unverified exact figure)** |
| Izumchenko et al. 2017, *Ann Oncol* 28:2595, doi:10.1093/annonc/mdx416 | PDX ("avatars"), 92 patients | Concordance about 87% **(unverified)** |

Implication for BRIDGE:
- Models agree with patients on direction about 70-80% of the time. That is good enough to act as
  an informative **prior**, but not as a substitute for patient labels.
- 5-FU/oxaliplatin in CRC PDO (Ooft 2019) is a documented failure case. The repo's own 5-FU
  difficulty is consistent with it.

### 1.5 Statistical design: partial pooling, borrowing, predictive vs prognostic

**Borrowing strength from pre-clinical data**
- Power priors: Ibrahim & Chen 2000, *Stat Sci* 15:46, doi:10.1214/ss/1009212673.
- Commensurate priors: Hobbs et al. 2011, *Biometrics* 67:1047, doi:10.1111/j.1541-0420.2011.01564.x.
- Both down-weight a source dataset by an amount learned from how well it agrees with the target.
  This is the principled version of the "GDSC prior strength τ" that TRIAD/ATLAS fix by hand.

**Partial pooling across cancers**
- EXNEX: Neuenschwander et al. 2016, *Pharm Stat* 15:123, doi:10.1002/pst.1730. From basket trials;
  pools only cancers that look exchangeable, which protects against forced pooling.
- Regularised horseshoe for sparse drug-specific effects: Piironen & Vehtari 2017, *Electron J Stat*
  11:5018, doi:10.1214/17-EJS1337SI.
- Local false sign rate (lfsr) for reporting effects: Stephens 2017, *Biostatistics* 18:275,
  doi:10.1093/biostatistics/kxw041.

**Predictive vs prognostic**
- Definitions and interaction test: Ballman 2015, *J Clin Oncol* 33:3968,
  doi:10.1200/JCO.2015.63.3651.
- Information-theoretic separation of prognostic and predictive markers: Sechidis et al. 2018,
  *Bioinformatics* 34:3365, doi:10.1093/bioinformatics/bty357.
- Treatment-effect heterogeneity:
  - virtual twins: Foster et al. 2011, *Stat Med* 30:2867, doi:10.1002/sim.4322;
  - causal forests: Wager & Athey 2018, *JASA* 113:1228, doi:10.1080/01621459.2017.1319839.
- Random-signature null: Venet et al. 2011, *PLoS Comput Biol* 7:e1002240,
  doi:10.1371/journal.pcbi.1002240. Already used in NIT.

**Benchmark lessons**
- Partin et al. 2023, *Front Med* 10:1086097, doi:10.3389/fmed.2023.1086097.
- The `drevalpy` benchmark library (it includes PharmaFormer).
- Lessons: naive baselines (tissue-only, mean-per-drug) are often competitive, and leakage through
  random splits inflates performance. BRIDGE must report a cancer-type-only baseline.

### 1.6 Take-aways that shape BRIDGE

1. **No method convincingly beats AUROC of about 0.65-0.70 on real patient chemotherapy labels.**
   ATLAS at 0.56 is low but not far off. Expect improvements of +0.02 to +0.06, not a leap.
2. **Only linear consensus or alignment methods keep pathway interpretability.** These are PRECISE,
   TRANSACT's linear part, the PDXGEM concordance filter and Celligner's cPCA. Deep models
   (CODE-AE, WISER, PharmaFormer) belong in BRIDGE only as *performance benchmarks*.
3. **Treat pre-clinical data as an informative prior with a learned discount**, not as the training
   set. The concordance meta-analyses (about 70-80%) and the repo's NIT failure both argue for this.
4. **Tissue matters more than networks** (TG-LASSO finding; TRIAD/ATLAS centring). Pool across
   cancers with partial pooling, not complete pooling.
5. **Predictive claims need randomised or comparator arms.** The repo already has some
   (BrighTNess, I-SPY2 arms, GSE20271, GSE41998, JBR.10/GSE14814, GSE42127).

---

## 2. Proposed strategy: BRIDGE

**One-line summary:** align pathway-score spaces across cell lines, PDX, PDO and patient tumours in
a *linear, cancer-aware consensus subspace*. Fit one *hierarchical Bayesian ordinal model* in which
pre-clinical coefficients act as commensurate priors for patient coefficients, with partial pooling
over drug, mechanism class and cancer. Read out per-(drug, cancer) pathway effects with lfsr, and
accept only effects that survive held-out-cancer, external-trial and treatment-interaction tests.

### 2.1 Inputs

**Domain D1: cell lines** (GDSC2, CTRP v2, PRISM; about 1,000 lines)
- Label: within-drug z-scored AUC, after the NIT deconfounding step (regress out proliferation
  and leave-drug-out general sensitivity).
- Role: prior only.
- Status: already curated.

**Domain D2: PDX** (Novartis PDXE: Gao et al. 2015, *Nat Med* 21:1318, doi:10.1038/nm.3954)
- Data:
  - RNA-seq (human-read FPKM) for about 400 models;
  - mRECIST category plus BestAvgResponse / time-to-doubling;
  - about 60 treatments in 6 indications (BRCA, CRC, CM, NSCLC, PDAC, GC);
  - 1x1x1 design.
- Label: ordinal mRECIST (CR/PR, SD, PD), plus continuous BestAvgResponse as a secondary label.
- Tooling: the Xeva R package (Mer et al. 2019, *Cancer Res* 79:4539,
  doi:10.1158/0008-5472.CAN-19-0349) parses it.

**Domain D3: PDO** (CRC Kong/LICOB, bladder Lee 2018, liver Broutier 2017/Ji 2023, pancreas
Tiriac 2018/Shi 2022, sarcoma Al Shihabi 2024, plus new sets)
- Label: within-drug z(AUC), deconfounded.
- Status: mostly curated.

**Domain T: patients (labelled)**
- TCGA RECIST, 2,769 patient x drug labels (training/LOCO).
- CGGA TMZ cohort: Zhao et al. 2021, *Genomics Proteomics Bioinformatics* 19:1,
  doi:10.1016/j.gpb.2020.10.005. Label: OS/PFS landmarked at TMZ start, with untreated or RT-only
  comparators.
- External trial cohorts (GEO about 25, ENLIGHT set, CTR-DB 2.0 additions): **frozen
  validation only**.

**Domain U: patients (unlabelled)**
- All TCGA primary tumours (about 10k) and CGGA, for alignment only.

**Features**
- 673 Reactome pathways scored per sample with a **rank-based single-sample score**: singscore
  (Foroutan et al. 2018, *BMC Bioinformatics* 19:404, doi:10.1186/s12859-018-2435-4), or the repo's
  rank score. Rank scores are platform-robust (RNA-seq vs microarray) and fast.

**Covariates**
- Proliferation (Hallmark E2F + G2M);
- tumour purity (ESTIMATE or consensus purity, patients only);
- treatment setting;
- domain;
- cancer type.

**Drug harmonisation**
- Map drugs across domains to ChEMBL/DrugBank IDs and the ATLAS mechanism classes.
- A (drug, cancer) cell enters the "transfer set" only if it has patients **and** at least one
  pre-clinical domain.
- Expected: about 15-25 drugs and about 40-60 drug x cancer cells.

### 2.2 Alignment step (linear, so pathways stay readable)

1. **Fidelity filter.** Run Celligner-style cPCA + MNN (or MOBER) on gene expression. Drop PDX/PDO
   models whose nearest-neighbour patient tumours are mostly from another cancer type. Report how
   many are dropped.
2. **Remove domain-specific axes.**
   - Run contrastive PCA in pathway space: patient tumours (background = models) give
     patient-only axes, which are typically stroma, immune and normal contamination.
   - Models (background = patients) give model-only axes: culture stress, mouse-stroma leakage,
     hypoxia/ECM.
   - Project the top k (k chosen by permutation, typically 3-6) out of *all* domains.
   - Immune/stroma pathways are therefore never learned from PDO/PDX. They can still enter
     through the patient-only layer of the model (2.3).
3. **Consensus subspace (PRECISE).**
   - For each cancer type c, compute PCA (d = 30-50) on patient pathway scores and on pooled model
     pathway scores.
   - Compute principal vectors and keep those with cosine similarity > 0.7, or above a
     permutation null.
   - Use the interpolated consensus representation (`precise` package, linear kernel = PRECISE;
     optionally TRANSACT with a mild RBF as a sensitivity analysis).
   - Result: K_c consensus factors per cancer, each with a pathway loading vector L_c.
4. **Optional OT variant (sensitivity analysis only).** Within each cancer, apply entropic OT
   (Sinkhorn, POT) barycentric mapping of model samples onto the patient distribution in the
   consensus space. Keep this only if it improves leave-one-cancer-out AUROC in step 2.5.
5. **No leakage.** Alignment uses only unlabelled expression. Patients in the held-out test
   cancer or test cohort are **excluded** from fitting alignment in strict mode. Transductive mode
   (alignment fitted on unlabelled test expression) is reported separately and labelled as such.

### 2.3 Model: hierarchical commensurate-prior ordinal regression

For sample i in domain m ∈ {CL, PDX, PDO, PAT}, drug d (class g(d)) and cancer c, with aligned
factors z_i ∈ R^K:

```
y_i ~ OrderedLogit(eta_i, cutpoints_{m,d})      # CL/PDO z(AUC) is binned to 3 levels or uses a Gaussian link
eta_i = a_{m,d,c} + s_i'gamma_m + z_i' theta_{m,d,c}
theta_{PAT,d,c} = beta_shared + beta_class[g] + beta_drug[d] + beta_dc[d,c]   (patient effects)
theta_{m,d,c}   = theta_{PAT,d,c} + delta_{m,d}   for m in {CL, PDX, PDO}        (model effects)
delta_{m,d}  ~ N(0, tau_m^2)        # commensurate discount, learned per domain
beta_drug    ~ regularised horseshoe
beta_dc      ~ N(0, sigma_dc^2)     # partial pooling across cancers (EXNEX-style mixture as a variant)
beta_class, beta_shared ~ N(0, sigma^2)
s_i = proliferation, purity (PAT only), setting (PAT only)
```

- **Learned discount.** τ_m is the data-driven "how much does PDX/PDO/cell-line biology transfer
  to patients". A large posterior τ_PDO means organoids are ignored automatically. This replaces
  ATLAS/TRIAD's fixed τ = 0.05.
- **Cell lines.** To keep compute down, D1 enters as a fixed prior mean: per-drug ridge
  coefficients in consensus space, with an extra variance term, instead of about 1,000 x 200 rows.
- **Fitting.**
  - Primary: MAP / empirical Bayes (L-BFGS, Laplace approximation for uncertainty).
  - Final per-drug models: NUTS (PyMC or numpyro on CPU), restricted to the transfer set.
- **Baselines** on the same splits:
  1. cancer-type-only;
  2. ATLAS (current);
  3. patient-only BRIDGE (τ_m → ∞);
  4. models-only zero-shot (no patient labels; PRECISE/TRANSACT-style);
  5. one deep benchmark (CODE-AE pretrained weights, or PharmaFormer through `drevalpy`) where
     feasible.

### 2.4 Pathway extraction

1. Back-project to pathways: w_{d,c} = L_c · (β_shared + β_class + β_drug + β_dc), computed per
   posterior draw. This gives the posterior mean, 90% CrI and **lfsr** for every
   (pathway, drug, cancer) cell.
2. **Second track (assumption-light).** For each (drug, cancer, pathway):
   - fit a per-domain marginal effect on aligned pathway scores, adjusted for proliferation,
     purity and setting;
   - combine across domains with a random-effects meta-analysis (DerSimonian-Laird or REML);
   - BH-FDR within drug.
   This gives classical q-values and an I² that shows whether models and patients agree.
3. **Acceptance rule for a reported biomarker** (pre-specified):
   - (a) lfsr < 0.05 in the patient layer **and** the same sign in at least 1 pre-clinical domain;
   - (b) stable in at least 80% of 200 bootstrap re-fits, including alignment re-fitting;
   - (c) passes the predictive test in 2.5.4 where a comparator exists, otherwise labelled
     "prognostic-or-unknown";
   - (d) the pathway is not a proliferation, purity or stroma proxy (|ρ| < 0.5 with those
     covariates).
4. **Output tables:**
   - `bridge_biomarkers_{drug}_{cancer}.tsv`: pathway, effect, CrI, lfsr, q, I², domain signs,
     predictive-test p, status tag;
   - one heat-map per drug: pathway x cancer.

### 2.5 Validation plan (all pre-registered in a `BRIDGE.md` before running)

1. **Leave-one-cancer-out (primary).**
   - Hold out all patients of cancer c. Keep the PDX/PDO/cell lines of c, since that is the
     intended use case ("new cancer, only models"). Re-fit alignment without c's patients.
   - Metric: mean AUROC over drug x cancer cells, the same 64 cells as ATLAS, plus the new
     PDXE-covered cells.
   - Paired bootstrap against ATLAS (0.558/0.561) and the cancer-only baseline.
   - Target: ≥ 0.60 and Δ ≥ +0.02 with a CI excluding 0.
2. **Leave-one-domain-out.**
   - Zero-shot models → patients quantifies pure transfer.
   - Patients-only quantifies what models add.
   - Report the posterior of τ_m.
3. **Frozen external validation.** Freeze the model on all of TCGA + models, then score:
   - the 21 GEO cohorts;
   - the ENLIGHT/trial cohorts;
   - CTR-DB 2.0 additions;
   - CGGA (TMZ).

   For each, report AUROC with DeLong CI, a 1,000-random-signature null (Venet), and adjustment for
   proliferation and subtype/receptor status. Pre-specify the list of (drug, cancer) tests; Holm
   correction across them.
4. **Predictive vs prognostic.**
   - (a) **Randomised comparators already in repo:**
     - BrighTNess (±carboplatin);
     - I-SPY2 arms vs control;
     - GSE20271 (±paclitaxel);
     - GSE41998 (paclitaxel vs ixabepilone);
     - JBR.10/GSE14814 (adjuvant cisplatin/vinorelbine vs observation);
     - GSE42127 (adjuvant chemo vs none, non-randomised).

     Test score x arm interaction (logistic for pCR, Cox for survival) as the primary
     predictive test.
   - (b) **Drug-swap specificity in TCGA:** the drug-d score must predict response to d better than
     to other drugs in the same cancer (difference in AUROC, permutation over drug labels).
   - (c) **Treated vs untreated Cox interaction**, landmarked at treatment start (as in NIT), for
     TCGA and CGGA (TMZ vs no TMZ, adjusted for MGMT, IDH, grade). Flagged as observational.
5. **Negative controls.**
   - whole-pipeline label permutation (≥ 200x on the MAP path);
   - random pathway sets of matched size;
   - proliferation-only and purity-only models;
   - a "shuffled-domain" control (pre-clinical labels permuted within drug), which should drive
     τ_m up.
6. **Calibration and utility.** Calibration slope and intercept on external cohorts; decision
   curves for the 2-3 cohorts large enough.

### 2.6 Expected pitfalls and mitigations

| Pitfall | Mitigation |
|---|---|
| PDO/PDX lack immune system and stroma (PDX has mouse stroma), so stroma and immune pathways cannot transfer | Project them out via cPCA; allow them only in the patient layer; tag biomarkers as "model-supported" vs "patient-only" |
| Viability AUC tracks proliferation and general sensitivity (shown in NIT for LICOB) | Keep NIT deconfounding; proliferation covariate in every domain |
| PDXE 1x1x1 design (one mouse per model x drug); mRECIST is noisy; single-agent in mice vs combination regimens in patients | Ordinal likelihood with domain-specific cut-points; map combinations to the mechanism-class layer, not the drug layer; secondary continuous label |
| Little drug overlap: PDXE is rich in targeted agents absent from TCGA RECIST | Restrict claims to the transfer set; other drugs get model-only biomarkers labelled "pre-clinical only" |
| Alignment removes real biology, or leaks test information | Strict vs transductive modes; alignment re-fitted inside every fold and bootstrap |
| TCGA RECIST labels are noisy, line of therapy is mixed, treatment is not randomised | Setting covariates; label-quality filter from `tcga_response_SUMMARY`; predictive claims only from randomised comparators |
| Small n per drug x cancer (15-50), so the posterior equals the prior | Report prior-to-posterior shrinkage; require patient-layer lfsr; do not call a biomarker from priors alone |
| Many tests (673 pathways x drugs x cancers) | lfsr plus hierarchical shrinkage; BH within drug for track 2; pre-registered external tests with Holm correction |
| Ceiling effect: the realistic AUROC ceiling for chemotherapy from bulk transcriptome is about 0.65-0.70 | Pre-specify success as a relative gain plus at least 1 externally validated predictive biomarker, not an absolute AUROC |
| Platform differences (microarray GEO vs RNA-seq) | Rank-based pathway scores; within-cohort centring at scoring time only |

### 2.7 Estimated compute (4 CPUs, no GPU, ≤ 16 GB RAM)

| Step | Size | Estimate |
|---|---|---|
| Rank-based pathway scoring (singscore-style) | about 20k samples x 673 sets | 10-30 min. gseapy ssGSEA would take hours, so avoid it. |
| Fidelity filter (cPCA + MNN on genes) | 1.5k models x 10k tumours | under 10 min |
| cPCA + PRECISE consensus per cancer | 673-dim; about 20 cancers | under 1 min per fit; about 1 h with bootstrap re-fits |
| OT variant (Sinkhorn, about 500 x 100 per cancer) | | seconds per cancer |
| MAP/EB hierarchical fit | about 5k labelled rows, K ≈ 30, about 5-10k parameters | 1-5 min per fit |
| LOCO (about 20 folds x 5 variants) | | 2-8 h |
| Bootstrap (200 x MAP, with alignment) | | 4-10 h (run overnight, 4 workers) |
| NUTS final per-drug models (about 20 drugs, 4 chains) | | 10-40 min each, about 4-12 h in total |
| External scoring and null signatures (1,000 x about 70 tests) | | under 1 h |
| **Total** | | **about 1-2 CPU-days wall time**; no GPU needed. CODE-AE/PharmaFormer benchmarks on CPU add about 0.5-1 day, or can be skipped by using published weights. |

### 2.8 Suggested implementation layout (for a later session; not created now)

- `robust/data_builders/build_pdxe.py`: PDXE expression + mRECIST into the catalog format.
- `robust/obd/bridge_align.py`: fidelity filter, cPCA, PRECISE consensus, OT variant.
- `robust/obd/bridge.py`: hierarchical model (MAP + NUTS), pathway back-projection, lfsr,
  track-2 meta-analysis.
- `robust/run_bridge.py`, `robust/run_bridge_external.py`, `robust/BRIDGE.md`.
- Dependencies: `precise`/`transact` (pip), `pymc` or `numpyro`, `POT`, `statsmodels`,
  `lifelines`.

---

## 3. Public datasets pairing PDO/PDX response with matched patient clinical response

Access tags:
- **Open**: downloadable without an application.
- **Controlled**: EGA or dbGaP data-access committee.
- **Verify**: the accession or access level was not confirmed against the paper's data statement
  during this survey.

| Dataset | Model / cancer | Matched patient response | Expression | Access |
|---|---|---|---|---|
| Vlachogiannis 2018, *Science* 359:920 | PDO, GI metastases (CRC, gastro-oesophageal) | 21 patients (RECIST/PFS) | Targeted sequencing; RNA-seq for a subset | **Controlled (EGA) - verify** |
| Ooft 2019, *Sci Transl Med* 11:eaay2574 (TUMOROID) | PDO, mCRC | Irinotecan, 5-FU/irinotecan, 5-FU/oxaliplatin; prospective | WES/RNA | **Controlled (EGA) - verify** |
| Yao 2020, *Cell Stem Cell* 26:17 | PDO, rectal, chemoradiation; 80 patients in a phase III trial | Tumour regression grade | Organoid responses: Mendeley Data doi:10.17632/yxvr85b69r.1 (**open**); sequencing not confirmed | Open (responses) / verify (expression) |
| Ganesh 2019, *Nat Med* 25:1607, doi:10.1038/s41591-019-0584-2 | PDO, rectal | Small matched set (chemo/RT) | RNA/WES | **Verify** |
| Tiriac 2018, *Cancer Discov* 8:1112, doi:10.1158/2159-8290.CD-18-0349 | PDO, PDAC | About 9 patients with clinical follow-up (supplement) | RNA-seq: dbGaP phs001611 (**controlled**); organoid AUC + expression via HCMI/CoderData already in repo | Partly open |
| Pasch 2019, *Cancer Discov* 9:852, doi:10.1158/2159-8290.CD-19-0289 | PDO, mixed (CRC, PDAC, breast and others) | Case-level concordance | Limited | **Verify** |
| Narasimhan 2020, *Clin Cancer Res* 26:3662, doi:10.1158/1078-0432.CCR-20-0073 | PDO, peritoneal metastases | Small matched set | Limited | **Verify** |
| de Witte 2020, *Cell Rep* 31:107762, doi:10.1016/j.celrep.2020.107762 | PDO, ovarian | Carboplatin/paclitaxel concordance in a subset | WGS/RNA | **Controlled (EGA) - verify** |
| Al Shihabi 2024, *Cell Stem Cell* 31:1524, doi:10.1016/j.stem.2024.08.010 | PDTO, sarcoma (24 subtypes) | Outcomes for a subset (supplement) | RNA-seq on Synapse syn61892224 (in repo) | Open (Synapse) - patient table: verify |
| Kim et al. 2024, *Cell Rep Med* (gastric PDO; PII S2666-3791(24)00331-8) | PDO, gastric | 11 of 12 patients concordant | **Verify** | Verify |
| Minoli 2023, *Nat Commun*, doi:10.1038/s41467-023-37696-2 | PDO, bladder (longitudinal) | Clinical history vs clonal evolution and drug response | WES/RNA | **Verify** |
| Lee 2018, *Cell* 173:515, doi:10.1016/j.cell.2018.03.017 | PDO, bladder | No matched clinical response (organoid-only) | GSE103990 | Open |
| Ji 2023, *Sci Transl Med* 15:eadg3358 (LICOB) | PDO, liver | Retrospective correlation for some drugs **(verify)** | figshare/CoderData (in repo) | Open (models) |
| Wong et al. 2025, *JHEP Rep* 7 | PDO, HCC (23 lines, 100 drugs) | Retrospective patient correlations | RNA-seq | **Verify** |
| Izumchenko 2017, *Ann Oncol* 28:2595 | PDX, 92 patients, mixed cancers | Matched patient response | Limited | **Verify** (likely supplement only) |
| NCI PDMR (pdmr.cancer.gov) | PDX/PDC/PDOrg, many cancers | Patient treatment history and best response for some donors | RNA-seq/WES | **Open** (registration) |
| Novartis PDXE (Gao 2015) | PDX, 6 indications | **No** matched patient response (PDX clinical trial only) | RNA-seq FPKM (supplement) | Open |
| HCMI (NCI/Wellcome) | PDO and other models + parental tumour | Clinical data (treatment/response sparse) | RNA-seq | Open (processed) / controlled (raw) |
| CTR-DB 2.0 (ctrdb.ncpsb.org.cn) | **Patients only** (no models) | 10,856 patients, 346 regimens | Pre-treatment transcriptome | Open |

Practical note:
- Truly *matched* model-plus-patient response datasets with open expression are rare and small
  (≤ 80 pairs; most ≤ 25).
- They suit a **sanity-check concordance analysis**: does the BRIDGE score computed on the PDO
  agree with the score on the parental tumour, and both with the patient? They are too small to
  train on.
- Yao 2020 (open responses) and NCI PDMR (open, with donor treatment history) are the best first
  targets.
- EGA sets (Vlachogiannis, Ooft, de Witte) fit the existing `docs/data_access/EGA_application_plan.md`
  route.

---

## 4. Top methods to borrow from (ranked for this repo)

1. **PRECISE / TRANSACT** (Mourragui 2019/2021): linear and kernel consensus subspaces between
   models and tumours. This is the only mature transfer family that keeps gene/pathway
   interpretability. Basis of BRIDGE step 2.2.
2. **PDXGEM** (Kim 2020): concordance filtering, keeping only features whose co-expression
   structure is conserved between PDX and patients. Cheap, interpretable, and portable to Reactome
   scores as a pre-filter.
3. **Hierarchical borrowing** (commensurate/power priors, EXNEX, horseshoe; TG-LASSO's tissue
   guidance). Learns how much to trust each pre-clinical domain and pools across cancers only
   where they look alike. This is BRIDGE's core.
4. **SELECT / ENLIGHT** (Lee 2021; Dinstag 2023): unsupervised, knowledge-based scoring validated
   across about 20-35 blinded trials. Use as a mandatory external baseline. The repo already holds
   the ENLIGHT cohorts.
5. **CODE-AE / PharmaFormer** (2022 / 2025): the strongest non-linear model-to-patient transfer
   results, including organoid fine-tuning. Use as performance ceilings or benchmarks, not as
   biomarker sources.
