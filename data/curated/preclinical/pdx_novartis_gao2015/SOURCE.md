# Novartis PDX Encyclopedia (Gao et al. 2015)

**Paper**: Gao H, Korn JM, Ferretti S, et al. *High-throughput screening using patient-derived tumor xenografts to
predict clinical trial drug response.* Nat Med 2015;21(11):1318-1325. doi:10.1038/nm.3954 (PMID 26479923).

**Accession**: Supplementary Table S1 (`nm.3954-S2.xlsx`, Springer static content), sheets `RNAseq_fpkm`,
`PCT curve metrics`, `PCT raw data`. Raw sequencing is not used.

| file | URL | sha256 |
|---|---|---|
| `Gao2015_nm3954_TableS1.xlsx` | https://static-content.springer.com/esm/art%3A10.1038%2Fnm.3954/MediaObjects/41591_2015_BFnm3954_MOESM10_ESM.xlsx | `c4b9a6903a4d1f76e3ddca4199039776d56bb99970aa5b7abe4f3abd732a0c6d` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `6e2f79000e79ecb430f2613b878d514302df7f81263b5db979010d032c3235a2` |

## Design
PDX clinical trial (PCT), 1 x 1 x 1 design: one mouse per PDX model x treatment, tumour volume followed for
~3 weeks or more. 280 models are treated; tumour types of the models with both data types:
CRC 51, BRCA 41, PDAC 38, CM 32, NSCLC 27 (BRCA breast, CRC colorectal, GC gastric, NSCLC lung,
PDAC pancreas, CM cutaneous melanoma).

## Derivation
* Response (`response.tsv`): **BestAvgResponse** from the paper (minimum over t >= 10 d of the running mean of the
  % tumour-volume change from day 0). **Lower = more sensitive** (negative = shrinkage). One row per model x
  treatment arm; the `untreated` arm and the low-dose `binimetinib-3.5mpk` arm are excluded (the latter is kept in
  `mrecist.tsv`). Combination arms are kept with the components joined by ` + ` (e.g. `ALPELISIB + ELGEMTUMAB`),
  `treatment_type` = single / combo. Novartis codes were mapped to INNs where one exists
  (BYL719 alpelisib, BKM120 buparlisib, LEE011 ribociclib, INC280 capmatinib, INC424 ruxolitinib, LDK378 ceritinib,
  LDE225 sonidegib, BGJ398 infigratinib, HDM201 siremadlin, LGH447 PIM447, LJM716 elgemtumab, WNT974,
  abraxane nab-paclitaxel); research codes without an INN are kept (CGM097, CLR457, CKX620, HSP990, LCL161,
  LFA102, LFW527, LGW813, LJC049, LKA136, LLM871, TAS266).
* `mrecist.tsv`: mRECIST class recomputed from BestResponse / BestAvgResponse with the thresholds of the paper's
  Online Methods: CR if BestResponse < -95% and BestAvgResponse < -40%; PR if BestResponse < -50% and
  BestAvgResponse < -20%; SD if BestResponse < 35% and BestAvgResponse < 30%; otherwise PD. `responder` = 1 for
  CR/PR. Counts: CR 143, PR 354, SD 1475, PD 2560. The paper's own
  `ResponseCategory` (which also tracks progression after the best response, e.g. `SD-->PD`) is kept as
  `paper_category`; its first class agrees with the recomputed call in 100.0% of arms.
* Expression: `RNAseq_fpkm` (FPKM, Novartis pipeline) -> log2(FPKM + 1); symbols re-mapped to current HGNC
  protein-coding symbols. Units: **log2(FPKM + 1)**. Mouse reads were removed by the authors.

## Caveats
* No matched patient clinical response is public for PDXE (there is no `patient_response.tsv`).
* Single animal per arm: noisy. Treat the binary responder call as the most robust endpoint.
* Response rates are low (CR+PR ~ 6% of arms), and most arms are targeted agents.

## Sample-ID matching
Model IDs (`X-1004` ...) shared by both sheets. Expression models: 399; treated models:
280; **overlap n = 192**. Genes: 18593; drugs/regimens: 61.
