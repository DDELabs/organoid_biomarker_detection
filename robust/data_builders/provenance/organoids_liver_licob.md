# liver_licob: LICOB primary liver cancer organoids (Ji et al. 2023)

**Paper:** Ji S, Feng L, Fu Z, et al. "Pharmaco-proteogenomic characterization of liver cancer organoids for precision oncology." *Sci Transl Med* 2023;15(706):eadg3358.

## Where the file came from

- https://raw.githubusercontent.com/wu-yc/iLICOB/master/data/data_ilicob_org: R serialized list from the iLICOB R package. The same file is stored in the repo as `data/external/LICOB_organoid.RData`.
  sha256 `5fb1460932286215d8368e1d0f00a481345ae1ac88fafb21e51068c72b43b135`
- The objects used are `[[4]][[1]]` (RNA, genes × 50 organoids) and `[[4]][[6]]` (drug AUC, 50 organoids × 76 drugs). This is the same parse as `robust.obd.cohorts.licob_organoids()`.

## Processing

- **Expression:** written exactly as distributed, with gene symbols as given. The values are on a **log2 scale**, apparently log2(TPM + small pseudocount), with a floor of about −5.6. They are not counts.
- **Response:** the metric is **AUC** as distributed (range about 0.03–2.0). **Lower = more sensitive.** I checked this against the data: ARV-771, methotrexate and panobinostat have the lowest medians, while FG-4592 and olaparib are about 1.
  - Drug names were converted to upper case common names (`_` → `-`, then DrugBank synonyms).
- **Overlap: 50 / 50 organoids** (HCCO*, ICCO*, CHCO*, HB*).

## Caveats

- **No IC50 is available.** The iLICOB release only has AUC. IC50 and clinical tables (for example, whether the patient received sorafenib) are in the STM supplementary data, which is blocked here (journal site and figshare). CoderData's liver build also reads that supplement through Synapse (syn64961953), which is also blocked.
- I inspected `data_ilicob_tissue`. It does **not** contain TCGA-LIHC data. It holds:
  - `[[1]]`: 76 caret elastic-net (`enet`) models, one per drug, used by `iLICOB_predict` for tissue input.
  - `[[2]]`: the 131 feature genes.
  - `[[3]]`: a 131-gene × 29-sample demo matrix that looks z-scored. Its columns are LICOB IDs (ICCO7, CHCO2, HCCO20 and so on), i.e., tissue profiles of LICOB patients.
  - It adds no response data.
- Tumour types are mixed (HCC, ICC, CHC, HB). Subset by prefix if needed.
