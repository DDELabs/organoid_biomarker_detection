# pancreas_tiriac2018: PDAC patient-derived organoids (Tiriac et al. 2018)

**Paper:** Tiriac H, et al. "Organoid Profiling Identifies Common Responders to Chemotherapy in Pancreatic Cancer." *Cancer Discov* 2018;8(9):1112-1129. doi:10.1158/2159-8290.CD-18-0349

## Where the files came from

The original sources (AACR supplementary xlsx and GDC project `ORGANOID-PANCREATIC`) are blocked from this environment. The data were taken from a public GitHub mirror that the authors of a later paper built from those sources:

- Repository: https://github.com/prassepaul/mlmed_transfer_learning (commit on `main`, cloned 2026-10-04). Paper: Prasse et al., *Cancers* 2022, 14(16):3950.
- `data/organoid_pancreas/organoid_pancreas_fpkm.txt`: GDC HTSeq-FPKM, 55 files × 60,483 Ensembl genes (GENCODE v22).
  sha256 `7e261ffb44f5ec0b622f85189d5e30a7c78a95ead535781bb1b3737a40bd66d2`
- `data/organoid_pancreas/organoid_value.tsv`: long table built from Tiriac Supplementary Table S4 (sheets 2–4: chemotherapy, targeted, and targeted-for-chemorefractory). sha256 `cfa21cc49c3d3447d1ea16a77fdae0cfb82aba8e3a051181b695c392b98c8eeb`
- `data/organoid_pancreas/file_name_to_organoid_mapping.tsv`: maps each GDC file to an organoid name, using GDC case → S1 '#' → organoid. sha256 `90d1c984c5120ac579ea667deb9c2e133c1d7fcec60c31da61c3e81bd5d7586c`
- The raw copies are in `scratchpad/org_raw/mlmed/data/organoid_pancreas/`.

## Processing (scratchpad `build_sets.py`)

- **Expression:** Ensembl IDs had their version removed and were mapped to symbols with `data/2017_07_31_biomart_protein_coding_genes.txt`, so only protein-coding genes are kept (19,669 genes). Values are **log2(FPKM + 1)**. Organoids with several GDC files were averaged. Duplicate symbols were also averaged. The result has 49 organoids.
- **Response:** the metric is **AUC** as published (Tiriac normalised dose-response AUC, about 0.23–0.98). **Lower = more sensitive.** I checked this against the data: bortezomib has the lowest median AUC (0.26), followed by SN-38, paclitaxel and gemcitabine, while celecoxib and ruxolitinib are about 0.86–0.87.
  - Drug names were converted to upper case DrugBank common names: 5-FU → FLUOROURACIL, "Disulfuram" → DISULFIRAM.
  - SN-38 is the active metabolite of irinotecan and keeps the name SN-38.
  - One sample ID had trailing whitespace ("hF74 "), which was stripped.
  - The result has 67 organoids and 25 drugs (5 chemotherapies plus 20 targeted agents).
  - mlmed kept only compounds that resolve to PubChem SMILES, so a few targeted agents from S4 may be missing.

## Matching

Organoid names such as hF2 and hM1A are used directly in both files.

- **Overlap (expression ∩ response): 45 organoids.**
- Expression only: hF34, hN30, hN31 (normal-derived), hT88.
- Response only: 22 organoids, including hM1E, hT102 and hT105.

## Caveats

- The GDC files were mapped to organoids through the **case** (patient), not the specimen. Three organoids have several files from one case: **hM1A** (4 files), **hF81** (3) and **hT64** (2). Their expression is a mean over files that may include sibling organoids from the same patient. For example, hM1E has a response entry but no expression of its own. Drop these three for a strict analysis.
- These are FPKM values from GDC GENCODE v22, not TPM.
- The table holds AUC only. The supplement also has no IC50 for most drugs.
