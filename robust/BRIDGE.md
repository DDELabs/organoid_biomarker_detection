# BRIDGE: patient-derived models → patient tumours

**Question:** do drug-response models trained on patient-derived tumour models carry over
to patients at the pathway level?
- PDX: Novartis PDX Encyclopedia, 399 models with RNA-seq, mRECIST-style response.
- Organoids: pancreas ×2, liver ×2, bladder, sarcoma.
- Cell lines: CTRPv2, PRISM.

Survey and design rationale: `docs/BRIDGE_strategy_research.md`.
Code: `obd/bridge.py`, `run_bridge.py`.
Results: `results/bridge/`, with these files:
- `summary.tsv`
- `pathways.tsv.gz`
- `concordance.tsv`
- `external.tsv`
- `domains.tsv`

## Method

1. Reactome rank scores for every model, standardised within tumour type or lineage.
2. **Alignment:** principal axes of the model data with less than one third of their
   variance in TCGA tumours (e.g. no stroma or immune cells, culture effects) are projected
   out. This is linear, so pathways stay interpretable.
3. **Per-domain sensitivity vector:** RidgeCV on rank-normalised sensitivity. Vectors are
   averaged within class (PDX / organoid / cell line), then across classes, to give the
   pre-clinical prior.
4. **Patient model:** TRIAD's prior-centred ridge logistic on TCGA RECIST labels. The
   **transfer strength s** (0–8 × τ) is learned per drug by leave-one-cancer-out (LOCO),
   **nested** for the reported AUROC.
5. **Benchmarks:**
   - **Zero-shot:** pre-clinical vector alone scores patients. Null = vector permuted
     across pathways.
   - **External:** frozen BRIDGE score in 30 trial cohorts, with a random-signature null.
6. **Pathways:** prior-free patient bootstrap z (a prior-centred bootstrap inflates z).
   "Consensus" = |z| ≥ 2 in patients and the same sign as the pre-clinical vector.

## Results (14 drugs with ≥ 2 evaluable cancers)

| | Mean within-cancer AUROC |
|---|---|
| Patient-only (LOCO) | 0.535 |
| BRIDGE, nested LOCO | 0.518 |
| Zero-shot pre-clinical vector | 0.537 |
| BRIDGE external (54 cohort × drug tests) | 0.529 (ATLAS on the same cohorts: 0.515) |

- **Transfer helps for a few drugs only:**

  | Drug | Patient-only | BRIDGE (nested) | Zero-shot |
  |---|---|---|---|
  | Docetaxel | 0.49 | 0.62 | 0.64 (null p = 0.02; driven by organoids) |
  | Epirubicin | – | – | 0.86 (null p = 0.01; cell lines; 2 cancers, n = 51) |
  | Vinorelbine | – | – | 0.59 (p = 0.10) |

  It hurts for doxorubicin and capecitabine. On average it does not help.
- **PDX:** the zero-shot paclitaxel AUROC is 0.575 (p = 0.20), the best of the three
  domains for that drug, but not significant. 5-FU and gemcitabine PDX vectors do not
  transfer.
- **External validation:**
  - Paclitaxel is significant in GSE41998 (0.64, p = 0.024), BrighTNess (0.60, p = 0.020),
    I-SPY2 (0.58, p = 0.036) and the ENLIGHT trastuzumab cohort (0.70, p = 0.048).
  - Doxorubicin is significant in BrighTNess (0.60, p = 0.024).
  - Epirubicin is borderline in GSE32646 (0.68, p = 0.056).

  The external signal is again concentrated in breast taxane/anthracycline cohorts.
- **Pathway concordance between patients and pre-clinical models is weak** (Spearman
  \|ρ\| < 0.17 for every drug). Sign agreement above chance (binomial; pathways overlap, so
  the p values are optimistic):

  | Drug | Sign agreement | p |
  |---|---|---|
  | Epirubicin | 72% | < 0.001 |
  | Vinorelbine | 75% | 0.001 |
  | Gemcitabine | 60% | 0.10 |

  Everything else is at chance.

## Cross-domain consensus pathways (patient z ≥ 2 and same sign pre-clinically)

These agree with the ATLAS findings:
- **Cisplatin:** cell cycle, DNA replication, M/G1 → response.
- **Epirubicin and doxorubicin:** Fanconi anaemia pathway and DNA repair → response.
  Epirubicin: z = 4.3, FDR q = 0.009.
- **Gemcitabine:** semaphorin interactions and TGF-β/EMT → resistance; E2F-mediated
  replication control → response.
- **Paclitaxel:** antigen-activated B-cell receptor → response.
- **5-FU:** Fanconi anaemia → response; BMP signalling and phase II conjugation →
  resistance.

## Conclusion

Patient-derived models are useful as **corroboration** for a handful of mechanistic
pathway biomarkers: DNA-repair/Fanconi → anthracyclines, cell cycle → platinum, semaphorin
→ gemcitabine. They are **not** a general predictor of patient response. Pathway effects
learned in PDX, organoids or cell lines are mostly uncorrelated with those in patients, and
learned transfer does not improve held-out AUROC on average.

The limiting factors are:
- PDX and organoids lack the immune system and human stroma;
- PDXE shares few drugs with TCGA;
- there is a single mouse per PDX arm;
- TCGA's RECIST labels are noisy.

Public PDO sets paired with the matched patient's response, which would test transfer
directly, are mostly in controlled-access EGA (see `docs/data_access`).
