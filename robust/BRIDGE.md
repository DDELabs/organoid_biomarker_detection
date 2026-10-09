# BRIDGE: patient-derived models → patient tumours

**Question:** do drug-response models trained on patient-derived tumour models carry over
to patients at the pathway level?
- PDX: Novartis PDX Encyclopedia (399 models with RNA-seq, mRECIST-style response) and Isella 2017
  (192 colorectal PDX, cetuximab).
- Organoids: pancreas ×2, liver ×2, bladder, sarcoma, colorectal.
- Cell lines: CTRPv2, PRISM.

Survey and design rationale: `docs/BRIDGE_strategy_research.md`.
Code: `obd/bridge.py`, `run_bridge.py`.
Results: `results/bridge/`, with these files:
- `summary.tsv`
- `pathways.tsv.gz`
- `concordance.tsv`
- `external.tsv`
- `domains.tsv`
- `cetuximab_transfer.tsv`, `cetuximab_signature.tsv`

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

Pre-clinical domains:
- PDX: Novartis PDXE and Isella 2017 colorectal.
- Organoids: 7 sets, including van de Wetering 2015 colorectal.
- Cell lines: CTRPv2 and PRISM.

| | Mean within-cancer AUROC |
|---|---|
| Patient-only (LOCO) | 0.535 |
| BRIDGE, nested LOCO | 0.499 |
| Zero-shot pre-clinical vector | 0.544 |
| BRIDGE external (54 cohort × drug tests) | 0.529 (ATLAS on the same cohorts: 0.515) |

- **Zero-shot transfer works for a few drug × domain combinations** (null = vector permuted
  across pathways):

  | Domain | Drug | AUROC | p |
  |---|---|---|---|
  | Organoids | Cisplatin | 0.62 | 0.04 |
  | Organoids | Docetaxel | 0.65 | 0.024 |
  | Cell lines | Epirubicin | 0.85 | 0.008 (2 cancers, n = 51) |
  | Cell lines | Vinorelbine | 0.59 | 0.096 |

  The gemcitabine organoid p = 0.016 is not credible (AUROC 0.49; the null is shifted).
- **Learned transfer does not beat patient-only on average.** The per-drug gains are
  fragile: docetaxel went from 0.49 to 0.62 with the first domain set, and back to 0.49 once
  the colorectal organoids were added. Choosing s inside the nested folds is noisy, because
  only 2–8 cancers are available per drug.
- **PDX:** zero-shot paclitaxel AUROC 0.575 (p = 0.20); 5-FU and gemcitabine PDX vectors do
  not transfer to TCGA.
- **External validation:**
  - Paclitaxel is significant in GSE41998 (0.64, p = 0.024), BrighTNess (0.60, p = 0.020),
    I-SPY2 (0.58, p = 0.036) and the ENLIGHT trastuzumab cohort (0.70, p = 0.048).
  - Doxorubicin is significant in BrighTNess (0.60, p = 0.024).
  - Epirubicin is borderline in GSE32646 (0.68, p = 0.056).

  The signal is again concentrated in breast taxane/anthracycline cohorts.
- **Pathway concordance between patients and pre-clinical models is weak** (Spearman
  \|ρ\| ≤ 0.17 for every drug). Sign agreement above chance (binomial; pathways overlap, so
  the p values are optimistic):

  | Drug | Sign agreement | p |
  |---|---|---|
  | Epirubicin | 67% | 0.001 |
  | Vinorelbine | 71% | 0.012 |

  Everything else is at chance (`concordance.tsv`).

### Cetuximab: PDX → PDX → patients (`run_cetuximab_transfer.py`, `cetuximab_transfer.tsv`)

Signature learned on 192 colorectal liver-metastasis PDX (Isella 2017; tumour volume change
after 3 weeks of cetuximab). Benchmark: AREG + EREG expression.

| Test | n | Signature | AREG + EREG |
|---|---|---|---|
| Isella PDX, 5-fold CV (Spearman) | 192 | **0.59** | 0.37 |
| Novartis PDXE colorectal (other lab; Spearman / AUROC benefit) | 41 | **0.40** (null p = 0.11) / **0.81** | 0.00 |
| Novartis PDXE lung | 25 | 0.00 | 0.07 |
| Head & neck patients, cetuximab + chemo (ENLIGHT; AUROC) | 40 | 0.42 | **0.79** |
| TCGA head & neck cetuximab (AUROC) | 20 | 0.16 | – |

The PDX signature transfers between PDX of the same tissue from different labs, where it
beats the ligand marker. It does not transfer across tissue (lung PDX), and it does not
transfer to patients. The colorectal monotherapy patient cohort (GSE5851) has no public
per-sample response.

## Cross-domain consensus pathways (patient z ≥ 2 and same sign pre-clinically)

These agree with the ATLAS findings:
- **Epirubicin:** DNA repair (z = 4.3, q = 0.01) and the Fanconi anaemia pathway (3.8,
  q = 0.02) → response. Doxorubicin points the same way (Fanconi, z = 2.6–2.7).
- **Cisplatin:** cell cycle and M/G1 → response (z ≈ 3.2, q = 0.19).
- **Gemcitabine:** semaphorin interactions → resistance; PERK/E2F replication control →
  response.
- **Paclitaxel:** antigen-activated B-cell receptor → response.
- **5-FU:** phase II conjugation (drug detoxification) → resistance (z = −3.9, q = 0.07).
- **Oxaliplatin:** IRAK2–TAK1 (innate immune signalling) → response (z = 4.1, q = 0.02).

## Conclusion

Patient-derived models are useful as **corroboration** for a handful of mechanistic
pathway biomarkers: DNA-repair/Fanconi → anthracyclines, cell cycle → platinum, semaphorin
→ gemcitabine. They are **not** a general predictor of patient response. Pathway effects
learned in PDX, organoids or cell lines are mostly uncorrelated with those in patients, and
learned transfer does not improve held-out AUROC on average.
Cetuximab shows where the gap opens: a PDX signature reproduces in another lab's PDX of the same
tissue (AUROC 0.81 for benefit), but not across tissue or in patients. Model-to-model transfer
works; model-to-patient transfer is the bottleneck.

The limiting factors are:
- PDX and organoids lack the immune system and human stroma;
- PDXE shares few drugs with TCGA;
- there is a single mouse per PDX arm;
- TCGA's RECIST labels are noisy.

Public PDO sets paired with the matched patient's response, which would test transfer
directly, are mostly in controlled-access EGA (see `docs/data_access`).
