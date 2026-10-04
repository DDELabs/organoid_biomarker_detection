# External validation of frozen 5-FU signatures (colorectal cancer)

Signatures were learned on the 19 van de Wetering organoids only and applied unchanged.
Pathway scores are per-sample rank scores (platform-independent), z-scored within each
cohort. HR is per SD of predicted resistance (> 1 = expected direction), adjusted for
stage, age and sex. "Random-null p" is the fraction of 1,000 random, equal-size,
random-sign Reactome signatures with |log HR| at least as large (Venet et al. 2011).

## GSE39582 (stage II-III; 141 patients on 5-FU-based adjuvant therapy vs 264 untreated)

| Signature | RFS HR treated (95% CI) | RFS interaction p | OS HR treated (95% CI) | OS p | OS interaction HR (p) | Random-null p, OS treated / interaction |
|---|---|---|---|---|---|---|
| paper_top7 (Kong et al.) | 0.95 (0.73-1.24) | 0.047* | 1.05 (0.79-1.39) | 0.74 | 1.07 (0.71) | 0.86 / 0.84 |
| robust_v1 | 1.00 (0.77-1.30) | 0.065* | 1.11 (0.84-1.48) | 0.45 | 1.12 (0.51) | 0.68 / 0.71 |
| v1_rank | 1.02 (0.79-1.32) | 0.075* | 1.12 (0.85-1.48) | 0.43 | 1.13 (0.50) | 0.68 / 0.71 |
| nit_soft | 1.11 (0.85-1.44) | 0.33 | 1.25 (0.94-1.67) | 0.13 | 1.31 (0.13) | 0.38 / 0.38 |
| **nit_deconf** | 1.19 (0.92-1.55) | 0.83 | **1.38 (1.03-1.86)** | **0.033** | **1.49 (0.030)** | 0.19 / 0.21 |

\* interaction HR < 1, the *opposite* of a predictive marker: the score is prognostic in
untreated patients (HR about 1.4) and flat in treated ones.

## Other cohorts

| Cohort | Endpoint | Best signature | Result |
|---|---|---|---|
| GSE14333 (87 presumed 5-FU, RFS) | RFS | nit_deconf | HR 1.27 (0.90-1.79), p = 0.17, same direction |
| GSE28702 (83 mFOLFOX6, response) | response AUC | paper_top7 | AUC 0.55 (p = 0.22); others 0.44-0.49 |
| TCGA-READ (66 treated, 9 deaths) | OS | none | HRs 0.53-0.78, opposite direction, underpowered |

## Conclusion

The organoid-derived NIT signature (deconfounded response, soft network prior) is the
only one that behaves like a 5-FU-*predictive* marker in the largest independent cohort:
nominal OS p = 0.03 and a treatment interaction in the expected direction, plus the same
direction in GSE14333. But it **does not beat random signatures** (random-null p of about
0.2), so it is **not validated**. The original paper signature shows no predictive value
externally. Next decisive data: the EGA-controlled paired organoid/patient cohorts
(see `docs/data_access/`) and GSE72970 once GEO is reachable.
