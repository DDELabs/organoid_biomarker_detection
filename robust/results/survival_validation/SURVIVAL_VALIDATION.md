# Survival validation of the top-10 pathway biomarkers (TCGA, held-out patients)

Signature score per cancer x drug = mean of its 5 response pathways - mean of its 5 resistance pathways. Held-out treated = received the drug but not used for discovery. HR per SD, Cox adjusted for age and stage (HR < 1 = higher score, longer survival). Interaction HR < 1 = the survival advantage of a high score is larger in treated than in untreated patients (predictive). Logrank = the original paper's median split.

## OS: independent pairs

- evaluable pairs: 20; HR < 1 (expected direction): 13/20 (binomial p = 0.132)
- nominal p < 0.05: 3 (2 in the expected direction)
- Stouffer z across pairs (expected direction positive): 2.49
- median-split log-rank p < 0.05 (paper style): 3
- untreated HR < 1: 13/18; interaction HR < 1: 12/20, of which p < 0.05: 0

Top 15 pairs by treated p:

| cancer | drug | tier | n_heldout | treated_events | treated_HR | treated_p | treated_logrank_p | untreated_HR | interaction_HR | interaction_p |
|---|---|---|---|---|---|---|---|---|---|---|
| LGG | TEMOZOLOMIDE | A20 | 128 | 45 | 0.495 | 1.07e-05 | 0.00143 | 0.987 | 0.583 | 0.0599 |
| COAD | FLUOROURACIL | A20 | 63 | 9 | 4.01 | 0.0132 | 0.049 | nan | 2.01 | 0.0558 |
| LUSC | CARBOPLATIN | A20 | 29 | 15 | 0.498 | 0.0371 | 0.0619 | 0.899 | 0.699 | 0.226 |
| HNSC | CISPLATIN | A20 | 49 | 22 | 0.608 | 0.0541 | 0.787 | 0.998 | 0.63 | 0.0655 |
| LUAD | PEMETREXED | A20 | 31 | 8 | 0.23 | 0.0572 | 0.161 | 0.917 | 0.428 | 0.0514 |
| LUAD | CARBOPLATIN | A20 | 34 | 16 | 0.542 | 0.0738 | 0.703 | 1.05 | 0.684 | 0.158 |
| OV | DOCETAXEL | B | 47 | 36 | 1.27 | 0.223 | 0.137 | 0.802 | 1.38 | 0.274 |
| OV | GEMCITABINE | B | 79 | 63 | 0.852 | 0.286 | 0.0195 | 0.728 | 1.97 | 0.0129 |
| LGG | LOMUSTINE | A | 20 | 13 | 0.699 | 0.329 | 0.988 | 0.882 | 0.769 | 0.531 |
| OV | PACLITAXEL | B | 61 | 39 | 0.87 | 0.418 | 0.278 | 0.986 | 0.911 | 0.725 |
| STAD | FLUOROURACIL | A20 | 36 | 16 | 1.19 | 0.61 | 0.669 | 0.835 | 1.29 | 0.389 |
| CESC | CISPLATIN | A20 | 57 | 18 | 0.863 | 0.65 | 0.466 | nan | 0.648 | 0.246 |
| OV | CARBOPLATIN | B | 66 | 42 | 0.931 | 0.707 | 0.889 | 1.01 | 0.857 | 0.576 |
| UCEC | CARBOPLATIN | A20 | 110 | 27 | 1.06 | 0.768 | 0.268 | 1.08 | 1.07 | 0.803 |
| LUSC | CISPLATIN | A20 | 41 | 12 | 1.09 | 0.776 | 0.701 | 0.896 | 1.01 | 0.972 |

Single pathways (OS, held-out treated): 128/222 in the expected direction (binomial p = 0.013; pathways overlap, so optimistic); p < 0.05 and concordant: 19, discordant: 3

## PFI: independent pairs

- evaluable pairs: 29; HR < 1 (expected direction): 17/29 (binomial p = 0.229)
- nominal p < 0.05: 4 (4 in the expected direction)
- Stouffer z across pairs (expected direction positive): 2.95
- median-split log-rank p < 0.05 (paper style): 3
- untreated HR < 1: 14/24; interaction HR < 1: 15/27, of which p < 0.05: 2

Top 15 pairs by treated p:

| cancer | drug | tier | n_heldout | treated_events | treated_HR | treated_p | treated_logrank_p | untreated_HR | interaction_HR | interaction_p |
|---|---|---|---|---|---|---|---|---|---|---|
| LGG | TEMOZOLOMIDE | A20 | 128 | 60 | 0.508 | 3.19e-06 | 3.64e-06 | 0.928 | 0.605 | 0.021 |
| HNSC | CISPLATIN | A20 | 49 | 23 | 0.563 | 0.0238 | 0.224 | 1.08 | 0.562 | 0.0191 |
| READ | OXALIPLATIN | A | 23 | 9 | 0.236 | 0.0373 | 0.0666 | nan | nan | nan |
| LUSC | CARBOPLATIN | A20 | 29 | 16 | 0.552 | 0.0425 | 0.0308 | 0.944 | 0.599 | 0.0734 |
| HNSC | CARBOPLATIN | A20 | 34 | 18 | 0.707 | 0.0515 | 0.602 | 0.933 | 0.817 | 0.283 |
| LUAD | CARBOPLATIN | A20 | 34 | 20 | 0.611 | 0.0522 | 0.648 | 0.956 | 0.702 | 0.126 |
| COAD | OXALIPLATIN | A20 | 59 | 23 | 0.703 | 0.054 | 0.0482 | nan | 0.644 | 0.051 |
| LGG | LOMUSTINE | A | 20 | 18 | 0.661 | 0.247 | 0.31 | 0.998 | 0.773 | 0.491 |
| LUSC | VINORELBINE | A | 23 | 8 | 1.73 | 0.283 | 0.604 | 1.03 | 1.51 | 0.39 |
| HNSC | PACLITAXEL | A20 | 25 | 10 | 0.678 | 0.303 | 0.926 | 0.864 | 0.763 | 0.477 |
| OV | CARBOPLATIN | B | 66 | 48 | 0.856 | 0.33 | 0.801 | 1.17 | 0.742 | 0.309 |
| UCEC | DOXORUBICIN | A | 21 | 16 | 0.729 | 0.337 | 0.461 | 0.815 | nan | nan |
| READ | FLUOROURACIL | A20 | 26 | 8 | 1.55 | 0.362 | 0.906 | nan | 1.38 | 0.523 |
| OV | CISPLATIN | B | 46 | 42 | 1.18 | 0.467 | 0.923 | 1.4 | 1.07 | 0.832 |
| STAD | FLUOROURACIL | A20 | 36 | 9 | 1.38 | 0.469 | 0.989 | 0.736 | 1.52 | 0.306 |

Single pathways (PFI, held-out treated): 165/282 in the expected direction (binomial p = 0.003; pathways overlap, so optimistic); p < 0.05 and concordant: 19, discordant: 8


## Notes

- **LGG temozolomide:** the result holds after adjusting for grade, histology and age.
  - OS HR 0.52 (p = 3e-5); PFI HR 0.52 (p = 1e-5); n = 128 held-out.
  - IDH and 1p/19q status are not in TCGA-CDR, so molecular-subtype confounding cannot be excluded.
- **Pairs with too few held-out patients:** most RECIST-labelled patients were used for
  discovery, so many pairs are left with too few others (e.g. UCS, STAD oxaliplatin).
- **Tier C pairs** (discovered on PFI) are not independent and are excluded from the summaries.
- **Treatment was not randomised.** "Untreated" means no recorded systemic therapy, which can
  also mean the therapy was not recorded.
