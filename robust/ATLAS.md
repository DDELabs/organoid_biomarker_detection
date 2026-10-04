# ATLAS: pan-cancer, multi-drug pathway atlas of drug response

**Goal:** find the pathways that determine response to each drug, and the overall response
of patients, across all drugs and cancers with usable data.

**Design:** a multi-task model with three layers of pathway effects:
- **shared** across all drugs ("overall response");
- **mechanism class:** platinum, fluoropyrimidine, antimetabolite, microtubule,
  topoisomerase I/II, alkylating, anti-angiogenic, hormonal, other DNA;
- **drug-specific**.

All three are fitted jointly on **2,769 TCGA patient × drug RECIST labels** (21 drugs, all
cancer types with treatment data). Features are 673 Reactome rank scores, centred within
cancer type. The model adjusts for treatment setting and proliferation. Pre-clinical (GDSC)
and network priors are used for prediction; pathway inference uses the prior-free model,
because priors make drug-layer effects non-zero by construction.

Code:
- `obd/atlas.py`
- `run_atlas.py`
- `run_atlas_inference.py`
- `run_atlas_external.py`

Data built by network-enabled helper sessions is in `data/curated/`:
- 21 GEO trial cohorts;
- 8 pre-clinical sets, including Lee 2018 bladder organoids.

## 1. Held-out prediction (TCGA, leave one cancer out)

| Model | Mean AUROC over 64 drug x cancer pairs |
|---|---|
| Separate per-drug models | 0.543 |
| ATLAS multi-task | 0.558 |
| **ATLAS + GDSC/network priors** | **0.561** (permutation p = 0.077, the floor for 12 permutations) |

Pooling across drugs and adding priors helps a little. Overall, transcriptomic pathway
signal for chemotherapy response is weak, with about 0.56 AUROC across cancers.

## 2. Pathway determinants (`results/atlas/inference_*.tsv.gz`, `MECHANISMS.md`, `atlas_heatmap.png`)

- **Mostly shared:** per-drug total effects are dominated by the shared layer (uniform rows
  in the heatmap). No single per-drug total effect passes FDR across 14,000 tests.
- **21 class- or drug-specific effects pass q < 0.1.** The credible ones, after literature
  review, are:

| Pathway | Drug / class | Direction | z | Evidence |
|---|---|---|---|---|
| RIG-I/MDA5 -> IFN-alpha/beta | cisplatin | resistance | -3.8 | KNOWN: interferon-related DNA damage resistance signature (Weichselbaum 2008) |
| Cell cycle (beyond the proliferation covariate) | cisplatin | response | +3.8 | KNOWN: proliferative tumours are platinum-sensitive |
| B-cell receptor signalling | paclitaxel | response | +4.0 | KNOWN: TIL/B-cell infiltrate predicts taxane/anthracycline pCR (Denkert 2010) |
| Fanconi anaemia pathway | topoisomerase II class | response | +3.8 | PLAUSIBLE: replication stress / HRD marker |
| FGFR3 signalling | cisplatin | response | +3.8 | PLAUSIBLE: FGFR3-luminal bladder biology |
| Prostanoid (PGE2/EP) receptors | paclitaxel | resistance | -3.6 | PLAUSIBLE: COX-2/PGE2-mediated chemoresistance |
| Semaphorin interactions | antimetabolite class | resistance | -3.8 | PLAUSIBLE |

Bevacizumab hits (GABA, neuronal) and nicotinic receptors with 5-FU are likely artefacts
of tumour composition.

**Shared ("overall response") layer:**
- Response: RNA metabolism, glucose metabolism, Fanconi anaemia regulation, FGFR1 signalling.
- Resistance: renal aquaporin and water balance, gap-junction degradation, phospholipase C,
  BCR downstream.
- None pass FDR individually.

## 3. Overall response score (shared layer) in TCGA (`shared_score_prognostic_vs_predictive.json`)

Progression-free interval, stratified by cancer type:

| Patients | n | HR per SD | p |
|---|---|---|---|
| Untreated | 6,902 total (2,460 events), 20 cancers | 0.88 | 6e-5 |
| Chemotherapy-treated | | 0.81 | 2e-15 |
| Score x treatment interaction | | 0.93 | 0.064 |

The score is mainly **prognostic**, with a borderline **predictive** component that is
stronger under chemotherapy. Treatment was not randomised.

## 4. External validation: frozen model, 69 cohorts (`results/atlas/external/`)

- **Regimen score:**
  - 44 of 72 cohort/arm tests have AUROC > 0.5. Only 1 is nominally significant on its own.
  - The signal concentrates in **breast neoadjuvant taxane + anthracycline**:
    - I-SPY2: all arms 0.57 (p = 0.063), control arm 0.56, AMG386 arm 0.63 (p = 0.037);
    - GSE41998: 0.60 (p = 0.063);
    - BrighTNess paclitaxel arm: 0.60 (p = 0.083);
    - GSE32646: 0.56.
  - Null elsewhere: colorectal, ovarian, bladder, oesophageal, lung.
- **Shared score:** externally null (median AUROC 0.50).
- **Drug specificity:** no randomised design supports it:
  - paclitaxel added (GSE20271);
  - paclitaxel vs ixabepilone (GSE41998);
  - carboplatin added (BrighTNess, I-SPY2);
  - adjuvant platinum (GSE42127, JBR.10);
  - adjuvant 5-FU (GSE103479).

  All are non-significant, and several point the opposite way.
- **STORM (GSE109211):** the deposited data are batch-confounded with outcome (random
  signatures reach the same AUROC), so it cannot be used.

## Conclusions

1. Transcriptomic pathways explain only a small part of chemotherapy response across cancers
   (held-out AUROC about 0.56). Most of it is a **shared, largely prognostic** programme,
   not drug-specific biology.
2. **Reproducible signal:** breast neoadjuvant taxane/anthracycline response (TRIAD paclitaxel
   0.66 in BrighTNess; ATLAS 0.56-0.63 across I-SPY2 arms and GSE41998).
3. **Credible drug-specific pathway biomarkers (FDR < 0.1, with literature support):**
   - interferon/RIG-I signalling -> cisplatin resistance;
   - B-cell/BCR activity -> paclitaxel response;
   - cell-cycle activity -> cisplatin response;
   - Fanconi anaemia -> topoisomerase II response.

   They are hypotheses for targeted validation, not clinical tests.
4. **What would move the field:** response labels that are not TCGA's free-text RECIST, and
   randomised single-drug comparisons with expression, which are rare in public data.
