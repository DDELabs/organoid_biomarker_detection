# Pipeline strategy: from the original paper to the pan-cancer atlas

This page explains the whole pipeline as it now stands and how it differs from the
original organoid → network → patient pipeline (Kong et al. 2020).

## 1. The original pipeline (paper, `python/`, `utilities/`)

| Step | What it does |
|---|---|
| Input | Colorectal (COAD) organoids: transcriptome plus 5-FU IC50; TCGA COAD patients: HTSeq FPKM-UQ expression, clinical data and drug records |
| Network | Genes close to the drug's targets on STRING (closest-distance proximity z, degree-matched null, one cut-off) |
| Pathways | ssGSEA on the network-proximal genes, z-scored within cohort |
| Model | Ridge / SVR / OLS fitted once on organoid IC50; top-k pathways by \|coef\| |
| Patient test | Score the 5-FU-treated TCGA patients; median split; log-rank on survival |
| Output | One drug, one cancer; best of 3 models × 9 k values reported (p = 0.076) |

Problems we found:
- **The data no longer exist.** HTSeq FPKM-UQ was replaced by STAR counts (GENCODE v36).
- **Selection inflates significance.** The reported best p = 0.076 has a whole-pipeline
  permutation p of 0.48.
- **Prognostic and predictive effects are confused.** No untreated comparison.
- **Organoid-only training does not transfer** to some drugs. Liver sorafenib was reversed
  by proliferation and general-sensitivity confounding.
- **Coverage is narrow.** One cancer, one drug.

## 2. What we built, layer by layer

| Layer | Module | Main idea | Status |
|---|---|---|---|
| **Robust** reproduction | `robust/obd`, `run_coad_5fu.py` | Exact reproduction of the paper; current GDC STAR FPKM-UQ; unit-invariant rank scores; stability selection; LOOCV; Cox, C-index and treatment interaction; whole-pipeline permutation; robustness scorecard | Paper result shown to be non-robust |
| **NIT** (next-gen) | `nextgen.py`, `run_nextgen.py` | Deconfound organoid response for proliferation (Hallmark E2F/G2M) and general drug sensitivity (PC1 of other drugs); soft network weights; GDSC transfer prior | Fixed the sorafenib reversal; still small-n |
| **TRIAD** | `triad.py`, `run_triad*.py` | Patient-anchored: logistic model on TCGA RECIST labels pooled pan-cancer, centred on a prior from organoid/GDSC/network | Paclitaxel validated in BrighTNess (AUC 0.66) and I-SPY2 |
| **ATLAS** | `atlas.py`, `run_atlas*.py` | Multi-task model over 21 drugs: pathway effects split into a shared layer ("overall response"), a mechanism-class layer and a drug layer; ~2,800 patient × drug labels; prior-free bootstrap inference with FDR; external validation in 69 trial cohorts | Held-out AUROC 0.56; credible FDR hits (IFN → cisplatin resistance, BCR → paclitaxel response, cell cycle → cisplatin response, Fanconi → topo-II response) |
| **Per-cancer lists** | `run_atlas_cancer_specific.py`, `run_atlas_relaxed.py` | Top 10 per cancer × drug: local within-cancer evidence plus the ATLAS pan-cancer prior; bootstrap stability; relaxed tiers (RECIST n ≥ 10, first-course outcome, PFI among treated with an untreated contrast) | 99 pairs, 22 cancers; few pass FDR locally |

## 3. Key differences from the paper

1. **Training signal:** patient response (TCGA RECIST, first-course outcome, PFI), with
   organoids/GDSC/network as *priors*. The paper used organoid IC50 alone.
2. **Scale:** many drugs and cancers jointly (multi-task, partial pooling), instead of one
   drug × one cancer.
3. **Features:** unit-invariant single-sample rank scores over all 673 Reactome pathways,
   centred within cancer type. The paper used ssGSEA on network-filtered genes only.
4. **Confounding:** treatment setting and proliferation as covariates; prognostic vs
   predictive split (untreated contrast, score × treatment interaction).
5. **Statistics:** permutation nulls of the whole pipeline, bootstrap z, BH-FDR,
   leave-one-cancer-out prediction, random-signature nulls. The paper reported a nominal
   best-of-27 p.
6. **External validation:** 25 GEO / trial cohorts with frozen models, plus
   drug-specificity tests in randomised designs.
7. **Data access:** Xena / GDC mirrors, curated GEO and pre-clinical sets in
   `data/curated` (reproducible builders), EGA application drafts in `docs/data_access`.

## 4. What the evidence supports (relaxed stringency)

- Transcriptomic pathways explain a small part of chemotherapy response (AUROC about 0.56).
  Most of it is a shared, largely prognostic programme.
- Reproducible: breast neoadjuvant taxane/anthracycline response.
- Credible drug-specific hypotheses:
  - interferon/RIG-I → cisplatin resistance;
  - B-cell/BCR → paclitaxel response;
  - cell cycle → cisplatin response;
  - Fanconi anaemia → anthracycline/topo-II response;
  - prostanoid receptors → taxane resistance;
  - FGFR3 → cisplatin response (bladder, head & neck).
- From the relaxed tiers, nominal hits that need validation:
  - NF-κB (RIP/TRAF6) activity → worse rituximab-regimen outcome in DLBCL (consistent with
    ABC-type biology);
  - circadian-clock genes → dacarbazine resistance in melanoma.

## 5. Next step

See `docs/BRIDGE_strategy_research.md` and `robust/BRIDGE.md`: train on patient-derived
models (PDX/PDO with transcriptome and response), align them to patient tumours, and
validate on patient cohorts.
