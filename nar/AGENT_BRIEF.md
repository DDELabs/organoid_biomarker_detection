# Shared brief for every agent working on this project (read fully before starting)

## The project in one paragraph
We start from Kong et al. 2020, Nature Communications, "Network-based machine learning in colorectal and bladder
organoid models predicts anti-cancer drug efficacy in patients" (DOI 10.1038/s41467-020-19313-8). Their method:
(1) find Reactome pathways close to a drug's targets in the STRING protein network (network proximity z-score);
(2) score those pathways in organoids with ssGSEA; (3) train Ridge / SVR / LinearRegression on organoid drug response;
(4) apply the model to TCGA patients treated with the drug, split at the median prediction, compare survival
(Kaplan-Meier + log-rank). They also did multi-pathway models, bootstrapping, isogenic resistant/sensitive cell
lines and a mutation-based concordance check. Goal now: a paper the user can submit to NAR Cancer.

## Where things are (branch claude/zealous-darwin-zeiopw, repo root /home/user/organoid_biomarker_detection)
- Original paper code: `python/`, `utilities/`; user's Python-3 notebooks for COAD/BLCA/LUAD/OVARY: `modified_codes_py3/`.
- Earlier work (READ FIRST): `docs/PROJECT_SUMMARY.md`, `robust/README.md`, `robust/results/SUMMARY.md`,
  `robust/results/external_validation/SUMMARY.md`, `docs/revamp/AUDIT_AND_POWER.md`, `docs/revamp/REVAMP_PLAN.md`,
  `robust/BRIDGE.md`, `robust/NEXTGEN.md`, `robust/TRIAD.md`, `robust/ATLAS.md`.
  Reusable code lives in `robust/` (loaders, STRING proximity, ssGSEA/rank scoring, survival tests). REUSE it; do not rewrite.
- Curated lab-model drug data: `data/curated/preclinical/` (CATALOG_PRECLINICAL.tsv: colorectal_vandewetering2015,
  bladder_lee2018, pancreas_tiriac2018, pancreas_shi2022, liver_broutier2017, liver_ji2023, sarcoma_alshihabi2024,
  colorectal_pdx_isella2017, pdx_novartis_gao2015, ctrpv2_ccle, prism_repurposing). Ovarian organoids: see
  `robust/data_builders/provenance/organoids_ovarian_vias2023.md`.
- Curated patient cohorts: `data/curated/trials/` (CATALOG_GEO.tsv, 29 cohorts), combined table
  `robust/results/training_table/` (20,766 samples), GSE39582 provenance in `robust/data_builders/provenance/`.
- Paper's own source data (STRING edges, proximity z-scores, organoid & TCGA NES for COAD and BLCA, bootstraps):
  download `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19313-8/MediaObjects/41467_2020_19313_MOESM3_ESM.xlsx`
  (also MOESM1_ESM.pdf = supplementary figures, MOESM2 = reporting summary). Cache under `nar/_cache/`.
- Reachable from this machine: GEO/NCBI FTP, GDC API (TCGA, open clinical + expression), cog.sanger.ac.uk (GDSC1/2,
  Cell Model Passports incl. rnaseq_all_20220624.zip, Project Score CRISPR v1 file
  Project_score_combined_Sanger_v1_Broad_20Q2_20210311.zip), ndownloader.figshare.com (DepMap files if you know the
  file id), STRING downloads, Reactome, Zenodo, Synapse, static-content.springer.com, PyPI.
  BLOCKED: cell.com/elsevier, AACR, nature.com pages, science.org, PMC/EuropePMC, PharmacoDB, DrugBank, GitHub,
  api.figshare.com, pytorch.org (torch from PyPI works, CPU only). No GPU. 4 CPUs, 15 GB RAM — be frugal: other
  agents run at the same time; use at most 2 worker processes.
- Python: pandas, numpy, scipy, sklearn, statsmodels, gseapy, networkx, matplotlib, xgboost installed.
  `lifelines` is NOT installable (use statsmodels: `statsmodels.duration.survfunc.survdiff` for log-rank,
  `statsmodels.duration.hazard_regression.PHReg` for Cox; KM curves via `SurvfuncRight`). Check whether `robust/`
  already has helpers first.

## Facts already established in this conversation (do not rediscover)
- TCGA treated counts: BLCA cisplatin 77, methotrexate 12; LUAD cisplatin 87, carboplatin 77, pemetrexed 60,
  gemcitabine 25; OV carboplatin 464, paclitaxel 448, cisplatin 156, olaparib 0; COAD: of 154 treated, most got
  FOLFOX (5-FU+oxaliplatin+leucovorin), only ~5 got 5-FU alone.
- User's LUAD attempt used GSE165611 = lung cancer patient-derived CELLS (not organoids), 102 samples, 30 drugs;
  its drug table (`supp_table_S2b.txt`) is NOT in the repo. User's OVARY attempt used GSE229119 = olaparib
  time-course in 3 organoids + cell lines (not a response panel). User's BLCA attempt used methotrexate (12 patients).
- GSE39582: 5FU 82, FUFOL 54, FOLFOX 23, FOLFIRI 12, untreated 326, with stage/MMR/KRAS/BRAF/OS/RFS.

## RULES (the user insisted on these)
1. Do EVERYTHING in your assignment completely: every figure, every supplementary figure, every table. No
   "partial, will do later". If something truly cannot be done (data blocked), say exactly what and why in one line
   and do the closest possible substitute.
2. Do NOT frame results as "negative", "failed", "not significant". Report the actual numbers (p-value, hazard ratio
   with 95% CI, n, number of events, correlation) and add a short comment "What I think:" — then let the user decide.
3. Write in plain, simple, human English. Short sentences. No jargon without a one-line explanation in brackets.
   No jumbled sentences. Explain each figure as if to a smart non-specialist.
4. Figures: matplotlib, 300 dpi PNG + PDF, readable fonts (>=8 pt), each with a caption file
   (`<figure>.caption.md`) written in plain language. Colour-blind-safe colours. Panel letters (a, b, c...).
5. Never fabricate data or numbers. Every number in a report must come from a file your code wrote.
6. Save code in your own folder (`nar/<your_folder>/code/`), outputs in `nar/<your_folder>/figures/`,
   `nar/<your_folder>/tables/`, plus `nar/<your_folder>/REPORT.md`. Make scripts re-runnable.
7. Do NOT run git commit/push (the coordinator does that). Do not edit files outside your own folder except
   reading. Do not delete anything.
8. Set random seeds. Record versions of data files used (file name + size or sha256) in your REPORT.
9. Your final message back to the coordinator: a concise summary (<= 400 words): what was done, key numbers,
   list of figures/tables produced, anything you could not do and why.
