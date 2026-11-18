# Project summary (state on 2026-11-16)

Goal: predict which cancer drug suits which patient and name the biological pathway responsible, starting from
Kong et al. 2020 (Nature Communications), "Network-based machine learning in colorectal and bladder organoid
models predicts anti-cancer drug efficacy in patients".

## 1. The paper's method

1. Protein network: find pathways that sit close to the drug's targets in the network (Barabasi-group "proximity"). This gives a shortlist.
2. Organoids: measure pathway activity in lab-grown tumours and train a model of how strongly the drug kills them.
3. Pick the single pathway the model weights most. This is the biomarker.
4. Patients: split treated TCGA patients at the median of that pathway and compare survival.

Headline results in the paper: colon cancer with 5-FU (pathway "activation of BH3-only proteins", p = 0.014, 114 patients)
and bladder cancer with cisplatin (pathway "amino acid synthesis and interconversion", p = 0.01, 77 patients).

## 2. What was done, in order

| Step | What | Result |
|---|---|---|
| 1 | Earlier session: five strategies (robust reproduction, NIT, TRIAD, ATLAS, BRIDGE) | all reported negative |
| 2 | Patient-trained survival model on 10,032 TCGA patients, all drugs, with doctors' response notes | strong general risk score; no drug-specific signal (best corrected p 0.135) |
| 3 | Integrated evidence tiers (network, cell lines, organoids, notes, survival) with relaxed criteria | 59 top-tier candidates against about 45 expected by chance |
| 4 | Reruns with all pathways, and with the paper's own shortlists | no gain either way |
| 5 | Exact reproduction of the paper's headline analysis | reproduces on the authors' data (p = 0.014; untreated 0.16) |
| 6 | Same analysis on today's TCGA data | colon 5-FU p = 0.866 (129 patients, 25 deaths; source `nar/02_organoid_only/tables/cross_set_summary.tsv`, row COAD current; the earlier `robust/` run gave 0.62 and is superseded, audit I11); bladder cisplatin p = 0.067 (79 patients, 28 deaths; `nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv`; earlier run 0.063), right direction |
| 7 | Paper's method on every cancer and drug with enough patients | 152 tests, 11 significant, about 8-9 expected by chance; 0 of 17 organoid-trained tests significant |
| 8 | Survey of public data and a combined training table | 20,766 tumour samples from TCGA and 41 other cohorts |
| 9 | Models trained on the combined table, tested on cohorts never seen | breast response predicted in 13 of 13 cohorts (0.67-0.85); other cancers 1 of 20 |
| 10 | Drug-specific tests in two randomised breast trials | no drug significant; carboplatin was the one lead |
| 11 | Carboplatin follow-up with the network shortlist | one strong transfer result between trials (0.80 against 0.58; 0 of 200 random shortlists as good) |
| 12 | Frozen carboplatin score tested in ovarian cancer | failed: no link to survival or response |

## 3. Verdict

- Holds: the paper's headline number reproduces on its own data; breast chemotherapy response can be predicted
  on unseen patient groups and adds to hormone-receptor and HER2 status.
- Does not hold: the paper's result on today's data; the method across other cancers and drugs; any
  drug-specific biomarker in patients, including the carboplatin lead.
- Main limit: patients with a known drug, a comparison group and an outcome. Breast has about 3,300 with a
  recorded response; no other cancer has more than about 720.

## 4. Errors made and corrected during the work

- The paper's secondary (multi-pathway) analysis was first tested as if it were the main one; corrected in step 5.
- Three trial datasets (VESPER, FOCUS, NeoTRIP) were announced as usable before their records were opened;
  they hold no public outcomes.
- CTR-DB was first described as blocked by the allow-list; the real cause is that the site moved and its
  secure connection is misconfigured. It was not obtained.

## 5. Where things are

| Path | Content |
|---|---|
| `robust/results/paper_single_pathway/` | exact reproduction of the paper's headline analysis |
| `robust/results/single_pathway_all/` | paper's method on all cancers and drugs |
| `robust/results/survlearn/`, `robust/results/integrated/` | patient-trained survival model; evidence tiers |
| `robust/results/training_table/` | combined training table (start with its README) |
| `robust/results/train_combined/` | held-out training results |
| `robust/results/carboplatin/` | carboplatin follow-up and the ovarian test |
| `docs/revamp/TRAINING_DATA_SURVEY.md` | data survey |
