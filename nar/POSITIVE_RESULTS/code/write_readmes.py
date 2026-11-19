"""Writes README.md per cancer folder, top README.md, INDEX.md and BORDERLINE.md.
Tables inside READMEs are built from the same source tables the index uses; narrative numbers were checked against
POSITIVE_INDEX.tsv (printed by collect.py). Run after collect.py and copy_files.py."""
import pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests
RULE = ('held-out/independent test, p < 0.05 in the pre-stated direction, or BH q < 0.10 in its own family; n and events reported')
ROOT = '/home/user/organoid_biomarker_detection'; OUT = f'{ROOT}/nar/POSITIVE_RESULTS'
T = lambda p: pd.read_csv(f'{ROOT}/{p}', sep='\t')
f3 = lambda x: f'{x:.3f}'
def fp(p):
    return f'{p:.1e}' if p < 0.001 else f'{p:.3f}'

# ---------- breast tables
bl = T('nar/19_model_benchmark/tables/breast_loco_all_models.tsv'); l2 = bl[bl.model == 'L2_fixed'].sort_values('test_cohort')
t1 = '| cohort | n | responders | AUC (95% CI) | p vs 0.5 |\n|---|---|---|---|---|\n' + '\n'.join(
    f'| {r.test_cohort} | {int(r.n)} | {int(r.responders)} | {f3(r.auc)} ({f3(r.ci_low)}-{f3(r.ci_high)}) | {fp(r.p_vs_chance)} |' for r in l2.itertuples())
v2 = T('nar/18_other_cohorts/tables/part1_v2.tsv'); v2 = v2[v2.table == 'v2_all_groups']
old = set(l2.test_cohort)
nb = v2[(v2.cancer == 'breast') & (~v2.test_cohort.isin(old))].sort_values('p').copy(); nb['q'] = multipletests(nb.p, method='fdr_bh')[1]
t2 = '| cohort | n | responders | AUC (95% CI) | p | q (BH over 15) |\n|---|---|---|---|---|---|\n' + '\n'.join(
    f'| {r.test_cohort} | {int(r.n)} | {int(r.responders)} | {f3(r.auc_pathways)} ({f3(r.auc_ci_low)}-{f3(r.auc_ci_high)}) | {fp(r.p)} | {fp(r.q)} |' for r in nb.itertuples())
mc = T('nar/19_model_benchmark/tables/model_comparison.tsv').sort_values('mean_auc', ascending=False)
t3 = '| model | mean AUC (95% CI), 13 cohorts | difference vs fixed L2 (p, q) |\n|---|---|---|\n' + '\n'.join(
    f'| {r.model} | {f3(r.mean_auc)} ({f3(r.mean_ci_low)}-{f3(r.mean_ci_high)}) | ' + ('reference' if r.model == 'L2_fixed' else f'{r.mean_diff_vs_L2:+.3f} (p {fp(r.p_mean_diff)}, q {fp(r.q_mean_diff_BH)})') + ' |' for r in mc.itertuples())
ak = T('nar/20_breast_pathways/tables/auc_vs_k.tsv')
def akv(rk, k): return ak[(ak.ranking == rk) & (ak.k == k)].mean_auc.iloc[0]
t4 = '| k pathways | chosen on training cohorts by one-at-a-time strength | chosen by L2 weight | random sets |\n|---|---|---|---|\n' + '\n'.join(
    f'| {k} | {f3(akv("univariate",k))} | {f3(akv("weight",k))} | {f3(akv("random",k))} |' for k in [1, 5, 20, 100, 200])

open(f'{OUT}/BREAST/README.md', 'w').write(f"""# Breast cancer: chemotherapy response is predicted in cohorts the model has never seen

## What was tested
Can a model that reads 673 Reactome pathway scores from a breast tumour predict whether the patient responded to chemotherapy (responder = 1 as coded in each cohort; the response definitions differ between cohorts)? The model is trained on all other breast cohorts and tested on one cohort it never saw (leave-one-cohort-out). All choices (penalty, pathway selection) are made on the training cohorts only.

## Data
- 13 original cohorts, **3,511 treated patients, 1,202 responders** (BRIGHTNESS, GSE123845, GSE163882, GSE20194, GSE20271, GSE22093, GSE23988, GSE240671, GSE25066, GSE32646, GSE41998, ISPY2, TCGA-BRCA).
- 15 new cohorts added later (W2): 1,327 patients, 416 responders.
- Pathway scores are z-scored inside each cohort. Model: L2-penalised logistic regression, penalty fixed at 100 (the original model). 17 other models were tried in W3.

## Result 1. Original 13 cohorts (rule met in 13 of 13)
Mean AUC 0.752 (95% CI 0.728-0.775). Per cohort (fixed L2, DeLong CI, p vs 0.5):

{t1}

The p column here is the DeLong test of the W3 re-run (identical AUCs). The original run in `robust/results/train_combined` uses a one-sided Mann-Whitney test; its largest p is 0.016 (TCGA-BRCA). BH q < 0.10 in 13 of 13 (`POSITIVE_INDEX.tsv` BR1). Sources: `tables/robust__part1_response_leave_one_cohort_out.tsv` (original run), `tables/19__breast_loco_all_models.tsv` (with CIs).

## Result 2. 15 new cohorts (rule met in 8 of 15 by p, 9 of 15 by BH q)
Trained on all breast cohorts of the v2 table except the held-out one:

{t2}

What I think: the new cohorts that did not reach p < 0.05 are small (21-114 patients, 7-17 responders) or have a weak AUC (GSE59515 0.471, GSE66305 0.613). The mean over all 15 is 0.682, lower than the original 13 (0.752). So the effect is real but not uniform. Source: `tables/18__part1_v2.tsv` (table `v2_all_groups`).

## Result 3. All-model benchmark (W3): the signal does not depend on the model type
Same 13 held-out cohorts, 17 model types, settings tuned only inside the training cohorts.

{t3}

- Baselines: ER+HER2 flags alone 0.693 (9 cohorts); proliferation score alone 0.646. The pathway model beats both (paired p < 0.001).
- 200 random sets of 183 pathways (the size of the network shortlist) average 0.742; the network shortlist gives 0.744. The shortlist adds nothing over random pathways.
- The best model (linear SVM, 0.764) is +0.012 over fixed L2 with p 0.023 but BH q 0.100, so it is a small and not clear gain. Source: `tables/19__model_comparison.tsv`, `tables/19__breast_loco_all_models.tsv`, `tables/19__random_pathway_sets_auc.tsv`; figures `19__fig_model_mean_auc`, `19__fig_model_heatmap`, `19__fig_roc_top3`.

## Result 4. Which pathways drive it (W4)
{t4}

- One pathway chosen on the training cohorts (EARLY_PHASE_OF_HIV_LIFE_CYCLE in 13 of 13 folds) already gives 0.680; five pathways chosen greedily plateau at about 0.72; about 100 or more are needed to be within 0.01 of the full model. Random sets of 100-200 pathways reach 0.730-0.743, so the information is spread over many correlated pathways.
- The top 30 weights are stable (28 of 30 selected in at least 80% of 100 bootstrap L1 fits) but are generic cell-machinery pathways (several are Reactome virus or drug pathways); weights of correlated pathways should not be read one by one.
- Not simply proliferation: proliferation alone gives 0.646; with proliferation regressed out of every pathway the model still gives 0.721 (0.694-0.747); dropping the 27 proliferation pathways gives 0.748. The model adds to proliferation + ER + HER2 in 11 of 13 cohorts (likelihood-ratio p < 0.05; Stouffer z 14.0).
- Subgroups (cohorts with enough responders and non-responders; receptor flags missing for GSE22093, GSE23988, TCGA-BRCA): ER+/HER2- AUC 0.670 (0.599-0.741), 1,111 patients, 8 cohorts; HER2+ 0.728 (0.669-0.788), 520 patients; triple-negative 0.720 (0.678-0.761), 1,470 patients. All three intervals exclude 0.5. In ER+/HER2- proliferation alone gives 0.627, so more of the signal there is proliferation.
Sources: `tables/20__auc_vs_k.tsv`, `20__top30_pathways.tsv`, `20__stability_selection.tsv`, `20__proliferation_tests.tsv`, `20__subgroup_auc_pooled.tsv`, `20__minimal_subset.tsv`, `20__family_weights.tsv`; figures `20__fig_*`.

## Caveats
- Cohorts mix drugs and response definitions (pCR, clinical response, RECIST-like). The test is "predict response in a new cohort", not "predict response to one named drug".
- TCGA-BRCA has 184 responders of 202, so its AUC (0.672) rests on only 18 non-responders.
- The network-proximity shortlist is from the paper's table, not a new STRING v12 run.
- Randomised HER2 trials (NOAH, CHER-LOB) do not show a drug-specific effect and are in `../BORDERLINE.md`, not here.

## What I think
This is the clearest positive in the project. A pathway model learned in other breast cohorts predicts chemotherapy response in 13 of 13 cohorts it never saw, and in 8 of 15 further cohorts, with AUC around 0.7-0.8 in the larger ones. It is a general chemotherapy-response signal that is not just proliferation or receptor status. It is not a drug-specific signal, and it does not need the network-proximity step: a random 183-pathway set does as well.

## Where things are
`figures/` (PNG + PDF + caption), `tables/` (copies of the source tables), `source_reports/` (copies of the agents' reports). Copy list with checksums: `../COPY_MANIFEST.tsv`.
""")

# ---------- lung
open(f'{OUT}/LUNG/README.md', 'w').write("""# Lung: a pathway score that predicts survival in people who got no chemotherapy

## What was tested
Lung cohorts were added in W1 (20 builders; 3,878 new samples; table v2). The question was whether pathway scores learned in lung chemo/no-chemo cohorts predict survival (a) in patients who did not get chemotherapy (general prognosis) and (b) as a chemo-benefit score. The prognosis score is trained only on the other cohorts, never on the cohort it is tested in.

## Data
14 tested cohorts of untreated or surgery-only patients (GSE14814 observation arm, GSE37745, GSE42127, GSE68465 no-chemo arm, GSE157011, GSE162520, GSE30219, GSE31210, GSE3141, GSE41271, GSE50081, GSE72094, GSE81089, GSE8894); **2,938 patients, 1,320 events**. GSE29013 (6 events) was not tested. Endpoint overall survival (GSE8894: relapse-free). Score: ridge Cox on 673 pathways, adjusted for age (and stage where known).

## Result 1 (strong). General prognostic score in untreated patients
| score | pathways | pooled HR per SD (95% CI) | p |
|---|---|---|---|
| general (main-effect) prognostic score | all 673 | **1.36 (1.29-1.44)** | 1.1e-28 |
| general prognostic score | 47 cisplatin shortlist | 1.21 (1.14-1.28) | 2.0e-11 |
| chemo-benefit score | all 673 | 0.98 (0.93-1.04) | 0.47 |
| chemo-benefit score | shortlist | 0.99 (0.94-1.05) | 0.75 |
- HR above 1 in 14 of 14 cohorts; p < 0.05 in 9 of 14 (GSE157011 1.20, GSE30219 1.49, GSE31210 1.80, GSE41271 1.72, GSE42127 1.37, GSE68465 1.38 (n=340, 162 deaths, 1.18-1.61), GSE72094 1.41, GSE81089 1.46, GSE8894 1.41). The observation arm of the randomised GSE14814 gives 1.06 (0.69-1.63), 30 deaths.
- Source: `tables/17__lung_untreated_prognosis.tsv`; figure `17__lung_untreated_prognosis`.

## Result 2 (strong). All-model benchmark on the same data (W3)
The original ridge Cox gives pooled C-index 0.601 (0.585-0.619). Network-smoothed L2 gives 0.622 (0.606-0.638), +0.021 (+0.009 to +0.033), paired p 0.010, q 0.020; plain ridge 0.620. The shortlist version is worse, and random shortlist-size sets reach the original in 44% of draws. Source: `tables/19__positives_all_models_lung.tsv`; figure `19__fig_bar_lung_untreated_prognosis`.

## Result 3 (nominal only). Two lab-model-trained tests in the 57-test survey
| model -> TCGA-LUAD patients | n / events | HR (95% CI) | log-rank p | BH q over 57 | shuffle p |
|---|---|---|---|---|---|
| lung adenocarcinoma cell lines, gemcitabine | 25 / 16 | 3.17 (1.06-9.44) | 0.030 | 0.339 | 0.138 |
| NSCLC PDX, erlotinib | 21 / 7 | 7.61 (0.81-71.0) | 0.040 | 0.354 | 0.045 |
Both meet the p < 0.05 arm of the rule but not the q arm, and the intervals are wide. Source: `tables/05__A3_master_results.tsv`; figure `05__A3_Fig1_master_forest`.

## Caveats
- The score is prognostic: it separates survival in untreated patients, but the chemo-specific part does not (HR 0.98), and in the randomised JBR.10-type trial GSE14814 the interaction is HR 1.02 (p 0.95). It is not evidence of chemo benefit prediction.
- Stage was approximated for GSE68465; GSE161537 and GSE162520 are targeted panels (185 of 673 pathways cannot be scored).
- Immunotherapy response (GSE218989 AUC 0.566 p 0.031, GSE161537 and GSE126044 with pan-cancer models) and the GSE29013 interaction are in `../BORDERLINE.md`.

## What I think
Lung gives one solid result: pathway scores carry prognostic information across 14 independent cohorts, and several models do slightly better than the original. It is a prognosis result, not a drug-response result. The two lab-model tests are weak leads.

## Where things are
`figures/`, `tables/`, `source_reports/17__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
""")

# ---------- bladder
open(f'{OUT}/BLADDER/README.md', 'w').write("""# Bladder: cisplatin organoid model, ERCC2 mutation check, isogenic cells and the joint model

## What was tested
The paper's bladder cisplatin result (Kong et al. 2020, Fig 2e-h, 4c, 5b) reproduced with the paper's own inputs and with current data; then extra tests. Method: pathways close to cisplatin targets in the STRING network; ssGSEA scores; ridge on 9 organoid responses (7 patients; two lines are recurrences of the same patients; audit I07); patients scored with the top pathway and split at the median; log-rank on overall survival.

## Data
Lee 2018 bladder organoids (9 with cisplatin response; response = refit AUC, the paper's IC50 is not public); TCGA-BLCA cisplatin-treated 79 patients (20 deaths on paper data; 28 on current GDC data); 288 untreated; GSE31684, GSE48276 etc. for the joint model.

## Results
| id | test | n / events | effect (95% CI) | p | tier |
|---|---|---|---|---|---|
| BL1 | paper-data rank-1 pathway, treated TCGA-BLCA | 79 / 20 | HR 2.63 (0.99-6.96) | 0.044 | nominal (BH q 0.354 over 57; shuffle p 0.30) |
| BL1 | paper's amino-acid pathway, same patients | 79 / 20 | HR 2.66 (1.06-6.69) | 0.031 | nominal |
| BL2 | current GDC data, response-shuffle test | 79 / 28 | log-rank p 0.067, HR 2.03 (0.94-4.38) | shuffle p **0.036** | nominal (the patient log-rank itself is not < 0.05) |
| BL3 | ERCC2-mutant vs wild type, paper data, amino-acid pathway | 36 / 367 | one-sided Mann-Whitney | **0.0030** | strong |
| BL3 | same, current data | 36 / 366 | | **0.00055** | strong |
| BL4 | isogenic T24 sensitive vs resistant, amino-acid pathway | 3 vs 3 | t 12.7 | **0.00022** | strong (tiny n) |
| BL5 | joint organoid+patient model, 4 held-out cohorts | 157 / 99 | HR 1.24 (1.03-1.50) per SD risk | 0.024 | nominal |
| BL5 | patient-only model, same cohorts | 157 / 99 | HR 1.45 (1.18-1.79) | 5.1e-04 | strong |
| BL6 | joint model, gemcitabine, TCGA-BLCA 5-fold | 64 / 30 | HR 1.53 (1.09-2.14) | 0.014 | nominal |
| BL7 | joint model, GSE48276 MVAC x score interaction | 53 (8 treated events) | HR 5.30 (1.69-16.66) | 0.0043 | nominal (1 of 32) |
| BL8 | current-data rank-1 signature, GSE31684 interaction | 35 treated / 25 events | HR 1.82 (1.07-3.09) | 0.028 | nominal (q 0.650 over 196) |
| BL9 | ERCC2, all-model benchmark, paper amino-acid pathway AUC | 397 / 36 | AUC 0.665 (0.571-0.752) | 5.4e-04 | strong |
Sources: `tables/05__A3_master_results.tsv`, `tables/02__ERCC2_mutation_concordance.tsv`, `tables/02__isogenic_T24_tests.tsv`, `tables/02__paper_biomarker_amino_acid_KM.tsv`, `tables/03__summary_main.tsv`, `tables/03__interaction_all.tsv`, `tables/05__B1_external_results.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/13__G8_ERCC2_null_and_adjustment.tsv`.

## How robust is it
- ERCC2 is the sturdy one. Share of all pathways (both signs) as strong as the amino-acid pathway: 0.017 (paper data), 0.002 (current data); share of random proximal-size sets: 0.063 and 0.000. It holds after tumour-mutation-burden adjustment (p 0.024 paper, 0.008 current; OR per SD 0.63 and 0.54). The paper's rank-1 pathway does not carry it (p 0.88 / 0.33).
- Survival: with 20 deaths a p of 0.04 is easy to reach by chance (29.6% of response-shuffled runs also reach p < 0.05 on paper data). One line per patient (7 lines) leaves the treated statistics unchanged (p 0.044) but the shuffle p moves to 0.263. The amino-acid pathway is rank 6 in our refit (rank 30 of 47 on current data, treated p 0.063, untreated p 0.021).
- Isogenic T24: only the paper's NES for 3 vs 3 samples exist, so p 0.00022 rests on the paper's processed values. Our rank-1 pathway (p 0.017) points the wrong way there.
- Joint model: the organoid-only model is not significant (HR 1.07, p 0.45; gemcitabine 0.90, p 0.59). The gain comes from the patient part (patient-only 1.45); joint minus patient-only is not above zero in any of 20 unit x outcome comparisons (nar/14 G16). Fold-aligned joint HR 1.27, p 0.013.
- All-model benchmark: no added model beats the paper's amino-acid pathway for ERCC2 (best MLP 0.637); a simple proliferation score gets 0.634, so this signal is not specific to the network method.

## Caveats
Nine organoids from seven patients; AUC instead of IC50; TCGA treatment records are "ever given"; untreated groups differ in age and stage; all survival tests have 20-28 events. The five bladder cohorts' paper-signature tests lean the right way (OR/HR 1.16-1.57) but none is below 0.05.

## What I think
ERCC2 concordance is the strongest bladder result and it is specific. The cisplatin survival split is at the edge: p 0.03-0.07 on 20-28 deaths, helped by the choice of pathway. The joint model is a real held-out signal in bladder, but it is carried by patient data.

## Where things are
`figures/`, `tables/`, `source_reports/03__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
""")

# ---------- colorectal
open(f'{OUT}/COLORECTAL/README.md', 'w').write("""# Colon: the paper's 5-FU result, isogenic cells, cetuximab vs BRAF V600E, and the response cohort GSE62080

## What was tested
The Kong et al. colon results re-run from the paper's own files and from current data: 5-FU model (19 van de Wetering organoids; rank-1 pathway "Activation of BH3-only proteins"), isogenic HCT-8 cells, cetuximab model vs BRAF V600E tumours; plus external colorectal cohorts.

## Results
| id | test | n | effect | p | tier |
|---|---|---|---|---|---|
| CR1 | 5-FU model, 5-FU-treated TCGA-COAD, authors' 2019 files, median split | 114 patients, 8 deaths (0 vs 8) | HR not estimable; continuous HR 1.52 (0.94-2.46) | **0.014** | nominal (BH q 0.273 over 57; shuffle p 0.034; bootstrap share 0.0072) |
| CR2 | isogenic HCT-8 sensitive vs 5-FU-resistant, BH3-only higher in sensitive | 3 vs 9 samples | t 3.27 | **0.0085** | strong direction, tiny n (recomputed from raw GSE81005: p 0.039) |
| CR3 | cetuximab model vs BRAF V600E, authors' data | 49 mutant / 369 wild type | one-sided Mann-Whitney | 0.031 | nominal |
| CR3 | same, current data (colon) | 50 / 377 | | **0.00072** | strong |
| CR4 | same, colon + rectum, current data | 51 / 522 | | **0.00048** | strong |
| CR5 | benchmark: L1 on shortlist vs original rank-1, AUC | 573 / 51 V600E | 0.704 vs 0.640, diff +0.064 (0.011-0.113) | paired 0.018, q 0.036 | q < 0.10 |
| CR6 | response, GSE62080 (metastatic FOLFIRI), trained on other colorectal cohorts | 21 / 9 responders | AUC 0.806 (0.605-1.000) | 0.021, q 0.056 | nominal, tiny n |
| CR7 | Gao CRC PDX 5-FU top-10 signature in GSE39582, OS, x 5-FU interaction | 139 / 40 | HR 1.46 (1.01-2.10) | 0.042 | nominal (q 0.650 over 196) |
| CR8 | organoid irinotecan top-10 signature, GSE72970 FOLFIRI vs FOLFOX PFS interaction | 84 / 80 | HR 1.65 (1.02-2.67) | 0.042 | nominal (q 0.650) |
Sources: `tables/02__cross_set_summary.tsv`, `tables/02__Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv`, `tables/02__Fig5a_cetuximab_BRAF_V600E.tsv`, `tables/10__Fig5_mutation_concordance_other_cancers.tsv`, `tables/13__G8_BRAF_null_and_adjustment.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/18__part1_v2.tsv`, `tables/05__B1_external_results.tsv`, `tables/05__A3_master_results.tsv`.

## How robust is it
- **5-FU (CR1) is fragile.** It reproduces exactly on the authors' 2019 inputs, but p moves to 0.866 on current GDC data (n 129, 25 deaths): the main drivers are the new organoid scores, the new patient expression pipeline and longer follow-up (nar/13 G2). With the 2019 setup but current survival the frozen BH3-only pathway gives p 0.002 (114 patients, 20 deaths), so the direction is not always lost. With the rank score instead of ssGSEA p becomes 0.32 (nar/14 G18). 94 of the 114 "5-FU" patients also received oxaliplatin; the split is almost entirely a FOLFOX group. In independent 5-FU cohorts the BH3-only pathway does not carry over (GSE39582 RFS HR 0.95, p 0.70). No cell-line support: BH3-only vs 5-FU in colorectal lines rho 0.05 (GDSC1, n 44), 0.02 (GDSC2), 0.26 (CTRPv2), 0.41 (PRISM).
- **Cetuximab/BRAF (CR3-CR5) is strong but not specific.** After adjusting for MSI, side, stage and TMB p is 0.026 (current) / 0.093 (authors); within MSS tumours p 0.007 (10 mutants). 17.7% of random pathway sets of the same size do as well (current data), so it is not special to the network. Refits leaving k organoids out still give p < 0.05 in 100% (k=1), 93%, 94%, 87% (k=4). KRAS mutations are not picked up (p 0.88).
- **Isogenic HCT-8 (CR2)**: same direction in both versions, 3 vs 9 samples.
- The paper's colon 5-FU result ranks q 0.27 across the 57 organoid-only analyses; none of the 57 passes BH.

## What I think
Two things hold up: the cetuximab/BRAF concordance (strong p, but partly a tumour-subtype signal) and the isogenic cell direction. The headline colon 5-FU survival split reproduces on the paper's data but disappears on current data, so I would present it as a reproduction with that limit stated. GSE62080 has 21 patients; treat it as a lead.

## Where things are
`figures/`, `tables/`. Checksums in `../COPY_MANIFEST.tsv`. Rectal cancer is in `../RECTAL/`.
""")

# ---------- rectal
open(f'{OUT}/RECTAL/README.md', 'w').write("""# Rectal cancer (GSE87211): 5-FU organoid signatures predict disease-free survival

## What was tested
Frozen pathway signatures built from organoid or PDX 5-FU data (top-10 pathways by absolute ridge coefficient, nothing refitted) were scored in an independent cohort: GSE87211, locally advanced rectal cancer given 5-FU-based chemoradiotherapy (5-FU + radiotherapy, with or without oxaliplatin or cetuximab). Outcome: disease-free survival. Direction fixed beforehand: higher score = predicted more resistant, so HR above 1.

## Data
189 patients with rectal tumours and DFS, 42 recurrences. All patients received 5-FU; there is no untreated comparator arm, so this tests prognosis in treated patients.

## Results
| id | signature | HR per SD (95% CI) | p | BH q (196 tests) | random-signature null p |
|---|---|---|---|---|---|
| RE1 | van de Wetering organoids, current data, top-10 | **1.91 (1.43-2.56)** | 1.1e-05 | **0.002** | 0.001 |
| RE2 | Gao colorectal PDX, top-10 | **1.65 (1.22-2.22)** | 0.001 | **0.097** | 0.026 |
| RE3 | paper (authors' 2019 files), rank-1 BH3-only | 1.38 (1.00-1.89) | 0.048 | 0.650 | 0.326 |
| RE4 | paper (authors' 2019 files), top-10 | 1.44 (1.04-1.98) | 0.026 | 0.650 | 0.214 |
| RE5 | all-model benchmark: frozen top-10 C-index | C 0.684 (0.597-0.767) | 1.1e-05 | | |
| RE6 | GSE87211 oxaliplatin-response x joint-score interaction (90 vs 111 patients, 13 vs 22 responders) | OR 2.24 (1.03-4.86) | 0.041 | not corrected (1 of 32) | |
(`POSITIVE_INDEX.tsv` gives RE1 as 1.91; the 05 report text rounds the same fit to 1.92.) Sources: `tables/05__B1_external_results.tsv`, `tables/05__B1_frozen_signatures.tsv`, `tables/19__positives_all_models_rest.tsv`, `tables/03__interaction_all.tsv`.

## How robust is it
- RE1 is the only one of 196 external tests to survive BH at q < 0.05; RE2 passes q < 0.10. Both beat random 10-pathway signatures (null p 0.001 and 0.026).
- All-model benchmark: no added model beats the frozen top-10 (best linear SVM 0.665, difference -0.019, p 0.51). Random shortlist-size sets reach 0.684 in 0.5% of draws. The rank-1 pathway alone is much weaker (C 0.573).
- The signal is consistent across signatures built from different lab models (organoids, PDX).

## Caveats
- Prognostic in treated patients only. With no untreated arm it cannot be called drug-specific, and the same kind of score separates survival in untreated patients in other cohorts (nar/05, nar/14).
- The W2 response model for the same cohort (n 201, 35 responders, response = no tumour at surgery) gives AUC 0.427 (p 0.18), i.e. no response signal; the DFS signal is not a response signal.
- 42 events.

## What I think
This is the best external validation of the organoid-derived 5-FU pathways, and it is stable across signatures. It most likely reflects a general recurrence-risk signature that the organoid pathways pick up, not proof that organoid response transfers.

## Where things are
`figures/`, `tables/`. Checksums in `../COPY_MANIFEST.tsv`.
""")

# ---------- ovarian
open(f'{OUT}/OVARIAN/README.md', 'w').write("""# Ovarian cancer: response cohorts, joint models, olaparib signature, A2780 isogenic, BRCA

## What was tested
Ovarian chemotherapy response prediction from the breast/pan-ovarian pathway model, the joint organoid+patient model, the user's olaparib data, and sanity checks.

## Results
| id | test | n / events | effect (95% CI) | p | tier |
|---|---|---|---|---|---|
| OV1 | GSE30161 response (carboplatin/cisplatin +/- paclitaxel; CR vs PR/PD), trained on other ovarian cohorts | 55 / 32 CR | AUC 0.728 (0.596-0.860) | **0.0043**, q 0.014 | strong |
| OV2 | GSE9891 response (RFS >= 12 months), same | 229 / 146 | AUC 0.624 (0.551-0.697) | **0.0018**, q 0.006 | strong |
| OV3 | olaparib correlation signature (6 ovarian lines, user's GSE229119) in platinum-treated GSE63885 | 75 / 66 deaths | HR 2.09 (1.24-3.50) | 0.0045 (exact permutation 0.014) | nominal |
| OV4 | BRCA1/2 somatic mutated vs wild type, platinum-treated TCGA-OV, PFS | 253 / 184 (16 mutated) | HR 0.42 (0.21-0.82) | 0.0087 | nominal |
| OV4 | same, OS | 253 / 149 | HR 0.47 (0.22-1.01) | 0.049 | nominal |
| OV5 | A2780 vs cisplatin-resistant (GSE15709), ovarian-organoid rank-1 pathway (polyamines) | 5 vs 5 | t 3.03 | 0.017 | nominal |
| OV5 | same pair, lung cell-line rank-1 pathway | 5 vs 5 | t -2.65 | 0.029 | nominal |
| OV6 | joint model, paclitaxel with platinum, 3 held-out cohorts | 564 / 426 | HR 1.14 (1.03-1.26) risk scale | 0.0099 | strong |
| OV6 | patient-only | 564 / 426 | HR 1.22 (1.10-1.34) | 1.1e-04 | strong |
| OV7 | joint model, platinum (cell-line proxy), 3 cohorts | 683 / 523 | HR 1.17 (1.07-1.29) | 6.6e-04 | strong |
| OV7 | patient-only | 683 / 523 | HR 1.19 (1.08-1.30) | 2.5e-04 | strong |
| OV8 | Vias organoid rank-1 signature, GSE9891 RFS | 240 / 179 | HR 1.23 (1.06-1.43) | 0.0059 | nominal (q 0.379 over 196) |
Sources: `tables/18__part1_v2.tsv`, `tables/22__ovary_signature_survival.tsv`, `tables/22__ovary_BRCA_somatic_survival_platinum.tsv`, `tables/10__Fig4bc_isogenic_other_cancers.tsv`, `tables/03__summary_main.tsv`, `tables/05__B1_external_results.tsv`, `tables/19__positives_all_models_rest.tsv`.

## Best ovarian result and whether it meets the rule
The best ovarian result by size and evidence is OV2 (GSE9891, 229 patients, AUC 0.624, p 0.0018) and the joint platinum survival model (OV7, C 0.583, p 6e-08 in the benchmark file). Both meet the rule. The best single lab-model-driven result is OV3 (olaparib signature, GSE63885) and it meets the p < 0.05 arm only nominally.

## How robust is it
- OV1/OV2: from the W2 run (pool of ovarian cohorts). GSE9891's AUC was 0.565 with the smaller v1 pool; the other three ovarian response cohorts (GSE15622 0.483, GSE51373 0.615, GSE63885 0.611) are not significant.
- OV3 is exploratory. The same signature gives p 0.65 (TCGA-OV) and 0.54 (GSE9891); with 6 lines a perfect correlation occurs by chance in 69% of orderings; ridge, SVM and L1 trained on the same lines give p 0.27-0.38; XGBoost has C 0.580 (not different from the original). Separately, the paper-method olaparib model points the other way in GSE9891 (HR 0.58, p 0.005) and TCGA-OV: see `../BORDERLINE.md`.
- OV6/OV7: carried by patient data (organoid- or cell-line-only HR 0.97 and 1.06, not significant). Platinum benchmark: joint C 0.583 and AUC 0.692 are not beaten by any added model.
- OV4 is a clinical sanity check on the cohort (germline BRCA is not public; 16 somatic mutants), not a result of the pathway model.

## What I think
Ovarian gives a real but modest signal in the larger response cohorts and in survival, mostly from patient-trained models. The olaparib signature in GSE63885 is intriguing, but with six lines and one cohort I would keep it exploratory.

## Where things are
`figures/`, `tables/`, `source_reports/22__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
""")

# ---------- pan-cancer DepMap
open(f'{OUT}/PAN_CANCER_DEPMAP/README.md', 'w').write("""# Cell lines (DepMap/GDSC/CTRPv2): the BH3-only pathway score tracks sensitivity to navitoclax and ABT-737

## What was tested
The colon 5-FU model's most stable pathway, "Activation of BH3-only proteins", was scored (single-sample rank score) in cancer cell lines and correlated with drug response, with cancer type removed (lineage-adjusted partial Spearman). Pre-stated reading: a higher score means the cell is closer to the death threshold, so more sensitive; negative rho for response metrics where higher = more resistant.

## Data
GDSC1 (378 drugs), GDSC2 (286), CTRPv2 (460), PRISM (1,442); about 500-960 lines per drug; CCLE/Cell Model Passports RNA-seq.

## Results (BH q over all drugs of the screen)
| screen | drug | n | partial rho | p | BH q | drug rank in screen | pathway rank of 672 (after BCL2, BCL2L1, MCL1, E2F removed) |
|---|---|---|---|---|---|---|---|
| GDSC2 | ABT-737 | 951 | -0.246 | 4.9e-14 | 1.4e-11 | 1 | 1 (added R2 0.0216, F-test p 7.6e-09) |
| GDSC2 | navitoclax | 961 | -0.202 | 6.2e-10 | 8.8e-08 | 2 | 1 (R2 0.0104, p 4.5e-05) |
| GDSC1 | navitoclax | 926 | -0.206 | 6.6e-10 | 2.5e-07 | 2 | 1 (R2 0.0098, p 2.9e-04) |
| CTRPv2 | ABT-737 | 706 | -0.238 | 3.3e-10 | 1.5e-07 | 2 | 1 (R2 0.0189, p 3.2e-05) |
| CTRPv2 | navitoclax | 772 | -0.211 | 5.8e-09 | 1.3e-06 | 5 | 1 (R2 0.0178, p 3.4e-05) |
Share of 1,000 random same-size gene sets as strong: 0.000 in all five. Other BH3 mimetics: venetoclax GDSC2 -0.104 (p 0.0017), CTRPv2 -0.226 (p 1.2e-06); WEHI-539 GDSC2 beta -0.212 (p 1.1e-06); obatoclax does not follow (GDSC1 +0.046, GDSC2 +0.036), which fits its off-target action. Sources: `tables/08__A_e_bh3_vs_all_drugs.tsv`, `tables/08__A_g_bh3_vs_bh3_mimetics.tsv`, `tables/13__G7_regression_added_value.tsv`, `tables/13__G7_partial_correlations.tsv`, `tables/08__S4_permutation_random_controls.tsv`.

## How robust is it
Seen in three screens (GDSC1, GDSC2, CTRPv2) and across five drug-screen pairs; survives removal of BCL2, BCL2L1, MCL1 and proliferation (partial r -0.202 -> -0.134 in GDSC2 navitoclax); ranks first of 672 pathways every time. PRISM is weaker (navitoclax -0.112, p 0.017, q 0.42; ABT-737 p 0.19).

## Caveats
- The effect size is modest (partial r about -0.2, added R2 about 1-2%).
- GDSC1, GDSC2 and CTRPv2 share many cell lines, so the three screens are not fully independent. The mimetic link was found by scanning all drugs (nar/08 part e) and then examined in more detail (nar/13 G7); q values correct for the drug scan within each screen.
- It is not a 5-FU result: 5-FU sits in the middle of the ranking (pan-cancer rho 0.02), and within colorectal lines the score does not follow 5-FU (rho 0.05, 0.02, 0.26, 0.41). Within colorectal lines alone navitoclax shows no link, except PRISM at the opposite sign (rho +0.405, n 27, p 0.036; `../BORDERLINE.md`). The 19 colorectal organoids cannot test it (navitoclax tested only to 2 uM; rho 0.05).

## What I think
A real, replicated biological readout: the score behaves like an apoptotic-priming marker for BCL-2/BCL-XL inhibitors. It supports the biology of the pathway, not the organoid-to-patient claim.

## Where things are
`figures/`, `tables/`, `source_reports/08__REPORT.md`. Checksums in `../COPY_MANIFEST.tsv`.
""")


# ---------- INDEX.md
pi = pd.read_csv(f'{OUT}/POSITIVE_INDEX.tsv', sep='\t').fillna('')
pi = pi[pi.tier != 'context']
order = ['BREAST', 'LUNG', 'BLADDER', 'COLORECTAL', 'RECTAL', 'OVARIAN', 'PAN_CANCER_DEPMAP']
pi['o'] = pi.cancer.map(order.index); pi = pi.sort_values(['o'], kind='stable')
def cell(x): return str(x).replace('|', '/')
lines = ['| cancer | id | result | tier | n | events / responders | effect | 95% CI | p | q | rule arm met | where |', '|' + '---|' * 12]
for r in pi.itertuples():
    where = f"`{r.cancer}/README.md`; source {cell(r.source)}"
    lines.append(f'| {r.cancer} | {r.id} | {cell(r.result)} | {r.tier} | {r.n} | {cell(r.events_or_responders)} | {cell(r.effect)} | {cell(r.ci95)} | {cell(r.p)} | {cell(r.q)} | {cell(r.rule_arm)} | {where} |')
cnt = pi.groupby(['cancer', 'tier']).size().unstack(fill_value=0).reindex(order)
open(f'{OUT}/INDEX.md', 'w').write(f"""# Index of positive results

Rule applied (also in `README.md`): **{RULE}**.
Tiers: **strong** = replicated, p < 0.01, or q < 0.10 with controls; **nominal** = meets p < 0.05 in the expected direction but is one test among many, is small, or fails family correction. Both tiers meet the rule; the tier tells you how much weight to give them. Results with p < 0.05 in the opposite direction, with no pre-stated direction, secondary analyses and unstable results are in `BORDERLINE.md`.
Every number below is read by `code/collect.py` from the source table in the last column (full table with robustness notes: `POSITIVE_INDEX.tsv`).

Counts per cancer and tier:

```
{cnt.to_string()}
```

""" + chr(10).join(lines) + chr(10))

# ---------- BORDERLINE.md
bi = pd.read_csv(f'{OUT}/BORDERLINE_INDEX.tsv', sep='\t').fillna('')
bl_lines = ['| cancer | id | result | why it is not a positive | n | events / responders | effect | 95% CI | p | q | source |', '|' + '---|' * 11]
for r in bi.itertuples():
    bl_lines.append(f'| {r.cancer} | {r.id} | {cell(r.result)} | {cell(r.why_not_positive)} | {r.n} | {cell(r.events_or_responders)} | {cell(r.effect)} | {cell(r.ci95)} | {cell(r.p)} | {cell(r.q)} | {cell(r.source)} |')
open(f'{OUT}/BORDERLINE.md', 'w').write(f"""# Borderline and opposite-direction results (not counted as positives)

These have a small p-value or look interesting but do not meet the rule (held-out test, p < 0.05 in the pre-stated direction, or q < 0.10), or they meet it only through a secondary or unstable analysis. Numbers are read from the source tables by `code/collect.py` (`BORDERLINE_INDEX.tsv`).

What I think about the groups:
- **HER2 randomised trials (NOAH, CHER-LOB).** The arm difference p 0.025 (NOAH) and 0.039 (CHER-LOB dual arm) is in the direction where the score separates responders better in the control arm (AUC 0.589 vs 0.328 in NOAH; 0.578 vs 0.242 in CHER-LOB), so it is not a drug-specific effect. No direction was stated beforehand. Interaction p 0.041 (NOAH), 0.128 and 0.379 (CHER-LOB).
- **Opposite direction.** Pancreas gemcitabine (Shi organoids, HR 0.52, p 0.023), Novartis melanoma PDX dacarbazine (HR 0.21, p 0.0011, BH q 0.063), lung cell-line pemetrexed (HR 0.19, p 0.0048), TCGA-PAAD joint-model interaction (HR 0.47), GSE63885 paclitaxel interactions, and the olaparib model in GSE9891 (HR 0.58, p 0.005). A small p in the wrong direction is not support for the method.
- **Not a model-vs-chance test.** Linear SVM beats fixed L2 in breast by +0.012 (p 0.023, BH q 0.100): a comparison of two models, and q is not below 0.10.
- **Secondary or post-hoc.** Registered gemcitabine score with proliferation added (HR 1.49, p 0.006; the primary test is p 0.199); GDSC2 models: only rank-1 pathways are the pre-stated test, so rank-2/3 picks such as DCC-mediated signalling are leads only.
- **Registered oxaliplatin test.** GSE72970 one-sided p 0.024 is just under the plan's 0.025 threshold but power was 30%, other oxaliplatin cohorts point away, and the gate (0 of 7 drugs) failed.
- **Unstable or not drug-specific.** Lung immunotherapy AUCs (GSE218989 0.566, GSE161537 0.681, GSE126044 0.836) depend on the model variant and training pool; STORM liver (sorafenib arm AUROC 0.965 but placebo 0.894-0.957); GDSC2 liver sorafenib rank-1 split (p 0.045, interval 0.95-12.84).
- **Opposite sign inside a subset.** PRISM navitoclax in colorectal lines (rho +0.405) goes against the pan-cancer finding; GSE87211 response AUC 0.427.

""" + chr(10).join(bl_lines) + chr(10))

# ---------- top README.md
open(f'{OUT}/README.md', 'w').write(f"""# POSITIVE_RESULTS

Everything in the project that meets the positive rule, one folder per cancer. Files were **copied** from `nar/` and `robust/results/`; no original was moved or edited (`COPY_MANIFEST.tsv` lists source, destination and checksum of every file).

## The rule
A result is positive if it is **{RULE}**.
- Held-out or independent means the test patients, cell lines or mutants were not used to train or tune the model (leave-one-cohort-out, a different cohort, TCGA patients scored by a lab-trained model, isogenic cell pairs).
- Pre-stated direction: the direction the method predicts (higher predicted resistance = worse survival, HR above 1; AUC above 0.5; mutant more resistant or sensitive as the drug biology says). Where a direction was not fixed beforehand, or the result goes the other way, it is not positive.
- Every positive is given a tier, **strong** or **nominal** (see `INDEX.md`). Nominal means it meets p < 0.05 but it is one test of many, small, or does not survive family correction. I did not hide those, because the rule asks only for p < 0.05 or q < 0.10; read the tier before you use one.
- Borderline and opposite-direction results are in `BORDERLINE.md`.

## Folders
`BREAST/`, `LUNG/`, `BLADDER/`, `COLORECTAL/`, `RECTAL/`, `OVARIAN/`, `PAN_CANCER_DEPMAP/`. Each has `README.md` (what was tested, data, model, numbers, robustness, caveats, "What I think", sources), `figures/` (PNG + PDF + caption; file names start with the number of the source folder), `tables/` (copies of the source tables), `source_reports/`.
No LIVER, PANCREAS, SARCOMA or melanoma folder: nothing there met the rule in the expected direction (see `BORDERLINE.md`).
Top-level: `INDEX.md` (one table of all positives), `POSITIVE_INDEX.tsv` and `BORDERLINE_INDEX.tsv` (machine-readable; numbers read from the source tables), `COPY_MANIFEST.tsv`, `code/` (`collect.py`, `copy_files.py`, `write_readmes.py`; run in that order, one process, seconds).

## Scan coverage
Reports and tables read: nar/02, 03, 05, 08, 09, 10, 13, 14, 17, 18, 19, 20, 22 and robust/results/train_combined/. nar/09 (registered test: 0 of 7 drugs pass its gates) and nar/14 (shuffle null, meta-analysis, joint-model grid; G9 shuffle p for the 57-test survey 0.030, pooled HR per SD 0.96) contribute only context and borderline rows. The G14 cell-line models, STORM liver and the other drugs and cancers (pancreas, sarcoma, melanoma, liver) produced no result in the expected direction that met the rule outside BORDERLINE.
""")
print('readmes, INDEX, BORDERLINE written')
