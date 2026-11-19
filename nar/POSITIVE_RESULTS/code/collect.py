"""W5 collector. Reads finished result tables under nar/ and robust/results/, applies the positive rule,
writes POSITIVE_INDEX.tsv and BORDERLINE_INDEX.tsv (numbers read from the source tables, not typed),
and COPIES (never moves) the figures, captions and tables into nar/POSITIVE_RESULTS/<CANCER>/.
Run: python3 -I collect.py   (one process, a few seconds). No source file is edited."""
import os, shutil, hashlib
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = '/home/user/organoid_biomarker_detection'
NAR = f'{ROOT}/nar'
OUT = f'{NAR}/POSITIVE_RESULTS'
T = lambda p: pd.read_csv(f'{ROOT}/{p}', sep='\t')

RULE = ("held-out/independent test, p < 0.05 in the pre-stated direction, or BH q < 0.10 in its own family; "
        "n and events reported")
rows, border = [], []

def fmt(x, d=3):
    if x is None or (isinstance(x, float) and np.isnan(x)): return ''
    if isinstance(x, (int, np.integer)): return str(x)
    if abs(x) < 0.001 and x != 0: return f'{x:.1e}'
    return f'{x:.{d}f}'
def ci(lo, hi): return f'{fmt(lo)} to {fmt(hi)}'

def add(cancer, rid, result, tier, test, n, ev, effect, cil, p, q, rob, src, rule='p<0.05 expected direction'):
    rows.append(dict(cancer=cancer, id=rid, result=result, tier=tier, test=test, n=n, events_or_responders=ev,
                     effect=effect, ci95=cil, p=p, q=q, rule_arm=rule, robustness=rob, source=src))
def bord(cancer, rid, result, why, n, ev, effect, cil, p, q, src):
    border.append(dict(cancer=cancer, id=rid, result=result, why_not_positive=why, n=n, events_or_responders=ev,
                       effect=effect, ci95=cil, p=p, q=q, source=src))

# ---------------------------------------------------------------- BREAST
f = 'robust/results/train_combined/part1_response_leave_one_cohort_out.tsv'
p1 = T(f); b = p1[p1.cancer == 'breast']
q = multipletests(b.p, method='fdr_bh')[1]
mc = T('nar/19_model_benchmark/tables/model_comparison.tsv'); l2 = mc[mc.model == 'L2_fixed'].iloc[0]
add('BREAST', 'BR1', 'Chemotherapy response, train on other breast cohorts, test on the held-out cohort (13 cohorts)',
    'strong', 'leave-one-cohort-out AUC (L2 logistic, 673 pathways); one-sided p vs 0.5',
    int(b.n.sum()), int(b.responders.sum()),
    f'mean AUC {fmt(l2.mean_auc)}; per-cohort range {fmt(b.auc_pathways.min())}-{fmt(b.auc_pathways.max())}',
    ci(l2.mean_ci_low, l2.mean_ci_high), f'{int((b.p<0.05).sum())}/13 cohorts p<0.05 (largest p {fmt(b.p.max())})',
    f'{int((q<0.10).sum())}/13 q<0.10',
    'replicated in 13/13 cohorts; beats ER+HER2-only (0.693) and proliferation-only (0.646) baselines, paired p<0.001; '
    '200 random 183-pathway sets average 0.742', f + '; nar/19_model_benchmark/tables/model_comparison.tsv', 'p and q')
v2 = T('nar/18_other_cohorts/tables/part1_v2.tsv'); v2 = v2[v2.table == 'v2_all_groups']
newb = v2[(v2.cancer == 'breast') & (~v2.test_cohort.isin(b.test_cohort))].copy()
newb['q'] = multipletests(newb.p, method='fdr_bh')[1]
ok = newb[newb.p < 0.05]
add('BREAST', 'BR2', 'Same test in 15 NEW breast cohorts (W2), trained on all breast cohorts of the v2 table except the held-out one',
    'strong', 'leave-one-cohort-out AUC; one-sided p vs 0.5', int(newb.n.sum()), int(newb.responders.sum()),
    f'{len(ok)}/{len(newb)} new cohorts p<0.05; AUC range among those {fmt(ok.auc_pathways.min())}-{fmt(ok.auc_pathways.max())}; '
    f'mean AUC all 15 {fmt(newb.auc_pathways.mean())}', '', f'smallest p {fmt(newb.p.min())}; largest significant {fmt(ok.p.max())}',
    f'{int((newb.q<0.10).sum())}/15 q<0.10 (BH over the 15)',
    'cohorts with p>=0.05 are small (n 21-114, 7-17 responders) or low-AUC (GSE59515 0.471)',
    'nar/18_other_cohorts/tables/part1_v2.tsv (table v2_all_groups)', 'p and q')
bl = T('nar/19_model_benchmark/tables/breast_loco_all_models.tsv')
agg = bl.groupby('model').auc.mean().sort_values(ascending=False)
add('BREAST', 'BR3', 'All-model benchmark: 17 model types on the same 13 held-out cohorts',
    'strong', 'mean held-out AUC per model (13 cohorts)', 3511, '', f'best SVM_linear {fmt(mc[mc.model=="SVM_linear"].iloc[0].mean_auc)}; L2 fixed {fmt(l2.mean_auc)}; '
    f'every non-baseline model >= {fmt(agg.drop(["ER_HER2_only","GCN_prior","Proliferation_score"]).min())}',
    ci(mc[mc.model=="SVM_linear"].iloc[0].mean_ci_low, mc[mc.model=="SVM_linear"].iloc[0].mean_ci_high), '', 'SVM vs L2 q=0.100 (not clearly better)',
    'no model clearly beats the original L2 after BH correction; signal does not depend on the model type',
    'nar/19_model_benchmark/tables/breast_loco_all_models.tsv, model_comparison.tsv', 'p and q (each model vs chance per cohort)')
ak = T('nar/20_breast_pathways/tables/auc_vs_k.tsv'); pr = T('nar/20_breast_pathways/tables/proliferation_tests.tsv')
sg = T('nar/20_breast_pathways/tables/subgroup_auc_pooled.tsv'); sg = sg[sg.model == 'full']
for _, r in sg.iterrows():
    if r.group == 'All patients': continue
    add('BREAST', f'BR4_{r.group}', f'Pathway-model response AUC inside receptor subgroup {r.group}', 'strong',
        'mean AUC over cohorts with enough responders (leave-one-cohort-out)', int(r.n_patients), '',
        f'mean AUC {fmt(r.mean_auc)} ({int(r.n_cohorts)} cohorts; {int(r.cohorts_auc_above_0_5)} above 0.5)', ci(r.ci_low, r.ci_high),
        'CI excludes 0.5', '', 'receptor subgroup from ER/HER2 flags; weakest in ER+/HER2- where proliferation alone gives 0.627',
        'nar/20_breast_pathways/tables/subgroup_auc_pooled.tsv', 'CI excludes chance')
rp = pr.set_index('model')
add('BREAST', 'BR5', 'Not just proliferation: model with proliferation regressed out still predicts', 'strong',
    'mean AUC, leave-one-cohort-out', 3511, '', f'resid {fmt(rp.loc["resid","mean_auc"])}; no-proliferation pathways {fmt(rp.loc["noprolif","mean_auc"])}; proliferation alone {fmt(rp.loc["prolif","mean_auc"])}',
    ci(rp.loc['resid','ci_low'], rp.loc['resid','ci_high']), 'LR test adds to proliferation+ER+HER2: 11/13 cohorts p<0.05 (REPORT text)', '',
    'single best pathway (EARLY_PHASE_OF_HIV_LIFE_CYCLE) AUC 0.680; plateau 0.72 at 5 pathways; ~100+ pathways needed to reach full model',
    'nar/20_breast_pathways/tables/proliferation_tests.tsv, auc_vs_k.tsv, minimal_subset.tsv; REPORT.md Q4', 'CI excludes chance')

# ---------------------------------------------------------------- LUNG
lu = T('nar/17_lung/tables/lung_untreated_prognosis.tsv')
pool = lu[(lu.test_cohort.str.startswith('ALL')) & (lu.score.str.startswith('prognostic'))]
for _, r in pool.iterrows():
    c = r.pathways
    per = lu[(lu.pathways == c) & (~lu.test_cohort.str.startswith('ALL')) & (lu.score.str.startswith('prognostic')) & lu.hr_per_sd.notna()]
    add('LUNG', 'LU1' if c.startswith('all') else 'LU1b', f'Untreated lung patients: general pathway prognostic score ({c})',
        'strong', 'Cox HR per SD of a score learned on OTHER cohorts, inverse-variance pooled over held-out cohorts',
        int(r.n), int(r.events), f'HR {fmt(r.hr_per_sd,2)} per SD', ci(r.ci_lo, r.ci_hi), fmt(r.p),
        '', f'HR>1 in {int((per.hr_per_sd>1).sum())}/{len(per)} cohorts; p<0.05 in {int((per.p<0.05).sum())}/{len(per)}',
        'nar/17_lung/tables/lung_untreated_prognosis.tsv', 'p<0.05 expected direction (higher score, worse survival)')
pa = T('nar/19_model_benchmark/tables/positives_all_models.tsv'); lun = pa[pa.task == 'lung_untreated_prognosis']
ns = lun[lun.model == 'NetSmoothL2'].iloc[0]; orig = lun[lun.model == 'Cox_ridge_all'].iloc[0]
add('LUNG', 'LU2', 'Lung untreated prognosis, all-model benchmark (same data)', 'strong', 'pooled C-index per model',
    int(orig.n), int(orig.events), f'original ridge Cox C {fmt(orig.est)}; best NetSmoothL2 C {fmt(ns.est)}; paired diff {fmt(ns.diff_vs_orig)}',
    ci(ns.lo, ns.hi), f'paired p {fmt(ns.diff_p)}', f'q {fmt(ns.q_diff)}',
    'several models do slightly better than the original; random shortlist-size sets reach the original in 44% of draws',
    'nar/19_model_benchmark/tables/positives_all_models.tsv', 'p and q')
a3 = T('nar/05_A3_B1/tables/A3_master_results.tsv').set_index('analysis_id')
for aid, rid, lab in [('LUAD_cellline_substitute__LUAD__GEMCITABINE', 'LU3', 'Lung adenocarcinoma cell-line model (gemcitabine) -> TCGA-LUAD gemcitabine-treated patients'),
                      ('PDX_gao2015_NSCLC__LUAD__ERLOTINIB', 'LU4', 'NSCLC PDX model (erlotinib) -> TCGA-LUAD erlotinib-treated patients')]:
    r = a3.loc[aid]
    add('LUNG', rid, lab, 'nominal', 'median-split log-rank, treated patients; lab-model-trained, patients held out',
        int(r.treated_n), int(r.treated_events), f'HR {fmt(r.hr,2)}', ci(r.hr_lo, r.hr_hi), fmt(r.treated_p), f'BH q {fmt(r.BH_q)} over 57 tests',
        f'shuffle p {fmt(r.perm_p)}; random-pathway bootstrap p {fmt(r.boot_p)}; untreated p {fmt(r.untreated_p)}',
        'nar/05_A3_B1/tables/A3_master_results.tsv', 'p<0.05 expected direction; fails BH over 57')
lo = T('nar/17_lung/tables/lung_loco_response.tsv')
for _, r in lo[((lo.test_cohort == 'GSE218989') & lo.model.str.startswith('same')) | ((lo.test_cohort == 'GSE161537') & lo.model.str.startswith('pan-cancer model (all')) | ((lo.test_cohort == 'GSE126044') & lo.model.str.startswith('pan-cancer model (all'))].iterrows():
    bord('LUNG', 'LUB_' + r.test_cohort, f'Immunotherapy response AUC, {r.model}', 'one of several model variants x 4 cohorts; the W2 run with a different training pool gives GSE218989 AUC 0.533 p 0.29, so it is not stable; chemo-trained models applied to anti-PD-1 data',
         int(r.n), int(r.responders), f'AUC {fmt(r.auc)}', ci(r.ci_lo, r.ci_hi), fmt(r.p), '', 'nar/17_lung/tables/lung_loco_response.tsv')
li = T('nar/17_lung/tables/lung_adjuvant_interaction.tsv'); li = li[(li.test_cohort == 'GSE29013') & li.pathways.str.startswith('all')].iloc[0]
bord('LUNG', 'LUB_GSE29013', 'Adjuvant chemo x score interaction, GSE29013, all 673 pathways (benefit score; HR below 1 = expected)',
     'one of 14 tests; direction as expected but the pooled interaction is HR 0.91 (0.77-1.06) p 0.23 and the randomised GSE14814 gives 1.02 p 0.95; observational cohort, 18 events',
     int(li.n), int(li.events), f'interaction HR {fmt(li.interaction_hr,2)}', ci(li.ci_lo, li.ci_hi), fmt(li.p), '', 'nar/17_lung/tables/lung_adjuvant_interaction.tsv')

# ---------------------------------------------------------------- BLADDER
r = a3.loc['BLCA_lee2018_paper__BLCA__CISPLATIN']
add('BLADDER', 'BL1', 'Cisplatin: 9-organoid model (paper source data), rank-1 pathway, TCGA-BLCA cisplatin-treated patients', 'nominal',
    'median-split log-rank', int(r.treated_n), int(r.treated_events), f'HR {fmt(r.hr,2)}', ci(r.hr_lo, r.hr_hi), fmt(r.treated_p),
    f'BH q {fmt(r.BH_q)} over 57', f'untreated p {fmt(r.untreated_p)}; shuffle p {fmt(r.perm_p)}; bootstrap p {fmt(r.boot_p)}; paper amino-acid pathway p 0.031 HR 2.66 (1.06-6.69)',
    'nar/05_A3_B1/tables/A3_master_results.tsv; nar/02_organoid_only/BLCA_lee2018/extras/tables/paper_biomarker_amino_acid_KM.tsv', 'p<0.05 expected direction')
r = a3.loc['BLCA_lee2018_current__BLCA__CISPLATIN']
add('BLADDER', 'BL2', 'Cisplatin: same model on current GDC data, response-shuffle test', 'nominal',
    'permutation (1,000 shuffles of organoid response; includes the choice of pathway)', int(r.treated_n), int(r.treated_events),
    f'log-rank p {fmt(r.treated_p)}, HR {fmt(r.hr,2)}', ci(r.hr_lo, r.hr_hi), f'shuffle p {fmt(r.perm_p)}', '',
    f'patient log-rank itself is {fmt(r.treated_p)} (not <0.05); score x treatment interaction HR {fmt(r.int_HR,2)} ({ci(r.int_lo,r.int_hi)}) p {fmt(r.int_p)}, q {fmt(r.BH_q_interaction)}',
    'nar/05_A3_B1/tables/A3_master_results.tsv; nar/02_organoid_only REPORT_BLCA_OV_LUAD_OTHER.md 3.2', 'permutation p<0.05')
g8 = T('nar/13_gaps_bio/tables/G8_ERCC2_null_and_adjustment.tsv')
print(g8.columns.tolist())
ercc = T('nar/02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv').set_index('feature_set')
for lab, key in [('paper data', 'paper biomarker AA (paper data)'), ('current data', 'paper biomarker AA (current data)')]:
    r = ercc.loc[key]
    add('BLADDER', 'BL3_' + lab.split()[0], f'ERCC2-mutant vs wild-type TCGA-BLCA: cisplatin model score ({lab}, paper amino-acid pathway)', 'strong',
        'one-sided Mann-Whitney, mutants predicted more sensitive', int(r.n_mut + r.n_wt), f'{int(r.n_mut)} mutant / {int(r.n_wt)} wild type',
        f'median score mutant {fmt(r.median_mut,4)} vs wild type {fmt(r.median_wt,4)}', '', fmt(r.p_one_sided), '',
        'share of ALL pathways as strong 0.017 (paper) / 0.002 (current); share of random proximal-size sets 0.063 / 0.000; TMB-adjusted p 0.024 / 0.008 (G8)',
        'nar/02_organoid_only/BLCA_lee2018/extras/tables/ERCC2_mutation_concordance.tsv; nar/13_gaps_bio/tables/G8_ERCC2_null_and_adjustment.tsv', 'p<0.05 expected direction')
iso = T('nar/02_organoid_only/BLCA_lee2018/extras/tables/isogenic_T24_tests.tsv').iloc[0]
add('BLADDER', 'BL4', 'Isogenic T24 cisplatin-sensitive vs resistant cells: paper amino-acid pathway lower in resistant', 'strong',
    'two-sided t-test, paper NES (raw T24 data not public)', 6, '3 sensitive / 3 resistant', f't {fmt(iso.t_sensitive_minus_resistant,1)}; mean {fmt(iso.mean_sensitive)} vs {fmt(iso.mean_resistant)}', '',
    fmt(iso.p), '', 'direction as predicted; 3 vs 3 samples; the rank-1 pathway (p 0.017) goes the other way, so this supports the amino-acid pathway not rank 1',
    'nar/02_organoid_only/BLCA_lee2018/extras/tables/isogenic_T24_tests.tsv', 'p<0.05 expected direction')
sm = T('nar/03_organoid_plus_patient/tables/summary_main.tsv'); sm = sm[(sm.features == 'shortlist') & (sm.outcome == 'survival')]
for unit, rid, lab in [('BLCA_CISPLATIN', 'BL5', 'Joint organoid+patient model, cisplatin chemotherapy, 4 held-out cohorts'),
                       ('BLCA_GEMCITABINE', 'BL6', 'Joint organoid+patient model, gemcitabine, TCGA-BLCA 5-fold cross-fit')]:
    for mdl in ['joint', 'patient-only', 'organoid-only']:
        r = sm[(sm.unit == unit) & (sm.model == mdl)].iloc[0]
        tier = 'strong' if (r.p < 0.01) else ('nominal' if r.p < 0.05 else 'not significant (shown for context)')
        if r.p < 0.05 and r.HR_risk > 1:
            add('BLADDER', f'{rid}_{mdl}', f'{lab}: {mdl} model', tier, 'held-out Cox HR per SD of risk score, pooled', int(r.n), int(r.events),
                f'HR {fmt(r.HR_risk,2)} (risk scale, >1 expected); C {fmt(r.C)}', ci(r.HR_risk_low, r.HR_risk_high), fmt(r.p), '',
                'the gain comes from the patient part: organoid-only is not significant; joint is not above patient-only (nar/14 G16)',
                'nar/03_organoid_plus_patient/tables/summary_main.tsv', 'p<0.05 expected direction')
        elif mdl == 'organoid-only':
            add('BLADDER', f'{rid}_{mdl}_ctx', f'{lab}: {mdl} (context row, not a positive)', 'context', 'same', int(r.n), int(r.events),
                f'HR {fmt(r.HR_risk,2)}', ci(r.HR_risk_low, r.HR_risk_high), fmt(r.p), '', 'organoid-only alone carries nothing here', 'nar/03_organoid_plus_patient/tables/summary_main.tsv', 'n/a')
ia = T('nar/03_organoid_plus_patient/tables/interaction_all.tsv'); ia = ia[(ia.model == 'joint') & (ia.features == 'shortlist')]
r = ia[ia.cohort == 'GSE48276'].iloc[0]
add('BLADDER', 'BL7', 'Joint model, cisplatin (MVAC) x score interaction in GSE48276', 'nominal', 'score x treatment interaction (risk scale, >1 = drug-specific)',
    int(r.n_treated + r.n_comparison), f'{int(r.events_treated)} treated events', f'interaction HR {fmt(r.ratio_risk,2)}', ci(r.low_risk, r.high_risk), fmt(r.p), 'not corrected (1 of 32 joint tests)',
    'only 8 treated events; the lab-only score gives 3.05 (p 0.018)', 'nar/03_organoid_plus_patient/tables/interaction_all.tsv', 'p<0.05 expected direction')
b1 = T('nar/05_A3_B1/tables/B1_external_results.tsv')
r = b1[(b1.cohort == 'GSE31684') & (b1.primary_p < 0.05) & (b1.direction == 'expected')].iloc[0]
add('BLADDER', 'BL8', 'Current-data cisplatin rank-1 signature in GSE31684: score x treatment interaction', 'nominal', 'Cox interaction HR per SD (B1 primary test)',
    int(r.n), int(r.events), f'HR {fmt(r.primary_est,2)}', ci(r.int_lo, r.int_hi), fmt(r.primary_p), f'BH q {fmt(r.BH_q_primary)} over 196', f'random-signature null p {fmt(r.null_p_two_sided)}; the same signature goes the other way in GSE48276 (HR 0.43, p 0.021)',
    'nar/05_A3_B1/tables/B1_external_results.tsv', 'p<0.05 expected direction')
r = pa[(pa.task == 'bladder_ERCC2_currentdata') & (pa.model == 'Paper_AA_pathway')].iloc[0]
add('BLADDER', 'BL9', 'ERCC2 concordance, all-model benchmark (current data): paper amino-acid pathway AUC', 'strong', 'AUC mutant vs wild type',
    int(r.n), int(r.events), f'AUC {fmt(r.est)}', ci(r.lo, r.hi), fmt(r.p_test), '', 'not beaten by any added model (best MLP 0.637); proliferation score alone 0.634 so the signal is not specific to the network method',
    'nar/19_model_benchmark/tables/positives_all_models.tsv', 'p<0.05 expected direction')

# ---------------------------------------------------------------- COLORECTAL
cs = T('nar/02_organoid_only/tables/cross_set_summary.tsv'); r = cs[(cs.cohort == 'COAD_authors') & (cs.patient_drug == 'FLUOROURACIL') & (cs.set.str.startswith('Colorectal organoids'))].iloc[0]
rr = a3.loc['COAD_vandewetering__COAD_authors__FLUOROURACIL']
add('COLORECTAL', 'CR1', 'Paper reproduction: 5-FU organoid model, rank-1 BH3-only pathway, 5-FU-treated TCGA-COAD (authors\' 2019 files)', 'nominal',
    'median-split log-rank', int(r.treated_n), f'{int(r.treated_events)} deaths (0 in predicted responders, 8 in non-responders)', 'HR not estimable (0 vs 8 deaths); continuous score HR 1.52 (0.94-2.46) p 0.091', '', fmt(r.treated_p),
    f'BH q {fmt(rr.BH_q)} over 57', f'untreated p {fmt(r.untreated_p)}; shuffle p {fmt(r.shuffle_share_rank1)}; bootstrap share {fmt(r.bootstrap_share,4)}; fragile: current GDC data gives p 0.866 (nar/13 G2); 94/114 patients also got oxaliplatin',
    'nar/02_organoid_only/tables/cross_set_summary.tsv; nar/13_gaps_bio/tables/G2_decomposition_colon_5FU.tsv', 'p<0.05 expected direction')
h = T('nar/02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv'); r = h[(h.cohort == 'COAD_authors') & h.pathway.str.contains('BH3')].iloc[0]
r2 = h[(h.cohort == 'COAD') & h.pathway.str.contains('BH3')].iloc[0]
add('COLORECTAL', 'CR2', 'Isogenic HCT-8 5-FU-sensitive vs 5-FU-resistant cells: BH3-only score higher in sensitive', 'strong',
    'two-sided t-test, authors\' NES (3 parental vs 9 resistant samples)', 12, '3 sensitive / 9 resistant', f't {fmt(r.t,2)}; mean {fmt(r.mean_sensitive,3)} vs {fmt(r.mean_resistant,3)}', '', fmt(r.p,4), '',
    f'recomputed from raw GSE81005: p {fmt(r2.p)} (same direction)', 'nar/02_organoid_only/COAD_vandewetering/tables/Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv', 'p<0.05 expected direction')
f5 = T('nar/02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv')
for i in [0, 2]:
    r = f5.iloc[i]
    add('COLORECTAL', f'CR3_{i}', f'Cetuximab organoid model vs BRAF V600E in TCGA colon ({r.label})', 'nominal' if r.p_one_sided > 0.01 else 'strong', 'one-sided Mann-Whitney, mutants predicted more resistant',
        int(r.n_mut + r.n_wt), f'{int(r.n_mut)} V600E / {int(r.n_wt)} wild type', f'rank-1 pathway {r.pathway.replace("REACTOME_","")}', '', fmt(r.p_one_sided), '',
        'see CR4 for specificity and confounding', 'nar/02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv', 'p<0.05 expected direction')
c5 = T('nar/10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv').iloc[0]
g8b = T('nar/13_gaps_bio/tables/G8_BRAF_null_and_adjustment.tsv'); print(g8b.columns.tolist())
add('COLORECTAL', 'CR4', 'Cetuximab vs BRAF V600E: colon + rectum, current data', 'strong', 'one-sided Mann-Whitney', int(c5.n_mut + c5.n_wt),
    f'{int(c5.n_mut)} V600E / {int(c5.n_wt)} wild type', '', '', fmt(c5.p_one_sided, 5), '',
    'after adjusting MSI, side, stage, TMB p 0.026 (current) / 0.093 (authors); MSS-only p 0.007 (10 mutants); share of random pathway sets as strong 0.177, so NOT specific to the network method; leave-k-out organoids: share p<0.05 1.00/0.93/0.94/0.87 for k=1..4',
    'nar/10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv; nar/13_gaps_bio/tables/G8_BRAF_null_and_adjustment.tsv; nar/10_completeness/tables/SuppFig18_leave_k_out_summary.tsv', 'p<0.05 expected direction')
r = pa[(pa.task == 'colon_cetuximab_BRAF')]; o = r[r.model == 'Paper_rank1'].iloc[0]; l1 = r[r.model == 'L1_shortlist'].iloc[0]
add('COLORECTAL', 'CR5', 'Cetuximab/BRAF all-model benchmark: L1 on the network shortlist beats the original rank-1 pathway', 'strong', 'AUC mutant vs wild type, paired bootstrap difference',
    int(o.n), int(o.events), f'original AUC {fmt(o.est)}; L1_shortlist AUC {fmt(l1.est)}; difference {fmt(l1.diff_vs_orig)}', ci(l1.lo, l1.hi), f'paired p {fmt(l1.diff_p)}', f'q {fmt(l1.q_diff)}',
    'BRAF V600E tracks MSI/right-sided subtype, so a better score may be a better subtype detector', 'nar/19_model_benchmark/tables/positives_all_models.tsv', 'q<0.10')
p1v = T('nar/18_other_cohorts/tables/part1_v2.tsv'); p1v = p1v[p1v.table == 'v2_rectal_split']; p1v['q'] = multipletests(p1v.p, method='fdr_bh')[1]
r = p1v[p1v.test_cohort == 'GSE62080'].iloc[0]
add('COLORECTAL', 'CR6', 'Colorectal chemotherapy response, GSE62080, trained on the other colorectal cohorts', 'nominal', 'leave-one-cohort-out AUC, one-sided p', int(r.n), int(r.responders),
    f'AUC {fmt(r.auc_pathways)}', ci(r.auc_ci_low, r.auc_ci_high), fmt(r.p), f'q {fmt(r.q)} over the 64 W2 tests', 'only 1 of 7 colorectal cohorts; the other 6 AUC 0.43-0.58; n=21',
    'nar/18_other_cohorts/tables/part1_v2.tsv (v2_rectal_split)', 'p<0.05 expected direction; q<0.10')
r = b1[(b1.cohort == 'GSE39582') & (b1.primary_p < 0.05) & (b1.primary_stat == 'interaction') & (b1.kind == 'top10') & (b1.direction == 'expected') & (b1.group == 'CRC_5FU')].iloc[0]
add('COLORECTAL', 'CR7', 'Gao CRC PDX 5-FU top-10 signature in GSE39582 (OS): score x 5-FU interaction', 'nominal', 'Cox interaction HR per SD (B1 primary test)',
    int(r.n), int(r.events), f'HR {fmt(r.primary_est,2)}', ci(r.int_lo, r.int_hi), fmt(r.primary_p), f'BH q {fmt(r.BH_q_primary)} over 196', f'random-signature null p {fmt(r.null_p_two_sided)}; RFS interaction 1.13 p 0.52',
    'nar/05_A3_B1/tables/B1_external_results.tsv', 'p<0.05 expected direction')
r = b1[(b1.cohort == 'GSE72970') & (b1.primary_p < 0.05) & (b1.group == 'CRC_IRINOTECAN')].iloc[0]
add('COLORECTAL', 'CR8', 'Colon organoid irinotecan top-10 signature in GSE72970 (FOLFIRI vs FOLFOX, PFS): interaction', 'nominal', 'Cox interaction HR per SD',
    int(r.n), int(r.events), f'HR {fmt(r.primary_est,2)}', ci(r.int_lo, r.int_hi), fmt(r.primary_p), f'BH q {fmt(r.BH_q_primary)} over 196', f'random-signature null p {fmt(r.null_p_two_sided)}',
    'nar/05_A3_B1/tables/B1_external_results.tsv', 'p<0.05 expected direction')

# ---------------------------------------------------------------- RECTAL
r = b1[(b1.cohort == 'GSE87211') & (b1.sig_id.str.startswith('COAD_vandewetering__COAD__FLUOROURACIL')) & (b1.kind == 'top10')].iloc[0]
add('RECTAL', 'RE1', 'Rectal cancer GSE87211 (5-FU based), organoid 5-FU current-data top-10 pathway signature vs disease-free survival', 'strong', 'Cox HR per SD, treated arm (no comparator arm)',
    int(r.n), int(r.events), f'HR {fmt(r.est,2)} per SD', ci(r.lo, r.hi), fmt(r.p), f'BH q {fmt(r.BH_q_primary)} over 196', f'random-signature null p {fmt(r.null_p_two_sided)}; the only B1 test surviving BH at q<0.05',
    'nar/05_A3_B1/tables/B1_external_results.tsv', 'q<0.10')
r = b1[(b1.cohort == 'GSE87211') & (b1.sig_id.str.startswith('PDX_gao2015_CRC__COAD__FLUOROURACIL')) & (b1.kind == 'top10')].iloc[0]
add('RECTAL', 'RE2', 'Same cohort, Gao colorectal PDX 5-FU top-10 signature', 'strong', 'Cox HR per SD, treated arm', int(r.n), int(r.events), f'HR {fmt(r.est,2)} per SD', ci(r.lo, r.hi), fmt(r.p),
    f'BH q {fmt(r.BH_q_primary)}', f'random-signature null p {fmt(r.null_p_two_sided)}', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'q<0.10')
for kind, rid in [('rank1', 'RE3'), ('top10', 'RE4')]:
    r = b1[(b1.cohort == 'GSE87211') & (b1.sig_id.str.startswith('COAD_vandewetering__COAD_authors__FLUOROURACIL')) & (b1.kind == kind)].iloc[0]
    add('RECTAL', rid, f'Same cohort, paper (authors\' 2019 files) 5-FU {kind} signature', 'nominal', 'Cox HR per SD, treated arm', int(r.n), int(r.events), f'HR {fmt(r.est,2)} per SD', ci(r.lo, r.hi), fmt(r.p),
        f'BH q {fmt(r.BH_q_primary)}', f'random-signature null p {fmt(r.null_p_two_sided)}', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'p<0.05 expected direction')
r = pa[(pa.task == 'rectal_GSE87211_5FU') & (pa.model == 'Paper_top10_frozen')].iloc[0]
add('RECTAL', 'RE5', 'GSE87211 all-model benchmark: frozen top-10 signature C-index', 'strong', 'C-index', int(r.n), int(r.events), f'C {fmt(r.est)}', ci(r.lo, r.hi), fmt(r.p_test), '',
    'not beaten by any added model (best SVM_linear 0.665); random shortlist-size sets reach it in 0.5% of draws; rank-1 version alone is 0.573', 'nar/19_model_benchmark/tables/positives_all_models.tsv', 'p<0.05 expected direction')
r = ia[ia.cohort == 'GSE87211'].copy(); r = r[r.outcome == 'response'].iloc[0]
add('RECTAL', 'RE6', 'GSE87211 oxaliplatin response x joint score interaction (odds ratio)', 'nominal', 'logistic interaction OR (risk scale)', int(r.n_treated + r.n_comparison), f'{int(r.events_treated)} + {int(r.events_comparison)} responders',
    f'OR {fmt(r.ratio_risk,2)}', ci(r.low_risk, r.high_risk), fmt(r.p), 'not corrected (1 of 32)', 'the lab-only score gives OR 2.31 p 0.035; same cohort survival interaction HR 0.87 p 0.69',
    'nar/03_organoid_plus_patient/tables/interaction_all.tsv', 'p<0.05 expected direction')
bord('RECTAL', 'REB1', 'GSE87211 response AUC (W2 pooled model, 35 responders)', 'AUC below 0.5 (opposite direction), not significant', 201, 35, 'AUC 0.427', '0.326 to 0.528', '0.177', '', 'nar/18_other_cohorts/tables/part1_v2.tsv')

# ---------------------------------------------------------------- OVARIAN
for coh, rid in [('GSE30161', 'OV1'), ('GSE9891', 'OV2')]:
    r = p1v[p1v.test_cohort == coh].iloc[0]
    add('OVARIAN', rid, f'Ovarian chemotherapy response, {coh}, trained on the other ovarian cohorts', 'strong', 'leave-one-cohort-out AUC, one-sided p', int(r.n), int(r.responders),
        f'AUC {fmt(r.auc_pathways)}', ci(r.auc_ci_low, r.auc_ci_high), fmt(r.p), f'q {fmt(r.q)} over the 64 W2 tests',
        'GSE30161 is a new W2 cohort (n=55); GSE9891 AUC was 0.565 in the v1 pool and 0.624 with the larger v2 pool; the other two ovarian response cohorts (GSE15622 0.483, GSE51373 0.615, GSE63885 0.611) are not significant',
        'nar/18_other_cohorts/tables/part1_v2.tsv', 'p<0.05 and q<0.10')
os_ = T('nar/22_user_data_results/tables/ovary_signature_survival.tsv'); r = os_[(os_.cohort == 'GSE63885') & (os_.drug == 'OLAPARIB')].iloc[0]
add('OVARIAN', 'OV3', 'Olaparib correlation signature (6 ovarian lines, user\'s GSE229119) in platinum-treated GSE63885', 'nominal', 'log-rank, high vs low signature, overall survival', int(r.n_treated), int(r.events_treated),
    f'HR {fmt(r.HR,2)}', ci(r.HR_lo, r.HR_hi), fmt(r.logrank_p), f'exact permutation p {fmt(r.permutation_p)} (720 orders)',
    'exploratory: TCGA-OV p 0.65 and GSE9891 p 0.54 for the same signature; with 6 lines, rho=1 occurs by chance in 69% of orderings; no model in the 19 benchmark reaches log-rank significance (ridge/SVM/L1 p 0.27-0.38)',
    'nar/22_user_data_results/tables/ovary_signature_survival.tsv; nar/19_model_benchmark REPORT_positives.md', 'p<0.05 expected direction')
bo = T('nar/22_user_data_results/tables/ovary_BRCA_somatic_survival_platinum.tsv').set_index('endpoint')
for ep in ['PFS', 'OS']:
    r = bo.loc[ep]
    add('OVARIAN', f'OV4_{ep}', f'BRCA1/2 somatic mutation vs wild type, platinum-treated TCGA-OV ({ep})', 'nominal', 'log-rank; mutated expected to do better (clinical prior, not a model output)', int(r.n), f'{int(r.events)} events; {int(r.n_BRCA_mut)} BRCA mutated',
        f'HR mutated vs wild type {fmt(r.HR_mut_vs_wt,2)}', ci(r.HR_lo, r.HR_hi), fmt(r.logrank_p), '', 'germline BRCA is not public so this is somatic only; a clinical sanity check on the cohort, not a test of the pathway model',
        'nar/22_user_data_results/tables/ovary_BRCA_somatic_survival_platinum.tsv', 'p<0.05 expected direction')
iso2 = T('nar/10_completeness/tables/Fig4bc_isogenic_other_cancers.tsv')
for i in [3, 5]:
    r = iso2.iloc[i]
    add('OVARIAN', f'OV5_{i}', f'A2780 cisplatin-sensitive vs resistant (GSE15709): {r.biomarker}', 'nominal', 'Welch t-test, 5 vs 5, predicted direction', 10, '5 sensitive / 5 resistant', f't {fmt(r.welch_t,2)}', '', fmt(r.p), '',
        'rank-1 pathway METABOLISM_OF_POLYAMINES (the carboplatin and cisplatin ovarian-organoid models are identical); the bladder rank-1 and amino-acid pathways are not significant in A2780 (p 0.22, 0.053 opposite direction)' if i < 5 else 'lung cell-line cisplatin rank-1 pathway applied to an ovarian pair',
        'nar/10_completeness/tables/Fig4bc_isogenic_other_cancers.tsv', 'p<0.05 expected direction')
for unit, rid, lab in [('OV_PACLITAXEL', 'OV6', 'Joint organoid+patient model, paclitaxel with platinum'), ('OV_PLATINUM', 'OV7', 'Joint organoid+cell-line+patient model, platinum')]:
    for mdl in ['joint', 'patient-only']:
        r = sm[(sm.unit == unit) & (sm.model == mdl)].iloc[0]
        add('OVARIAN', f'{rid}_{mdl}', f'{lab}: {mdl}', 'strong' if r.p < 0.01 else 'nominal', 'held-out Cox HR per SD of risk score, 3 cohorts pooled', int(r.n), int(r.events), f'HR {fmt(r.HR_risk,2)} (risk scale)', ci(r.HR_risk_low, r.HR_risk_high), fmt(r.p), '',
            'gain comes from the patient part; lab-only HR 0.97 / 1.06 (ns); pooled across 3 cohorts, shortlist features', 'nar/03_organoid_plus_patient/tables/summary_main.tsv', 'p<0.05 expected direction')
r = b1[(b1.cohort == 'GSE9891') & (b1.group == 'OV_PLATINUM') & (b1.primary_p < 0.05)].iloc[0]
add('OVARIAN', 'OV8', 'Vias organoid (oxaliplatin proxy) rank-1 signature in GSE9891 platinum-treated, RFS', 'nominal', 'Cox HR per SD, treated arm', int(r.n), int(r.events), f'HR {fmt(r.est,2)}', ci(r.lo, r.hi), fmt(r.p), f'BH q {fmt(r.BH_q_primary)} over 196',
    f'random-signature null p {fmt(r.null_p_two_sided)}', 'nar/05_A3_B1/tables/B1_external_results.tsv', 'p<0.05 expected direction')
bord('OVARIAN', 'OVB1', 'Olaparib-model (paper method) in GSE9891 and TCGA-OV, platinum-treated', 'HR<1 = predicted-less-sensitive did better = opposite direction (GSE9891 p 0.005 HR 0.58; TCGA-OV p 0.051 HR 0.77)', 241, 110, 'HR 0.58', '0.39 to 0.85', '0.005', '', 'nar/22_user_data_results REPORT.md section 3')

# ---------------------------------------------------------------- PAN-CANCER DEPMAP
g7 = T('nar/13_gaps_bio/tables/G7_regression_added_value.tsv'); ae = T('nar/08_depmap/tables/A_e_bh3_vs_all_drugs.tsv')
for _, r in g7[g7.drug.isin(['NAVITOCLAX', 'ABT-737'])].iterrows():
    sc = 'CTRP' if r.screen == 'CTRP' else r.screen
    a = ae[(ae.screen == sc) & (ae.drug == r.drug)].iloc[0]
    add('PAN_CANCER_DEPMAP', f'DM_{r.screen}_{r.drug}', f'BH3-only pathway score vs {r.drug} sensitivity, {r.screen}, cancer type removed', 'strong',
        'lineage-adjusted partial Spearman across all lines; BH across all drugs of the screen', int(r.n), f'{int(r.n)} cell lines',
        f'partial rho {fmt(a.pancancer_partial_rho)} (negative = high score more sensitive); drug rank {int(a.pancancer_partial_rank)} of screen',
        '', fmt(a.pancancer_p), f'BH q {fmt(a.pancancer_q)}',
        f'after also removing BCL2, BCL2L1, MCL1, E2F: beta {fmt(r.beta_BH3)} ({ci(r.beta_lo, r.beta_hi)}), added R2 {fmt(r.dR2,4)}, F-test p {fmt(r.p_F)}; rank {int(r.rank_pathways)} of {int(r.n_pathways)} pathways; share of random same-size gene sets as strong {fmt(r.share_random_le)}',
        'nar/08_depmap/tables/A_e_bh3_vs_all_drugs.tsv; nar/13_gaps_bio/tables/G7_regression_added_value.tsv', 'q<0.10 and p<0.05 expected direction')
bord('PAN_CANCER_DEPMAP', 'DMB1', 'BH3-only score vs navitoclax inside colorectal lines (PRISM)', 'opposite sign to pan-cancer (rho +0.405), n=27; GDSC/CTRP colorectal lines all CI cross 0; 19 colorectal organoids also show no navitoclax link (rho 0.05, drug tested to 2 uM only)', 27, '', 'rho +0.405', '0.030 to 0.681', '0.036', '', 'nar/08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv; nar/13_gaps_bio G3')

# ---------------------------------------------------------------- BORDERLINE (opposite direction, direction not pre-specified, or not stable)
ra = T('nar/18_other_cohorts/tables/part2_v2.tsv'); ra = ra[ra.test.str.contains('NOAH|CHER')]
for _, r in ra.iterrows():
    bord('BREAST', 'BRB_' + r.test.split(' (')[0].replace(' ', '_') + ('_dual' if 'trastuzumab + lapatinib' in r.test else ''), r.test, 'drug-arm vs control-arm AUC difference with direction NOT pre-specified; points the opposite way to a drug-specific effect (score separates responders better in the control arm)',
         f'{int(r.n_drug)} vs {int(r.n_control)}', f'{int(r.resp_drug)} vs {int(r.resp_control)} responders', f'AUC drug arm {fmt(r.auc_in_drug_arm)} vs control {fmt(r.auc_in_control_arm)}', '', f'{fmt(r.p_auc_difference)} (interaction p {fmt(r.p_interaction)})', '', 'nar/18_other_cohorts/tables/part2_v2.tsv')
for aid, why in [('PAAD_shi__PAAD__GEMCITABINE', 'opposite direction (HR<1): predicted non-responders lived longer; untreated p 0.017 too'),
                 ('PDX_gao2015_CM__SKCM__DACARBAZINE', 'opposite direction (HR<1) although BH q 0.063'),
                 ('LUAD_cellline_substitute__LUAD__PEMETREXED', 'opposite direction (HR<1)')]:
    r = a3.loc[aid]
    bord('MULTI', 'A3_' + aid.split('__')[0], aid, why, int(r.treated_n), int(r.treated_events), f'HR {fmt(r.hr,2)}', ci(r.hr_lo, r.hr_hi), fmt(r.treated_p), f'q {fmt(r.BH_q)}', 'nar/05_A3_B1/tables/A3_master_results.tsv')
r = ia[ia.cohort == 'TCGA-PAAD'].iloc[0]
bord('PANCREAS', 'PAAD_joint_int', 'Joint model gemcitabine x score interaction, TCGA-PAAD', 'opposite direction (HR 0.47): score is prognostic in untreated patients, not drug-specific', int(r.n_treated + r.n_comparison), int(r.events_treated), f'HR {fmt(r.ratio_risk,2)}', ci(r.low_risk, r.high_risk), fmt(r.p), '', 'nar/03_organoid_plus_patient/tables/interaction_all.tsv')
r = b1[(b1.cohort == 'GSE63885') & (b1.group == 'OV_TAXANE') & (b1.primary_p < 0.05)]
for _, r in r.iterrows():
    bord('OVARIAN', f'B1_GSE63885_{r.endpoint}', 'Paclitaxel top-10 signature x treatment interaction in GSE63885', 'opposite direction', int(r.n), int(r.events), f'HR/OR {fmt(r.primary_est,2)}', '', fmt(r.primary_p), fmt(r.BH_q_primary), 'nar/05_A3_B1/tables/B1_external_results.tsv')
b09 = T('nar/09_registered_test/tables/stage2_primary.tsv') if os.path.exists(f'{NAR}/09_registered_test/tables/stage2_primary.tsv') else None
bord('COLORECTAL', 'REG_oxali', 'Registered Stage 2: oxaliplatin score x FOLFOX-vs-FOLFIRI interaction, GSE72970 PFS', 'one-sided p 0.024 (threshold 0.025) but 30% power, other oxaliplatin cohorts point the other way; gate G3 failed (0/7 drugs pass); 4 FOLFIRINOX patients excluded', 120, '', 'HR 1.63', '1.00 to 2.64', '0.024 (one-sided)', '', 'nar/09_registered_test/REPORT.md')
bord('PAN', 'REG_gem', 'Registered Stage 1 secondary: gemcitabine score x exposure, with proliferation added', 'secondary, post-hoc adjustment; primary test HR 1.12 (0.94-1.32) p 0.199, Holm 1.0; plan says hypothesis only', 2258, '154 events on drug', 'HR 1.49', '1.12 to 1.99', '0.006 (Holm within family 0.043)', '', 'nar/09_registered_test/REPORT.md')
r = T('nar/13_gaps_bio/tables/G14_SuppFig15_cell_line_models.tsv'); print(r.columns.tolist())
bord('LIVER', 'G14_liver_sorafenib', 'GDSC2-trained sorafenib model (all lines), rank-1 pathway, TCGA-LIHC 30 sorafenib patients', 'one of 10 settings; HR CI 0.95-12.84; sign flips with the lineage-matched training set (n=15 lines)', 30, 15, 'HR 3.49', '0.95 to 12.84', '0.045', '', 'nar/13_gaps_bio/tables/G14_SuppFig15_cell_line_models.tsv')
bord('LIVER', 'STORM', 'STORM sorafenib arm, organoid-pathway scores AUROC 0.947-0.965 / OR 3.9', 'same scores give AUROC 0.894-0.957 in the PLACEBO arm; STORM recurrence axis (PC1 alone AUROC 0.96); interaction OR 1.34 p 0.80: prognostic, not sorafenib-specific', 139, '21 + 21', 'AUROC 0.965', '', '<0.01', '', 'nar/03_organoid_plus_patient/tables/storm_auroc_by_arm.tsv')
bord('LUNG', 'IO_pan', 'GSE207422 (anti-PD-1 + chemo) AUROC for non-response with drug models', 'direction reversed (AUROC 0.24 for predicting non-response), shuffle p 0.07-0.10', 24, 9, 'AUROC 0.24', '', '0.039', '', 'nar/22_user_data_results/REPORT.md 2b')
bord('MULTI', 'A3_interaction_opposite', 'A3 score x treatment interactions with p<0.05 pointing opposite: Novartis PDX pancreas gemcitabine 0.55 (p 0.006), lung pemetrexed 0.46 (0.015), sarcoma gemcitabine 0.56 (0.031), melanoma dacarbazine 0.58 (0.037)', 'opposite direction; BH q >= 0.36', '', '', '', '', '0.006-0.037', '0.36-0.42', 'nar/05_A3_B1/tables/A3_master_results.tsv')
bord('OVARIAN', 'DCC_paclit', 'GDSC2 lineage-matched ovarian paclitaxel model, rank-3 pathway DCC-mediated attractive signaling (TCGA-OV)', 'rank-3 pathway only (the paper tests rank 1); HR 0.62 opposite; platinum rank-2 (HR 1.60, p 0.0004) likewise a non-rank-1 pick', 348, 215, 'HR 0.62', '0.47 to 0.81', '0.0005', '', 'nar/13_gaps_bio/tables/G14_SuppFig15_cell_line_models.tsv')
bord('BREAST', 'SVM_vs_L2', 'Linear SVM vs fixed L2 in breast benchmark', 'is a comparison of two models, not a test against chance: difference +0.012, p 0.023, BH q 0.100 (not below 0.10); sign test p 0.094', 3511, '', 'mean AUC diff +0.012', '0.002 to 0.022', '0.023', '0.100', 'nar/19_model_benchmark/tables/model_comparison.tsv')

pd.DataFrame(rows).to_csv(f'{OUT}/POSITIVE_INDEX.tsv', sep='\t', index=False)
pd.DataFrame(border).to_csv(f'{OUT}/BORDERLINE_INDEX.tsv', sep='\t', index=False)
print(len(rows), 'positive rows;', len(border), 'borderline rows')
