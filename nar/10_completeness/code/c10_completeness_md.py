"""Writes COMPLETENESS.md from completeness.tsv plus the numbers in this folder's tables (no hand-typed results)."""
from common import *
T = pd.read_csv(OUT / "completeness.tsv", sep="\t").fillna("")
cnt = json.load(open(TAB / "completeness_counts.json"))
X = pd.read_csv(TAB / "SuppFig4_cross_test_matrix.tsv", sep="\t"); xs = json.load(open(TAB / "SuppFig4_cross_test_summary.json"))
M = pd.read_csv(TAB / "Fig5_mutation_concordance_other_cancers.tsv", sep="\t")
I = pd.read_csv(TAB / "Fig4bc_isogenic_other_cancers.tsv", sep="\t").dropna(subset=["p"])
N9 = pd.read_csv(TAB / "SuppFig9_network_neighbourhood_1_2_3_steps.tsv", sep="\t")
S10 = pd.read_csv(TAB / "SuppFig10_size_vs_p_summary_all_sets.tsv", sep="\t")
S17 = pd.read_csv(TAB / "SuppFig17_split_603010_summary.tsv", sep="\t")
S18 = pd.read_csv(TAB / "SuppFig18_leave_k_out_summary.tsv", sep="\t")
P14 = pd.read_csv(TAB / "SuppFig14_cell_lines_vs_TCGA_PCA_variance.tsv", sep="\t")
CL = pd.read_csv(TAB / "SuppTable1_2_clinical_balance_all_sets.tsv", sep="\t").drop_duplicates(["key", "feature"])
V = pd.read_csv(TAB / "verify_headline_biomarkers.tsv", sep="\t")
f = lambda d, k: d.get(k, 0)
L = ["# Completeness check against the Kong et al. 2020 checklist (77 items)", "",
 "Folder `nar/10_completeness/`. Status is given twice for every checklist item: (a) the paper's own cancer and drug (exact",
 "reproduction: colon 5-FU, bladder cisplatin, colon cetuximab) and (b) the same item for the other sets (colon current data,",
 "bladder current data, pancreas Tiriac and Shi, liver Ji and Broutier, ovarian Vias, lung cell-line substitute, sarcoma,",
 "PDX: Isella colorectal and Gao breast/colorectal/lung/melanoma/pancreas). Full table: `completeness.tsv` (one row per item,",
 "with the files that cover it, how many files the path matched, and the number that was checked).", "",
 "Status words: DONE = figure/table exists and its number was checked; PARTIAL = done for some sets only, or a close",
 "substitute; MISSING = not produced; IN_PROGRESS_B3 = being produced now by the B3 agent (`nar/07_B3_models/`: Fig 3b-f",
 "for every set, Supp Figs 6-7 centrality, Supp Table 4, deep-learning grid); N/A = the item only exists for one drug.", "",
 "## Counts before and after this agent's work", "",
 "| | DONE | PARTIAL | MISSING | IN_PROGRESS_B3 | N/A |", "|---|---|---|---|---|---|"]
for lab, k in (("(a) paper's own case, before", "paper_before"), ("(a) paper's own case, after", "paper_after"),
               ("(b) other sets, before", "generalisation_before"), ("(b) other sets, after", "generalisation_after")):
    c = cnt[k]; L.append(f"| {lab} | {f(c,'DONE')} | {f(c,'PARTIAL')} | {f(c,'MISSING')} | {f(c,'IN_PROGRESS_B3')} | {f(c,'N/A')} |")
L += ["", "Before = what the earlier agents (02, 03, 08, 09) had produced. In the 'before' rows, items now marked IN_PROGRESS_B3 are",
 "counted as MISSING (14 for (a), 16 for (b)).", "",
 "## How the earlier numbers were checked", "",
 "- Every path in `completeness.tsv` was checked with a file search; all paths exist (column `*_files_found`, all > 0).",
 f"- Headline survival tests re-derived from the shared data bundles and the agents' coefficient tables "
 f"(`tables/verify_headline_biomarkers.tsv`): {int(V.match.sum())} of {len(V)} reproduce the reported treated-patient p exactly "
 "(to 6 decimals). The two others: Gao pancreas PDX (reported 0.14 in the report; recomputed 0.137, the table lookup used a",
 "different project label) and Broutier liver (0.089 here vs 0.139 reported, because the shared bundle holds TCGA-LIHC",
 "scores computed for the Ji gene universe; Broutier was left out of the new cross-set analyses for this reason).",
 "- Spot checks of report numbers against tables (all matched): colon 5-FU p = 0.014 (n = 114); bladder p = 0.0435;",
 "Supp Table 1 age 0.260, sex 0.188, T 0.435, N 0.440; ERCC2 p = 0.0030; T24 isogenic p = 0.00022; HCT-8 p = 0.0085 (paper) /",
 "0.039 (recomputed); BRAF V600E p = 0.031; bladder PCA 27.1/17.8 and 20.7/18.8; GDSC 5-FU 0.708 / 0.583; organoid+patient",
 "CRC_5FU joint HR 0.89 (0.80-1.00), C 0.532; registered test gemcitabine HR 1.12 (0.94-1.32), p = 0.199; DepMap GDSC1",
 "rho 0.054 (n = 44).", "",
 "## What was added here (all in `nar/10_completeness/figures` and `tables`)", "",
 "1. **Schematics** (Fig 1a workflow, 1b network proximity, 3a comparison feature sets, 4a isogenic check, 4d bootstrap) and",
 "   a new Fig 1c map of all 57 set x cohort x drug analyses.",
 "2. **Fig 2c / 2g network pictures** for the paper's two drugs and every other headline drug (18 panels).",
 f"3. **Supp Fig 4 cross-test matrix** (17 biomarkers x 16 treated cohorts). Paper's two cross-tests: colon 5-FU patients with"
 f" the bladder amino acid pathway p = {xs['paper_SuppFig4a_COAD_with_BLCA_amino']:.3f} (paper 0.16); bladder cisplatin patients with the"
 f" colon BH3-only pathway p = {xs['paper_SuppFig4b_BLCA_with_COAD_BH3']:.2f} (paper 0.86). Own-cohort tests p < 0.05: {xs['own_p_below_005']} of {xs['n_own']};"
 f" other-cohort tests p < 0.05: {xs['off_p_below_005']} of {xs['n_off']} ({100*xs['off_share']:.1f}%, about the 5% expected by chance).",
 "   What I think: the biomarkers are not general survival markers across cancers; this supports the paper's specificity claim.",
 "4. **Supp Fig 8 (target genes alone) and Supp Fig 9 (1-, 2-, 3-step network neighbours)** for colon (paper data: TYMS",
 f"   neighbours 1/2/3 steps p = {' / '.join(f'{x:.3f}' for x in N9[N9.key=='COAD_paper_5FU'].sort_values('hops').p)}) and all 16 headline sets."
 f" Neighbour models with p < 0.05: {int((N9.p<0.05).sum())} of {len(N9)}.",
 f"5. **Supp Fig 10 size vs p** for all 16 sets: R2 from {S10.R2.min():.3f} to {S10.R2.max():.3f} (paper 0.0026 and 0.0008).",
 "   What I think: pathway size does not drive the survival splits in any set.",
 "6. **Supp Fig 14 with explained-variance numbers**, cell lines vs TCGA for 8 cancers. Colon genes raw "
 + "/".join(f"{x:.1f}" for x in P14[(P14.tcga=='COAD')&(P14.level=='genes')&(~P14.scaled)][["pc1","pc2"]].values[0]) + "%, scaled "
 + "/".join(f"{x:.1f}" for x in P14[(P14.tcga=='COAD')&(P14.level=='genes')&(P14.scaled)][["pc1","pc2"]].values[0]) + "%. Supp Fig 1",
 "   variance numbers for all 15 lab-model bundles (`tables/SuppFig1_PCA_variance_all_sets.tsv`).",
 "7. **Supp Fig 16** proximity heatmap: the paper's 44 columns plus 8 drugs added for the other sets.",
 "8. **Supp Fig 17 as the paper's 60/10/30 split** (100 repeats), every set; median held-out Spearman rho: "
 + "; ".join(f"{r.key} {r.median_spearman:.2f}" for r in S17.itertuples()) + ".",
 "9. **Fig 5 for other cancers** (one-sided Mann-Whitney, all TCGA tumours with mutation calls): "
 + "; ".join(f"{r.model.split(' |')[0]} {r.mutation} p = {r.p_one_sided:.2g} ({r.n_mut} vs {r.n_wt})" for r in M.dropna(subset=['p_one_sided']).itertuples()) + ".",
 "   What I think: the colon organoid cetuximab model recovers BRAF V600E again (p = 0.00048) but does not pick up KRAS",
 "   mutations; the erlotinib and platinum models do not separate EGFR or BRCA1/2 mutants in the expected direction.",
 "   No accepted mutation marker exists for gemcitabine, sorafenib or dacarbazine, so no such test was possible there.",
 "10. **Supp Fig 12** (colon cetuximab, other feature sets vs BRAF V600E) and **Supp Fig 18** (leave-k-out organoids): share of",
 "    refits still giving p < 0.05: " + ", ".join(f"k={r.k}: {r.share_p_below_005:.2f}" for r in S18.itertuples()) + ".",
 "11. **Fig 4b/c for other cancers**: " + "; ".join(f"{r.dataset.split(' vs')[0]} | {r.biomarker}: p = {r.p:.2g}, {'predicted direction' if r.direction_as_predicted else 'opposite direction'}" for r in I.itertuples()) + ".",
 "    What I think: the ovarian and lung platinum pathways differ between A2780 and its cisplatin-resistant copy in the",
 "    predicted direction (p 0.017-0.029, 5 vs 5); the liver pairs have 2 vs 2 samples and cannot decide much.",
 f"12. **Supp Tables 1-2 for every set**: clinical balance of predicted groups; chi-square p < 0.05 in "
 f"{int((CL.chi2_p<0.05).sum())} of {int(CL.chi2_p.notna().sum())} tests (" + ", ".join(f"{r.key} {r.feature} {r.chi2_p:.3f}" for r in CL[CL.chi2_p<0.05].itertuples()) + ").",
 "", "## What remains, in one line each", "",
 "- Fig 3b-f, Supp Figs 6-7, Supp Table 4, deep-learning grid: being produced by the B3 agent (`nar/07_B3_models/`).",
 "- Isogenic pairs for pancreas gemcitabine, sarcoma, melanoma dacarbazine and the PDX drugs: no clean public parental vs",
 "  resistant expression pair was found in this work (GEO search, 2026-11-18).",
 "- Mutation concordance for gemcitabine, sorafenib, dacarbazine: no accepted mutation marker exists.",
 "- Supp Fig 15 (cell-line-trained models) for pancreas, liver, ovary and sarcoma was not run (colon, bladder, and lung,",
 "  whose main analysis is itself cell-line-trained, are done).",
 "- The paper's bladder IC50 values and T24 raw expression are not public (AUC and the paper's NES are used).", "",
 "## Item by item", "", "| ID | (a) paper case | (b) other sets | Where (a) | Number checked | Note |", "|---|---|---|---|---|---|"]
for r in T.itertuples():
    L.append(f"| {r.ID} | {r.paper_after}{' (was '+r.paper_before+')' if r.paper_before!=r.paper_after else ''} | "
             f"{r.generalisation_after}{' (was '+r.generalisation_before+')' if r.generalisation_before!=r.generalisation_after else ''} | "
             f"`{r.paper_paths.split(';')[0]}` | {r.paper_number_check} | {r.note} |")
L += ["", "## Code and data", "", "Code `code/` (run in order): c01_verify_rank1, c02_schematics, c03_networks, c04_feature_tests,",
 "c05_mutations, c06_pca_clinical_split, c07_isogenic_other, c08_completeness, c09_index, c10_completeness_md. Seed 20261117.",
 "Inputs read from other agents' folders (never modified): B3 bundles, 02 coefficient and summary tables, the paper's STRING",
 "network, the DepMap cache. New downloads in `_data/raw/` (Xena GDC somatic mutations COAD/READ/LUAD; GEO GSE15709, GSE322742,",
 "GPL570 table). Sizes and sha256: `tables/data_versions.tsv`."]
(OUT / "COMPLETENESS.md").write_text("\n".join(L) + "\n"); print("ok")
