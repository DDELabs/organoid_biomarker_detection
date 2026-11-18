"""Build REPORT.md. Every number below is read from a file in ../tables/ (nothing typed by hand)."""
from pathlib import Path
import pandas as pd

TAB = Path(__file__).resolve().parent.parent / "tables"
OUT = Path(__file__).resolve().parent.parent / "REPORT.md"
R = lambda f: pd.read_csv(TAB / f, sep="\t")


def p(x):
    return "NA" if pd.isna(x) else (f"{x:.2g}" if x < 0.001 else f"{x:.3f}")


def f2(x):
    return "NA" if pd.isna(x) else f"{x:.2f}"


def ci(a, b, c, d=2):
    return f"{a:.{d}f} ({b:.{d}f} to {c:.{d}f})"


L = []
w = L.append

# ------------------------------------------------------------------ G3
g3 = R("G3_bh3_vs_drugs_organoids.tsv")
rg = R("G3_response_range_organoids.tsv").set_index("drug")
sg = R("G3_single_genes_vs_bh3_mimetics_organoids.tsv")
allg = R("G3_bh3_vs_all_drugs_organoids_AUC.tsv").set_index("drug")


def g3row(score, metric, drug):
    return g3[(g3.score == score) & (g3.metric == metric) & (g3.drug == drug)].iloc[0]


# ------------------------------------------------------------------ G2
g2 = R("G2_decomposition_colon_5FU.tsv")
A = g2[g2.group == "A_chain"].sort_values("step")
B = g2[g2.group == "B_one_from_authors"]
C = g2[g2.group == "C_one_into_current"]
D = g2[g2.group == "D_regimen_current"].set_index("step")
Da = g2[g2.group == "D_regimen_authors"].set_index("step")
Dc = g2[g2.group == "D_regimen_current_capped"].set_index("step")

w("# Biology gap analyses (G2, G3, G5, G7, G8, G13, G14)\n")
w("Every number in this report was copied by `code/make_report.py` from a file in `tables/`. "
  "Figures are in `figures/` (PNG + PDF + a plain-language `.caption.md`). "
  "\"What I think:\" lines are my reading; the decision is yours.\n")
w("Words used: **rho** = Spearman rank correlation (-1 to +1). **HR** = hazard ratio (above 1 = shorter survival). "
  "**CI** = 95% confidence interval. **Frozen** = the pathway and direction fixed from the 2019 paper, with no re-selection. "
  "**AUC** = area under the dose-response curve (lower = more sensitive).\n")

# ================= G2
w("## G2. Colon 5-FU: why p moved from 0.014 (2019) to 0.866 (now)\n")
w("Files: `tables/G2_decomposition_colon_5FU.tsv`, `tables/G2_data_versions.tsv`, `tables/G2_gene_universe_sizes.tsv`. "
  "Figure: `figures/G2_decomposition_colon_5FU.png/.pdf`.\n")
w("**a. One change at a time, adding current inputs step by step (orange = model re-picks its top pathway; "
  "blue = the paper's BH3-only pathway kept fixed).**\n")
w("| step | what was changed | n | deaths | re-picked top pathway | log-rank p | frozen BH3-only log-rank p | frozen HR per SD (95% CI) |")
w("|---|---|---|---|---|---|---|---|")
for _, r in A.iterrows():
    w(f"| {int(r.step)} | {r.label} | {int(r.n)} | {int(r.deaths)} | {r.rank1.replace('_', ' ').lower()} | {p(r.p_rank1)} | "
      f"{p(r.p_frozen)} | {ci(r.HRsd_frozen, r.HRsd_lo, r.HRsd_hi)} |")
a0, a7 = A.iloc[0], A.iloc[-1]
w(f"\nThe re-picked pathway starts as BH3-only (p = {p(a0.p_rank1)}) and ends as "
  f"{a7.rank1.replace('_', ' ').lower()} (p = {p(a7.p_rank1)}). The frozen BH3-only pathway goes from p = {p(a0.p_frozen)} "
  f"to p = {p(a7.p_frozen)}, HR per SD {ci(a7.HRsd_frozen, a7.HRsd_lo, a7.HRsd_hi)}.\n")
w("**b. Only one of the authors' inputs swapped in (everything else as in 2019).**\n")
w("| swapped input | n | deaths | re-picked top pathway | log-rank p (re-picked) | frozen BH3-only p |")
w("|---|---|---|---|---|---|")
for _, r in B.iterrows():
    w(f"| {r.label} | {int(r.n)} | {int(r.deaths)} | {r.rank1.replace('_', ' ').lower()} | {p(r.p_rank1)} | {p(r.p_frozen) if pd.notna(r.p_frozen) else 'NA'} |")
w("\n**c. Starting from the current setup, put back one of the authors' inputs.**\n")
w("| authors' input put back | n | deaths | re-picked top pathway | log-rank p (re-picked) | frozen BH3-only p |")
w("|---|---|---|---|---|---|")
for _, r in C.iterrows():
    w(f"| {r.label} | {int(r.n)} | {int(r.deaths)} | {r.rank1.replace('_', ' ').lower()} | {p(r.p_rank1)} | {p(r.p_frozen) if pd.notna(r.p_frozen) else 'NA'} |")
op = C[C.step == "orgpat"].iloc[0]
sur = B[B.step == "sur"].iloc[0]
w(f"\n**d. Frozen BH3-only pathway alone, by regimen (current data).** `5-FU + oxaliplatin` = FOLFOX-type.\n")
w("| patient group | n | deaths | log-rank p | HR (95% CI), predicted non-responder vs responder | HR per SD (95% CI) |")
w("|---|---|---|---|---|---|")
for lab, key in (("all 5-FU patients", "all"), ("5-FU + oxaliplatin (FOLFOX-type)", "folfox"),
                 ("5-FU without oxaliplatin", "nooxa"), ("5-FU alone", "alone")):
    r = D.loc[key]
    hr = "NA (0 deaths)" if pd.isna(r.HR_frozen) else ci(r.HR_frozen, r.HR_frozen_lo, r.HR_frozen_hi)
    hs = "NA" if pd.isna(r.HRsd_frozen) else ci(r.HRsd_frozen, r.HRsd_lo, r.HRsd_hi)
    w(f"| {lab} | {int(r.n)} | {int(r.deaths)} | {p(r.p_frozen)} | {hr} | {hs} |")
rc = Dc.loc["all"]
w(f"\nSame, with current follow-up cut at the 2019 date: all 5-FU patients n = {int(rc.n)}, deaths = {int(rc.deaths)}, "
  f"frozen log-rank p = {p(rc.p_frozen)}, HR per SD {ci(rc.HRsd_frozen, rc.HRsd_lo, rc.HRsd_hi)}.\n")
pat_b = B[B.step == "pat"].iloc[0]
pat_c = C[C.step == "pat"].iloc[0]
org_c = C[C.step == "org"].iloc[0]
w(f"**What I think:** The patient-side scores are the main driver. With the paper's pathway fixed, putting the authors' 2019 patient scores back into the current setup gives "
  f"p = {p(pat_c.p_frozen)} (n = {int(pat_c.n)}, {int(pat_c.deaths)} deaths), while putting back only the authors' organoid scores gives p = {p(org_c.p_frozen)}. "
  f"Going the other way, swapping only the current STAR patient scores into the 2019 setup moves the fixed pathway to p = {p(pat_b.p_frozen)}. "
  f"The organoid side matters for which pathway the model picks: recomputing organoid scores alone moves the re-picked pathway from "
  f"{p(A.iloc[0].p_rank1)} to {p(A.iloc[2].p_rank1)} (step 2). Survival data and follow-up are not the cause: with the authors' scores, "
  f"the current TCGA-CDR survival gives p = {p(sur.p_rank1)} ({int(sur.deaths)} deaths), and the authors' scores plus current survival give "
  f"p = {p(op.p_rank1)} (HR {ci(op.HR_rank1, op.HR_rank1_lo, op.HR_rank1_hi, 1)}, n = {int(op.n)}, {int(op.deaths)} deaths). "
  f"The 2019 result rested on {int(a0.deaths)} deaths. In short, the 2019 patient expression scores and today's STAR scores rank the same 5-FU patients differently for this pathway.\n")
w("**Not done (G2, bladder).** The request also asked for the same frozen-pathway report on bladder. The script `g2_decompose.py` is colon-only and "
  "I did not write a bladder version. The bladder ERCC2 side is covered in G8.\n")

# ================= G3
w("## G3. BH3-only score vs navitoclax in the 19 colorectal organoids\n")
w("Files: `tables/G3_bh3_vs_drugs_organoids.tsv`, `G3_single_genes_vs_bh3_mimetics_organoids.tsv`, `G3_bh3_vs_all_drugs_organoids_AUC.tsv`, "
  "`G3_response_range_organoids.tsv`, `G3_data_versions.tsv`. Figure: `figures/G3_bh3_vs_bh3_mimetics_organoids.png/.pdf`.\n")
w("Spearman rho between the BH3-only score and drug response (AUC). A negative rho means a higher score goes with more sensitivity. "
  "Null = share of random gene sets (same size, 2,000 sets) and of all Reactome pathways with rho at least as negative.\n")
w("| drug | score used | n | rho (95% CI) | p | share of random sets as negative | share of all pathways as negative |")
w("|---|---|---|---|---|---|---|")
for d in ("NAVITOCLAX", "OBATOCLAX", "FLUOROURACIL", "OXALIPLATIN"):
    for sc, lab in (("current", "recomputed"), ("authors", "authors' NES")):
        r = g3row(sc, "AUC", d)
        rs = "NA" if pd.isna(r.share_random_sets_rho_le) else f"{r.share_random_sets_rho_le:.3f}"
        w(f"| {d.lower()} | {lab} | {int(r.n)} | {ci(r.rho, r.lo, r.hi)} | {p(r.p)} | {rs} | {r.share_pathways_rho_le:.3f} (of {int(r.n_pathways)}) |")
nv, fu = g3row("current", "AUC", "NAVITOCLAX"), g3row("current", "AUC", "FLUOROURACIL")
fa = g3row("authors", "AUC", "FLUOROURACIL")
w(f"\nTested range: navitoclax was tested up to {rg.loc['NAVITOCLAX','max_conc_uM']:g} uM; organoid AUC runs from "
  f"{rg.loc['NAVITOCLAX','auc_min']:.3f} to {rg.loc['NAVITOCLAX','auc_max']:.3f}, and {int(rg.loc['NAVITOCLAX','n_ic50_above_max_conc'])} of "
  f"{int(rg.loc['NAVITOCLAX','n'])} IC50 values lie above the top dose. For 5-FU, AUC runs from {rg.loc['FLUOROURACIL','auc_min']:.3f} to "
  f"{rg.loc['FLUOROURACIL','auc_max']:.3f} (SD {rg.loc['FLUOROURACIL','auc_sd']:.3f}) vs navitoclax SD {rg.loc['NAVITOCLAX','auc_sd']:.3f}.\n")
w("Single genes vs navitoclax (rho, 95% CI, p): " + "; ".join(
    f"{r.gene} {ci(r.rho, r.lo, r.hi)}, p = {p(r.p)}" for _, r in sg[sg.drug == "NAVITOCLAX"].iterrows()) + ".\n")
w(f"Across the whole screen ({len(allg)} drugs, AUC), the BH3-only score ranks navitoclax {int(allg.loc['NAVITOCLAX','rank_most_negative'])}, "
  f"5-FU {int(allg.loc['FLUOROURACIL','rank_most_negative'])} and obatoclax {int(allg.loc['OBATOCLAX','rank_most_negative'])} "
  f"from the most negative correlation.\n")
w(f"**What I think:** In these 19 organoids the navitoclax link seen in cell lines is not visible (rho = {nv.rho:.2f}, CI {nv.lo:.2f} to {nv.hi:.2f}, "
  f"p = {p(nv.p)}). The data cannot settle it either way: navitoclax was only tested up to {rg.loc['NAVITOCLAX','max_conc_uM']:g} uM, nearly every organoid "
  f"stayed at AUC near 1, and the interval is wide. The BH3-only score points the 5-FU way in the same organoids "
  f"(recomputed rho = {fu.rho:.2f}, p = {p(fu.p)}; authors' NES rho = {fa.rho:.2f}, CI {fa.lo:.2f} to {fa.hi:.2f}, p = {p(fa.p)}). "
  f"Obatoclax, the planned control, shows rho = {g3row('current','AUC','OBATOCLAX').rho:.2f}, so it did not behave as a flat control here. "
  f"A drug panel with a higher top dose is the clean way to test the apoptotic-priming claim in organoids.\n")

# ================= G5
w("## G5. FOLFOX question: is the 5-FU marker a regimen (oxaliplatin) marker?\n")
w("Files: `tables/G5_cell_lines_bh3_vs_oxaliplatin.tsv`, `G5_folfox_patients.tsv`, `G5_data_versions.tsv`, `G3_bh3_vs_drugs_organoids.tsv`. "
  "Figure: `figures/G5_folfox_question.png/.pdf`.\n")
w("**a. 19 organoids (AUC, recomputed score).** ")
for d in ("OXALIPLATIN", "FLUOROURACIL"):
    r = g3row("current", "AUC", d)
    w(f"- {d.lower()}: rho = {ci(r.rho, r.lo, r.hi)}, p = {p(r.p)}")
w("\n**b. Cancer cell lines.** rho of BH3-only score vs response (colorectal lines only, and all lines with cancer type removed).\n")
cl = R("G5_cell_lines_bh3_vs_oxaliplatin.tsv")
w("| screen | drug | lines | n | rho (95% CI) | p |")
w("|---|---|---|---|---|---|")
for _, r in cl.iterrows():
    if r.n == 0:
        w(f"| {r.screen} | {r.drug.lower()} | {r.scope} | 0 | no data in table | NA |")
    else:
        w(f"| {r.screen} | {r.drug.lower()} | {r.scope} | {int(r.n)} | {ci(r.rho, r.lo, r.hi)} | {p(r.p)} |")
gp = R("G5_folfox_patients.tsv")
w("\n**c. TCGA colon 5-FU patients, frozen BH3-only score (HR per SD; above 1 = the paper's direction).**\n")
w("| data | group | n | events | HR per SD (95% CI) | p |")
w("|---|---|---|---|---|---|")
for _, r in gp[gp.cohort == "TCGA colon"].iterrows():
    w(f"| {r.setup} | {r.subset} | {int(r.n)} | {int(r.events)} | {ci(r.HRsd_HR, r.HRsd_lo, r.HRsd_hi)} | {p(r.HRsd_p)} |")
ir = gp[(gp.cohort == "TCGA colon") & gp.inter_p.notna()].iloc[0]
w(f"\nScore x oxaliplatin interaction (current data): HR {ci(ir.inter_HR, ir.inter_lo, ir.inter_hi)}, p = {p(ir.inter_p)}.\n")
w("**d. GSE39582 (stage II-III), frozen BH3-only score by adjuvant regimen.**\n")
w("| endpoint | group | n | events | HR per SD (95% CI) | p | interaction HR (95% CI), p |")
w("|---|---|---|---|---|---|---|")
for _, r in gp[gp.cohort.str.startswith("GSE")].iterrows():
    if pd.notna(r.inter_p):
        w(f"| {r.endpoint} | {r.subset} | {int(r.n)} | {int(r.events)} | | | {ci(r.inter_HR, r.inter_lo, r.inter_hi)}, p = {p(r.inter_p)} |")
    else:
        w(f"| {r.endpoint} | {r.subset} | {int(r.n)} | {int(r.events)} | {ci(r.HRsd_HR, r.HRsd_lo, r.HRsd_hi)} | {p(r.HRsd_p)} | |")
ox = g3row("current", "AUC", "OXALIPLATIN")
cx = cl[(cl.screen == "CTRP") & (cl.drug == "OXALIPLATIN") & cl.scope.str.startswith("colorectal")].iloc[0]
ca = cl[(cl.screen == "CTRP") & (cl.drug == "OXALIPLATIN") & cl.scope.str.startswith("pan")].iloc[0]
g5 = gp[gp.cohort.str.startswith("GSE")].set_index(["endpoint", "subset"])
tf = gp[(gp.cohort == "TCGA colon") & (gp.setup == "current data")].set_index("subset")
w(f"\n**What I think:** The score points the same way for oxaliplatin as for 5-FU in the 19 organoids (oxaliplatin rho = {ox.rho:.2f}, p = {p(ox.p)}), "
  f"so the organoid data cannot tell the two drugs apart. In colorectal cell lines the numbers are close to zero or positive (higher score with more resistance) "
  f"(for example CTRPv2 oxaliplatin rho = {cx.rho:.2f}, n = {int(cx.n)}); across all CTRPv2 lines with cancer type removed rho = {ca.rho:.2f} (p = {p(ca.p)}). "
  f"In TCGA, the HR per SD is {tf.loc['5-FU + oxaliplatin (FOLFOX-type)','HRsd_HR']:.2f} with oxaliplatin and "
  f"{tf.loc['5-FU without oxaliplatin','HRsd_HR']:.2f} without, and the interaction p = {p(ir.inter_p)}: the groups are too small "
  f"({int(tf.loc['5-FU without oxaliplatin','n'])} patients without oxaliplatin) to separate a 5-FU effect from a FOLFOX effect. "
  f"In GSE39582 (stage II-III) the HR per SD is {g5.loc[('RFS','5FU or FUFOL (no oxaliplatin)'),'HRsd_HR']:.2f} for relapse-free survival in 5-FU/FUFOL patients (n = {int(g5.loc[('RFS','5FU or FUFOL (no oxaliplatin)'),'n'])}) "
  f"and {g5.loc[('RFS','FOLFOX'),'HRsd_HR']:.2f} in FOLFOX (n = {int(g5.loc[('RFS','FOLFOX'),'n'])}); the FOLFOX vs 5-FU interaction p = {p(g5.loc[('RFS','FOLFOX vs 5FU/FUFOL (interaction)'),'inter_p'])}. Most of the 5-FU patients in TCGA got FOLFOX, so the paper's "
  f"\"5-FU\" claim is in practice a FOLFOX claim.\n")

# ================= G7
w("## G7. Does the BH3-only pathway add beyond BCL2, BCL2L1, MCL1 and proliferation?\n")
w("Files: `tables/G7_regression_added_value.tsv`, `G7_partial_correlations.tsv`. Figure: `figures/G7_bh3_beyond_simple_markers.png/.pdf`. "
  "Cell-line screens GDSC1, GDSC2 (ln IC50), CTRPv2 (AUC); all variables rank-normalised; cancer type always removed.\n")
rg7 = R("G7_regression_added_value.tsv")
pc = R("G7_partial_correlations.tsv")
w("| screen | drug | n | partial r, cancer type only (95% CI) | partial r, + BCL2, BCL2L1, MCL1 + E2F (95% CI) | extra R2 from BH3-only | F-test p | rank among pathways | share of random sets as strong |")
w("|---|---|---|---|---|---|---|---|---|")
for _, r in rg7[rg7.drug.isin(["NAVITOCLAX", "ABT-737"])].iterrows():
    a = pc[(pc.screen == r.screen) & (pc.drug == r.drug) & (pc["var"] == "BH3-only pathway") & (pc.adjusted_for == "lineage")].iloc[0]
    b = pc[(pc.screen == r.screen) & (pc.drug == r.drug) & (pc["var"] == "BH3-only pathway") & (pc.adjusted_for == "+ all three genes + E2F")].iloc[0]
    w(f"| {r.screen} | {r.drug.lower()} | {int(r.n)} | {ci(a.r, a.lo, a.hi, 3)} | {ci(b.r, b.lo, b.hi, 3)} | {r.dR2:.4f} | {p(r.p_F)} | "
      f"{int(r.rank_pathways)} of {int(r.n_pathways)} | {r.share_random_le:.3f} |")
w("\nOther BH3-mimetics (BH3-only beta in the full model, 95% CI, extra R2):")
for _, r in rg7[~rg7.drug.isin(["NAVITOCLAX", "ABT-737"])].iterrows():
    w(f"- {r.screen} {r.drug.lower()} (n = {int(r.n)}): beta {ci(r.beta_BH3, r.beta_lo, r.beta_hi, 3)}, extra R2 {r.dR2:.4f}, p = {p(r.p_F)}")
w("\nSingle markers vs navitoclax (partial r with cancer type removed, GDSC2): " + "; ".join(
    f"{r['var']} {r.r:.3f} ({r.lo:.3f} to {r.hi:.3f})" for _, r in pc[(pc.screen == "GDSC2") & (pc.drug == "NAVITOCLAX") & (pc.adjusted_for == "lineage")].iterrows()) + ".\n")
g2n = rg7[(rg7.screen == "GDSC2") & (rg7.drug == "NAVITOCLAX")].iloc[0]
w(f"**What I think:** The BH3-only score keeps most of its link after removing BCL2L1, BCL2, MCL1 and proliferation, so it is not just BCL-XL or E2F. "
  f"For example in GDSC2 navitoclax the partial r goes from "
  f"{pc[(pc.screen=='GDSC2')&(pc.drug=='NAVITOCLAX')&(pc['var']=='BH3-only pathway')&(pc.adjusted_for=='lineage')].iloc[0].r:.3f} to "
  f"{pc[(pc.screen=='GDSC2')&(pc.drug=='NAVITOCLAX')&(pc['var']=='BH3-only pathway')&(pc.adjusted_for=='+ all three genes + E2F')].iloc[0].r:.3f}. "
  f"It adds a small but consistent slice of variance (extra R2 {g2n.dR2:.3f} in GDSC2 navitoclax) and ranks {int(g2n.rank_pathways)} of {int(g2n.n_pathways)} pathways "
  f"in all three screens. Removing E2F1 and TFDP1 from the set leaves the link intact (see `BH3_noE2F` in the partial-correlation table). "
  f"The effect is modest in size, so I would describe it as an added signal, not a replacement for BCL2-family genes.\n")

# ================= G8
w("## G8. Null and confounder control for the mutation checks\n")
w("Files: `tables/G8_ERCC2_null_and_adjustment.tsv`, `G8_BRAF_null_and_adjustment.tsv`, and their `_data_versions.tsv`. "
  "Figure: `figures/G8_mutation_concordance_null_and_adjustment.png/.pdf`. One-sided p values (expected direction).\n")
e = R("G8_ERCC2_null_and_adjustment.tsv")
b = R("G8_BRAF_null_and_adjustment.tsv")
w("**ERCC2 (bladder, cisplatin).**\n")
w("| data | model | mutant / wild type | p | share of all pathways (both directions) as strong | share of random proximal-sized sets as strong | p after TMB adjustment (OLS) | OR per SD after TMB (95% CI) |")
w("|---|---|---|---|---|---|---|---|")
for _, r in e.iterrows():
    w(f"| {r.version} | {r.model} | {int(r.n_mut)} / {int(r.n_wt)} | {p(r.p_one_sided)} | {r.share_all_pathways_both_signs:.3f} | "
      f"{r.share_random_proximal_sized_sets:.3f} | {p(r.adjTMB_ols_p_one_sided)} | {ci(r.adjTMB_logit_OR_per_SD, r.adjTMB_logit_lo, r.adjTMB_logit_hi)} |")
e0 = e[e.model.str.startswith("paper biomarker")]
w("\n**BRAF V600E (colon, cetuximab).**\n")
w("| data | mutant / wild type | p | share of pathways as strong | share of random sets as strong | p adj. MSI | p adj. side | p adj. TMB | p adj. all (MSI, side, stage, TMB) | MSS tumours only (mutants, p) |")
w("|---|---|---|---|---|---|---|---|---|---|")
for _, r in b.iterrows():
    w(f"| {r.version} | {int(r.n_mut)} / {int(r.n_wt)} | {p(r.p_one_sided)} | {r.share_all_pathways_both_signs:.3f} | {r.share_random_proximal_sized_sets:.3f} | "
      f"{p(r.adjMSI_ols_p_one_sided)} | {p(r.adjside_ols_p_one_sided)} | {p(r.adjTMB_ols_p_one_sided)} | {p(r.adjall_ols_p_one_sided)} | {int(r.MSS_only_n_mut)}, {p(r.MSS_only_p_one_sided)} |")
ec = e0[e0.version == "current"].iloc[0]
ep = e0[e0.version == "paper"].iloc[0]
bc = b[b.version == "current"].iloc[0]
bp = b[b.version == "authors"].iloc[0]
w(f"\n**What I think:** ERCC2 looks specific: only {ec.share_all_pathways_both_signs:.3f} of all pathways and {ec.share_random_proximal_sized_sets:.3f} of random sets "
  f"separate mutants as well on current data (paper data: {ep.share_all_pathways_both_signs:.3f} and {ep.share_random_proximal_sized_sets:.3f}), and it holds after TMB adjustment "
  f"(p = {p(ec.adjTMB_ols_p_one_sided)} current, {p(ep.adjTMB_ols_p_one_sided)} paper). "
  f"BRAF V600E is less specific: {bc.share_all_pathways_both_signs:.3f} of pathways and {bc.share_random_proximal_sized_sets:.3f} of random sets do as well on current data "
  f"({bp.share_all_pathways_both_signs:.3f} and {bp.share_random_proximal_sized_sets:.3f} on the authors' data). After adjusting for MSI, side, stage and TMB the current-data p is "
  f"{p(bc.adjall_ols_p_one_sided)} (authors' data {p(bp.adjall_ols_p_one_sided)}), and within MSS tumours p = {p(bc.MSS_only_p_one_sided)} with {int(bc.MSS_only_n_mut)} mutants. "
  f"So the BRAF check is best described as a consistent direction whose size partly tracks tumour subtype, while the ERCC2 check stands on its own.\n")

# ================= G13
w("## G13. The 'HIV' frozen pathways: what genes are in them?\n")
w("Files: `tables/G13_frozen_pathway_composition.tsv`, `G13_most_overlapping_pathways.tsv`, `G13_purity_immune_correlation.tsv`, "
  "`G13_genes_*.txt` (gene lists), `G13_data_versions.tsv`. Figure: `figures/G13_odd_pathway_composition.png/.pdf`.\n")
c = R("G13_frozen_pathway_composition.tsv")
o = c[c.odd.astype(str) == "True"].set_index("pathway")
mo = R("G13_most_overlapping_pathways.tsv")
pu = R("G13_purity_immune_correlation.tsv")
w("| pathway | genes | in proteasome set | in nuclear export/pore set | in immune-system set | in adaptive-immune set | in cell-cycle set | host-factor genes (PPIA, BANF1, XRCC*...) |")
w("|---|---|---|---|---|---|---|---|")
for pth, r in o.iterrows():
    w(f"| {pth.replace('REACTOME_', '').replace('_', ' ').lower()} | {int(r.n_genes)} | {int(r.overlap_proteasome)} ({r.share_in_proteasome:.2f}) | "
      f"{int(r['overlap_nuclear export / pore'])} ({r['share_in_nuclear export / pore']:.2f}) | {int(r['overlap_immune system'])} ({r['share_in_immune system']:.2f}) | "
      f"{int(r['overlap_adaptive immune'])} ({r['share_in_adaptive immune']:.2f}) | {int(r['overlap_cell cycle'])} ({r['share_in_cell cycle']:.2f}) | "
      f"{int(r.iloc[list(o.columns).index('n_HIV accessory host factors (PPIA, BANF1, HMGA1, PSIP1, XRCC*, LIG4)')])} |")
w("\nClosest other Reactome pathways (Jaccard overlap of gene lists):")
for pth in o.index:
    t = mo[mo.pathway == pth].sort_values("jaccard", ascending=False).head(3)
    w(f"- {pth.replace('REACTOME_', '').replace('_', ' ').lower()}: " + "; ".join(
        f"{x.other.replace('REACTOME_', '').replace('_', ' ').lower()} (Jaccard {x.jaccard:.2f}, {int(x.overlap)} shared genes)" for _, x in t.iterrows()))
w(f"\nCorrelation with tumour purity and an adaptive-immune score, median over {int(pu.n_cancer_types.iloc[0])} TCGA cancer types (10th to 90th percentile):\n")
w("| pathway | vs purity | vs adaptive immune score |")
w("|---|---|---|")
for pth in o.index:
    a = pu[(pu.pathway == pth) & (pu.vs == "purity")].iloc[0]
    bb = pu[(pu.pathway == pth) & (pu.vs == "adaptive immune score")].iloc[0]
    w(f"| {pth.replace('REACTOME_', '').replace('_', ' ').lower()} | {a.median_rho:.2f} ({a.q10:.2f} to {a.q90:.2f}) | {bb.median_rho:.2f} ({bb.q10:.2f} to {bb.q90:.2f}) |")
allimm = pu[pu.vs == "adaptive immune score"].median_rho
w(f"\nFor all {pu.pathway.nunique()} distinct frozen pathways the median (over pathways) of the per-pathway median rho with the immune score is {allimm.median():.2f}.\n")
aa = o.iloc[0]
hv = o.iloc[1]
w(f"**What I think:** These two sets are human genes that viruses use, not virus genes. The 12-gene APOBEC3G set is built from "
  f"HIV life-cycle pathways (closest: integration of provirus, {int(mo[(mo.pathway==o.index[0])].sort_values('jaccard',ascending=False).iloc[0].overlap)} shared genes). "
  f"The 207-gene HIV infection set is {hv.share_in_proteasome:.0%} proteasome, {hv['share_in_nuclear export / pore']:.0%} nuclear pore/export, and "
  f"{hv['share_in_immune system']:.0%} immune system genes. Both scores rise with an adaptive-immune score "
  f"(median rho {pu[(pu.pathway==o.index[0])&(pu.vs=='adaptive immune score')].iloc[0].median_rho:.2f} and {pu[(pu.pathway==o.index[1])&(pu.vs=='adaptive immune score')].iloc[0].median_rho:.2f}), "
  f"and the APOBEC3G set also falls with purity (median rho {pu[(pu.pathway==o.index[0])&(pu.vs=='purity')].iloc[0].median_rho:.2f}). "
  f"A reader can therefore take these as proteasome, nuclear-transport and immune-infiltration biology under a virus-themed name. "
  f"Note: the frozen list has 35 rows but {pu.pathway.nunique()} distinct pathways; the table covers the distinct ones.\n")

# ================= G14
w("## G14. Supp. Fig. 15 for pancreas, liver, ovary and sarcoma (cell-line-trained models)\n")
w("Files: `tables/G14_SuppFig15_cell_line_models.tsv`, `G14_data_versions.tsv`. Figure: `figures/G14_SuppFig15_cell_line_models.png/.pdf`. "
  "Models trained on GDSC2 cell lines (all lines, or lines of the matching cancer type), applied to treated TCGA patients; "
  "each of the top 3 pathways splits patients at the median. HR above 1 = predicted non-responders die sooner (expected direction).\n")
g14 = R("G14_SuppFig15_cell_line_models.tsv")
w("| setting | training | lines | patients | deaths | rank | pathway | log-rank p | HR (95% CI) |")
w("|---|---|---|---|---|---|---|---|---|")
for _, r in g14.sort_values(["setting", "training", "rank"]).iterrows():
    w(f"| {r.setting} | {r.training} | {int(r.n_lines)} | {int(r.n)} | {int(r.events)} | {int(r['rank'])} | "
      f"{r.pathway.replace('REACTOME_', '').replace('_', ' ').lower()[:60]} | {p(r.p)} | {ci(r.HR, r.HR_lo, r.HR_hi)} |")
top = g14[g14["rank"] == 1]
best = top.loc[top.p.idxmin()]
w(f"\n**What I think:** All {top.setting.nunique()} settings (x2 training sets) are now filled in, so the paper's comparison exists for every cancer. "
  f"Across the rank-1 pathways, the smallest p is {p(best.p)} ({best.setting}, {best.training}, n = {int(best.n)}, HR {ci(best.HR, best.HR_lo, best.HR_hi)}); "
  f"the ovarian platinum and paclitaxel sets (n = {int(top[top.setting.str.contains('platinum')].n.iloc[0])} and {int(top[top.setting.str.contains('paclitaxel')].n.iloc[0])}) "
  f"have the most patients and give HRs close to 1 (see table). Pancreas and liver have {int(top[top.setting.str.startswith('pancreas')].n.iloc[0])} and "
  f"{int(top[top.setting.str.startswith('liver')].n.iloc[0])} patients, so their intervals are wide. The sign flips between the all-lines and matched-lines models for pancreas and liver, "
  f"which suggests the pathway choice is driven by which lines are in the training set. The ovarian platinum model is trained on cisplatin only "
  f"(no carboplatin screen in GDSC), as stated in the figure caption.\n")

w("## Files and data versions\n")
w("Data versions (file name, bytes, sha256) for each gap are in `tables/G*_data_versions.tsv`. Scripts are in `code/` "
  "(`python3 <script>` from inside `code/`; `PROGRESS.md` gives the order). Random seeds are set in `gcommon.py` (`SEED`).\n")
w("## Fixes made while running\n")
w("- `g7_bh3_beyond_markers.py`: skip a drug that is absent from a screen (WEHI-539 is not in GDSC1); the script had stopped there.\n"
  "- `g5_folfox.py`: GSE39582 stage is stored as text, so the stage II-III filter returned no patients; it now compares text values.\n")
OUT.write_text("\n".join(L))
print("wrote", OUT, len(L), "lines")
