"""G11: append reconciliation rows (numbers that differ between agents) to MASTER_NUMBERS.tsv. Values are read
from the named source tables; the reason column is the auditor's finding (see AUDIT.md)."""
from pathlib import Path
import pandas as pd
NAR = Path(__file__).resolve().parents[2]
M = pd.read_csv(NAR / "12_audit/tables/MASTER_NUMBERS.tsv", sep="\t")
M = M[~M.id.str.startswith("G11_")]
t = lambda p: pd.read_csv(NAR / p, sep="\t")
braf = t("02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv")
b02 = braf[braf.label == "organoids, current data, all"].iloc[0]
f10 = t("10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv")
b10 = f10[f10.model.str.contains("vdW") & (f10.mutation == "BRAF V600E")].iloc[0]
gi = t("02_organoid_only/tables/cross_set_summary.tsv")
c5 = gi[(gi.cohort == "COAD") & (gi.model_drug == "FLUOROURACIL")].iloc[0]
br = gi[gi.folder == "LIHC_broutier"].iloc[0]
vb = t("10_completeness/tables/verify_headline_biomarkers.tsv").set_index("key").loc["LIHC_broutier_SOR"]
ot = t("02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv")
bl = ot[(ot.set == "BLCA_lee2018_paper") & (ot.drug == "CISPLATIN")].iloc[0]
cx = ot[(ot.set == "PDX_gao2015_CRC") & (ot.drug == "CETUXIMAB")].iloc[0]
cg = gi[(gi.cohort == "COAD") & (gi.model_drug == "CETUXIMAB")].iloc[0]
s2 = t("09_registered_test/tables/stage2_controlled.tsv")
g72 = s2[(s2.drug == "OXALIPLATIN") & (s2.cohort == "GSE72970") & (s2.endpoint == "PFS")].iloc[0]
s39 = s2[(s2.drug == "FLUOROURACIL") & (s2.cohort == "GSE39582") & (s2.endpoint == "RFS")].iloc[0]
# after the audit fixes (nar/15_fixes): read the re-run 03 / 05 / 09-sensitivity numbers
I3 = t("03_organoid_plus_patient/tables/interaction_all.tsv")
I3 = I3[(I3.features == "shortlist") & (I3.model == "joint")]
i3_39 = I3[(I3.unit == "CRC_5FU") & (I3.cohort == "GSE39582")].iloc[0]
i3_ox = I3[(I3.unit == "CRC_OXALIPLATIN") & (I3.cohort == "GSE72970") & I3.outcome.str.startswith("survival")].iloc[0]
i3_ir = I3[(I3.unit == "CRC_IRINOTECAN") & (I3.cohort == "GSE72970") & I3.outcome.str.startswith("survival")].iloc[0]
B1 = t("05_A3_B1/tables/B1_external_results.tsv")
b39 = B1[(B1.cohort == "GSE39582") & B1.comparison.str.startswith("5-FU-based") & (B1.endpoint == "RFS")].iloc[0]
b72 = B1[(B1.cohort == "GSE72970") & (B1.comparison == "FOLFOX vs FOLFIRI") & (B1.endpoint == "PFS")].iloc[0]
h39 = t("09_registered_test/tables/stage2_sensitivity_GSE39582_harmonised.tsv")
h39 = h39[h39.definition.str.startswith("shared") & (h39.endpoint == "RFS")].iloc[0]
rows = [
 ("G11_BRAF", "Colon cetuximab model vs BRAF V600E, current data (one-sided MWU)",
  f"COAD only: p={b02.p_one_sided:.5f} ({b02.n_mut} vs {b02.n_wt}) [use in manuscript, matches paper design]; COAD+READ: p={b10.p_one_sided:.5f} ({b10.n_mut} vs {b10.n_wt})",
  "nar/02_organoid_only/COAD_vandewetering/tables/Fig5a_cetuximab_BRAF_V600E.tsv; nar/10_completeness/tables/Fig5_mutation_concordance_other_cancers.tsv",
  "02 = TCGA-COAD only, GDC masked MAF, 02 ssGSEA scores; 10 = COAD+READ pooled, Xena wxs calls, B3 bundle scores. Both correct; different cohort. Quote 0.00072 and label 0.00048 as COAD+READ."),
 ("G11_COAD5FU_current", "Colon 5-FU, paper method re-run on current GDC data",
  f"p={c5.treated_p:.3f}, n={c5.treated_n}, deaths={c5.treated_events}, rank-1 = purine ribonucleoside monophosphate biosynthesis",
  "nar/02_organoid_only/tables/cross_set_summary.tsv",
  "FIXED (I11, nar/15_fixes): docs/PROJECT_SUMMARY.md now quotes 0.866 with this source; the earlier robust/ value 0.62 is marked superseded."),
 ("G11_Broutier", "Liver Broutier organoids, sorafenib, TCGA-LIHC treated",
  f"p={br.treated_p:.4f} (authoritative); 10_completeness recompute p={vb.treated_p:.4f}",
  "nar/02_organoid_only/tables/cross_set_summary.tsv; nar/10_completeness/tables/verify_headline_biomarkers.tsv",
  "10 used the LIHC_ji bundle (TCGA-LIHC ssGSEA on the Ji gene universe) because B3 built no Broutier bundle. 0.139 uses Broutier's own gene universe and was rebuilt identically by 05 (nar/05_A3_B1/_data/a3_check_gi.tsv). Use 0.139; 07/10 should either build a Broutier bundle or state Broutier is excluded."),
 ("G11_BLCA_n", "Bladder cisplatin-treated TCGA patients, paper reproduction",
  f"ours n={bl.treated_n} (deaths {bl.treated_events}), untreated n={bl.untreated_n}; paper n=77, untreated 294",
  "nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv",
  "Paper patient lists are not public; 80 cisplatin records parsed from the biotab, 79 with RNA-seq and survival. Report 79 and state the paper's 77; 03 (TCGA-BLCA cisplatin vs none) uses 68 treated with a different survival source."),
 ("G11_cetuximab_n", "TCGA-COAD cetuximab-treated patients",
  f"02 GI n={cg.treated_n}/{cg.treated_events}; 02 OTHER (Gao PDX) n={cx.treated_n}/{cx.treated_events}",
  "nar/02_organoid_only/tables/cross_set_summary.tsv; nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv",
  "Two survival sources (TCGA-CDR vs Xena GDC survival) and two parsers. Harmonise on TCGA-CDR."),
 ("G11_GSE72970", "GSE72970 FOLFOX vs FOLFIRI arm sizes",
  f"09: FOLFOX {g72.n_drug_arm} vs FOLFIRI {g72.n_comparator}; 03 (after fix): oxaliplatin unit {i3_ox.n_treated} vs {i3_ox.n_comparison}, irinotecan unit {i3_ir.n_treated} vs {i3_ir.n_comparison}; 05 B1 (after fix): {int(b72.n)} vs {int(b72.comp_n)}",
  "nar/09_registered_test/tables/stage2_controlled.tsv; nar/03_organoid_plus_patient/tables/interaction_all.tsv; nar/05_A3_B1/tables/B1_external_results.tsv",
  "FIXED (I02, nar/15_fixes): FOLFIRINOX (4 patients) excluded from both arms in 03 and 05, as in 09 (before: 03/05 40 vs 84, 03 irinotecan 88 vs 36). 09's comparator has 1 fewer patient (XELIRI+bevacizumab not in its FOLFIRI list)."),
 ("G11_GSE39582", "GSE39582 '5-FU-treated' group",
  f"shared definition (nar/15_fixes/definitions/GSE39582_treatment.tsv): 5FU|FUFOL|FOLFOX, stage II-III = 139 treated (138 with RFS) vs 264 untreated stage II-III; 03 (after fix) {i3_39.n_treated} vs {i3_39.n_comparison} with RFS; 05 B1 (after fix) {int(b39.n)} vs {int(b39.comp_n)}; 09 primary (frozen) {s39.n_drug_arm} vs {s39.n_comparator}, 09 harmonised sensitivity {h39.n_drug_arm} vs {h39.n_comparator} (HR_int {h39.HR_int:.2f}, one-sided p={h39.p_one_sided:.3f})",
  "nar/15_fixes/definitions/GSE39582_treatment.tsv; nar/03_organoid_plus_patient/tables/interaction_all.tsv; nar/05_A3_B1/tables/B1_external_results.tsv; nar/09_registered_test/tables/stage2_sensitivity_GSE39582_harmonised.tsv",
  "FIXED (I01, nar/15_fixes): 03 and 05 re-run with the shared definition (before: 03 200 vs 288, 05 202 vs 294). 09 primary is frozen and keeps its definition; the shared definition is its sensitivity analysis. Small n differences between folders come from patients without RFS or with missing expression."),
]
add = pd.DataFrame([dict(id=a, description=b, authoritative_value=c, source_table=d, note=e, expected_text="", reports_quoting="", status="reconciled (G11)") for a, b, c, d, e in rows])
pd.concat([M, add]).to_csv(NAR / "12_audit/tables/MASTER_NUMBERS.tsv", sep="\t", index=False)
print(add[["id", "authoritative_value"]].to_string())
