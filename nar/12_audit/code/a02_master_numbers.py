"""Build nar/12_audit/tables/MASTER_NUMBERS.tsv: one authoritative value per headline number, read directly
from the source table that the owning agent's code wrote, plus where the number is quoted in the reports and
whether the quoted text matches. Read-only on every other folder. Run: python3 a02_master_numbers.py"""
import re
from pathlib import Path

import pandas as pd

NAR = Path(__file__).resolve().parents[2]
OUT = NAR / "12_audit" / "tables"
REPORTS = {
    "02_GI": NAR / "02_organoid_only/REPORT_GI.md",
    "02_OTHER": NAR / "02_organoid_only/REPORT_BLCA_OV_LUAD_OTHER.md",
    "03": NAR / "03_organoid_plus_patient/REPORT.md",
    "05": NAR / "05_A3_B1/REPORT.md",
    "08": NAR / "08_depmap/REPORT.md",
    "09": NAR / "09_registered_test/REPORT.md",
    "10": NAR / "10_completeness/COMPLETENESS.md",
}
TXT = {k: v.read_text() for k, v in REPORTS.items()}


def t(rel):
    return pd.read_csv(NAR / rel, sep="\t")


rows = []


def add(nid, desc, value, src, fmt, quoted_as=None):
    """fmt: list of strings as they should appear in reports (any one counts as a match)."""
    found = [k for k, s in TXT.items() if any(f in s for f in fmt)]
    rows.append(dict(id=nid, description=desc, authoritative_value=value, source_table=src,
                     expected_text=" | ".join(fmt), reports_quoting=";".join(found) or "-",
                     note=quoted_as or ""))


gi = t("02_organoid_only/tables/cross_set_summary.tsv")
ot = t("02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv")
G = lambda c, d: gi[(gi.cohort == c) & (gi.model_drug == d)].iloc[0]
O = lambda s, d: ot[(ot.set == s) & (ot.drug == d)].iloc[0]
S_GI, S_OT = "nar/02_organoid_only/tables/cross_set_summary.tsv", "nar/02_organoid_only/BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv"

r = G("COAD_authors", "FLUOROURACIL")
add("COAD5FU_paper_p", "Colon 5-FU, authors' 2019 files: treated log-rank p (rank-1 BH3-only)", f"{r.treated_p:.5f}", S_GI, ["0.014"])
add("COAD5FU_paper_n", "Colon 5-FU, authors' files: treated n / deaths (resp vs non-resp)",
    f"{r.treated_n} / {r.treated_events} ({r.treated_events_resp} vs {r.treated_events_nonresp})", S_GI, ["114"])
add("COAD5FU_paper_untreated", "Colon 5-FU, authors' files: untreated n, p, HR",
    f"n={r.untreated_n}, p={r.untreated_p:.4f}, HR {r.untreated_HR:.2f} ({r.untreated_HR_lo:.2f}-{r.untreated_HR_hi:.2f})", S_GI, ["0.161"])
add("COAD5FU_paper_shuffle", "Colon 5-FU: response-shuffle share (rank-1 / best-of-30)",
    f"{r.shuffle_share_rank1:.4f} / {r.shuffle_share_best:.4f}", S_GI, ["0.034"])
add("COAD5FU_paper_boot", "Colon 5-FU: random-pathway bootstrap share (rank-1 p statistic)", f"{r.bootstrap_share:.4f}", S_GI,
    ["0.007"], "02 report table also quotes 0.0074 = paper's top-rank statistic (74/10,000); label the two separately")
r = G("COAD", "FLUOROURACIL")
add("COAD5FU_current", "Colon 5-FU, current GDC: n/deaths, p, HR",
    f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.3f}, HR {r.treated_HR:.2f} ({r.treated_HR_lo:.2f}-{r.treated_HR_hi:.2f})", S_GI, ["0.866"])
r = G("COAD", "CETUXIMAB")
add("COAD_cetux_current", "Colon cetuximab (vdW, current): n/deaths, p", f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.3f}", S_GI, ["0.329"],
    "same cohort+drug in 02_OTHER PDX_gao2015_CRC has n=14/6 (different survival source)")
r = G("PAAD", "GEMCITABINE") if False else gi[(gi.folder == "PAAD_shi") & (gi.model_drug == "GEMCITABINE")].iloc[0]
add("PAAD_shi_gem", "Pancreas Shi gemcitabine: n/deaths, p, HR, untreated p",
    f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.4f}, HR {r.treated_HR:.2f}, untreated p={r.untreated_p:.4f}", S_GI, ["0.023"])
for f, lab in [("LIHC_ji", "Ji"), ("LIHC_broutier", "Broutier")]:
    r = gi[gi.folder == f].iloc[0]
    add(f"{f}_sor", f"Liver {lab} sorafenib: n/deaths, p, HR",
        f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.4f}, HR {r.treated_HR:.2f} ({r.treated_HR_lo:.2f}-{r.treated_HR_hi:.2f})",
        S_GI, [f"{r.treated_p:.3f}"],
        "10_completeness recomputed 0.089 for Broutier using the Ji gene-universe bundle; 0.139 is authoritative (02 pipeline, rebuilt identically by 05 A3)" if lab == "Broutier" else "")

r = O("BLCA_lee2018_paper", "CISPLATIN")
add("BLCA_paper_rank1", "Bladder cisplatin, paper NES: n/deaths, p, HR (rank-1 transport of glucose...)",
    f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.4f}, HR {r.treated_hr:.2f} ({r.treated_hr_lo:.2f}-{r.treated_hr_hi:.2f})", S_OT,
    ["0.044", "0.0435"], "paper: n=77 treated, 294 untreated; ours 79 / 288")
add("BLCA_paper_untreated", "Bladder cisplatin, paper NES: untreated n, p", f"n={r.untreated_n}, p={r.untreated_p:.3f}", S_OT, ["0.93"])
add("BLCA_paper_perm", "Bladder cisplatin, paper NES: permutation p", f"{r.perm_empirical_p:.3f}", S_OT, ["0.30", "0.297"])
add("BLCA_paper_boot", "Bladder cisplatin, paper NES: bootstrap p (top-rank stat / rank-1 p stat)",
    f"{r.boot_empirical_p_rank:.4f} / {r.boot_empirical_p_rank1p:.4f}", S_OT, ["0.044", "0.066"],
    "02 report quotes 0.044 (top-rank); 05 A3 quotes 0.066 (rank-1 p). Both correct, different statistics")
r = O("BLCA_lee2018_current", "CISPLATIN")
add("BLCA_current", "Bladder cisplatin, current GDC: n/deaths, p, HR, perm p",
    f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.4f}, HR {r.treated_hr:.2f}, perm={r.perm_empirical_p:.3f}", S_OT, ["0.067", "0.036"])
for s, d in [("PDX_gao2015_CM", "DACARBAZINE"), ("LUAD_cellline_substitute", "PEMETREXED"), ("LUAD_cellline_substitute", "GEMCITABINE"),
             ("PDX_gao2015_NSCLC", "ERLOTINIB")]:
    r = O(s, d)
    add(f"{s}_{d}", f"{s} {d}: n/deaths, p, HR, perm p",
        f"{r.treated_n}/{r.treated_events}, p={r.treated_p:.4f}, HR {r.treated_hr:.2f}, perm={r.perm_empirical_p:.3f}", S_OT,
        [f"{r.treated_p:.4f}"[:5] if r.treated_p > 0.01 else f"{r.treated_p:.4f}".rstrip("0")])

ex = t("05_A3_B1/tables/A3_expected_by_chance.tsv").set_index("family")
SA = "nar/05_A3_B1/tables/A3_expected_by_chance.tsv"
add("A3_treated", "A3: treated log-rank p<0.05 out of 57 (binomial p)",
    f"{int(ex.loc['treated log-rank p (rank-1)', 'observed_p_below_0.05'])}/57 (binom p={ex.loc['treated log-rank p (rank-1)', 'binomial_p_at_least_observed']:.3f})", SA, ["7 of 57"])
add("A3_untreated", "A3: untreated log-rank p<0.05 out of 57 (binomial p, assumes independence)",
    f"{int(ex.loc['untreated log-rank p', 'observed_p_below_0.05'])}/57 (binom p={ex.loc['untreated log-rank p', 'binomial_p_at_least_observed']:.1e})", SA, ["12 of 57"],
    "the 57 untreated tests share untreated patients (e.g. 11 LUAD tests on the same 327); binomial p overstated")
m = t("05_A3_B1/tables/A3_master_results.tsv")
add("A3_minBH", "A3: smallest BH q / Holm p across 57", f"{m.BH_q.min():.3f} / {m.Holm_p.min():.3f}", "nar/05_A3_B1/tables/A3_master_results.tsv", ["0.063"])
b1 = t("05_A3_B1/tables/B1_summary_counts.tsv").set_index("group").loc["all"]
add("B1_counts", "B1: tests / treated p<0.05 (expected dir) / interaction p<0.05 of with-comparator / null p<0.05 / BH<0.05",
    f"{b1.tests} / {b1['treated_p<0.05']} ({b1['treated_p<0.05_expected_dir']}) / {b1['interaction_p<0.05']} of {b1.with_comparator} / {b1['null_p<0.05']} / {b1['BH_q_primary<0.05']}",
    "nar/05_A3_B1/tables/B1_summary_counts.tsv", ["196"])

ag = t("08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv")
for sc in ["GDSC1", "GDSC2", "CTRP", "PRISM"]:
    x = ag[(ag.screen == sc) & (ag.drug == "NAVITOCLAX") & (ag.scope.str.startswith("pan"))].iloc[0]
    add(f"DepMap_navitoclax_{sc}", f"DepMap BH3-only score vs navitoclax, {sc}, pan-cancer lineage-adjusted",
        f"rho={x.rho:.3f} (n={x.n}, p={x.p:.1e})", "nar/08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv", [f"{x.rho:.3f}".replace("-", "−"), f"{x.rho:.3f}"])
x = ag[(ag.screen == "PRISM") & (ag.drug == "NAVITOCLAX") & (ag.scope == "colorectal")].iloc[0]
add("DepMap_navitoclax_PRISM_CRC", "DepMap navitoclax PRISM colorectal lines",
    f"rho={x.rho:.3f} (95% CI {x.ci_low:.3f}-{x.ci_high:.3f}), n={x.n}, p={x.p:.3f}", "nar/08_depmap/tables/A_g_bh3_vs_bh3_mimetics.tsv",
    ["+0.405"], "08 REPORT says 'every confidence interval crosses 0' within colorectal lines: false for this row")
bb = t("08_depmap/tables/A_b_bh3_vs_5fu_colon.tsv")
add("DepMap_BH3_5FU_GDSC1", "DepMap BH3-only vs 5-FU, GDSC1 colorectal", f"rho={bb.iloc[0].rho:.3f}, n={int(bb.iloc[0].n)}" if "rho" in bb else str(bb.iloc[0].to_dict()),
    "nar/08_depmap/tables/A_b_bh3_vs_5fu_colon.tsv", ["0.054"])

s1 = t("09_registered_test/tables/stage1_primary.tsv")
for _, x in s1.iterrows():
    add(f"Reg_S1_{x.drug}", f"Registered Stage 1 interaction HR (PFI), {x.drug}",
        f"{x.PFI_HR_int:.2f} ({x.PFI_lo:.2f}-{x.PFI_hi:.2f}), p={x.PFI_p:.3f}, Holm={x.PFI_holm_p:.2f}, null p={x.random_null_p_one_sided:.3f}",
        "nar/09_registered_test/tables/stage1_primary.tsv", [f"{x.PFI_HR_int:.2f} ({x.PFI_lo:.2f}–{x.PFI_hi:.2f})", f"{x.PFI_HR_int:.2f} ({x.PFI_lo:.2f}-{x.PFI_hi:.2f})"])
s2 = t("09_registered_test/tables/stage2_controlled.tsv")
x = s2[(s2.drug == "OXALIPLATIN") & (s2.cohort == "GSE72970") & (s2.endpoint == "PFS")].iloc[0]
add("Reg_S2_ox_GSE72970", "Registered Stage 2: oxaliplatin GSE72970 PFS interaction (one-sided)",
    f"HR {x.estimate:.2f} ({x.lo:.2f}-{x.hi:.2f}), p1={x.p_one_sided:.3f}, Holm={x.holm_all_primary_stage2:.3f}; FOLFOX {x.n_drug_arm} vs FOLFIRI {x.n_comparator}",
    "nar/09_registered_test/tables/stage2_controlled.tsv", ["1.63 (1.00–2.64)", "1.63 (1.00-2.64)"],
    "after audit fix I02 (nar/15_fixes) 03 and 05 also exclude FOLFIRINOX: 36 vs 84 (09's 83 also drops 1 XELIRI+bevacizumab patient)")

sm = t("03_organoid_plus_patient/tables/summary_main.tsv")
sm = sm[(sm.features == "shortlist") & (sm.outcome == "survival")]
for u in ["CRC_5FU", "BLCA_CISPLATIN", "OV_PACLITAXEL", "OV_PLATINUM", "BLCA_GEMCITABINE"]:
    x = sm[(sm.unit == u) & (sm.model == "joint")].iloc[0]
    # after audit fix I05 (nar/15_fixes): 03 reports HR per SD of the RISK score, HR > 1 = expected, like all folders
    add(f"Joint_{u}", f"03 joint model held-out HR per SD of risk score (HR>1 = expected, harmonised I05), {u}",
        f"C={x.C:.3f}, HR {x.HR_risk:.2f} ({x.HR_risk_low:.2f}-{x.HR_risk_high:.2f}), p={x.p:.3f}, n={int(x.n)}/{int(x.events)}",
        "nar/03_organoid_plus_patient/tables/summary_main.tsv", [f"{x.HR_risk:.2f} ({x.HR_risk_low:.2f}-{x.HR_risk_high:.2f})"])

df = pd.DataFrame(rows)
df["status"] = ["quoted" if r != "-" else "not found as text" for r in df.reports_quoting]
df.to_csv(OUT / "MASTER_NUMBERS.tsv", sep="\t", index=False)
print(df[["id", "authoritative_value", "reports_quoting"]].to_string())
