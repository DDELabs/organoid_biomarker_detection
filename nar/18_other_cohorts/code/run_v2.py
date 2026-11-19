"""W2 analyses: run_train_combined parts 1-3 on training_table_v2 (outputs in nar/18_other_cohorts/tables, never in robust/results/train_combined)."""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "robust"))
warnings.filterwarnings("ignore")
import run_train_combined as RTC
from obd import survlearn as SL
OUT = HERE.parent / "tables"; OUT.mkdir(exist_ok=True)
RTC.OUT = OUT
V2 = RTC.RESULTS / "training_table_v2"
LUNG_ONLY = ROOT / "data" / "external" / "training_table_v2_lung_only"
RECTAL = ["GSE35452", "GSE119409", "GSE150082", "GSE45404_GPL4133", "GSE45404_GPL570", "GSE145037", "GSE133057"]
CHEMO_BREAST = "CHEMOTHERAPY|DOXORUBICIN|EPIRUBICIN|CYCLOPHOSPHAMIDE|PACLITAXEL|DOCETAXEL|FLUOROURACIL|CAPECITABINE|ANTHRACYCLINE|TAXANE|CARBOPLATIN|CISPLATIN|METHOTREXATE"

def hm_ci(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
    return max(0, a - 1.96 * se), min(1, a + 1.96 * se)

def load(path):
    RTC.T = path
    return RTC.load()

def part1_tables(P, Z, tag):
    R = RTC.part1(P, Z)
    ci = [hm_ci(r.auc_pathways, r.responders, r.n - r.responders) for r in R.itertuples()]
    R["auc_ci_low"], R["auc_ci_high"] = [c[0] for c in ci], [c[1] for c in ci]
    R.insert(0, "table", tag)
    return R

def cox_ci(r, name):
    b, se = r.loc[name, "coef"], r.loc[name, "se"]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), r.loc[name, "p"]

def part3_ci(P, Z, name, cohorts, drug_pat, use_stage=True):
    """Copy of run_train_combined.part3 with a 95% CI added (same maths)."""
    d = RTC.surv_frame(P, cohorts, use_stage)
    d = d[d["treated"].eq(0) | d["drugs"].fillna("").str.contains(drug_pat)]
    rows, base = [], ["age_n", "stage_n"]
    for c in sorted(d["cohort"].unique()):
        tr, te = d[d["cohort"] != c], d[d["cohort"] == c]
        if len(tr) == 0 or te["event"].sum() < 15 or te["treated"].sum() < 15 or (1 - te["treated"]).sum() < 15:
            continue
        _, gam = SL.fit_interaction(tr, Z.loc[tr.index], base, 300.0, 300.0)
        te = te.assign(score=Z.loc[te.index].values @ gam.values)
        te["score"] = (te["score"] - te["score"].mean()) / te["score"].std()
        te["tx_score"] = te["treated"] * te["score"]
        r = SL.cox_test(te, base + ["treated", "score", "tx_score"])
        hr, lo, hi, p = cox_ci(r, "tx_score")
        rows.append(dict(setting=name, test_cohort=c, n=len(te), treated=int(te["treated"].sum()), events=int(te["event"].sum()),
                         train_n=len(tr), train_cohorts=tr["cohort"].nunique(), interaction_hazard_ratio=hr, hr_ci_low=lo, hr_ci_high=hi, p_interaction=p))
    return rows

def crossfit_survival(P, Z, cohort, tx, label, endpoints, reps=5):
    """Randomised trial, drug vs control: cross-fitted treatment x score interaction for survival (model never sees the scored patient)."""
    rows = []
    for ep in endpoints:
        d = P[P["cohort"].eq(cohort) & P[f"{ep}_event"].notna() & P[f"{ep}_months"].gt(0)].copy()
        d["time"], d["event"], d["project"], d["treated"] = d[f"{ep}_months"], d[f"{ep}_event"], d["cohort"], tx.loc[d.index].astype(float)
        sc = pd.Series(0.0, index=d.index)
        for seed in range(reps):
            f = SL.folds_by_patient(d.index, 5, seed)
            for k in range(5):
                tr, te = d.index[f != k], d.index[f == k]
                _, g = SL.fit_interaction(d.loc[tr], Z.loc[tr], [], 300.0, 300.0)
                sc.loc[te] += Z.loc[te].values @ g.values / reps
        d["score"] = (sc - sc.mean()) / sc.std(); d["tx_score"] = d["treated"] * d["score"]
        r = SL.cox_test(d, ["treated", "score", "tx_score"])
        hr, lo, hi, p = cox_ci(r, "tx_score")
        rows.append(dict(test=label, endpoint=ep.upper(), n=len(d), treated=int(d["treated"].sum()), events=int(d["event"].sum()),
                         interaction_hazard_ratio=hr, hr_ci_low=lo, hr_ci_high=hi, p_interaction=p))
    return rows

def safe(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception as e:
        print("FAILED", a[2:3], repr(e)[:150]); return []

def main():
    P, Z = load(V2)
    # ---- acceptance: old breast AUCs reproduced from the lung-only v2 copy
    P0, Z0 = load(LUNG_ONLY)
    old = pd.read_csv(ROOT / "robust/results/train_combined/part1_response_leave_one_cohort_out.tsv", sep="\t")
    r0 = part1_tables(P0, Z0, "lung_only_v2")
    mrg = old.merge(r0, on=["cancer", "test_cohort"], suffixes=("_v1", "_v2lung"))
    mrg["same_auc"] = np.isclose(mrg["auc_pathways_v1"], mrg["auc_pathways_v2lung"], atol=1e-9)
    mrg[["cancer", "test_cohort", "n_v1", "auc_pathways_v1", "auc_pathways_v2lung", "same_auc", "auc_trained_on_other_cancers_v1", "auc_trained_on_other_cancers_v2lung"]].to_csv(OUT / "reproduce_v1_part1.tsv", sep="\t", index=False)
    print("reproduced:", int(mrg.same_auc.sum()), "of", len(mrg))
    # ---- part 1, full v2 (as originally defined) and with rectal chemoradiotherapy split out
    R1 = part1_tables(P, Z, "v2_all_groups")
    P2 = P.copy(); P2.loc[P2["cohort"].isin(RECTAL), "cancer_group"] = "rectal chemoradiotherapy"
    R1b = part1_tables(P2, Z, "v2_rectal_split")
    pd.concat([R1, R1b]).to_csv(OUT / "part1_v2.tsv", sep="\t", index=False)
    # ---- part 2
    R2 = RTC.part2(P, Z)
    rows = []
    L = P[P["has_response"] & P["responder"].notna()]
    noah = L[L["cohort"].eq("GSE50948")]
    rows.append(safe(lambda *_: [RTC.drug_specific(P, Z, noah.index[noah["arm"].str.contains("trastuzumab")], noah.index[noah["arm"].str.contains("alone \\(HER2-positive")], "NOAH (randomised): trastuzumab + chemo vs chemo alone", ["er"])], 0, 0, "noah"))
    cl = L[L["cohort"].eq("GSE66305")]
    a = cl.index[cl["arm"].eq("chemotherapy + trastuzumab")]
    for arm in ("chemotherapy + lapatinib", "chemotherapy + trastuzumab + lapatinib"):
        rows.append(safe(lambda *_: [RTC.drug_specific(P, Z, cl.index[cl["arm"].eq(arm)], a, f"CHER-LOB (randomised): {arm} vs chemotherapy + trastuzumab", [])], 0, 0, arm))
    ap = L[L["cohort"].eq("GSE9782") & L["randomised"].eq(True)]
    rows.append(safe(lambda *_: [RTC.drug_specific(P, Z, ap.index[ap["arm"].eq("bortezomib")], ap.index[ap["arm"].eq("dexamethasone")], "APEX (randomised): bortezomib vs dexamethasone, response", [])], 0, 0, "apex"))
    # transfer: anti-HER2 trained elsewhere, tested in NOAH
    tr = L.index[L["drugs"].fillna("").str.contains("TRASTUZUMAB") & ~L["cohort"].isin(["GSE50948"]) & L["cancer_group"].eq("breast") & L["possible_duplicate_of"].isna()]
    rows.append(safe(lambda *_: [RTC.transfer(P, Z, tr, noah.index[noah["arm"].str.contains("trastuzumab")], noah.index[noah["arm"].str.contains("alone \\(HER2-positive")], "TRANSFER trastuzumab: trained in other trastuzumab-treated breast cohorts, tested in NOAH", ["er"])], 0, 0, "xfer"))
    R2n = pd.DataFrame([x for r in rows for x in r])
    # ICON7 bevacizumab and APEX survival (cross-fitted interaction)
    ic = P[P["cohort"].eq("GSE140082")]
    s = safe(crossfit_survival, P, Z, "GSE140082", ic["drugs"].fillna("").str.contains("BEVACIZUMAB"), "ICON7 (randomised): bevacizumab x score", ["pfs", "os"])
    ap2 = P[P["cohort"].eq("GSE9782") & P["randomised"].eq(True)]
    s += safe(crossfit_survival, P.loc[ap2.index], Z, "GSE9782", ap2["arm"].eq("bortezomib"), "APEX (randomised): bortezomib x score", ["pfs", "os"])
    pd.concat([R2.assign(table="original_part2_rerun"), R2n.assign(table="new_drug_specific_response")]).to_csv(OUT / "part2_v2.tsv", sep="\t", index=False)
    pd.DataFrame(s).to_csv(OUT / "new_randomised_survival_tests.tsv", sep="\t", index=False)
    pd.concat([R2n.assign(endpoint="response"), pd.DataFrame(s)]).to_csv(OUT / "new_randomised_tests.tsv", sep="\t", index=False)
    # ---- part 3
    r = RTC_p3 = []
    r += part3_ci(P, Z, "colorectal 5-FU after surgery", ["GSE39582", "GSE14333", "GSE103479", "TCGA-COAD", "TCGA-READ"], "FLUOROURACIL|CHEMOTHERAPY|CAPECITABINE")
    r += part3_ci(P, Z, "bladder chemotherapy", ["GSE48276", "GSE31684", "GSE13507", "TCGA-BLCA"], "CISPLATIN|CHEMOTHERAPY|GEMCITABINE|CARBOPLATIN", use_stage=False)
    br = ["GSE6532_GPL96", "GSE6532_GPL570", "GSE9195", "GSE17705", "GSE12093", "GSE7390", "GSE11121", "GSE45255", "GSE58644", "TCGA-BRCA"]
    r += part3_ci(P, Z, "breast tamoxifen vs none", br, "TAMOXIFEN", use_stage=False)
    r += part3_ci(P, Z, "breast chemotherapy vs none", br + ["GSE31519"], CHEMO_BREAST, use_stage=False)
    r += part3_ci(P, Z, "ovarian platinum vs none", ["TCGA-OV", "GSE32062", "GSE26712", "GSE49997", "GSE17260", "GSE30161", "GSE9891", "GSE63885"], "PLATINUM|CARBOPLATIN|CISPLATIN", use_stage=False)
    pd.DataFrame(r).to_csv(OUT / "part3_v2.tsv", sep="\t", index=False)
    print(pd.DataFrame(r).round(3).to_string())
    print(R1b.round(3).to_string()); print(R2n.round(3).to_string()); print(pd.DataFrame(s).round(3).to_string())

if __name__ == "__main__":
    main()
