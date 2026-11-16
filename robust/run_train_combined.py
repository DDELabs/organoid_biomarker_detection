"""Train on the combined training table (robust/results/training_table) and report honest, held-out results.

Part 1  response model per cancer: train on all other cohorts of that cancer, test on a cohort never seen
        (leave-one-cohort-out). Compared with a clinical-only baseline where one exists (breast: ER, HER2).
Part 2  drug-specific test in randomised trials: a model trained on patients who got drug X must predict response
        in X-treated patients better than in control-arm patients (score x arm interaction). Includes a transfer test
        between two independent trials of the same drug (carboplatin: BrighTNess <-> I-SPY2).
Part 3  drug-specific survival test for colorectal 5-FU: treated x score interaction learned on some cohorts and
        tested on a cohort never seen.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
from obd import RESULTS  # noqa: E402
from obd import survlearn as SL  # noqa: E402

T = RESULTS / "training_table"
OUT = RESULTS / "train_combined"
LAM = 100.0
STAGE = {"1": 1, "2": 2, "3": 3, "4": 4, "I": 1, "II": 2, "III": 3, "IV": 4}


def load():
    P = pd.read_csv(T / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    S = pd.concat([pd.read_csv(T / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")]).loc[P.index]
    P = P[P["possible_duplicate_of"].isna()].copy()
    yes = {"1", "1.0", "positive", "pos", "p", "+", "yes"}; no = {"0", "0.0", "negative", "neg", "n", "-", "no"}
    for c in ("er", "her2"):
        v = P[c].astype(str).str.strip().str.lower()
        P[c] = np.where(v.isin(yes), 1.0, np.where(v.isin(no), 0.0, np.nan))
    Z = SL.zscore_within(S.loc[P.index], P["cohort"])       # standardise inside each cohort (platforms differ)
    return P, Z


def fit(X, y, lam=LAM):
    return LogisticRegression(C=1.0 / lam, max_iter=3000).fit(X, y)


def auc(y, s):
    y, s = np.asarray(y), np.asarray(s)
    if len(set(y)) < 2:
        return np.nan, np.nan
    return roc_auc_score(y, s), mannwhitneyu(s[y == 1], s[y == 0], alternative="two-sided").pvalue


def auc_se(y, s):
    """Hanley-McNeil standard error of an AUROC."""
    a = roc_auc_score(y, s); n1, n0 = int(np.sum(y)), int(len(y) - np.sum(y))
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return a, np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))


# ------------------------------------------------------------------------------------------- part 1
def part1(P, Z):
    L = P[P["has_response"] & P["treated"].eq(1)]
    rows = []
    for grp, d in L.groupby("cancer_group"):
        ok = [c for c, x in d.groupby("cohort") if x["responder"].sum() >= 8 and (1 - x["responder"]).sum() >= 8]
        if len(ok) < 2:
            continue
        d = d[d["cohort"].isin(ok)]
        for c in ok:
            tr, te = d.index[d["cohort"] != c], d.index[d["cohort"] == c]
            ytr, yte = d.loc[tr, "responder"].astype(int), d.loc[te, "responder"].astype(int)
            a, p = auc(yte, fit(Z.loc[tr].values, ytr).decision_function(Z.loc[te].values))
            # same test cohort, model trained on every OTHER cancer: is the signal specific to this cancer?
            oth = L.index[L["cancer_group"] != grp]
            a_o, _ = auc(yte, fit(Z.loc[oth].values, L.loc[oth, "responder"].astype(int)).decision_function(Z.loc[te].values))
            r = dict(cancer=grp, test_cohort=c, n=len(te), responders=int(yte.sum()), train_n=len(tr), auc_pathways=a, p=p,
                     auc_trained_on_other_cancers=a_o)
            if grp == "breast":
                k = ["er", "her2"]
                dtr, dte = d.loc[tr, k].dropna(), d.loc[te, k].dropna()
                if len(dte) > 20 and dte.nunique().max() > 1:
                    cm = LogisticRegression(C=1e3).fit(dtr.values, ytr.loc[dtr.index])
                    cs = pd.Series(cm.decision_function(dte.values), dte.index)
                    r["auc_er_her2_only"] = auc(yte.loc[dte.index], cs)[0]
                    strata = (dte["er"].astype(int).astype(str) + dte["her2"].astype(int).astype(str)).values
                    sc = pd.Series(fit(Z.loc[tr].values, ytr).decision_function(Z.loc[dte.index].values), dte.index)
                    r["auc_pathways_within_er_her2_groups"] = SL.stratified_auc(sc.values, yte.loc[dte.index].values, strata, 5)[0]
            rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "part1_response_leave_one_cohort_out.tsv", sep="\t", index=False)
    return R


# ------------------------------------------------------------------------------------------- part 2
def interaction_p(y, score, arm, cov):
    """Logistic response ~ arm + score + arm*score (+ covariates); returns interaction z and p."""
    import statsmodels.api as sm
    s = (score - score.mean()) / score.std()
    X = pd.DataFrame({"arm": arm, "score": s, "arm_x_score": arm * s}, index=score.index)
    for c in cov.columns:
        if cov[c].nunique() > 1:
            X[c] = cov[c].values
    r = sm.Logit(y.astype(float), sm.add_constant(X.astype(float))).fit(disp=0)
    return r.params["arm_x_score"], r.pvalues["arm_x_score"]


def drug_specific(P, Z, treated_idx, control_idx, label, cov_cols, seeds=10):
    """Cross-fitted model of response in the treated arm; applied to held-out treated and to all control patients."""
    y = P["responder"].fillna(-1).astype(int)   # only labelled rows are indexed below
    sc_t = pd.Series(0.0, index=treated_idx); sc_c = pd.Series(0.0, index=control_idx)
    for seed in range(seeds):
        f = SL.folds_by_patient(treated_idx, 5, seed)
        for k in range(5):
            tr, te = f.index[f != k], f.index[f == k]
            m = fit(Z.loc[tr].values, y.loc[tr])
            sc_t.loc[te] += m.decision_function(Z.loc[te].values) / seeds
            sc_c += m.decision_function(Z.loc[control_idx].values) / (5 * seeds)
    at, se_t = auc_se(y.loc[treated_idx].values, sc_t.values)
    ac, se_c = auc_se(y.loc[control_idx].values, sc_c.values)
    idx = list(treated_idx) + list(control_idx)
    arm = pd.Series([1.0] * len(treated_idx) + [0.0] * len(control_idx), index=idx)
    sc = pd.concat([(sc_t - sc_t.mean()) / sc_t.std(), (sc_c - sc_c.mean()) / sc_c.std()])
    b, p = interaction_p(y.loc[idx], sc, arm, P.loc[idx, cov_cols].fillna(0))
    zdiff = (at - ac) / np.sqrt(se_t ** 2 + se_c ** 2)
    return dict(test=label, n_drug=len(treated_idx), resp_drug=int(y.loc[treated_idx].sum()), n_control=len(control_idx),
                resp_control=int(y.loc[control_idx].sum()), auc_in_drug_arm=at, auc_in_control_arm=ac,
                p_auc_difference=2 * norm.sf(abs(zdiff)), interaction_coef=b, p_interaction=p)


def transfer(P, Z, train_idx, drug_idx, control_idx, label, cov_cols):
    """Model trained in one trial's drug arm, tested in another trial's drug arm vs its control arm."""
    y = P["responder"].fillna(-1).astype(int)   # only labelled rows are indexed below
    m = fit(Z.loc[train_idx].values, y.loc[train_idx])
    st, scn = pd.Series(m.decision_function(Z.loc[drug_idx].values), drug_idx), pd.Series(m.decision_function(Z.loc[control_idx].values), control_idx)
    at, se_t = auc_se(y.loc[drug_idx].values, st.values); ac, se_c = auc_se(y.loc[control_idx].values, scn.values)
    idx = list(drug_idx) + list(control_idx)
    arm = pd.Series([1.0] * len(drug_idx) + [0.0] * len(control_idx), index=idx)
    b, p = interaction_p(y.loc[idx], pd.concat([st, scn]), arm, P.loc[idx, cov_cols].fillna(0))
    return dict(test=label, n_drug=len(drug_idx), resp_drug=int(y.loc[drug_idx].sum()), n_control=len(control_idx),
                resp_control=int(y.loc[control_idx].sum()), auc_in_drug_arm=at, auc_in_control_arm=ac,
                p_auc_difference=2 * norm.sf(abs((at - ac) / np.sqrt(se_t ** 2 + se_c ** 2))), interaction_coef=b, p_interaction=p)


def part2(P, Z):
    rows = []
    I = P[P["cohort"] == "ISPY2"]
    neg = I[I["her2"] == 0]
    ctrl = neg.index[neg["arm"] == "Paclitaxel"]
    for arm, d in neg.groupby("arm"):
        if arm == "Paclitaxel" or len(d) < 50 or d["responder"].sum() < 10:
            continue
        rows.append(drug_specific(P, Z, d.index, ctrl, f"I-SPY2 (HER2-negative): {arm} vs paclitaxel alone", ["er"]))
    B = P[P["cohort"] == "BRIGHTNESS"]
    bc, bcarb = B.index[B["arm"] == "PACLITAXEL"], B.index[B["arm"] != "PACLITAXEL"]
    rows.append(drug_specific(P, Z, bcarb, bc, "BrighTNess: paclitaxel + carboplatin (with or without veliparib) vs paclitaxel alone", []))
    tn = neg[neg["er"] == 0]
    vc, tnc = tn.index[tn["arm"].str.contains("ABT 888")], tn.index[tn["arm"] == "Paclitaxel"]
    rows.append(transfer(P, Z, bcarb, vc, tnc, "TRANSFER carboplatin: trained in BrighTNess, tested in I-SPY2 triple-negative", []))
    allvc = neg.index[neg["arm"].str.contains("ABT 888")]
    rows.append(transfer(P, Z, allvc, bcarb, bc, "TRANSFER carboplatin: trained in I-SPY2, tested in BrighTNess", []))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "part2_drug_specific_randomised.tsv", sep="\t", index=False)
    return R


# ------------------------------------------------------------------------------------------- part 3
def surv_frame(P, cohorts, use_stage=True):
    d = P[P["cohort"].isin(cohorts) & P["treated"].notna()].copy()
    for a in ("rfs", "dfs", "pfs", "os"):
        use = d.get("time", pd.Series(np.nan, index=d.index)).isna() & d[f"{a}_months"].notna() & d[f"{a}_event"].notna()
        d.loc[use, "time"], d.loc[use, "event"] = d.loc[use, f"{a}_months"], d.loc[use, f"{a}_event"]
    d = d[d["time"] > 0]
    st = d["stage"].astype(str).str.upper().str.replace("STAGE", "").str.strip().str.extract(r"^(IV|III|II|I|[1-4])")[0].map(STAGE)
    if use_stage:
        d = d[st.notna() & (st < 4)]                      # stage I-III: relapse after surgery
        d["stage_n"] = st.loc[d.index] - 2
    else:
        d["stage_n"] = 0.0
    d["age_n"] = ((pd.to_numeric(d["age"], errors="coerce") - 65) / 10).fillna(0)
    d["project"] = d["cohort"]
    return d


def part3(P, Z, name, cohorts, drug_pat, use_stage=True):
    d = surv_frame(P, cohorts, use_stage)
    d = d[d["treated"].eq(0) | d["drugs"].fillna("").str.contains(drug_pat)]
    rows = []
    base = ["age_n", "stage_n"]
    for c in sorted(d["cohort"].unique()):
        tr, te = d[d["cohort"] != c], d[d["cohort"] == c]
        if len(tr) == 0 or te["event"].sum() < 15 or te["treated"].sum() < 15 or (1 - te["treated"]).sum() < 15:
            continue
        _, gam = SL.fit_interaction(tr, Z.loc[tr.index], base, 300.0, 300.0)
        te = te.assign(score=Z.loc[te.index].values @ gam.values)
        te["score"] = (te["score"] - te["score"].mean()) / te["score"].std()
        te["tx_score"] = te["treated"] * te["score"]
        r = SL.cox_test(te, base + ["treated", "score", "tx_score"])
        rows.append(dict(setting=name, test_cohort=c, n=len(te), treated=int(te["treated"].sum()), events=int(te["event"].sum()),
                         train_n=len(tr), interaction_hazard_ratio=np.exp(r.loc["tx_score", "coef"]), p_interaction=r.loc["tx_score", "p"]))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P, Z = load()
    pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 80)
    R1 = part1(P, Z); print(R1.round(3).to_string())
    R2 = part2(P, Z); print(R2.round(3).to_string())
    rows = part3(P, Z, "colorectal 5-FU after surgery", ["GSE39582", "GSE14333", "GSE103479", "TCGA-COAD", "TCGA-READ"], "FLUOROURACIL|CHEMOTHERAPY|CAPECITABINE")
    rows += part3(P, Z, "bladder chemotherapy", ["GSE48276", "GSE31684", "TCGA-BLCA"], "CISPLATIN|CHEMOTHERAPY|GEMCITABINE|CARBOPLATIN", use_stage=False)
    R3 = pd.DataFrame(rows); R3.to_csv(OUT / "part3_survival_interaction.tsv", sep="\t", index=False); print(R3.round(3).to_string())


if __name__ == "__main__":
    main()
