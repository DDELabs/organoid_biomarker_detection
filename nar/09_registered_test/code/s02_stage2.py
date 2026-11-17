"""Stage 2: registered one-sided tests of each frozen score in external cohorts (ANALYSIS_PLAN.md §3)."""
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import common as C
from obd import scoring
from obd.cohorts import to_canonical
from obd.curated import load_curated_trial
from obd.reference import gene_sets

warnings.filterwarnings("ignore")
LOG_EFF = np.log(1 / 0.7)

FOLFOX_72970 = ["FOLFOX", "FOLFOX+BEVACIZUMAB"]
FOLFIRI_72970 = ["FOLFIRI", "FOLFIRI+BEVACIZUMAB", "FOLFIRI+ERBITUX"]
OX_104645 = ["FOLFOX", "FOLFOX+Bev", "SOX", "SOX+Bev", "XELOX", "XELOX+Bev"]
IRI_104645 = ["FOLFIRI", "FOLFIRI+Bev", "FOLFIRI+Pmab", "IRIS", "IRIS+Bev"]
FU_39582 = ["5FU", "FUFOL", "FOLFOX", "FOLFIRI"]
ISPY_HER2_T = ["Paclitaxel + Trastuzumab", "Paclitaxel + Pertuzumab + Trastuzumab", "Paclitaxel + AMG 386 + Trastuzumab",
               "Paclitaxel + MK-2206 + Trastuzumab"]

# (drug, cohort, drug-arm list, comparator-arm list, endpoint type, endpoint, role, randomised, extra filter)
CONTROLLED = [
    ("FLUOROURACIL", "GSE39582", FU_39582, ["no chemotherapy"], "tte", "rfs", "primary", False, None),
    ("FLUOROURACIL", "GSE39582", FU_39582, ["no chemotherapy"], "tte", "os", "secondary", False, None),
    ("FLUOROURACIL", "GSE103479", ["adjuvant chemo"], ["surgery only"], "tte", "pfs", "primary", False, None),
    ("FLUOROURACIL", "GSE103479", ["adjuvant chemo"], ["surgery only"], "tte", "os", "secondary", False, None),
    ("FLUOROURACIL", "GSE14333", ["adjuvant chemotherapy (type not given)"], ["no chemotherapy"], "tte", "dfs", "primary", False, None),
    ("OXALIPLATIN", "GSE72970", FOLFOX_72970, FOLFIRI_72970, "tte", "pfs", "primary", False, None),
    ("OXALIPLATIN", "GSE72970", FOLFOX_72970, FOLFIRI_72970, "tte", "os", "secondary", False, None),
    ("OXALIPLATIN", "GSE72970", FOLFOX_72970, FOLFIRI_72970, "bin", "responder", "secondary", False, None),
    ("OXALIPLATIN", "GSE39582", ["FOLFOX"], ["5FU", "FUFOL"], "tte", "rfs", "primary", False, None),
    ("OXALIPLATIN", "GSE39582", ["FOLFOX"], ["5FU", "FUFOL"], "tte", "os", "secondary", False, None),
    ("OXALIPLATIN", "GSE104645", OX_104645, IRI_104645, "bin", "responder", "primary", False, None),
    ("OXALIPLATIN", "GSE87211", ["5-FU + Oxaliplatin + RT"], ["5-FU + RT"], "tte", "dfs", "primary", False, None),
    ("OXALIPLATIN", "GSE87211", ["5-FU + Oxaliplatin + RT"], ["5-FU + RT"], "bin", "responder", "secondary", False, None),
    ("IRINOTECAN", "GSE72970", FOLFIRI_72970, FOLFOX_72970, "tte", "pfs", "primary", False, None),
    ("IRINOTECAN", "GSE72970", FOLFIRI_72970, FOLFOX_72970, "tte", "os", "secondary", False, None),
    ("IRINOTECAN", "GSE72970", FOLFIRI_72970, FOLFOX_72970, "bin", "responder", "secondary", False, None),
    ("IRINOTECAN", "GSE104645", IRI_104645, OX_104645, "bin", "responder", "primary", False, None),
    ("IRINOTECAN", "GSE39582", ["FOLFIRI"], ["5FU", "FUFOL"], "tte", "rfs", "primary", False, None),
    ("PACLITAXEL", "GSE20271", ["T/FAC"], ["FAC"], "bin", "responder", "primary", True, None),
    ("PACLITAXEL", "GSE41998", ["AC->paclitaxel"], ["AC->ixabepilone"], "bin", "responder", "primary", True, None),
    ("PACLITAXEL", "ISPY2", ISPY_HER2_T, ["T-DM1 + Pertuzumab"], "bin", "responder", "primary", True, "her2pos"),
    ("PACLITAXEL", "GSE9891", ["platinum+taxane"], ["platinum"], "tte", "rfs", "primary", False, None),
    ("PACLITAXEL", "GSE9891", ["platinum+taxane"], ["platinum"], "tte", "os", "secondary", False, None),
    ("PACLITAXEL", "GSE63885", ["taxane/platinum"], ["platinum/cyclophosphamide"], "tte", "os", "primary", False, None),
    ("PACLITAXEL", "GSE15622", ["paclitaxel"], ["carboplatin"], "bin", "responder", "primary", False, None),
    ("SORAFENIB", "GSE109211", ["sorafenib"], ["placebo"], "bin", "responder", "primary", True, None),
    ("GEMCITABINE", "GSE42127", ["adjuvant chemo"], ["observation"], "tte", "os", "exploratory", False, None),
    ("PEMETREXED", "GSE42127", ["adjuvant chemo"], ["observation"], "tte", "os", "exploratory", False, None),
]
SINGLE = [
    ("FLUOROURACIL", "GSE87211", ["5-FU + RT"], "tte", "dfs"),
    ("OXALIPLATIN", "GSE28702", ["mFOLFOX6"], "bin", "responder"),
    ("OXALIPLATIN", "GSE19860", ["mFOLFOX6"], "bin", "responder"),
    ("IRINOTECAN", "GSE62080", ["FOLFIRI"], "bin", "responder"),
    ("PACLITAXEL", "BRIGHTNESS", ["PACLITAXEL"], "bin", "responder"),
    ("PACLITAXEL", "ISPY2", ["Paclitaxel"], "bin", "responder"),
    ("PACLITAXEL", "GSE32646", ["P->FEC"], "bin", "responder"),
    ("PACLITAXEL", "GSE51373", ["platinum-taxane"], "bin", "responder"),
    ("PACLITAXEL", "GSE25066", None, "tte", "rfs"),
]


def storm():
    expr, clin = load_curated_trial("GSE109211", canonical=False)
    e = to_canonical(expr)
    gs = {k: [g for g in v if g in e.index] for k, v in gene_sets(("REACTOME",)).items()}
    s = scoring.rank_score(e, {k: v for k, v in gs.items() if len(v) >= 5}).T
    c = pd.DataFrame(index=clin.index)
    c["cohort"], c["arm"], c["responder"] = "GSE109211", clin["arm"], clin["responder"].astype(float)
    for k in ("age", "sex", "stage", "her2"):
        c[k] = np.nan
    s.index = [f"GSE109211:{i}" for i in s.index]
    c.index = s.index
    return c, s


def cohort(P, S, cid, H):
    p = P[P["cohort"] == cid].copy()
    s = S.loc[p.index]
    Zm = (s - s.mean()) / s.std()
    for drug, w in H.items():
        x = (Zm[list(w)] * pd.Series(w)).mean(axis=1)
        p[f"score_{drug}"] = (x - x.mean()) / x.std()
    return p


def covars(p):
    X = pd.DataFrame(index=p.index)
    a = pd.to_numeric(p["age"], errors="coerce")
    if a.notna().sum() > 0.5 * len(p) and a.std() > 0:
        X["age10"] = ((a - a.mean()) / 10).fillna(0)
    m = p["sex"].astype(str).str.upper().str[0].map({"M": 1.0, "F": 0.0})
    if m.notna().sum() > 0.5 * len(p) and m.std() > 0:
        X["male"] = (m - m.mean()).fillna(0)
    st = p["stage"].map(C.stage_num)
    return X, st


def fit_tte(p, score, arm, ep):
    t, e = pd.to_numeric(p[f"{ep}_months"], errors="coerce"), pd.to_numeric(p[f"{ep}_event"], errors="coerce")
    ok = t.notna() & e.notna() & (t > 0) & score.notna()
    p, t, e, score = p[ok], t[ok], e[ok], score[ok]
    arm = None if arm is None else arm[ok]
    X, st = covars(p)
    strata = st.fillna(-1).values if st.notna().sum() > 0.5 * len(p) else None
    if arm is None:
        D = float(e.sum())
        M = np.column_stack([score.values, X.values]) if X.shape[1] else score.values[:, None]
        b, se, _ = C.cox_fit(t.values, e.values, M, strata=strata)
        return {"n": int(ok.sum()), "events": D, "b": b[0], "se": se[0], "expected_se": 1 / np.sqrt(D)}
    D1, D0 = float(e[arm == 1].sum()), float(e[arm == 0].sum())
    out = {"n": int(ok.sum()), "n_drug_arm": int((arm == 1).sum()), "n_comparator": int((arm == 0).sum()),
           "events_drug_arm": D1, "events_comparator": D0}
    if D1 < 5 or D0 < 5:
        return {**out, "fitted": False}
    M = np.column_stack([score.values, arm.values, score.values * arm.values] + ([X.values] if X.shape[1] else []))
    b, se, cov = C.cox_fit(t.values, e.values, M, strata=strata)
    se_arm1 = np.sqrt(se[0] ** 2 + se[2] ** 2 + 2 * cov[0, 2])
    return {**out, "fitted": True, "b": b[2], "se": se[2], "expected_se": np.sqrt(1 / D1 + 1 / D0),
            "eff_comparator": np.exp(b[0]), "eff_drug_arm": np.exp(b[0] + b[2]),
            "eff_drug_arm_lo": np.exp(b[0] + b[2] - 1.96 * se_arm1), "eff_drug_arm_hi": np.exp(b[0] + b[2] + 1.96 * se_arm1),
            "eff_comparator_lo": np.exp(b[0] - 1.96 * se[0]), "eff_comparator_hi": np.exp(b[0] + 1.96 * se[0]),
            "covariates": ";".join(list(X.columns) + (["strata:stage"] if strata is not None else []))}


def fit_bin(p, score, arm):
    y = pd.to_numeric(p["responder"], errors="coerce")
    ok = y.notna() & score.notna()
    p, y, score = p[ok], y[ok], score[ok]
    arm = None if arm is None else arm[ok]
    X, st = covars(p)
    if st.notna().sum() > 0.5 * len(p) and st.nunique() > 1:
        X = pd.concat([X, pd.get_dummies(st.fillna(-1), prefix="stage", drop_first=True).astype(float)], axis=1)
    if arm is None:
        n1, n0 = int(y.sum()), int((1 - y).sum())
        M = sm.add_constant(np.column_stack([score.values] + ([X.values] if X.shape[1] else [])))
        r = sm.Logit(y.values, M).fit(disp=0)
        return {"n": int(ok.sum()), "responders": n1, "nonresponders": n0, "b": r.params[1], "se": r.bse[1],
                "expected_se": np.sqrt(1 / (len(y) * y.mean() * (1 - y.mean())))}
    a = arm.values
    r1, r0 = y[a == 1], y[a == 0]
    out = {"n": int(ok.sum()), "n_drug_arm": int(len(r1)), "n_comparator": int(len(r0)),
           "responders_drug_arm": int(r1.sum()), "nonresponders_drug_arm": int((1 - r1).sum()),
           "responders_comparator": int(r0.sum()), "nonresponders_comparator": int((1 - r0).sum())}
    if min(out["responders_drug_arm"], out["nonresponders_drug_arm"], out["responders_comparator"], out["nonresponders_comparator"]) < 5:
        return {**out, "fitted": False}
    used = list(X.columns)
    for attempt in range(2):
        M = sm.add_constant(np.column_stack([score.values, a, score.values * a] + ([X[used].values] if used else [])))
        try:
            r = sm.Logit(y.values, M).fit(disp=0, maxiter=200)
            if np.all(np.isfinite(r.bse)):
                break
        except Exception:
            pass
        used = []
    cov = r.cov_params()
    se_arm1 = np.sqrt(cov[1, 1] + cov[3, 3] + 2 * cov[1, 3])
    exp_se = np.sqrt(1 / (len(r1) * r1.mean() * (1 - r1.mean())) + 1 / (len(r0) * r0.mean() * (1 - r0.mean())))
    return {**out, "fitted": True, "b": r.params[3], "se": r.bse[3], "expected_se": exp_se,
            "eff_comparator": np.exp(r.params[1]), "eff_drug_arm": np.exp(r.params[1] + r.params[3]),
            "eff_drug_arm_lo": np.exp(r.params[1] + r.params[3] - 1.96 * se_arm1),
            "eff_drug_arm_hi": np.exp(r.params[1] + r.params[3] + 1.96 * se_arm1),
            "eff_comparator_lo": np.exp(r.params[1] - 1.96 * r.bse[1]), "eff_comparator_hi": np.exp(r.params[1] + 1.96 * r.bse[1]),
            "covariates": ";".join(used)}


def one_sided(b, se, kind):
    z = b / se
    return float(stats.norm.sf(z)) if kind == "tte" else float(stats.norm.cdf(z))  # tte: HR>1; bin: OR<1


def main():
    H = C.frozen()
    P, S = C.load_table()
    P = P[P["source"] != "TCGA"]
    S = S.loc[P.index]
    c, s = storm()
    s = s.reindex(columns=S.columns)
    P, S = pd.concat([P, c]), pd.concat([S, s])
    coh = {cid: cohort(P, S, cid, H) for cid in sorted({x[1] for x in CONTROLLED} | {x[1] for x in SINGLE})}
    rows = []
    for drug, cid, arm_d, arm_c, kind, ep, role, rand, filt in CONTROLLED:
        p = coh[cid]
        if filt == "her2pos":
            p = p[p["her2"].astype(str).isin(["1", "1.0", "positive"])]
        p = p[p["arm"].isin(arm_d + arm_c)]
        arm = p["arm"].isin(arm_d).astype(float)
        sc = p[f"score_{drug}"]
        r = fit_tte(p, sc, arm, ep) if kind == "tte" else fit_bin(p, sc, arm)
        row = {"drug": drug, "cohort": cid, "drug_arm": "|".join(arm_d), "comparator": "|".join(arm_c),
               "endpoint": ep.upper() if kind == "tte" else "response", "kind": kind, "role": role, "randomised": rand, **r}
        if r.get("fitted"):
            row.update({"effect": "HR_int" if kind == "tte" else "OR_int", "estimate": np.exp(r["b"]),
                        "lo": np.exp(r["b"] - 1.96 * r["se"]), "hi": np.exp(r["b"] + 1.96 * r["se"]),
                        "p_one_sided": one_sided(r["b"], r["se"], kind), "p_two_sided": 2 * stats.norm.sf(abs(r["b"] / r["se"])),
                        "direction_as_expected": (r["b"] > 0) if kind == "tte" else (r["b"] < 0),
                        "power_obs": C.power_one_sided(LOG_EFF, r["se"]), "power_expected": C.power_one_sided(LOG_EFF, r["expected_se"]),
                        "min_detectable_obs": np.exp((stats.norm.ppf(0.975) + stats.norm.ppf(0.8)) * r["se"])})
        rows.append(row)
    T = pd.DataFrame(rows)
    m = T["role"].eq("primary") & T["fitted"].eq(True)
    T.loc[m, "holm_all_primary_stage2"] = C.holm(T.loc[m, "p_one_sided"].values)
    T.to_csv(C.TAB / "stage2_controlled.tsv", sep="\t", index=False)

    srows = []
    for drug, cid, arms, kind, ep in SINGLE:
        p = coh[cid]
        if arms is not None:
            p = p[p["arm"].isin(arms)]
        sc = p[f"score_{drug}"]
        r = fit_tte(p, sc, None, ep) if kind == "tte" else fit_bin(p, sc, None)
        srows.append({"drug": drug, "cohort": cid, "arm": "all" if arms is None else "|".join(arms),
                      "endpoint": ep.upper() if kind == "tte" else "response", "kind": kind, **r,
                      "effect": "HR" if kind == "tte" else "OR", "estimate": np.exp(r["b"]),
                      "lo": np.exp(r["b"] - 1.96 * r["se"]), "hi": np.exp(r["b"] + 1.96 * r["se"]),
                      "p_one_sided": one_sided(r["b"], r["se"], kind),
                      "direction_as_expected": (r["b"] > 0) if kind == "tte" else (r["b"] < 0),
                      "power_obs": C.power_one_sided(LOG_EFF, r["se"])})
    pd.DataFrame(srows).to_csv(C.TAB / "stage2_single_arm_supportive.tsv", sep="\t", index=False)

    # pooling and gate G3 (primary, fitted, non-exploratory)
    pool, gate = [], []
    for drug in C.DRUGS:
        g = T[(T["drug"] == drug) & T["role"].eq("primary") & T["fitted"].eq(True)]
        for kind, gk in g.groupby("kind"):
            sign = 1 if kind == "tte" else -1   # put both on "positive = expected direction"
            m = C.dersimonian_laird(gk["b"].values, gk["se"].values)
            pool.append({"drug": drug, "kind": kind, "k": m["k"], "cohorts": ";".join(gk["cohort"]),
                         "FE": np.exp(m["fe"]), "FE_lo": np.exp(m["fe"] - 1.96 * m["se_fe"]), "FE_hi": np.exp(m["fe"] + 1.96 * m["se_fe"]),
                         "FE_p_one_sided": float(stats.norm.sf(sign * m["fe"] / m["se_fe"])),
                         "RE": np.exp(m["re"]), "RE_lo": np.exp(m["re"] - 1.96 * m["se_re"]), "RE_hi": np.exp(m["re"] + 1.96 * m["se_re"]),
                         "RE_p_one_sided": float(stats.norm.sf(sign * m["re"] / m["se_re"])), "I2": m["I2"], "tau2": m["tau2"]})
        if len(g) == 0:
            gate.append({"drug": drug, "n_controlled_primary": 0, "G3_pass": False, "note": "no controlled external cohort"})
            continue
        powered_hit = g[(g["p_one_sided"] < 0.025) & (g["power_obs"] >= 0.8)]
        same_dir = bool(g["direction_as_expected"].all())
        gate.append({"drug": drug, "n_controlled_primary": len(g), "n_adequately_powered": int((g["power_obs"] >= 0.8).sum()),
                     "n_one_sided_p_lt_0.025": int((g["p_one_sided"] < 0.025).sum()),
                     "n_same_direction": int(g["direction_as_expected"].sum()),
                     "G3_powered_hit": len(powered_hit) > 0, "G3_all_same_direction": same_dir,
                     "G3_pass": (len(powered_hit) > 0) and same_dir, "note": ""})
    pd.DataFrame(pool).to_csv(C.TAB / "stage2_pooled.tsv", sep="\t", index=False)
    pd.DataFrame(gate).to_csv(C.TAB / "stage2_gate_G3.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(T[["drug", "cohort", "endpoint", "role", "n", "fitted", "estimate", "lo", "hi", "p_one_sided", "power_obs"]].to_string())
    print(pd.DataFrame(gate).to_string())


if __name__ == "__main__":
    main()
