"""A3: statistical controls across all 57 organoid-only analyses.

Inputs: per-patient tables from a3_build_gi.py / a3_build_other.py (exact reproductions of the organoid-only
results), the two cross-set summary tables (permutation and random-pathway bootstrap p), TCGA-CDR (age, sex, stage),
ABSOLUTE purity, MSI (Xena COADREAD clinical matrix + 2019 BCR biotab), and drug records (regimen).
Outputs: tables/A3_master_results.tsv and supporting tables.
"""
import glob

import numpy as np
import pandas as pd
from scipy.stats import binom

from common import (DATA, ORG, ORG2, REPO, TAB, bh, cox, holm, logrank, mdhr, power_logrank, retrodesign,
                    stage_num)

np.random.seed(20261117)
EXCL_CODRUG = {"LEUCOVORIN", "__ANY_RECORD__"}


# ------------------------------------------------------------------------------------------ inputs
def summaries():
    gi = pd.read_csv(ORG / "tables" / "cross_set_summary.tsv", sep="\t")
    gi["analysis_id"] = gi["folder"] + "__" + gi["cohort"] + "__" + gi["model_drug"]
    gi = gi.rename(columns={"treated_HR": "hr", "treated_HR_lo": "hr_lo", "treated_HR_hi": "hr_hi",
                            "shuffle_share_rank1": "perm_p", "bootstrap_share": "boot_p", "rank1_pathway": "pathway"})
    gi["set_label"] = gi["set"]
    gi["drug"] = gi["patient_drug"]
    gi["n_boot"] = 10000
    gi["project"] = gi["cohort"]
    ot = pd.read_csv(ORG2 / "SUMMARY" / "tables" / "cross_set_summary.tsv", sep="\t")
    ot["analysis_id"] = ot["set"] + "__" + ot["project"] + "__" + ot["drug"]
    ot = ot.rename(columns={"treated_hr": "hr", "treated_hr_lo": "hr_lo", "treated_hr_hi": "hr_hi",
                            "perm_empirical_p": "perm_p", "boot_empirical_p_rank1p": "boot_p", "biomarker": "pathway",
                            "untreated_hr": "untreated_HR", "untreated_hr_lo": "untreated_HR_lo",
                            "untreated_hr_hi": "untreated_HR_hi"})
    ot["n_boot"] = np.where(ot["drug"].eq("CISPLATIN") & ot["set"].str.startswith("BLCA"), 10000, 2000)
    cols = ["analysis_id", "set_label", "project", "drug", "model_drug", "n_models", "pathway", "treated_n",
            "treated_events", "treated_p", "hr", "hr_lo", "hr_hi", "untreated_n", "untreated_events", "untreated_p",
            "untreated_HR", "untreated_HR_lo", "untreated_HR_hi", "perm_p", "boot_p", "n_boot"]
    return pd.concat([gi[cols], ot[cols]], ignore_index=True)


def clinical():
    cdr = pd.read_csv(ORG / "_data" / "raw" / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv", sep="\t",
                      dtype=str).drop_duplicates("_PATIENT").set_index("_PATIENT")
    st = cdr["ajcc_pathologic_tumor_stage"].map(stage_num)
    st = st.fillna(cdr["clinical_stage"].map(stage_num))
    cl = pd.DataFrame({"age": pd.to_numeric(cdr["age_at_initial_pathologic_diagnosis"], errors="coerce"),
                       "male": cdr["gender"].str.upper().map({"MALE": 1.0, "FEMALE": 0.0}), "stage": st})
    ab = pd.read_csv(REPO / "robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    ab = ab[ab["array"].str[13:15] == "01"].copy()
    ab["patient"] = ab["array"].str[:12]
    cl["purity"] = ab.groupby("patient")["purity"].mean().reindex(cl.index)
    # MSI: Xena COADREAD matrix (MSI_updated_Oct62011, then BCR field), then 2019 COAD biotab
    x = pd.read_csv(DATA / "raw" / "COADREAD_clinicalMatrix", sep="\t")
    x["patient"] = x["sampleID"].str[:12]
    m1 = x["MSI_updated_Oct62011"].map({"MSI-H": 1.0, "MSI-L": 0.0, "MSS": 0.0})
    m2 = x["microsatellite_instability"].map({"YES": 1.0, "NO": 0.0})
    x["msi"] = m1.fillna(m2)
    msi = x.dropna(subset=["msi"]).groupby("patient")["msi"].max()
    b = pd.read_csv(glob.glob(str(REPO / "data/TCGA_COAD/clinical_patient/*/nationwidechildrens.org_clinical_patient_coad.txt"))[0],
                    sep="\t", skiprows=[1, 2]).set_index("bcr_patient_barcode")
    msi = msi.combine_first(b["microsatellite_instability"].map({"YES": 1.0, "NO": 0.0}).dropna())
    cl["msi"] = msi.reindex(cl.index)
    return cl


def add_cov(df, name, v, min_known=0.5):
    """Add covariate v (Series aligned to df) with median imputation + missing indicator. Returns used names."""
    known = v.notna().mean()
    if known < min_known or v.nunique(dropna=True) < 2:
        return []
    out = [name]
    if v.isna().any():
        df[name + "_missing"] = v.isna().astype(float)
        if df[name + "_missing"].sum() >= 3:
            out.append(name + "_missing")
        else:
            df.drop(columns=name + "_missing", inplace=True)
    df[name] = v.fillna(v.median()).astype(float)
    return out


def design(d, cl, project, drug, arm_specific=True):
    """Covariate matrix for patients d (treated arm): age/10, sex, stage, purity, MSI (colorectal), regimen."""
    X = pd.DataFrame(index=d.index)
    c = cl.reindex(d["patient"].values)
    c.index = d.index
    used = []
    used += add_cov(X, "age10", c["age"] / 10)
    if c["male"].nunique() == 2 and min((c["male"] == 1).sum(), (c["male"] == 0).sum()) >= 5:
        used += add_cov(X, "male", c["male"])
    used += add_cov(X, "stage", c["stage"])
    used += add_cov(X, "purity", c["purity"])
    if project in ("COAD", "READ", "COADREAD", "COAD_authors"):
        used += add_cov(X, "msi", c["msi"], min_known=0.3)
    reg = []
    if arm_specific:
        dr = d["codrugs"].fillna("").str.split(";")
        counts = pd.Series([x for L in dr for x in set(L) if x and x not in EXCL_CODRUG and x != drug]).value_counts()
        n = len(d)
        for k, v in counts.items():
            if 0.10 * n <= v <= 0.90 * n and len(reg) < 2:
                X["reg_" + k] = dr.map(lambda L: float(k in L))
                reg.append("reg_" + k)
    X = X[used + reg]
    # drop binary covariates that cannot be estimated (a level with < 3 patients or no events)
    E = d["E"].values
    for c in list(X.columns):
        v = X[c].values
        if set(np.unique(v)) <= {0.0, 1.0}:
            if min(v.sum(), (1 - v).sum()) < 3 or E[v == 1].sum() == 0 or E[v == 0].sum() == 0:
                X = X.drop(columns=c)
    used = [c for c in used if c in X.columns]
    reg = [c for c in reg if c in X.columns]
    return X, used, reg


# ------------------------------------------------------------------------------------------ main
def main():
    S = summaries()
    pat = pd.concat([pd.read_csv(DATA / f"a3_patients_{k}.tsv.gz", sep="\t") for k in ("gi", "other")])
    cl = clinical()
    rows, covrows = [], []
    for _, s in S.iterrows():
        aid = s["analysis_id"]
        A = pat[pat["analysis_id"] == aid].reset_index(drop=True)
        t = A[A.group == "treated"].reset_index(drop=True)
        u = A[A.group == "untreated"].reset_index(drop=True)
        r = s.to_dict()
        D = int(t.E.sum())
        r["events_nonresp"] = int(t.loc[t.split == 1, "E"].sum())
        r["events_resp"] = int(t.loc[t.split == 0, "E"].sum())
        r["logrank_p_check"] = logrank(t["T"], t["E"], t["split"])
        if np.isfinite(r["hr"]) and 1e-3 < r["hr"] < 1e3:
            r["direction"] = "expected (HR>1)" if r["hr"] > 1 else "opposite (HR<1)"
        else:
            r["direction"] = ("expected (HR>1)" if r["events_nonresp"] > r["events_resp"] else "opposite (HR<1)") + \
                             " [HR not estimable]"
        r["multi_k10_logrank_p"] = logrank(t["T"], t["E"], (t["score10"] > t["score10"].median()).astype(int))
        # continuous score, per SD (treated)
        z = (t["score1"] - t["score1"].mean()) / t["score1"].std()
        c0 = cox(t["T"], t["E"], z.values, ["score"])
        if c0 is not None:
            r.update(cont_HR=c0.loc["score", "HR"], cont_lo=c0.loc["score", "lo"], cont_hi=c0.loc["score", "hi"],
                     cont_p=c0.loc["score", "p"])
        # adjusted Cox
        X, used, reg = design(t, cl, s["project"], s["drug"])
        r["covariates"] = ";".join(used + reg)
        r["n_covariates"] = X.shape[1]
        r["events_per_variable"] = D / (X.shape[1] + 1)
        keep = X.columns[X.std() > 0]
        X = X[keep]
        ca = cox(t["T"], t["E"], np.column_stack([t["split"].values, X.values]), ["split"] + list(X.columns))
        if ca is not None:
            r.update(adj_HR=ca.loc["split", "HR"], adj_lo=ca.loc["split", "lo"], adj_hi=ca.loc["split", "hi"],
                     adj_p=ca.loc["split", "p"])
            for cv in ca.index:
                covrows.append({"analysis_id": aid, "term": cv, **ca.loc[cv].to_dict()})
        cac = cox(t["T"], t["E"], np.column_stack([z.values, X.values]), ["score"] + list(X.columns))
        if cac is not None:
            r.update(adj_cont_HR=cac.loc["score", "HR"], adj_cont_lo=cac.loc["score", "lo"],
                     adj_cont_hi=cac.loc["score", "hi"], adj_cont_p=cac.loc["score", "p"])
        # clinical-only (no regimen) adjusted, as a sensitivity
        Xc = X[[c for c in X.columns if not c.startswith("reg_")]]
        cb = cox(t["T"], t["E"], np.column_stack([t["split"].values, Xc.values]), ["split"] + list(Xc.columns))
        if cb is not None:
            r.update(adjclin_HR=cb.loc["split", "HR"], adjclin_lo=cb.loc["split", "lo"],
                     adjclin_hi=cb.loc["split", "hi"], adjclin_p=cb.loc["split", "p"])
        # interaction score x treatment (treated vs untreated = no drug record)
        if len(u) >= 6:
            B = pd.concat([t.assign(arm=1.0), u.assign(arm=0.0)], ignore_index=True)
            zz = (B["score1"] - B["score1"].mean()) / B["score1"].std()
            ci = cox(B["T"], B["E"], np.column_stack([zz, B["arm"], zz * B["arm"]]), ["score", "arm", "inter"], check=3)
            if ci is not None:
                r.update(int_HR=ci.loc["inter", "HR"], int_lo=ci.loc["inter", "lo"], int_hi=ci.loc["inter", "hi"],
                         int_p=ci.loc["inter", "p"], untreated_cont_HR=ci.loc["score", "HR"])
            Xb, ub, _ = design(B, cl, s["project"], s["drug"], arm_specific=False)
            Xb = Xb.loc[:, Xb.std() > 0]
            cia = cox(B["T"], B["E"], np.column_stack([zz, B["arm"], zz * B["arm"], Xb.values]),
                      ["score", "arm", "inter"] + list(Xb.columns), check=3)
            if cia is not None:
                r.update(int_adj_HR=cia.loc["inter", "HR"], int_adj_lo=cia.loc["inter", "lo"],
                         int_adj_hi=cia.loc["inter", "hi"], int_adj_p=cia.loc["inter", "p"])
            cs = cox(B["T"], B["E"], np.column_stack([B["split"], B["arm"], B["split"] * B["arm"]]),
                     ["split", "arm", "inter"])
            if cs is not None:
                r.update(int_split_HR=cs.loc["inter", "HR"], int_split_p=cs.loc["inter", "p"])
            r["interaction_n"] = len(B)
            r["interaction_events"] = int(B.E.sum())
        # power (Schoenfeld, 1:1 median split, two-sided 0.05)
        r["power_HR1.5"] = power_logrank(D, 1.5)
        r["power_HR2.0"] = power_logrank(D, 2.0)
        r["MDHR_80"] = mdhr(D)
        r["se_logHR_schoenfeld"] = 2 / np.sqrt(D) if D else np.nan
        rows.append(r)
    M = pd.DataFrame(rows)
    M["BH_q"] = bh(M["treated_p"])
    M["Holm_p"] = holm(M["treated_p"])
    M["BH_q_perm"] = bh(M["perm_p"])
    M["BH_q_interaction"] = bh(M["int_p"])
    M["Holm_interaction"] = holm(M["int_p"])
    M.to_csv(TAB / "A3_master_results.tsv", sep="\t", index=False)
    pd.DataFrame(covrows).to_csv(TAB / "A3_adjusted_cox_all_terms.tsv", sep="\t", index=False)

    # expected-by-chance
    n = M["treated_p"].notna().sum()
    obs = int((M["treated_p"] < 0.05).sum())
    exp = 0.05 * n
    ch = pd.DataFrame([{"family": "treated log-rank p (rank-1)", "n_tests": n, "observed_p_below_0.05": obs,
                        "expected_by_chance": exp, "binomial_p_at_least_observed": float(binom.sf(obs - 1, n, 0.05)),
                        "observed_p_below_0.01": int((M["treated_p"] < 0.01).sum()), "expected_0.01": 0.01 * n,
                        "BH_q_below_0.05": int((M["BH_q"] < 0.05).sum()), "BH_q_below_0.10": int((M["BH_q"] < 0.10).sum()),
                        "Holm_below_0.05": int((M["Holm_p"] < 0.05).sum()),
                        "expected_direction_among_p_below_0.05": int(((M["treated_p"] < 0.05) & M["direction"].str.startswith("expected")).sum())}])
    for fam, col in (("permutation p (response shuffle)", "perm_p"), ("random-pathway bootstrap p", "boot_p"),
                     ("untreated log-rank p", "untreated_p"), ("interaction p (score x treated)", "int_p"),
                     ("multi-pathway k=10 log-rank p", "multi_k10_logrank_p")):
        v = M[col].dropna()
        o = int((v < 0.05).sum())
        ch = pd.concat([ch, pd.DataFrame([{"family": fam, "n_tests": len(v), "observed_p_below_0.05": o,
                                           "expected_by_chance": 0.05 * len(v),
                                           "binomial_p_at_least_observed": float(binom.sf(o - 1, len(v), 0.05)),
                                           "observed_p_below_0.01": int((v < 0.01).sum()), "expected_0.01": 0.01 * len(v),
                                           "BH_q_below_0.05": int((bh(v) < 0.05).sum()),
                                           "BH_q_below_0.10": int((bh(v) < 0.10).sum()),
                                           "Holm_below_0.05": int((holm(v) < 0.05).sum())}])], ignore_index=True)
    ch.to_csv(TAB / "A3_expected_by_chance.tsv", sep="\t", index=False)

    # winner's curse for nominal hits
    wc = []
    for _, r in M[M["treated_p"] < 0.05].iterrows():
        s_sch = r["se_logHR_schoenfeld"]
        s_cox = np.log(r["hr_hi"] / r["hr_lo"]) / (2 * 1.96) if np.isfinite(r["hr"]) else np.nan
        for true_hr in (1.25, 1.5, 2.0):
            pw, ts, ex = retrodesign(np.log(true_hr), s_sch)
            wc.append({"analysis_id": r["analysis_id"], "treated_n": r["treated_n"], "events": r["treated_events"],
                       "treated_p": r["treated_p"], "observed_HR": r["hr"], "observed_HR_lo": r["hr_lo"],
                       "observed_HR_hi": r["hr_hi"], "direction": r["direction"], "se_logHR_used": s_sch,
                       "se_logHR_cox": s_cox, "assumed_true_HR": true_hr, "power": pw, "type_S_error": ts,
                       "exaggeration_ratio_typeM": ex,
                       "expected_significant_HR": float(np.exp(np.log(true_hr) * ex))})
    pd.DataFrame(wc).to_csv(TAB / "A3_winners_curse.tsv", sep="\t", index=False)
    print(ch.to_string())
    print(M[["analysis_id", "treated_n", "treated_events", "treated_p", "hr", "adj_HR", "adj_p", "int_HR", "int_p",
             "BH_q", "MDHR_80", "covariates"]].to_string())


if __name__ == "__main__":
    main()
