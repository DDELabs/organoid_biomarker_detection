"""Audit 2: events per cancer x drug in TCGA, power for main effects and treatment x biomarker
interactions, and testability of every row of FINAL_BIOMARKER_SURVIVAL.

Targets: per-SD HR 0.7 (main effect, treated patients) and interaction HR 0.7 (treated vs untreated),
two-sided alpha 0.05, 80% power; also at an FDR-controlling per-test alpha for ~40 rows
(Bonferroni 0.05/40 and Jung-2005 BH-equivalent alpha with 10% true effects, q = 0.1).

Outputs (robust/revamp_audit/out/):
  power_targets.tsv            required events / patients for each target
  tcga_events_by_pair.tsv      every cancer x drug with >= 10 treated patients (with expression)
  tcga_events_by_drug_pooled.tsv  pooled across cancers (stratified analysis), per drug
  final_rows_testability.tsv   the 43 FINAL_BIOMARKER_SURVIVAL rows: events, power, MDHR, testability
  followup_by_cancer.tsv       median potential follow-up (reverse KM), OS/PFI event rates
  power_simulation.tsv         Monte-Carlo check of the interaction formula (statsmodels PHReg)
"""
import json

import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg

from common import (OUT, REPO, LOG_HR, R, events_interaction, events_main, expr_patients, jung_alpha, mdhr, n_for_auc,
                    power_auc, power_interaction, power_main, reverse_km_median, systemic_long)

M_ROWS = 40
A_BONF = 0.05 / M_ROWS
A_JUNG = jung_alpha(M_ROWS, 4, q=0.10, power=0.8)
A_JUNG05 = jung_alpha(M_ROWS, 4, q=0.05, power=0.8)


def targets():
    rows = []
    for label, a in [("alpha 0.05", 0.05), (f"BH q=0.10, 4/40 true (alpha*={A_JUNG:.4f})", A_JUNG),
                     (f"BH q=0.05, 4/40 true (alpha*={A_JUNG05:.4f})", A_JUNG05), ("Bonferroni 0.05/40", A_BONF)]:
        rows.append({"target": "main effect HR 0.7/SD", "alpha": label, "events": events_main(0.7, a)})
        rows.append({"target": "main effect HR 0.7/SD, marker R2=0.3 with covariates", "alpha": label, "events": events_main(0.7, a, r2=0.3)})
        rows.append({"target": "main effect HR 0.8/SD", "alpha": label, "events": events_main(0.8, a)})
        rows.append({"target": "interaction HR 0.7/SD (events 50:50 treated:untreated)", "alpha": label, "events": events_interaction(0.7, a)})
        rows.append({"target": "interaction HR 0.7/SD (events 20:80)", "alpha": label, "events": events_interaction(0.7, a, frac_events_treated=0.2)})
        rows.append({"target": "interaction HR 0.8/SD (50:50)", "alpha": label, "events": events_interaction(0.8, a)})
        for auc in (0.56, 0.60, 0.65, 0.70):
            rows.append({"target": f"binary response AUROC {auc} (prevalence 0.5)", "alpha": label, "patients": n_for_auc(auc, 0.5, a)})
            rows.append({"target": f"binary response AUROC {auc} (prevalence 0.2)", "alpha": label, "patients": n_for_auc(auc, 0.2, a)})
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "power_targets.tsv", sep="\t", index=False)
    return T


def simulate_interaction(n_rep=300, n=1000, frac_treated=0.5, hr_int=0.7, seed=1):
    """Exponential survival, marker N(0,1), prognostic marker effect HR 0.85 in both arms,
    treatment HR 0.8, interaction hr_int; uniform censoring tuned so ~50% events."""
    rng = np.random.default_rng(seed)
    hits, evs = 0, []
    for _ in range(n_rep):
        x = rng.standard_normal(n)
        t = (rng.random(n) < frac_treated).astype(float)
        lp = np.log(0.85) * x + np.log(0.8) * t + np.log(hr_int) * x * t
        T = rng.exponential(1 / np.exp(lp))
        C = rng.uniform(0, 1.6, n)
        time, ev = np.minimum(T, C), (T <= C).astype(float)
        X = np.column_stack([x * t, x, t])
        r = PHReg(time, X, status=ev).fit(disp=0)
        hits += r.pvalues[0] < 0.05
        evs.append(ev.sum())
    D = float(np.mean(evs))
    return {"n": n, "mean_events": D, "sim_power": hits / n_rep, "formula_power": power_main(D * 0.25, hr_int),
            "formula_events_for_80pct": events_interaction(hr_int)}


def main():
    T = targets()
    print(T.round(1).to_string())
    long = systemic_long()
    cdr = R.load_cdr()
    systemic = set(long["patient"])
    pairs, fu = [], []
    for proj in sorted(long["project"].unique()):
        pts = expr_patients(proj)
        if not pts:
            continue
        c = cdr.reindex(sorted(pts))
        c = c[c["OS.time"].notna()]
        untreated = c.index.difference(list(systemic))
        fu.append({"project": proj, "n_expr": len(c), "OS_events": int(c["OS"].sum()), "PFI_events": int(c["PFI"].sum()),
                   "OS_event_rate": float(c["OS"].mean()), "median_followup_months": reverse_km_median(c["OS.time"], c["OS"]) / 30.44,
                   "n_untreated": len(untreated), "untreated_OS_events": int(c.loc[untreated, "OS"].sum()),
                   "untreated_PFI_events": int(c.loc[untreated, "PFI"].sum())})
        Du_os, Du_pfi = c.loc[untreated, "OS"].sum(), c.loc[untreated, "PFI"].sum()
        g = long[long["project"] == proj]
        for drug, gd in g.groupby("drug"):
            tr = sorted(set(gd["patient"]) & set(c.index))
            if len(tr) < 10:
                continue
            labelled = set(gd.loc[gd["recist"].notna(), "patient"]) & set(tr)
            held = sorted(set(tr) - labelled)
            lab = gd[gd["recist"].notna()].drop_duplicates("patient")
            lab = lab[lab["patient"].isin(tr)]
            n1, n0 = int(lab["responder"].sum()), int((lab["responder"] == 0).sum())
            ct, ch = c.loc[tr], c.loc[held]
            Dt_os, Dt_pfi = ct["OS"].sum(), ct["PFI"].sum()
            Dh_os, Dh_pfi = ch["OS"].sum(), ch["PFI"].sum()
            pairs.append({
                "project": proj, "drug": drug, "n_treated": len(tr), "OS_events_treated": int(Dt_os), "PFI_events_treated": int(Dt_pfi),
                "n_heldout": len(held), "OS_events_heldout": int(Dh_os), "PFI_events_heldout": int(Dh_pfi),
                "n_untreated": len(untreated), "OS_events_untreated": int(Du_os), "PFI_events_untreated": int(Du_pfi),
                "recist_n": n1 + n0, "recist_responders": n1, "recist_nonresponders": n0,
                "power_main_PFI_treated": power_main(Dt_pfi), "power_main_OS_treated": power_main(Dt_os),
                "power_main_PFI_treated_FDR": power_main(Dt_pfi, alpha=A_JUNG),
                "power_int_PFI": power_interaction(Dt_pfi, Du_pfi), "power_int_PFI_FDR": power_interaction(Dt_pfi, Du_pfi, alpha=A_JUNG),
                "MDHR_PFI_treated": mdhr(Dt_pfi), "MDHR_OS_heldout": mdhr(Dh_os),
                "power_recist_auc065": power_auc(n1, n0, 0.65), "power_recist_auc060": power_auc(n1, n0, 0.60),
            })
    P = pd.DataFrame(pairs)

    def testable(r):
        best = max(r["PFI_events_treated"], r["OS_events_treated"])
        tags = []
        if best >= events_main(0.7):
            tags.append("main@0.05")
        if best >= events_main(0.7, A_JUNG):
            tags.append("main@FDR")
        Dt, Du = r["PFI_events_treated"], r["PFI_events_untreated"]
        if power_interaction(Dt, Du) >= 0.8:
            tags.append("interaction@0.05")
        if power_interaction(Dt, Du, alpha=A_JUNG) >= 0.8:
            tags.append("interaction@FDR")
        if r["power_recist_auc065"] >= 0.8:
            tags.append("RECIST AUROC0.65")
        return ";".join(tags) or "none"

    P["testable"] = P.apply(testable, axis=1)
    P = P.sort_values("PFI_events_treated", ascending=False)
    P.to_csv(OUT / "tcga_events_by_pair.tsv", sep="\t", index=False)
    F = pd.DataFrame(fu)
    F.to_csv(OUT / "followup_by_cancer.tsv", sep="\t", index=False)

    # pooled per drug across cancers (stratified Cox: information adds across strata)
    pooled = P.groupby("drug").agg(cancers=("project", "nunique"), n_treated=("n_treated", "sum"),
                                   OS_events_treated=("OS_events_treated", "sum"), PFI_events_treated=("PFI_events_treated", "sum"),
                                   PFI_events_untreated=("PFI_events_untreated", "sum"), recist_n=("recist_n", "sum"),
                                   recist_responders=("recist_responders", "sum")).reset_index()
    pooled["power_main_PFI"] = pooled["PFI_events_treated"].map(power_main)
    pooled["power_int_PFI"] = [power_interaction(a, b) for a, b in zip(pooled["PFI_events_treated"], pooled["PFI_events_untreated"])]
    pooled["MDHR_PFI"] = pooled["PFI_events_treated"].map(mdhr)
    pooled = pooled.sort_values("PFI_events_treated", ascending=False)
    pooled.to_csv(OUT / "tcga_events_by_drug_pooled.tsv", sep="\t", index=False)

    # ---------------- the FINAL table rows
    W = pd.read_csv(REPO / "robust/results/survival_validation/FINAL_BIOMARKER_SURVIVAL.tsv", sep="\t")
    rows = []
    for _, r in W.iterrows():
        d = {"cancer": r["cancer"], "drug": r["drug"], "verdict": r["verdict"]}
        for ep in ("OS", "PFI"):
            Dh = r.get(f"{ep}_heldout_events")
            Du = r.get(f"{ep}_untreated_events")
            d[f"{ep}_events_heldout"] = Dh
            d[f"{ep}_events_untreated"] = Du
            d[f"{ep}_power_main_HR0.7"] = power_main(Dh) if pd.notna(Dh) else np.nan
            d[f"{ep}_power_int_HR0.7"] = power_interaction(Dh, Du) if pd.notna(Dh) and pd.notna(Du) else np.nan
            d[f"{ep}_MDHR_main"] = mdhr(Dh) if pd.notna(Dh) else np.nan
            if pd.notna(Dh) and pd.notna(Du) and Dh > 0 and Du > 0:
                d[f"{ep}_MD_interaction_HR"] = float(np.exp(-2.8016 * np.sqrt(1 / Dh + 1 / Du)))
        best = np.nanmax([d.get("OS_power_main_HR0.7", np.nan), d.get("PFI_power_main_HR0.7", np.nan)]) if any(
            pd.notna([d.get("OS_power_main_HR0.7"), d.get("PFI_power_main_HR0.7")])) else np.nan
        d["best_power_main"] = best
        d["power_class"] = ("adequate (>=80%)" if best >= 0.8 else "marginal (50-80%)" if best >= 0.5 else
                            "underpowered (<50%)" if best >= 0 else "no data")
        rows.append(d)
    FR = pd.DataFrame(rows)
    FR.to_csv(OUT / "final_rows_testability.tsv", sep="\t", index=False)

    sim = [simulate_interaction(n=n) for n in (500, 1000)]
    pd.DataFrame(sim).to_csv(OUT / "power_simulation.tsv", sep="\t", index=False)

    s = {
        "events_main_HR0.7_alpha0.05": events_main(0.7), "events_main_HR0.7_FDR": events_main(0.7, A_JUNG),
        "events_main_HR0.7_bonf": events_main(0.7, A_BONF), "events_int_HR0.7_alpha0.05": events_interaction(0.7),
        "events_int_HR0.7_FDR": events_interaction(0.7, A_JUNG), "alpha_jung_q0.1": A_JUNG,
        "n_pairs_ge10_treated": int(len(P)),
        "n_pairs_main_testable_0.05": int(P["testable"].str.contains("main@0.05").sum()),
        "n_pairs_main_testable_FDR": int(P["testable"].str.contains("main@FDR").sum()),
        "n_pairs_int_testable_0.05": int(P["testable"].str.contains("interaction@0.05").sum()),
        "n_pairs_int_testable_FDR": int(P["testable"].str.contains("interaction@FDR").sum()),
        "n_pairs_recist_auc065": int(P["testable"].str.contains("RECIST").sum()),
        "final_rows_power_class": FR["power_class"].value_counts().to_dict(),
        "final_rows_median_heldout_OS_events": float(FR["OS_events_heldout"].median()),
        "final_rows_median_heldout_PFI_events": float(FR["PFI_events_heldout"].median()),
        "final_rows_max_heldout_OS_events": float(FR["OS_events_heldout"].max()),
        "final_rows_median_MDHR_OS": float(FR["OS_MDHR_main"].median()),
        "final_rows_n_power_int_ge0.8": int((FR[["OS_power_int_HR0.7", "PFI_power_int_HR0.7"]].max(axis=1) >= 0.8).sum()),
        "final_rows_median_power_int": float(FR[["OS_power_int_HR0.7", "PFI_power_int_HR0.7"]].max(axis=1).median()),
        "simulation": sim,
    }
    json.dump(s, open(OUT / "power_summary.json", "w"), indent=1, default=float)
    print(json.dumps(s, indent=1, default=float))
    print(P.head(40)[["project", "drug", "n_treated", "OS_events_treated", "PFI_events_treated", "PFI_events_untreated", "recist_n",
                      "recist_responders", "power_main_PFI_treated", "power_int_PFI", "testable"]].round(2).to_string())
    print(pooled.head(20).round(2).to_string())
    print(FR.round(2).to_string())
    print(F.round(1).to_string())


if __name__ == "__main__":
    main()
