"""Audit 3: purity / composition confounding, prognostic confounding (LGG), multiple testing,
discovery instability, LOCO AUROC noise, and external-cohort dominant-axis checks.

Tumour purity = TCGA PanCanAtlas ABSOLUTE (TCGA_mastercalls.abs_tables_JSedit.fixed.txt, GDC file
4f277128-f793-4354-a13d-30cc7fe9f6b5), copied to robust/revamp_audit/data/.

Outputs (robust/revamp_audit/out/):
  purity_pathway_rho.tsv          per pathway: median within-cancer Spearman rho with purity
  purity_biomarker_rows.tsv       per FINAL row: score-purity rho, held-out Cox HR without / with purity
  purity_response.tsv             purity vs RECIST responder (within-cancer AUROC) and PFI among treated
  lgg_tmz_confounding.tsv         LGG TMZ score vs grade/histology; HR before/after grade adjustment
  multiple_testing.tsv            all p-values in FINAL_BIOMARKER_SURVIVAL with BH q
  loco_noise.tsv                  observed vs sampling variance of per-pair LOCO AUROCs
  external_dominant_axis.tsv      PC1 of pathway scores vs response in external cohorts used as evidence
  confounding_summary.json
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from common import HERE, OUT, REPO, R, ROBUST, bh, systemic_long

sys.path.insert(0, str(ROBUST))
import run_final_biomarker_survival as F  # noqa: E402
import run_triad as RT  # noqa: E402
from obd import triad as T  # noqa: E402

RES = REPO / "robust" / "results"


def load_purity():
    a = pd.read_csv(HERE / "data" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    a = a[a["call status"].eq("called") & a["array"].str[13:15].isin(["01", "03"])]
    a["patient"] = a["array"].str[:12]
    return a.drop_duplicates("patient").set_index("patient")["purity"].astype(float)


def main():
    summ = {}
    pur = load_purity()
    long = systemic_long()
    cdr = R.load_cdr()
    W = pd.read_csv(RES / "survival_validation/FINAL_BIOMARKER_SURVIVAL.tsv", sep="\t")
    projects = sorted(set(long["project"]) & {p.name.split("_")[1] for p in (RES / "cache").glob("TCGA_*_reactome_rank_all.tsv.gz")})
    M, proj = RT.patient_matrix(projects)
    Xc = T.center_within(M, proj)
    pp = pur.reindex(Xc.index)
    summ["patients_with_expr"] = int(len(Xc))
    summ["patients_with_purity"] = int(pp.notna().sum())

    # ---------------- 1. pathway-purity correlation (within cancer)
    rhos = {}
    for p in sorted(set(proj)):
        idx = proj.index[(proj == p) & pp.notna().values]
        if len(idx) < 30:
            continue
        r = Xc.loc[idx].rank().corrwith(pp[idx].rank())
        rhos[p] = r
    RH = pd.DataFrame(rhos)
    out = pd.DataFrame({"median_rho": RH.median(1), "q10": RH.quantile(0.1, axis=1), "q90": RH.quantile(0.9, axis=1),
                        "frac_cancers_abs_gt_0.3": (RH.abs() > 0.3).mean(1)}).sort_values("median_rho")
    out.to_csv(OUT / "purity_pathway_rho.tsv", sep="\t")
    summ["n_cancers_purity"] = int(RH.shape[1])
    summ["frac_pathways_abs_median_rho_gt_0.3"] = float((out["median_rho"].abs() > 0.3).mean())
    summ["frac_pathways_abs_median_rho_gt_0.5"] = float((out["median_rho"].abs() > 0.5).mean())
    P = "REACTOME_"
    hits = {"BCR signalling (paclitaxel response)": P + "SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR",
            "RIG-I/MDA5 IFN (cisplatin resistance)": P + "RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS",
            "Cell cycle (cisplatin response)": P + "CELL_CYCLE", "Fanconi anaemia (topo-II response)": P + "FANCONI_ANEMIA_PATHWAY",
            "FGFR3 mutants (cisplatin response)": P + "SIGNALING_BY_FGFR3_MUTANTS", "Prostanoid receptors (taxane resistance)": P + "PROSTANOID_LIGAND_RECEPTORS",
            "Semaphorin interactions (antimetabolite resistance)": P + "SEMAPHORIN_INTERACTIONS",
            "CTLA4 inhibitory signalling": P + "CTLA4_INHIBITORY_SIGNALING",
            "Immunoregulatory lymphoid/non-lymphoid": P + "IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
            "Collagen formation (stroma)": P + "COLLAGEN_FORMATION", "IFN-gamma signalling": P + "INTERFERON_GAMMA_SIGNALING",
            "Gap junction degradation": P + "GAP_JUNCTION_DEGRADATION", "Renal aquaporins": P + "REGULATION_OF_WATER_BALANCE_BY_RENAL_AQUAPORINS",
            "IRAK2-TAK1 (oxaliplatin)": P + "IRAK2_MEDIATED_ACTIVATION_OF_TAK1_COMPLEX_UPON_TLR7_8_OR_9_STIMULATION"}
    summ["atlas_hit_purity_rho"] = {k: round(float(out.loc[v, "median_rho"]), 3) for k, v in hits.items() if v in out.index}

    # ---------------- 2. purity vs response / PFI
    lab = long[long["recist"].notna()].sort_values("start_days").drop_duplicates(["patient", "drug"])
    lab = lab[lab["patient"].isin(Xc.index)]
    lab = lab.assign(purity=pur.reindex(lab["patient"]).values).dropna(subset=["purity"])
    pr = []
    for (p, d), g in lab.groupby(["project", "drug"]):
        if len(g) >= 20 and 5 <= g["responder"].sum() <= len(g) - 5:
            pr.append({"project": p, "drug": d, "n": len(g), "auc_purity_responder": roc_auc_score(g["responder"], g["purity"])})
    PR = pd.DataFrame(pr)
    PR.to_csv(OUT / "purity_response.tsv", sep="\t", index=False)
    summ["purity_auc_pairs"] = int(len(PR))
    summ["purity_auc_median"] = float(PR["auc_purity_responder"].median())
    summ["purity_auc_abs_dev_median"] = float((PR["auc_purity_responder"] - 0.5).abs().median())
    # pooled within-cancer purity (rank within cancer) vs responder, patient-level
    lab["pur_rank"] = lab.groupby("project")["purity"].rank(pct=True)
    summ["purity_auc_pooled_within_cancer_rank"] = float(roc_auc_score(lab["responder"], lab["pur_rank"]))

    # ---------------- 3. FINAL rows: score-purity rho, Cox with/without purity on held-out patients
    top = pd.concat([pd.read_csv(RES / "atlas/cancer_specific/top10_by_cancer_drug.tsv", sep="\t"),
                     pd.read_csv(RES / "atlas/cancer_specific/relaxed_top10.tsv", sep="\t")], ignore_index=True)
    systemic = set(long["patient"])
    rows = []
    for row in F.ROWS:
        cancer, drug, _, resp, res, _, rule, _ = row
        if cancer not in set(proj) or (not resp and not res) or (resp and str(resp[0])[:4] in ("SIG:", "GENE")):
            continue
        if resp == "TOP10":
            g = top[(top["cancer"] == cancer) & (top["drug"] == drug)]
            resp, res = list(g[g.direction == "response"].pathway), list(g[g.direction == "resistance"].pathway)
        inC = list(proj.index[proj == cancer])
        score = F.make_score(Xc.loc[inC], resp, res)
        if score is None:
            continue
        pu = pur.reindex(inC)
        ok = pu.notna()
        rho = stats.spearmanr(score[ok], pu[ok])[0] if ok.sum() > 20 else np.nan
        gl = long[(long["project"] == cancer) & (long["drug"] == drug)]
        treated = sorted(set(gl["patient"]) & set(inC))
        if rule == "recist":
            disc = set(gl[gl["recist"].notna()]["patient"])
        elif rule == "firstcourse":
            first = gl[gl["setting_proxy"].eq("early")]["patient"].unique()
            oc = cdr.reindex(first)["treatment_outcome_first_course"]
            disc = set(oc[oc.isin(F.CR_LIKE | F.NR_LIKE)].index)
        else:
            disc = set()
        held = [p for p in treated if p not in disc and pd.notna(pu.get(p))]
        C = F.tcga_covars(cdr, inC)
        d = {"cancer": cancer, "drug": drug, "rho_score_purity": rho, "n_heldout_with_purity": len(held)}
        for ep in ("OS", "PFI"):
            c = cdr.reindex(inC)
            E = pd.DataFrame({"time": c[f"{ep}.time"].astype(float) / 30.44, "event": c[ep].astype(float)}, index=inC)
            base = pd.concat([E, C], axis=1).assign(score=score, purity=(pu - pu.mean()) / pu.std())
            a = F.cox(base.loc[held], ["score", "age", "stage_II", "stage_III", "stage_IV"])
            b = F.cox(base.loc[held], ["score", "purity", "age", "stage_II", "stage_III", "stage_IV"])
            d.update({f"{ep}_events": a["events"], f"{ep}_HR": a["HR"], f"{ep}_p": a["p"], f"{ep}_HR_adj_purity": b["HR"], f"{ep}_p_adj_purity": b["p"]})
        rows.append(d)
    PB = pd.DataFrame(rows)
    PB.to_csv(OUT / "purity_biomarker_rows.tsv", sep="\t", index=False)
    summ["final_rows_abs_rho_score_purity_median"] = float(PB["rho_score_purity"].abs().median())
    summ["final_rows_abs_rho_gt_0.3"] = int((PB["rho_score_purity"].abs() > 0.3).sum())
    summ["final_rows_scored"] = int(len(PB))
    lhr = np.log(PB[["OS_HR", "OS_HR_adj_purity"]].dropna())
    summ["final_rows_median_abs_change_logHR_OS_purity_adj"] = float((lhr.iloc[:, 0] - lhr.iloc[:, 1]).abs().median()) if len(lhr) else np.nan

    # ---------------- 4. LGG TMZ: prognostic confounding by grade / histology
    g = top[(top["cancer"] == "LGG") & (top["drug"] == "TEMOZOLOMIDE")]
    inC = list(proj.index[proj == "LGG"])
    s = F.make_score(Xc.loc[inC], list(g[g.direction == "response"].pathway), list(g[g.direction == "resistance"].pathway))
    c = cdr.reindex(inC)
    gl = long[(long["project"] == "LGG") & (long["drug"] == "TEMOZOLOMIDE")]
    treated = set(gl["patient"]) & set(inC)
    disc = set(gl[gl["recist"].notna()]["patient"])
    grp = pd.Series("untreated", index=inC)
    grp[[p for p in inC if p in long["patient"].values]] = "other systemic"
    grp[[p for p in inC if p in treated]] = "TMZ labelled (discovery)"
    grp[[p for p in inC if p in treated and p not in disc]] = "TMZ held-out"
    G3 = (c["grade"] == "G3").astype(float)
    lg = []
    for name in ["TMZ labelled (discovery)", "TMZ held-out", "untreated", "other systemic"]:
        m = grp == name
        lg.append({"group": name, "n": int(m.sum()), "frac_G3": float(G3[m].mean()), "OS_events": int(c.loc[m, "OS"].sum()),
                   "median_score": float(s[m].median()), "frac_astrocytoma": float(c.loc[m, "histology"].astype(str).str.contains("Astro").mean())})
    rho_grade = stats.spearmanr(s, G3)[0]
    auc_grade = roc_auc_score(G3, s)
    hist = pd.get_dummies(c["histology"].astype(str), prefix="h", drop_first=True).astype(float)
    base = pd.DataFrame({"time": c["OS.time"].astype(float) / 30.44, "event": c["OS"].astype(float), "score": s, "age": c["age"].astype(float),
                         "G3": G3}, index=inC).join(hist)
    held = [p for p in inC if grp[p] == "TMZ held-out"]
    hcols = list(hist.columns)
    r0 = F.cox(base.loc[held], ["score", "age"])
    r1 = F.cox(base.loc[held], ["score", "age", "G3"])
    r2 = F.cox(base.loc[held], ["score", "age", "G3"] + hcols)
    unt = [p for p in inC if grp[p] == "untreated"]
    u0 = F.cox(base.loc[unt], ["score", "age"])
    allp = F.cox(base, ["score", "age"])
    allg = F.cox(base, ["score", "age", "G3"] + hcols)
    LG = pd.DataFrame(lg)
    LG.to_csv(OUT / "lgg_tmz_confounding.tsv", sep="\t", index=False)
    summ["lgg"] = {"rho_score_G3": float(rho_grade), "auc_score_for_G3": float(auc_grade),
                   "heldout_HR_age": r0, "heldout_HR_age_grade": r1, "heldout_HR_age_grade_hist": r2,
                   "untreated_HR_age": u0, "all_LGG_HR_age": allp, "all_LGG_HR_age_grade_hist": allg,
                   "groups": LG.round(3).to_dict(orient="records")}

    # ---------------- 5. multiple testing over the FINAL table
    pcols = ["OS_heldout_p", "OS_heldout_logrank_p", "OS_interaction_p", "PFI_heldout_p", "PFI_heldout_logrank_p",
             "PFI_interaction_p", "OS_alltreated_p", "PFI_alltreated_p", "OS_untreated_p", "PFI_untreated_p", "ext_p",
             "ext_interaction_p", "ext_best_interaction_p"]
    L = W.melt(id_vars=["cancer", "drug", "verdict"], value_vars=[c for c in pcols if c in W], var_name="test", value_name="p").dropna()
    L = L.drop_duplicates(subset=["cancer", "drug", "test", "p"])
    L["q_BH_all"] = bh(L["p"])
    L = L.sort_values("p")
    L.to_csv(OUT / "multiple_testing.tsv", sep="\t", index=False)
    summ["n_pvalues_final_table"] = int(len(L))
    summ["n_nominal_p_lt_0.05"] = int((L["p"] < 0.05).sum())
    summ["expected_false_pos_if_all_null"] = float(0.05 * len(L))
    summ["n_q_lt_0.1"] = int((L["q_BH_all"] < 0.1).sum())
    summ["tests_q_lt_0.1"] = L[L["q_BH_all"] < 0.1][["cancer", "drug", "test", "p", "q_BH_all"]].round(5).to_dict(orient="records")
    pred = L[L["test"].isin(["OS_heldout_p", "PFI_heldout_p", "OS_interaction_p", "PFI_interaction_p", "ext_p", "ext_interaction_p", "ext_best_interaction_p"])].copy()
    pred["q_BH_primary"] = bh(pred["p"])
    summ["primary_tests_n"] = int(len(pred))
    summ["primary_q_lt_0.1"] = pred[pred["q_BH_primary"] < 0.1][["cancer", "drug", "test", "p", "q_BH_primary"]].round(5).to_dict(orient="records")
    summ["COAD_oxaliplatin_external_q"] = pred[(pred.cancer == "COAD") & (pred.drug == "OXALIPLATIN") & pred.test.str.startswith("ext")][["test", "p", "q_BH_primary"]].round(4).to_dict(orient="records")
    # rows evaluated only on all-treated (includes discovery patients)
    summ["rows_alltreated_only"] = int(W["OS_alltreated_HR"].notna().sum() | 0)
    summ["rows_alltreated_only_list"] = W.loc[W["OS_alltreated_HR"].notna() | W["PFI_alltreated_HR"].notna(), ["cancer", "drug"]].apply(" ".join, axis=1).tolist()

    # ---------------- 6. discovery instability of the top-10 lists
    t1 = pd.read_csv(RES / "atlas/cancer_specific/top10_by_cancer_drug.tsv", sep="\t")
    t2 = pd.read_csv(RES / "atlas/cancer_specific/relaxed_top10.tsv", sep="\t")
    summ["top10_standard"] = {"rows": int(len(t1)), "median_stability": float(t1["top10_stability"].median()),
                              "frac_stability_ge_0.5": float((t1["top10_stability"] >= 0.5).mean()),
                              "frac_local_q_lt_0.1": float((t1["local_q"] < 0.1).mean()),
                              "frac_abs_atlas_z_gt_abs_local_z": float((t1["atlas_z"].abs() > t1["local_z"].abs()).mean())}
    summ["top10_relaxed"] = {"rows": int(len(t2)), "median_stability": float(t2["top10_stability"].median()),
                             "frac_stability_ge_0.5": float((t2["top10_stability"] >= 0.5).mean()),
                             "frac_local_q_lt_0.1": float((t2["local_q"] < 0.1).mean())}

    # ---------------- 7. LOCO AUROC: how much of the spread is sampling noise?
    lp = pd.read_csv(RES / "atlas/loco_pairs.tsv", sep="\t")
    labs = long[long["recist"].notna()].sort_values("start_days").drop_duplicates(["patient", "drug"])
    cnt = labs.groupby(["project", "drug"])["responder"].agg(["size", "sum"]).rename(columns={"size": "nl", "sum": "n1"})
    lp = lp.join(cnt, on=["project", "drug"])
    lp["n1"] = lp["n1"].fillna(lp["n"] * 0.6)
    lp["n0"] = (lp["n"] - lp["n1"]).clip(lower=1)
    lp["n1"] = lp["n1"].clip(lower=1)
    lp["var_null"] = (lp["n1"] + lp["n0"] + 1) / (12 * lp["n1"] * lp["n0"])
    rows = []
    for mdl, g in lp.groupby("model"):
        w = 1 / g["var_null"]
        mu = float(np.average(g["auc"], weights=w))
        Q = float((w * (g["auc"] - mu) ** 2).sum())
        k = len(g)
        tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
        z = (mu - 0.5) / np.sqrt(1 / w.sum())
        rows.append({"model": mdl, "pairs": k, "mean_auc": float(g["auc"].mean()), "precision_weighted_auc": mu,
                     "z_vs_0.5_ignoring_dependence": float(z), "obs_sd": float(g["auc"].std()), "mean_null_sd": float(np.sqrt(g["var_null"]).mean()),
                     "Q": Q, "df": k - 1, "p_heterogeneity": float(1 - stats.chi2.cdf(Q, k - 1)), "tau_between_pairs": float(np.sqrt(tau2)),
                     "frac_pairs_auc_gt_0.5": float((g["auc"] > 0.5).mean()), "median_n": float(g["n"].median())})
    LN = pd.DataFrame(rows)
    LN.to_csv(OUT / "loco_noise.tsv", sep="\t", index=False)
    summ["loco_noise"] = LN.round(4).to_dict(orient="records")

    # ---------------- 8. external cohorts used as evidence: does PC1 of pathway scores track response?
    ext = []
    cache = RES / "cache"
    cat = pd.read_csv(REPO / "data/curated/trials/CATALOG_GEO.tsv", sep="\t")
    for cid in ["GSE104580", "GSE72970", "GSE109211", "EMTAB3267", "GSE63885", "GSE41998", "GSE20271", "GSE32646"]:
        f = cache / f"ATLASEXT_geo_{cid}_reactome_rank.tsv.gz"
        cf = REPO / f"data/curated/trials/{cid}/clinical.tsv"
        if not f.exists() or not cf.exists():
            continue
        S = pd.read_csv(f, sep="\t", index_col=0)
        cl = pd.read_csv(cf, sep="\t").set_index("sample")
        if "responder" not in cl:
            continue
        sm = [x for x in S.columns if x in cl.index and pd.notna(cl.loc[x, "responder"])]
        if len(sm) < 20:
            continue
        X = S[sm].T
        X = (X - X.mean()) / X.std().replace(0, 1)
        X = X.fillna(0)
        u, sv, vt = np.linalg.svd(X.values, full_matrices=False)
        y = cl.loc[sm, "responder"].astype(int).values
        aucs = [roc_auc_score(y, u[:, k]) for k in range(3)]
        aucs = [max(a, 1 - a) for a in aucs]
        per = np.array([max(roc_auc_score(y, X.iloc[:, j]), 1 - roc_auc_score(y, X.iloc[:, j])) for j in range(X.shape[1])])
        ext.append({"cohort": cid, "n": len(sm), "responders": int(y.sum()), "var_PC1": float(sv[0] ** 2 / (sv ** 2).sum()),
                    "AUC_PC1": aucs[0], "AUC_PC2": aucs[1], "AUC_PC3": aucs[2],
                    "frac_pathways_AUC_gt_0.65": float((per > 0.65).mean()), "max_pathway_AUC": float(per.max()),
                    "catalog_qc_flag": bool(cat.set_index("cohort_id").loc[cid, "qc_flag"]) if cid in set(cat["cohort_id"]) else None})
    EX = pd.DataFrame(ext)
    EX.to_csv(OUT / "external_dominant_axis.tsv", sep="\t", index=False)
    summ["external_dominant_axis"] = EX.round(3).to_dict(orient="records")

    json.dump(summ, open(OUT / "confounding_summary.json", "w"), indent=1, default=lambda o: float(o) if np.isscalar(o) else str(o))
    print(json.dumps(summ, indent=1, default=lambda o: float(o) if np.isscalar(o) else str(o)))
    print(PB.round(3).to_string())
    print(out.head(10).round(3).to_string())
    print(out.tail(10).round(3).to_string())


if __name__ == "__main__":
    main()
