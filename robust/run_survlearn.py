"""Patient-trained drug-specific pathway biomarkers: survival + treatment (+ response notes), all drugs.

Usage: python robust/run_survlearn.py [--min-treated 30] [--random 0] [--drugs A B ...]
Writes robust/results/survlearn/{summary.tsv, coefficients.tsv.gz, oof_scores.tsv.gz}
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from run_atlas_relaxed import SUPPORTIVE, THERAPY  # noqa: E402

OUT = RESULTS / "survlearn"
PROXY = {"CARBOPLATIN": "CISPLATIN", "OXALIPLATIN": "CISPLATIN"}   # no protein target of their own: platinum class
BASE = ["age", "stage", "stage_missing", "purity", "purity_missing", "other_systemic", "g"]
LAMS = (30.0, 100.0, 300.0, 1000.0)
STAGE = {"I": 1, "II": 2, "III": 3, "IV": 4}


def stage_num(s):
    if not isinstance(s, str):
        return np.nan
    t = s.upper().replace("STAGE", "").strip()
    for k in ("IV", "III", "II", "I"):
        if t.startswith(k):
            return STAGE[k]
    return np.nan


def load():
    M, proj = RT.patient_matrix(R.PROJECTS)
    cdr = R.load_cdr()
    pats = [p for p in M.index if p in cdr.index]
    cdr = cdr.loc[pats].copy()
    cdr["project"] = proj.loc[pats]
    cdr["time"], cdr["event"] = cdr["PFI.time"], cdr["PFI"]
    cdr = cdr[(cdr["time"] > 0) & cdr["event"].notna()]
    Z = SL.zscore_within(M.loc[cdr.index], cdr["project"])
    by = cdr.groupby("project")
    cdr["age"] = ((cdr["age"] - by["age"].transform("mean")) / 10).fillna(0)
    st = cdr["stage"].map(stage_num)
    cdr["stage_missing"] = st.isna().astype(float)
    cdr["stage"] = (st - st.groupby(cdr["project"]).transform("mean")).fillna(0)
    pur = pd.read_csv(Path(__file__).parent / "revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    pur = pur.assign(patient=pur["array"].str[:12]).drop_duplicates("patient").set_index("patient")["purity"]
    p = pur.reindex(cdr.index)
    cdr["purity_missing"] = p.isna().astype(float)
    cdr["purity"] = (p - p.groupby(cdr["project"]).transform("mean")).fillna(0)
    long = R.load_response_labels()
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    long = long[long["patient"].isin(cdr.index)].copy()
    # first-line exposure = drug started before the progression-free interval ended
    pfi = long["patient"].map(cdr["time"])
    long["before_event"] = (long["start_days"] < pfi) | (long["start_days"].isna() & long["setting_proxy"].eq("early"))
    return Z, cdr, long


PAPER = False    # True: use the authors' own proximity table (data/coad_blca_organoid_drugs_zscore_result_reactome.txt)


def proximity(drug):
    if PAPER:
        from obd import network as N
        old = N.precomputed_proximity()
        return old[drug].dropna() if drug in old.columns else None
    f = RESULTS / "cache" / f"proximity_STRINGv12_{PROXY.get(drug, drug)}.tsv.gz"
    return pd.read_csv(f, sep="\t", index_col=0)["z"] if f.exists() else None


def choose_lambda(D, Zs, k_inner=3, seed=1):
    """Cross-validated partial likelihood for the interaction penalty (inside the training data only)."""
    f = SL.folds_by_patient(D.index, k_inner, seed)
    T = D["treated"].values[:, None]
    X = np.hstack([D[BASE].values, D[["treated"]].values, Zs.values, Zs.values * T])
    k, nb = Zs.shape[1], len(BASE) + 1
    full = SL.StratCox(D["time"], D["event"], D["project"])
    best, best_l = -np.inf, LAMS[-1]
    for lam in LAMS:
        pen = np.r_[np.full(nb, 1e-3), np.full(k, 100.0), np.full(k, lam)]
        cv = 0.0
        for j in range(k_inner):
            tr = (f != j).values
            m = SL.StratCox(D["time"][tr], D["event"][tr], D["project"][tr])
            b = m.fit(X[tr], pen)
            cv += full.loglik(b, X) - m.loglik(b, X[tr])
        if cv > best:
            best, best_l = cv, lam
    return best_l


def one_drug(drug, Z, cdr, long, fold, paths=None, tag="proximal", min_treated=30):
    rec = long[long["drug"] == drug]
    first = set(rec.loc[rec["before_event"], "patient"])
    per = cdr.loc[sorted(first), "project"].value_counts()
    cancers = list(per[per >= 10].index)
    lab = rec[rec["responder"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")
    if paths is None:
        paths, how = SL.shortlist(proximity(drug), Z.columns)
    else:
        how = tag
    row = {"drug": drug, "pathway_source": how, "n_pathways": len(paths), "cancers": ",".join(cancers),
           "n_any_line": rec["patient"].nunique(), "n_labelled": len(lab)}
    D = cdr[cdr["project"].isin(cancers)].copy()
    D["treated"] = D.index.isin(first).astype(float)
    other = set(long.loc[long["before_event"] & (long["drug"] != drug), "patient"])
    D["other_systemic"] = D.index.isin(other).astype(float)
    row.update(n=len(D), n_treated=int(D["treated"].sum()), events_treated=int(D.loc[D.treated == 1, "event"].sum()),
               events_untreated=int(D.loc[D.treated == 0, "event"].sum()))
    surv_ok = row["n_treated"] >= min_treated and row["events_treated"] >= 15 and row["events_untreated"] >= 15
    labp = lab[lab.index.isin(Z.index)]
    y_all = labp["responder"].astype(int)
    notes_ok = len(y_all) >= 30 and y_all.value_counts().min() >= 8 if len(y_all) else False
    scores = pd.DataFrame(index=sorted(set(D.index) | set(labp.index)))
    coefs = {}
    if surv_ok:
        Zs = Z.loc[D.index, paths]
        s = pd.Series(np.nan, index=D.index)
        lams = []
        for f in sorted(fold.unique()):
            tr = D.index[(fold.loc[D.index] != f).values]
            te = D.index[(fold.loc[D.index] == f).values]
            lam = choose_lambda(D.loc[tr], Zs.loc[tr])
            lams.append(lam)
            _, d = SL.fit_interaction(D.loc[tr], Zs.loc[tr], BASE, 100.0, lam)
            s.loc[te] = -(Zs.loc[te].values @ d.values)
        scores["surv"] = s
        lam = choose_lambda(D, Zs)
        a, d = SL.fit_interaction(D, Zs, BASE, 100.0, lam)
        coefs["prognostic"], coefs["drug_specific"] = a, d
        row["lambda_int"] = lam
    if notes_ok:
        Zl = Z.loc[y_all.index, paths]
        allp = scores.index
        s = pd.Series(np.nan, index=allp)
        fl = fold.reindex(allp)
        for f in sorted(fold.unique()):
            tr = y_all.index[(fold.loc[y_all.index] != f).values]
            c = SL.response_model(Zl.loc[tr], y_all.loc[tr], lam=30.0)
            te = allp[(fl == f).values]
            if c is not None and len(te):
                s.loc[te] = Z.loc[te, paths].values @ c.values
        scores["notes"] = s
        c = SL.response_model(Zl, y_all, lam=30.0)
        if c is not None:
            coefs["notes"] = c
    if "surv" in scores and "notes" in scores:
        zz = lambda v: (v - v.mean()) / v.std()  # noqa: E731
        scores["combined"] = (zz(scores["surv"]) + zz(scores["notes"])) / 2
    # ---- honest evaluation on out-of-fold scores
    for name in [c for c in ("surv", "notes", "combined") if c in scores]:
        sc = scores[name]
        if surv_ok:
            E = D.join(sc.rename("score"), how="inner").dropna(subset=["score"])
            E["score"] = (E["score"] - E["score"].mean()) / E["score"].std()
            E["tx_score"] = E["treated"] * E["score"]
            try:
                r = SL.cox_test(E, BASE + ["treated", "score", "tx_score"])
                row[f"{name}_int_HR"], row[f"{name}_int_p"], row[f"{name}_int_z"] = np.exp(r.loc["tx_score", "coef"]), r.loc["tx_score", "p"], r.loc["tx_score", "z"]
                row[f"{name}_untreated_HR"] = np.exp(r.loc["score", "coef"])
                rt = SL.cox_test(E[E.treated == 1], BASE[:-2] + ["g", "score"])
                row[f"{name}_treated_HR"], row[f"{name}_treated_p"] = np.exp(rt.loc["score", "coef"]), rt.loc["score", "p"]
            except Exception as exc:  # noqa: BLE001
                row[f"{name}_error"] = str(exc)[:80]
        if len(y_all) >= 20:
            ok = sc.reindex(y_all.index).notna()
            if ok.sum() >= 20:
                yy = y_all[ok]
                auc, pairs = SL.stratified_auc(sc.reindex(yy.index).values, yy.values, labp.loc[yy.index, "project"].values)
                row[f"{name}_resp_AUC"], row[f"{name}_resp_pairs"] = auc, pairs
    row["surv_model"], row["notes_model"] = bool(surv_ok), bool(notes_ok)
    return row, coefs, scores.assign(drug=drug)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-treated", type=int, default=30)
    ap.add_argument("--random", type=int, default=0, help="random pathway shortlists per drug (negative control)")
    ap.add_argument("--drugs", nargs="*")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--paper-proximity", action="store_true", help="use the paper's own shortlists; only drugs in its table")
    a = ap.parse_args()
    global PAPER, OUT
    if a.paper_proximity:
        PAPER, OUT = True, RESULTS / "survlearn_paperprox"
    OUT.mkdir(parents=True, exist_ok=True)
    Z, cdr, long = load()
    fold = SL.folds_by_patient(cdr.index, 5, seed=0)
    cdr["g"] = SL.general_risk(Z, cdr, "time", "event", fold)
    cdr["g"] = cdr["g"] / cdr["g"].std()
    gt = SL.cox_test(cdr, ["age", "stage", "stage_missing", "g"])
    print(f"patients {len(cdr)}, PFI events {int(cdr.event.sum())}; general risk score (out-of-fold) HR per SD = "
          f"{np.exp(gt.loc['g', 'coef']):.2f}, p = {gt.loc['g', 'p']:.1e}", flush=True)
    n = long.drop_duplicates(["patient", "drug"]).groupby("drug").size().sort_values(ascending=False)
    drugs = a.drugs or list(n[n >= a.min_treated].index)
    if PAPER:
        drugs = [d for d in drugs if proximity(d) is not None]
    fixed = {}
    if PAPER:        # the paper's exact rule: z <= -1.2816 in the authors' table (resolved here, not in the workers)
        fixed = {d: [p for p in proximity(d).index[proximity(d) <= -1.2816] if p in Z.columns] for d in drugs}
    res = Parallel(n_jobs=a.jobs)(delayed(one_drug)(d, Z, cdr, long, fold, paths=fixed.get(d), tag="proximal",
                                                    min_treated=a.min_treated) for d in drugs)
    S = pd.DataFrame([r[0] for r in res])
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    co = [pd.DataFrame(c).assign(drug=r["drug"]).rename_axis("pathway").reset_index() for r, c, _ in res if c]
    if co:
        pd.concat(co).to_csv(OUT / "coefficients.tsv.gz", sep="\t", index=False)
    pd.concat([s.rename_axis("patient").reset_index() for _, _, s in res]).to_csv(OUT / "oof_scores.tsv.gz", sep="\t", index=False)
    if a.random:
        rng = np.random.default_rng(7)
        jobs = []
        for r in S.itertuples():
            if r.surv_model and r.pathway_source == "proximal":
                for i in range(a.random):
                    jobs.append((r.drug, list(rng.choice(Z.columns, r.n_pathways, replace=False)), i))
        rr = Parallel(n_jobs=a.jobs)(delayed(one_drug)(d, Z, cdr, long, fold, paths=p, tag=f"random{i}") for d, p, i in jobs)
        pd.DataFrame([x[0] for x in rr]).to_csv(OUT / "random_shortlists.tsv", sep="\t", index=False)
    cols = [c for c in ["drug", "n_pathways", "n_treated", "events_treated", "events_untreated", "n_labelled", "surv_int_HR", "surv_int_p",
                        "surv_resp_AUC", "notes_int_HR", "notes_int_p", "notes_resp_AUC", "combined_int_HR", "combined_int_p", "combined_resp_AUC"] if c in S]
    print(S[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
