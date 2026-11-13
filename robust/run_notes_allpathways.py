"""Response-note model per drug: network shortlist vs all 673 pathways (same nested tuning for both).

Out-of-fold within-cancer AUROC on TCGA, and frozen AUROC in the independent trial cohorts.
Usage: python robust/run_notes_allpathways.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_integrated as RI  # noqa: E402
import run_survlearn as RS  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.curated import load_curated_trial, trial_catalog  # noqa: E402

LAMS = (10.0, 30.0, 100.0, 300.0, 1000.0)


def pick_lambda(X, y, groups, seed=1):
    f = SL.folds_by_patient(X.index, 3, seed)
    best, bl = -1, LAMS[-1]
    for lam in LAMS:
        s = pd.Series(np.nan, index=X.index)
        for j in range(3):
            tr, te = X.index[(f != j).values], X.index[(f == j).values]
            c = SL.response_model(X.loc[tr], y.loc[tr], lam)
            if c is not None:
                s.loc[te] = X.loc[te].values @ c.values
        ok = s.notna()
        auc = SL.stratified_auc(s[ok].values, y[ok].values, groups[ok].values)[0] if ok.sum() > 10 else np.nan
        if np.isfinite(auc) and auc > best:
            best, bl = auc, lam
    return bl


def one(drug, Z, cdr, long, fold, trials, paths_all):
    rec = long[long["drug"] == drug]
    lab = rec[rec["responder"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")
    y, grp = lab["responder"].astype(int), cdr.loc[lab.index, "project"]
    if len(y) < 30 or y.value_counts().min() < 8:
        return []
    short, how = SL.shortlist(RS.proximity(drug), paths_all)
    out = []
    for name, paths in (("network shortlist", short), ("all pathways", paths_all)):
        if name == "network shortlist" and how.startswith("no targets"):
            continue
        X = Z.loc[y.index, paths]
        s = pd.Series(np.nan, index=y.index)
        for f in sorted(fold.unique()):
            tr, te = y.index[(fold.loc[y.index] != f).values], y.index[(fold.loc[y.index] == f).values]
            lam = pick_lambda(X.loc[tr], y.loc[tr], grp.loc[tr])
            c = SL.response_model(X.loc[tr], y.loc[tr], lam)
            if c is not None:
                s.loc[te] = X.loc[te].values @ c.values
        ok = s.notna()
        auc, pairs = SL.stratified_auc(s[ok].values, y[ok].values, grp[ok].values)
        rng = np.random.default_rng(5)
        null = []
        for _ in range(500):
            yp = y[ok].values.copy()
            for g in grp[ok].unique():
                i = np.flatnonzero((grp[ok] == g).values)
                yp[i] = yp[rng.permutation(i)]
            null.append(SL.stratified_auc(s[ok].values, yp, grp[ok].values)[0])
        row = {"drug": drug, "pathways": name, "n_pathways": len(paths), "n": int(ok.sum()), "responders": int(y[ok].sum()),
               "tcga_auc": auc, "tcga_p": (np.sum(np.array(null) >= auc) + 1) / 501}
        lam = pick_lambda(X, y, grp)
        c = SL.response_model(X, y, lam)
        for key, want in (("ext", True), ("ext_other", False)):
            S, Y, G = [], [], []
            for cid, (sc, clin) in trials.items():
                got = clin["drugs"].fillna("").astype(str).str.upper().str.split(";").map(lambda x: drug in x)
                if not got.any():
                    continue
                sel = [x for x in clin.index[(got == want) & clin["responder"].notna()] if x in sc.index]
                yy = clin.loc[sel, "responder"].astype(int)
                if len(sel) < 15 or yy.value_counts().min() < 4:
                    continue
                S.append(sc.loc[sel, paths].values @ c.values)
                Y.append(yy.values)
                G.append(np.full(len(sel), cid))
            if Y:
                s2, y2, g2 = np.concatenate(S), np.concatenate(Y), np.concatenate(G)
                row[f"{key}_auc"], row[f"{key}_n"], row[f"{key}_cohorts"] = SL.stratified_auc(s2, y2, g2)[0], len(y2), len(np.unique(g2))
        out.append(row)
    return out


def main():
    Z, cdr, long = RS.load()
    paths = list(Z.columns)
    fold = SL.folds_by_patient(cdr.index, 5, seed=0)
    cat = trial_catalog()
    trials = {}
    for cid in cat.loc[cat["usable"].astype(str).str.lower().eq("true"), "cohort_id"]:
        try:
            _, clin = load_curated_trial(cid, canonical=False)
            if "responder" in clin and "drugs" in clin:
                trials[cid] = (RI.cohort_scores(cid).reindex(columns=paths).fillna(0), clin)
        except Exception:  # noqa: BLE001
            pass
    n = long[long["responder"].notna()].drop_duplicates(["patient", "drug"]).groupby("drug").size()
    drugs = list(n[n >= 30].sort_values(ascending=False).index)
    res = Parallel(n_jobs=4)(delayed(one)(d, Z, cdr, long, fold, trials, paths) for d in drugs)
    R = pd.DataFrame([r for rr in res for r in rr])
    out = RESULTS / "survlearn" / "notes_shortlist_vs_all.tsv"
    R.to_csv(out, sep="\t", index=False)
    W = R.pivot(index="drug", columns="pathways", values=["tcga_auc", "ext_auc"]).reindex(drugs)
    W["n"] = R.drop_duplicates("drug").set_index("drug")["n"]
    print(W.round(2).to_string())
    for k in ("tcga_auc", "ext_auc"):
        both = W[k].dropna()
        print(k, "median  shortlist %.3f  all %.3f  | all better in %d of %d drugs" % (both["network shortlist"].median(), both["all pathways"].median(),
              (both["all pathways"] > both["network shortlist"]).sum(), len(both)))
    a = R[R.pathways == "all pathways"]
    print(a[["drug", "n", "tcga_auc", "tcga_p", "ext_auc", "ext_n", "ext_other_auc", "ext_other_n"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
