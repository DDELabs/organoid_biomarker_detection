"""Cell-line tasks (GDSC + CTRPv2) for the B3 benchmark.

1. Pair table: (cell line, drug) for every drug that has a network-proximity fingerprint (paper's 44-column z table +
   8 drugs recomputed on the paper network by the organoid-only agent; carboplatin/oxaliplatin use cisplatin's).
   GDSC: ln IC50 (GDSC2 if >= 300 lines, else GDSC1). CTRPv2: AUC. Higher = more resistant in both.
   Features: Reactome rank scores already computed by the DepMap agent (nar/_cache/depmap/derived/reactome_*.pkl),
   z-scored inside each screen; response z-scored inside each (screen, drug).
2. Cell-blind CV (5 folds grouped by cell line, GDSC and CTRP copies of a line kept together) and drug-blind CV
   (5 folds grouped by fingerprint, so drugs with identical targets are held out together). Same folds for every model.
   Metric: Spearman (predicted vs observed) across held-out cells, per (screen, drug), averaged.
3. Learning curve: fraction of training cell lines (cell-blind fold split).
4. Final fingerprint MLP pretrained on all pairs -> _data/fpnet.pt (used by bench.py for zero-shot and fine-tuning).
Usage: python cells.py cv | lc | gcn | pretrain
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import ElasticNet, ElasticNetCV, RidgeCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

import b3common as B
import nets

OUT = B.DATA / "cells"
OUT.mkdir(exist_ok=True)


def build():
    f = OUT / "pairs.pkl"
    if f.exists():
        return pickle.load(open(f, "rb"))
    F = B.fingerprints()
    sc_g, sc_c = pd.read_pickle(B.DEPMAP / "reactome_cmp.pkl"), pd.read_pickle(B.DEPMAP / "reactome_ctrp.pkl")
    P = sorted(set(sc_g.index) & set(sc_c.index))
    g = pd.read_pickle(B.DEPMAP / "gdsc_long.pkl")
    g = g[g["ln_ic50"].notna() & g["sanger_model_id"].notna()]
    g["drug"] = g["drug"].str.upper()
    c = pd.read_pickle(B.DEPMAP / "ctrp_long.pkl")
    c["drug"] = c["drug"].str.upper()
    samp = pd.read_pickle(B.DEPMAP / "ctrp_samples.pkl").set_index("sample")
    rows = []
    for d in F.columns:
        gd = g[g.drug == d]
        n2 = gd[gd.dataset == "GDSC2"].sanger_model_id.nunique()
        ds = "GDSC2" if n2 >= 300 else "GDSC1"
        gd = gd[gd.dataset == ds].groupby("sanger_model_id")["ln_ic50"].median()
        gd = gd[gd.index.isin(sc_g.columns)]
        if len(gd) >= 50:
            rows.append(pd.DataFrame({"screen": "GDSC", "dataset": ds, "cell": gd.index, "group": gd.index, "drug": d,
                                      "y_raw": gd.values}))
        cd = c[c.drug == d].groupby("sample")["response"].median()
        cd = cd[cd.index.isin(sc_c.columns)]
        if len(cd) >= 50:
            grp = [samp["Sanger_Model_ID"].get(s) if pd.notna(samp["Sanger_Model_ID"].get(s)) else s for s in cd.index]
            rows.append(pd.DataFrame({"screen": "CTRP", "dataset": "CTRPv2", "cell": cd.index, "group": grp, "drug": d,
                                      "y_raw": cd.values}))
    T = pd.concat(rows, ignore_index=True)
    T["y"] = T.groupby(["screen", "drug"])["y_raw"].transform(lambda v: (v - v.mean()) / v.std())
    Xg = B.zcols(sc_g.loc[P].T.values)[0]
    Xc = B.zcols(sc_c.loc[P].T.values)[0]
    feat = {"GDSC": pd.DataFrame(Xg, index=sc_g.columns, columns=P), "CTRP": pd.DataFrame(Xc, index=sc_c.columns, columns=P)}
    fp = fp_matrix(P)
    grp = B.fp_group(F)
    T["fpgroup"] = T.drug.map(grp)
    D = dict(T=T, feat=feat, fp=fp, P=P)
    pickle.dump(D, open(f, "wb"), protocol=4)
    return D


def fp_matrix(P):
    F = B.fingerprints().reindex(P)
    return (-F.clip(-4, 4) / 2).fillna(0.0)          # positive = close to the drug targets


def design(D, idx):
    T = D["T"].loc[idx]
    Xc = np.vstack([D["feat"][s].loc[c].values for s, c in zip(T.screen, T.cell)])
    Fd = D["fp"][T.drug].T.values
    return Xc.astype(np.float32), Fd.astype(np.float32), T.y.values.astype(np.float32)


def folds(D, kind, k=5):
    T = D["T"]
    groups = T.group.values if kind == "cell" else T.fpgroup.values
    return list(GroupKFold(n_splits=k).split(T, groups=groups))


def per_unit(T, pred):
    out = []
    for (s, d), sub in T.assign(pred=pred).groupby(["screen", "drug"]):
        if len(sub) >= 10 and sub.pred.std() > 0:
            r = B.spearman_ci(sub.pred, sub.y)
            out.append(dict(screen=s, drug=d, n=len(sub), rho=r["rho"], p=r["p"]))
    return out


POOLED = ["Fingerprint MLP (cell lines)", "XGBoost (pathways + fingerprint)", "Elastic net (pathways + fingerprint)",
          "Random forest (pathways + fingerprint)"]
PERDRUG = ["Paper: single network pathway", "Paper: top-2 network pathways", "Ridge: all network pathways",
           "Elastic net (all pathways)"]


def fit_pooled(name, Xc, Fd, y, seed=0):
    X = np.hstack([Xc, Fd])
    if name.startswith("Fingerprint"):
        net = nets.fp_train(Xc, Fd, y, epochs=10, seed=seed)
        return lambda Xc2, Fd2: nets.fp_predict(net, Xc2, Fd2)
    if name.startswith("XGBoost"):
        m = XGBRegressor(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
                         n_jobs=2, random_state=seed, verbosity=0).fit(X, y)
    elif name.startswith("Elastic"):
        m = ElasticNet(alpha=1e-3, l1_ratio=0.5, max_iter=3000).fit(X, y)
    else:
        m = RandomForestRegressor(n_estimators=100, max_features=0.1, min_samples_leaf=20, max_samples=0.3, n_jobs=2,
                                  random_state=seed).fit(X, y)
    return lambda Xc2, Fd2: m.predict(np.hstack([Xc2, Fd2]))


def perdrug_predict(name, Xtr, ytr, Xte, F):
    """Paper models trained inside one (screen, drug); F = proximal pathway column indices."""
    if name.startswith("Elastic"):
        m = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], alphas=20, cv=3, max_iter=3000, n_jobs=1).fit(Xtr, ytr)
        return m.predict(Xte)
    if len(F) == 0:
        return np.full(len(Xte), np.nan)
    m = RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(Xtr[:, F], ytr)
    c = m.coef_
    if name.startswith("Ridge"):
        return Xte[:, F] @ c
    k = 1 if "single" in name else 2
    o = np.argsort(-np.abs(c))[:k]
    return Xte[:, F][:, o] @ c[o]


def run_cv():
    D = build()
    T = D["T"]
    print("pairs", len(T), "drugs", T.drug.nunique(), "units", T.groupby(["screen", "drug"]).ngroups, flush=True)
    Z = B.fingerprints().reindex(D["P"])
    res = []
    for kind in ("cell", "drug"):
        for fi, (tr, te) in enumerate(folds(D, kind)):
            Xc, Fd, y = design(D, T.index[tr])
            Xc2, Fd2, _ = design(D, T.index[te])
            for name in POOLED:
                t0 = time.time()
                pr = fit_pooled(name, Xc, Fd, y, seed=fi)(Xc2, Fd2)
                for u in per_unit(T.iloc[te], pr):
                    res.append(dict(task=f"{kind}-blind", fold=fi, model=name, **u))
                print(kind, fi, name, round(time.time() - t0), "s", flush=True)
            if kind == "cell":
                Ttr, Tte = T.iloc[tr], T.iloc[te]
                Xall = {s: D["feat"][s] for s in D["feat"]}
                for (s, d), sub in Tte.groupby(["screen", "drug"]):
                    trs = Ttr[(Ttr.screen == s) & (Ttr.drug == d)]
                    if len(sub) < 10 or len(trs) < 20:
                        continue
                    Xtr, Xte = Xall[s].loc[trs.cell].values, Xall[s].loc[sub.cell].values
                    zz = Z[d].values
                    Fi = np.where(np.nan_to_num(zz, nan=0) <= B.CUTOFF)[0]
                    for name in PERDRUG:
                        pr = perdrug_predict(name, Xtr, trs.y.values, Xte, Fi)
                        if np.all(np.isfinite(pr)) and np.std(pr) > 0:
                            r = B.spearman_ci(pr, sub.y)
                            res.append(dict(task="cell-blind", fold=fi, model=name, screen=s, drug=d, n=len(sub),
                                            rho=r["rho"], p=r["p"]))
                print(kind, fi, "per-drug models done", flush=True)
            pd.DataFrame(res).to_csv(B.TABLES / "cellline_cv_per_unit.tsv", sep="\t", index=False)
    summarise()


def gcn_cells():
    """GCN per drug, cell-blind folds, GDSC lines only, for the drugs used in the organoid/PDX sets."""
    D = build()
    T = D["T"]
    E = pd.read_pickle(B.DEPMAP / "expr_cmp.pkl")
    drugs = sorted({d for _, d in B.combos()} & set(T.drug))
    fl = folds(D, "cell")
    res = []
    for d in drugs:
        sub = T[(T.screen == "GDSC") & (T.drug == d)]
        cells = [c for c in sub.cell if c in E.columns]
        sub = sub[sub.cell.isin(cells)]
        Ed = E[sub.cell.values]
        var_order = list(Ed.var(axis=1).sort_values(ascending=False).index)
        keep, A, prior = nets.gcn_graph(d, list(Ed.index), var_order)
        X = B.zcols(Ed.loc[keep].T.values)[0]
        pos = {i: k for k, i in enumerate(sub.index)}
        for fi, (tr, te) in enumerate(fl):
            itr = [pos[i] for i in T.index[tr] if i in pos]
            ite = [pos[i] for i in T.index[te] if i in pos]
            if len(ite) < 10:
                continue
            m = nets.gcn_fit(A, prior, X[itr], sub.y.values[itr], seed=fi)
            r = B.spearman_ci(nets.gcn_predict(m, X[ite]), sub.y.values[ite])
            res.append(dict(task="cell-blind", fold=fi, model="GCN on STRING (per drug)", screen="GDSC", drug=d,
                            n=len(ite), rho=r["rho"], p=r["p"]))
        print("gcn", d, len(sub), len(keep), flush=True)
        pd.DataFrame(res).to_csv(B.TABLES / "cellline_cv_per_unit_gcn.tsv", sep="\t", index=False)
    summarise()


def summarise():
    parts = [pd.read_csv(f, sep="\t") for f in (B.TABLES / "cellline_cv_per_unit.tsv", B.TABLES / "cellline_cv_per_unit_gcn.tsv")
             if f.exists()]
    R = pd.concat(parts)
    rows = []
    for (task, model), sub in R.groupby(["task", "model"]):
        u = sub.groupby(["screen", "drug"])["rho"].mean()            # average over folds per unit
        m, s, n = u.mean(), u.std(ddof=1), len(u)
        rows.append(dict(task=task, model=model, n_units=n, mean_rho=m, lo=m - 1.96 * s / np.sqrt(n),
                         hi=m + 1.96 * s / np.sqrt(n), share_units_rho_pos=(u > 0).mean()))
    pd.DataFrame(rows).to_csv(B.TABLES / "cellline_cv_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())


def learning_curve():
    D = build()
    T = D["T"]
    tr, te = folds(D, "cell")[0]
    rng = np.random.default_rng(B.SEED)
    cells = np.unique(T.group.values[tr])
    Xc2, Fd2, _ = design(D, T.index[te])
    res = []
    for frac in (0.05, 0.1, 0.25, 0.5, 1.0):
        for rep in range(3 if frac < 1 else 1):
            keep = set(rng.choice(cells, max(5, int(frac * len(cells))), replace=False))
            sel = [i for i in tr if T.group.values[i] in keep]
            Xc, Fd, y = design(D, T.index[sel])
            for name in POOLED[:3]:
                pr = fit_pooled(name, Xc, Fd, y, seed=rep)(Xc2, Fd2)
                u = pd.DataFrame(per_unit(T.iloc[te], pr))
                res.append(dict(frac=frac, rep=rep, n_cells=len(keep), model=name, mean_rho=u.rho.mean()))
                print(frac, rep, name, round(u.rho.mean(), 3), flush=True)
    pd.DataFrame(res).to_csv(B.TABLES / "learning_curve_cell_lines.tsv", sep="\t", index=False)


def pretrain():
    D = build()
    Xc, Fd, y = design(D, D["T"].index)
    net = nets.fp_train(Xc, Fd, y, epochs=10, seed=0)
    torch.save({"state": net.state_dict(), "P": D["P"]}, B.DATA / "fpnet.pt")
    print("saved fpnet", len(y), flush=True)


if __name__ == "__main__":
    {"cv": run_cv, "lc": learning_curve, "gcn": gcn_cells, "pretrain": pretrain, "summary": summarise}[sys.argv[1]]()
