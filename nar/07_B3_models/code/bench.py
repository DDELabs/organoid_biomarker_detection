"""B3 organoid / PDX / lab-model benchmark: every model on the same splits for one (set, drug) pair.

For each pair (>= 15 treated TCGA patients):
  * organoid task: out-of-fold prediction of the lab-model drug response (leave-one-out if <= 40 models, else
    10-fold; identical folds for every model), Spearman rho with 95% CI (Fisher z);
  * patient task: model trained on all lab models, applied to TCGA patients given the drug: C-index (bootstrap CI),
    HR per SD of predicted resistance (Cox), log-rank of the median split, treated-vs-untreated interaction (Cox).
Risk = predicted resistance (higher = expected shorter survival); HR > 1 / C > 0.5 = the expected direction.
Lab-model features are z-scored inside each training fold; patient features inside the patient cohort (as the paper).
Usage: python bench.py <worker> <n_workers>
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold
from xgboost import XGBRegressor

import b3common as B
import nets

OUT = B.DATA / "bench"
OUT.mkdir(exist_ok=True)
N_RANDOM = 100
ALPHAS = np.arange(0.1, 1, 0.1)          # paper: RidgeCV alpha 0.1 .. 0.9, 3-fold

MODELS = [
    # name, family
    ("Paper: single network pathway", "paper"),
    ("Paper: top-2 network pathways", "paper"),
    ("Ridge: all network pathways", "paper"),
    ("Whole transcriptome (ridge)", "paper comparison"),
    ("All pathways (ridge)", "paper comparison"),
    ("First-neighbour genes (ridge)", "paper comparison"),
    ("Bolis correlation genes (ridge)", "paper comparison"),
    ("Deep learning, MOLI-style (genes)", "paper comparison"),
    ("Random pathway sets (median of 100)", "paper comparison"),
    ("Degree-central pathways (rank 1)", "paper comparison"),
    ("Betweenness-central pathways (rank 1)", "paper comparison"),
    ("Closeness-central pathways (rank 1)", "paper comparison"),
    ("Elastic net (all pathways)", "modern"),
    ("Random forest (all pathways)", "modern"),
    ("XGBoost (all pathways)", "modern"),
    ("Fingerprint MLP: cell lines only (zero-shot)", "modern"),
    ("Fingerprint MLP: fine-tuned on lab models", "modern"),
    ("GCN on STRING (genes)", "modern"),
    ("Network propagation genes (ridge)", "modern"),
    ("Paper: single network pathway + alignment", "alignment"),
    ("Ridge: all network pathways + alignment", "alignment"),
    ("Elastic net + alignment", "alignment"),
    ("Random forest + alignment", "alignment"),
    ("XGBoost + alignment", "alignment"),
]


# ------------------------------------------------------------------------------------------------ context
class Ctx:
    def __init__(self, bname, drug):
        b = B.bundle(bname)
        d = b["drugs"][drug]
        self.bname, self.drug, self.meta, self.d = bname, drug, b["meta"], d
        y = d["y"]
        S = [s for s in y.index if s in b["osc"].columns and s in b["gm"].columns]
        self.S, self.y = S, y.loc[S].values.astype(float)
        self.U = sorted(set(b["osc"].index) & set(b["tsc"].index))
        pats = [p for p in b["tsc"].columns if p in b["gp"].columns and p in b["surv"].index]
        self.pats = pats
        self.Xo_pw = b["osc"].loc[self.U, S].T.values.astype(float)
        self.Xp_pw = B.zcols(b["tsc"].loc[self.U, pats].T.values)[0]
        self.genes = list(b["gm"].index)
        self.Xo_g = b["gm"][S].T.values.astype(float)
        self.Xp_g = B.zcols(b["gp"][pats].T.values)[0].astype(np.float32)
        tr = set(d["treated"])
        un = set(b["untreated"])
        self.ti = np.array([i for i, p in enumerate(pats) if p in tr])
        self.ui = np.array([i for i, p in enumerate(pats) if p in un and p not in tr])
        sv = b["surv"].loc[pats]
        self.T, self.E = np.clip(sv["T"].values, 0, None), sv["E"].values.astype(int)
        z = d["z"].reindex(self.U)
        self.z = z
        self.F = [i for i, p in enumerate(self.U) if pd.notna(z[p]) and z[p] <= B.CUTOFF]
        self.var_order = list(pd.Series(self.Xo_g.var(0), index=self.genes).sort_values(ascending=False).index)
        self.gidx = {g: i for i, g in enumerate(self.genes)}
        # alignment (nar/03): organoid pathway matrix aligned to the patient matrix (unsupervised, no response used)
        Xm = B.zcols(self.Xo_pw)[0]
        self.Xo_pwA, self.n_drop = B.align_models(Xm, self.Xp_pw)
        # rank-score features for the fingerprint MLP (same scoring as the cell lines)
        ck = torch.load(B.DATA / "fpnet.pt", weights_only=False)
        self.P = ck["P"]
        rs_o = B.rank_scores(b["gm"][S]).reindex(self.P)
        rs_p = B.rank_scores(b["gp"][pats]).reindex(self.P)
        self.Ro = np.nan_to_num(B.zcols(rs_o.T.values)[0])
        self.Rp = np.nan_to_num(B.zcols(rs_p.T.values)[0])
        net = nets.FPNet(len(self.P))
        net.load_state_dict(ck["state"])
        net.eval()
        self.fpnet = net
        from cells import fp_matrix
        fpm = fp_matrix(self.P)
        self.fp = fpm[drug].values.astype(np.float32) if drug in fpm.columns else np.zeros(len(self.P), np.float32)
        # centrality-selected pathway sets (Supp Fig 6/7): top |F| pathways by mean gene centrality
        cen = pd.read_csv(B.DATA / "centrality_pathways.tsv", sep="\t", index_col=0).reindex(self.U)
        k = max(len(self.F), 2)
        self.central = {c: [self.U.index(p) for p in cen[c].dropna().sort_values(ascending=False).index[:k]]
                        for c in ("degree", "betweenness", "closeness")}
        n = len(S)
        self.folds = ([(np.setdiff1d(np.arange(n), [i]), np.array([i])) for i in range(n)] if n <= 40 else
                      list(KFold(10, shuffle=True, random_state=0).split(np.arange(n))))
        self.cv = "leave-one-out" if n <= 40 else "10-fold"


# ------------------------------------------------------------------------------------------------ helpers
def ridge_paper(X, y):
    return RidgeCV(cv=min(3, len(y)), alphas=ALPHAS).fit(X, y).coef_


def zy(y):
    return (y - y.mean()) / (y.std() or 1)


def rank_model(Xo, Xp, y, tr, te, cols, k, full):
    Xtr, mu, sd = B.zcols(Xo[tr][:, cols])
    c = ridge_paper(Xtr, zy(y[tr]))
    o = np.argsort(-np.abs(c), kind="stable")[: (len(c) if k is None else k)]
    po = B.zcols(Xo[te][:, cols], mu, sd)[0][:, o] @ c[o]
    pp = Xp[:, cols][:, o] @ c[o] if full else None
    return po, pp, [cols[i] for i in o], c[o]


def bolis(Xtr, ytr, rng, n_rep=10, p_cut=0.05):
    """Bolis et al. 2018-style selection: genes with Spearman p < 0.05 vs response in every one of 10 random halves."""
    from scipy.stats import t as tdist
    n = len(ytr)
    R = np.apply_along_axis(lambda v: pd.Series(v).rank().values, 0, Xtr)
    count = np.zeros(Xtr.shape[1])
    for _ in range(n_rep):
        h = rng.choice(n, max(4, n // 2), replace=False)
        a = R[h] - R[h].mean(0)
        b = pd.Series(ytr[h]).rank().values
        b = b - b.mean()
        r = (a * b[:, None]).sum(0) / np.sqrt((a ** 2).sum(0) * (b ** 2).sum(0) + 1e-12)
        m = len(h)
        tt = r * np.sqrt((m - 2) / np.maximum(1 - r ** 2, 1e-12))
        p = 2 * tdist.sf(np.abs(tt), m - 2)
        count += p < p_cut
    sel = np.where(count == n_rep)[0]
    if len(sel) < 5:                                    # fallback: the 10 most often selected genes
        sel = np.argsort(-count, kind="stable")[:10]
    return sel


# ------------------------------------------------------------------------------------------------ one model, one fold
def fit_predict(name, C, tr, te, full=False, extra=None):
    """Returns (pred on lab models te, pred on all patients or None, info)."""
    info = {}
    if name.startswith("Paper: single") or name.startswith("Paper: top-2") or name.startswith("Ridge: all network"):
        Xo = C.Xo_pwA if "alignment" in name else C.Xo_pw
        k = 1 if "single" in name else (2 if "top-2" in name else None)
        po, pp, cols, c = rank_model(Xo, C.Xp_pw, C.y, tr, te, C.F, k, full)
        info = {"pathways": [C.U[i] for i in cols], "coef": list(c)}
        return po, pp, info
    if name.startswith("All pathways"):
        po, pp, _, _ = rank_model(C.Xo_pw, C.Xp_pw, C.y, tr, te, list(range(len(C.U))), None, full)
        return po, pp, info
    if "central" in name:
        key = name.split("-")[0].lower()
        po, pp, cols, c = rank_model(C.Xo_pw, C.Xp_pw, C.y, tr, te, C.central[key], 1, full)
        return po, pp, {"pathways": [C.U[i] for i in cols]}
    if name.startswith("Whole transcriptome") or name.startswith("First-neighbour") or name.startswith("Network propagation") \
            or name.startswith("Bolis"):
        if name.startswith("Whole"):
            cols = np.arange(len(C.genes))
        elif name.startswith("First"):
            g = B.paper_network()
            t = [x for x in B.targets(C.drug) if x in g]
            nb = set(t) | {n for x in t for n in g.neighbors(x)}
            cols = np.array([C.gidx[x] for x in nb if x in C.gidx])
        elif name.startswith("Network"):
            pr = nets.rwr(C.drug)
            cols = np.array([C.gidx[g] for g, _ in sorted(pr.items(), key=lambda kv: -kv[1]) if g in C.gidx][:300])
        if len(np.atleast_1d(cols)) == 0 if not name.startswith("Bolis") else False:
            return np.full(len(te), np.nan), (np.full(len(C.pats), np.nan) if full else None), {"n_genes": 0}
        Xtr_all = C.Xo_g[tr]
        if name.startswith("Bolis"):
            cols = bolis(Xtr_all, C.y[tr], np.random.default_rng(B.SEED))
        Xtr, mu, sd = B.zcols(Xtr_all[:, cols])
        c = ridge_paper(Xtr, zy(C.y[tr]))
        po = B.zcols(C.Xo_g[te][:, cols], mu, sd)[0] @ c
        pp = C.Xp_g[:, cols] @ c if full else None
        return po, pp, {"n_genes": int(len(cols)), "genes": [C.genes[i] for i in cols][:200]}
    if name.startswith("Elastic") or name.startswith("Random forest") or name.startswith("XGBoost"):
        Xo = C.Xo_pwA if "alignment" in name else C.Xo_pw
        Xtr, mu, sd = B.zcols(Xo[tr])
        ytr = zy(C.y[tr])
        if name.startswith("Elastic"):
            m = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], n_alphas=30, cv=min(3, len(tr)), max_iter=5000,
                             random_state=0).fit(Xtr, ytr)
            info = {"alpha": float(m.alpha_), "l1_ratio": float(m.l1_ratio_), "n_nonzero": int((m.coef_ != 0).sum())}
        elif name.startswith("Random"):
            m = RandomForestRegressor(n_estimators=300, max_features="sqrt", min_samples_leaf=2, random_state=0,
                                      n_jobs=1).fit(Xtr, ytr)
        else:
            m = XGBRegressor(n_estimators=300, max_depth=2, learning_rate=0.05, subsample=0.8, colsample_bytree=0.3,
                             random_state=0, n_jobs=1, verbosity=0).fit(Xtr, ytr)
        po = m.predict(B.zcols(Xo[te], mu, sd)[0])
        pp = m.predict(C.Xp_pw) if full else None
        return po, pp, info
    if name.startswith("Fingerprint"):
        if "zero-shot" in name:
            po = nets.fp_predict(C.fpnet, C.Ro[te], C.fp)
            pp = nets.fp_predict(C.fpnet, C.Rp, C.fp) if full else None
            return po, pp, info
        net = nets.fp_finetune(C.fpnet, C.Ro[tr], C.fp, zy(C.y[tr]).astype(np.float32), seed=0)
        po = nets.fp_predict(net, C.Ro[te], C.fp)
        pp = nets.fp_predict(net, C.Rp, C.fp) if full else None
        return po, pp, info
    if name.startswith("GCN"):
        if not hasattr(C, "gcn_graph"):
            C.gcn_graph = nets.gcn_graph(C.drug, C.genes, C.var_order)
        keep, A, prior = C.gcn_graph
        cols = [C.gidx[g] for g in keep]
        Xtr, mu, sd = B.zcols(C.Xo_g[tr][:, cols])
        m = nets.gcn_fit(A, prior, Xtr, zy(C.y[tr]), seed=0)
        po = nets.gcn_predict(m, B.zcols(C.Xo_g[te][:, cols], mu, sd)[0])
        pp = nets.gcn_predict(m, C.Xp_g[:, cols]) if full else None
        return po, pp, {"n_nodes": len(keep)}
    if name.startswith("Deep learning"):
        cols = [C.gidx[g] for g in C.var_order[:2000]]
        Xtr, mu, sd = B.zcols(C.Xo_g[tr][:, cols])
        yb = (C.y[tr] > np.median(C.y[tr])).astype(float)
        m = nets.moli_fit(Xtr, yb, extra["hp"], seed=0)
        po = nets.moli_predict(m, B.zcols(C.Xo_g[te][:, cols], mu, sd)[0])
        pp = nets.moli_predict(m, C.Xp_g[:, cols]) if full else None
        return po, pp, info
    raise KeyError(name)


def moli_grid(C):
    """Hyper-parameter grid with mean 3-fold CV AUC on the lab models (paper MOESM3 'Figure 3f' layout)."""
    cols = [C.gidx[g] for g in C.var_order[:2000]]
    yb = (C.y > np.median(C.y)).astype(int)
    rows = []
    skf = StratifiedKFold(3, shuffle=True, random_state=0)
    for gi, hp in enumerate(nets.MOLI_GRID):
        aucs = []
        for tr, te in skf.split(C.Xo_g, yb):
            Xtr, mu, sd = B.zcols(C.Xo_g[tr][:, cols])
            m = nets.moli_fit(Xtr, yb[tr].astype(float), hp, seed=gi)
            p = nets.moli_predict(m, B.zcols(C.Xo_g[te][:, cols], mu, sd)[0])
            if len(set(yb[te])) == 2:
                aucs.append(roc_auc_score(yb[te], p))
        rows.append(dict(hp, n_folds=3, average_AUC=float(np.mean(aucs)) if aucs else np.nan))
    g = pd.DataFrame(rows)
    best = int(g["average_AUC"].fillna(-1).values.argmax())
    g["selected"] = False
    g.loc[best, "selected"] = True
    return g, nets.MOLI_GRID[best]


# ------------------------------------------------------------------------------------------------ evaluation
def patient_metrics(C, pp, boot=True):
    r = pp[C.ti]
    out = {"n_treated": len(C.ti), "events_treated": int(C.E[C.ti].sum())}
    if not np.all(np.isfinite(r)) or np.std(r) == 0:
        return out
    if boot:
        c, lo, hi = B.cindex_ci(C.T[C.ti], C.E[C.ti], r)
    else:
        c, lo, hi = B.cindex(C.T[C.ti], C.E[C.ti], r), np.nan, np.nan
    out.update(C_index=c, C_lo=lo, C_hi=hi)
    out.update(B.cox_sd(C.T[C.ti], C.E[C.ti], r))
    ms, hi_g = B.median_split(C.T[C.ti], C.E[C.ti], r)
    out.update(ms)
    if len(C.ui):
        out.update(B.interaction(C.T[C.ti], C.E[C.ti], r, C.T[C.ui], C.E[C.ui], pp[C.ui]))
        out["n_untreated"] = len(C.ui)
        if boot:
            out["C_untreated"] = B.cindex(C.T[C.ui], C.E[C.ui], pp[C.ui])
    return out


def run_pair(bname, drug):
    f = OUT / f"{bname}__{drug}.pkl"
    if f.exists():
        return
    t0 = time.time()
    C = Ctx(bname, drug)
    res = dict(bundle=bname, drug=drug, meta=C.meta, n_models=len(C.S), n_proximal=len(C.F), cv=C.cv,
               n_treated=len(C.ti), n_untreated=len(C.ui), n_dropped_align=C.n_drop, y=C.y, S=C.S,
               F=[C.U[i] for i in C.F], T=C.T[C.ti], E=C.E[C.ti], rows=[], oof={}, risk={}, info={})
    if len(C.F) < 2:
        res["skip"] = f"{len(C.F)} proximal pathways"
    grid, hp = moli_grid(C)
    res["moli_grid"] = grid
    for name, fam in MODELS:
        if name.startswith("Random pathway"):
            continue
        if len(C.F) < 2 and ("network pathway" in name or "central" in name):
            continue
        try:
            oof = np.full(len(C.S), np.nan)
            for tr, te in C.folds:
                oof[te] = fit_predict(name, C, tr, te, extra={"hp": hp})[0]
            _, pp, info = fit_predict(name, C, np.arange(len(C.S)), np.arange(len(C.S)), full=True, extra={"hp": hp})
        except Exception as e:                       # noqa: BLE001
            print("ERR", bname, drug, name, repr(e), flush=True)
            continue
        row = dict(model=name, family=fam, **{f"loo_{k}": v for k, v in B.spearman_ci(oof, C.y).items()})
        row.update(patient_metrics(C, pp))
        res["rows"].append(row)
        res["oof"][name], res["risk"][name], res["info"][name] = oof, pp[C.ti], info
        res["risk_untreated_" + name] = pp[C.ui] if len(C.ui) else None
    # random pathway sets of the same size as the proximal set
    rng = np.random.default_rng(B.SEED)
    rr = []
    k = max(len(C.F), 2)
    for i in range(N_RANDOM):
        cols = list(rng.choice(len(C.U), k, replace=False))
        oof = np.full(len(C.S), np.nan)
        for tr, te in C.folds:
            oof[te] = rank_model(C.Xo_pw, C.Xp_pw, C.y, tr, te, cols, 1, False)[0]
        _, pp, sel, _ = rank_model(C.Xo_pw, C.Xp_pw, C.y, np.arange(len(C.S)), np.arange(len(C.S)), cols, 1, True)
        m = patient_metrics(C, pp, boot=False)
        m.update({f"loo_{a}": b for a, b in B.spearman_ci(oof, C.y).items()}, pathway=C.U[sel[0]])
        rr.append(m)
    rr = pd.DataFrame(rr)
    res["random_sets"] = rr
    med = {c: rr[c].median() for c in rr.columns if c != "pathway" and rr[c].dtype != object}
    res["rows"].append(dict(model="Random pathway sets (median of 100)", family="paper comparison", **med))
    pickle.dump(res, open(f, "wb"), protocol=4)
    print("done", bname, drug, len(C.S), "models", round(time.time() - t0), "s", flush=True)


if __name__ == "__main__":
    w, n = int(sys.argv[1]), int(sys.argv[2])
    todo = B.combos()
    only = sys.argv[3:] or None
    for i, (b, d) in enumerate(todo):
        if i % n == w and (only is None or b in only):
            run_pair(b, d)
