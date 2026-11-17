"""Learning curves: held-out performance as more lab models / more training patients are added.

For each unit (network shortlist features) the hyper-parameters are fixed at the median of the values the
main run chose by inner cross-validation (tables/<unit>__shortlist__folds.tsv), so each point costs one fit per fold.
Metric: pooled held-out concordance (C-index pairs for survival cohorts and AUROC pairs for response cohorts,
counted within cohorts only). Writes tables/learning_curves.tsv.
Usage: python learning_curves.py [UNIT ...] [--jobs 2]
"""
import argparse

import numpy as np
import pandas as pd

import common as C
from run_joint import blocks_for, lab_data, shortlist, trainable, unit_data
from units import UNITS

ORG_FRAC = (0.25, 0.5, 0.75, 1.0)
PAT_FRAC = (0.1, 0.25, 0.5, 0.75, 1.0)
REPS = 20


def folds(u, feats):
    D, ZD, Zall = unit_data(u)
    train_c = trainable(D)
    out = []
    if len(train_c) >= 2:
        for c in train_c:
            blocks, _ = blocks_for(D, [x for x in train_c if x != c], feats, ZD)
            te = D[(D["unit_cohort"] == c) & (D["group"] == "treated")]
            out.append((blocks, te))
    else:
        rng = np.random.default_rng(C.SEED)
        tr_idx = D.index[D["group"] == "treated"]
        fold = pd.Series(rng.permutation(len(tr_idx)) % 5, index=tr_idx)
        for f in range(5):
            blocks, _ = blocks_for(D.drop(fold.index[fold == f]), [train_c[0]], feats, ZD)
            out.append((blocks, D.loc[fold.index[fold == f]]))
    return out, ZD, Zall


def pooled(score_fn, fl, ZD, feats):
    num = den = 0.0
    for blocks, te, beta in score_fn(fl):
        s = ZD.loc[te.index, feats].values @ beta
        ts = te["has_surv"].values & (te["event"].sum() > 0)
        if te["has_surv"].sum() >= 10 and te.loc[te["has_surv"], "event"].sum() >= 5:
            m = te["has_surv"].values
            _, n_, d_ = C.concord(s[m], te.loc[m, "time"], te.loc[m, "event"])
        elif te["has_resp"].sum() >= 10:
            m = te["has_resp"].values
            a = C.auc_test(s[m], te.loc[m, "responder"])
            n_, d_ = a["num"], a["den"]
        else:
            continue
        num += n_
        den += d_
    return num / den if den else np.nan


def run(u):
    feats, _ = shortlist(u["short"], None or list(C.patients()[1].columns))
    fi = pd.read_csv(C.TABLES / f"{u['id']}__shortlist__folds.tsv", sep="\t")
    hp = fi.groupby("model")[["w", "lam"]].median()
    wj, lj = float(hp.loc["joint", "w"]), float(hp.loc["joint", "lam"])
    lp = float(hp.loc["patient-only", "lam"])
    fl, ZD, Zall = folds(u, feats)
    Xo_full, yo, _ = lab_data(u, Zall)
    Xo, yv = Xo_full[feats].values, yo.values
    rng = np.random.default_rng(C.SEED)
    rows = []
    n_lab = len(yv)
    for fr in ORG_FRAC:
        for rep in range(REPS if fr < 1 else 1):
            k = max(4, int(round(fr * n_lab)))
            sel = np.sort(rng.choice(n_lab, k, replace=False)) if fr < 1 else np.arange(n_lab)
            b_org, _ = C.organoid_only(Xo[sel], yv[sel])
            c_org = pooled(lambda fl: [(b, t, b_org) for b, t in fl], fl, ZD, feats)
            c_j = pooled(lambda fl: [(b, t, C.fit_joint(Xo[sel], yv[sel], b, wj, lj)[0]) for b, t in fl], fl, ZD, feats)
            rows += [dict(unit=u["id"], axis="lab models", fraction=fr, n=k, rep=rep, model="organoid-only", concordance=c_org),
                     dict(unit=u["id"], axis="lab models", fraction=fr, n=k, rep=rep, model="joint", concordance=c_j)]
    for fr in PAT_FRAC:
        for rep in range(10 if fr < 1 else 1):
            sub = []
            ntot = 0
            for b, t in fl:
                bs = []
                for bl in b:
                    n = bl.X.shape[0]
                    k = max(10, int(round(fr * n))) if fr < 1 else n
                    idx = np.sort(rng.choice(n, min(k, n), replace=False))
                    s = bl.subset(idx)
                    if s.kind == "cox" and s.event.sum() < 2:
                        continue
                    if s.kind == "logit" and len(np.unique(s.y)) < 2:
                        continue
                    bs.append(s)
                sub.append((bs, t))
                ntot = max(ntot, sum(x.X.shape[0] for x in bs))
            if any(len(bs) == 0 for bs, _ in sub):
                continue
            c_p = pooled(lambda fl_: [(b, t, C.fit_joint(None, None, b, 0.0, lp)[0]) for b, t in fl_], sub, ZD, feats)
            c_j = pooled(lambda fl_: [(b, t, C.fit_joint(Xo, yv, b, wj, lj)[0]) for b, t in fl_], sub, ZD, feats)
            rows += [dict(unit=u["id"], axis="patients", fraction=fr, n=ntot, rep=rep, model="patient-only", concordance=c_p),
                     dict(unit=u["id"], axis="patients", fraction=fr, n=ntot, rep=rep, model="joint", concordance=c_j)]
    out = pd.DataFrame(rows)
    out.to_csv(C.TABLES / f"{u['id']}__learning_curve.tsv", sep="\t", index=False)
    print(u["id"], "learning curve done", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("units", nargs="*")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    us = [u for u in UNITS if not a.units or u["id"] in a.units]
    from joblib import Parallel, delayed
    res = Parallel(n_jobs=a.jobs)(delayed(run)(u) for u in us)
    pd.concat(res).to_csv(C.TABLES / "learning_curves.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
