"""W3 positives, ovarian joint organoid/cell-line + patient model (nar/03_organoid_plus_patient), platinum and paclitaxel.

Original: shared pathway weights fitted to lab models AND treated patients (Cox blocks for survival cohorts, logistic blocks
for response cohorts), leave-one-cohort-out over the ovarian cohorts, 673 pathways ('all' feature set).  Its held-out scores
are read from nar/03_organoid_plus_patient/cache/<unit>__all__scores.tsv.gz (written by run_joint.py).
Zoo here: same training cohorts, same lab rows (aligned lab features, target = minus the rank-normalised lab sensitivity) and
same held-out patients.  Patient target (higher = worse): per-cohort standardised martingale residual of a null Cox model
(survival cohorts) or non-response (response cohorts).  Hyper-parameters: cohort-grouped CV on the training patients
(lab rows always in the training part), never the held-out cohort.
Usage: python pos_ovarian.py OV_PLATINUM OV_PACLITAXEL
"""
import os, sys, pickle, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
import pzoo as Z, pos_stats as S
ROOT = Z.ROOT; NAR = ROOT / "nar"
sys.path.insert(0, str(NAR / "03_organoid_plus_patient/code"))
import run_joint as RJ                                    # unit_data, lab_data, blocks_for, trainable, shortlist
from units import UNITS
OUT = Z.HERE / "tables" / "preds_pos"; OUT.mkdir(parents=True, exist_ok=True)
KW = ["CELL_CYCLE", "MITOTIC", "G1_S", "G2_M", "E2F", "DNA_REPLICATION", "PRE_REPLICATIVE", "CHECKPOINT"]
N_RANDOM = 100


def null_resid(T, E):
    T, E = np.asarray(T, float), np.asarray(E, float)
    o = np.argsort(T); t, e = T[o], E[o]
    risk = len(t) - np.arange(len(t))
    H = np.cumsum(e / risk)
    pos = np.searchsorted(t, T, side="right") - 1
    m = E - H[pos]
    return (m - m.mean()) / (m.std() or 1.0)


def run(uid):
    t0 = time.time()
    u = [x for x in UNITS if x["id"] == uid][0]
    D, ZD, Zall = RJ.unit_data(u)
    feats = list(Zall.columns)
    Xo, yo, labinfo = RJ.lab_data(u, Zall)
    Xo = Xo[feats].values; ylab = -np.asarray(yo.values, float)           # risk orientation: high = more resistant
    ylab = (ylab - ylab.mean()) / (ylab.std() or 1)
    short_names, how = RJ.shortlist(u["short"], feats)
    short = [feats.index(s) for s in short_names]
    prior = np.zeros(len(feats)); prior[short] = 1.0            # GCN node prior: member of the network-proximity shortlist
    train_c = RJ.trainable(D)
    orig = pd.read_csv(NAR / f"03_organoid_plus_patient/cache/{uid}__all__scores.tsv.gz", sep="\t", index_col=0)
    prol = [i for i, f in enumerate(feats) if any(k in f for k in KW)]
    names = ["Joint_original", "PatientOnly_original", "Proliferation_score"] + Z.MODEL_NAMES
    held = D.index[(D["group"] == "treated")]
    scores = pd.DataFrame(np.nan, index=held, columns=names)
    hps, rnd = [], {}
    rng = np.random.default_rng(Z.SEED)
    for c in train_c:
        blocks, used = RJ.blocks_for(D, [x for x in train_c if x != c], feats, ZD)
        rows, ys, gs = [], [], []
        for coh, kind, idx in used:
            x = D.loc[idx]
            if kind == "survival":
                y = null_resid(x["time"].values, x["event"].values)
            else:
                nr = 1.0 - x["responder"].values.astype(float); y = (nr - nr.mean()) / (nr.std() or 1)
            rows.append(ZD.loc[idx, feats].values); ys.append(y); gs += [coh] * len(idx)
        Xp = np.vstack(rows); yp = np.concatenate(ys); gp = np.array(gs)
        Xtr = np.vstack([Xp, Xo]); ytr = np.concatenate([yp, ylab]); n_p = len(yp)
        lab_ix = np.arange(n_p, len(ytr))
        folds = [(np.concatenate([a, lab_ix]), b) for a, b in GroupKFold(min(4, len(np.unique(gp)))).split(Xp, groups=gp)]
        gfull = np.concatenate([gp, np.array(["lab"] * len(ylab))])
        te = held[(D.loc[held, "unit_cohort"] == c).values]
        Zt = ZD.loc[te, feats].values
        print(f"--- {uid} held-out {c}: train patients {n_p} (cohorts {sorted(set(gp))}) + {len(ylab)} lab rows, test {len(te)} ({time.time()-t0:.0f}s)", flush=True)
        scores.loc[te, "Joint_original"] = -orig.loc[te, "joint"].values
        scores.loc[te, "PatientOnly_original"] = -orig.loc[te, "patient-only"].values
        from scipy.stats import spearmanr
        r = spearmanr(Xp[:, prol].mean(1), yp)[0]; sg = -1.0 if r < 0 else 1.0
        scores.loc[te, "Proliferation_score"] = sg * Zt[:, prol].mean(1)
        for n in Z.MODEL_NAMES:
            f, hp, best = Z.fit_model(n, Xtr, ytr, feats, short, groups=gfull, prior=prior, folds=folds)
            scores.loc[te, n] = f(Zt)
            hps.append(dict(held_out=c, model=n, inner_cv_spearman=best, **hp))
            print(f"   {n} {hp} inner={best:.3f} ({time.time()-t0:.0f}s)", flush=True)
        rs = []
        for i in range(N_RANDOM):
            ix = np.sort(rng.choice(len(feats), len(short), replace=False))
            f, hp, _ = Z.fit_model("L2_shortlist", Xtr, ytr, feats, ix, groups=gfull, prior=prior, folds=folds)
            rs.append(f(Zt))
        rnd[c] = pd.DataFrame(np.column_stack(rs), index=te)
    ev = D.loc[held, ["unit_cohort", "time", "event", "endpoint", "responder", "has_surv", "has_resp"]].rename(columns={"unit_cohort": "cohort", "time": "T", "event": "E"})
    res = dict(unit=uid, scores=scores, eval=ev, hps=pd.DataFrame(hps), random=rnd, n_short=len(short), n_feat=len(feats), train_c=train_c, how=how)
    pickle.dump(res, open(OUT / f"ovarian_joint_{uid}.pkl", "wb"))
    print("saved", uid, f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    for uid in sys.argv[1:]:
        run(uid)
