"""W3 positives, lung untreated-prognosis score (nar/17_lung, analysis c).

Original: ridge Cox (lambda 300) with a chemo x pathway interaction learned on the chemo/no-chemo lung cohorts
(training_table_v2); its main-effect part = 'prognostic score'; tested in patients who got no chemo, one cohort at a time,
the tested cohort never in the training data (leave-one-cohort-out); HR per SD (Cox, adjusted for age/stage) pooled
by inverse variance.  Here the same training cohorts / test patients / pooling are used for every model of the zoo.
Zoo target = martingale residual of a stratified Cox model (age, stage, chemo) in the training cohorts (higher = worse
than expected), standardised inside each cohort.  Hyper-parameters: cohort-grouped 4-fold CV inside the training cohorts.
Usage: python pos_lung.py [--models M1,M2,...]
"""
import os, sys, pickle, time, argparse
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd
import pzoo as Z, pos_stats as S
ROOT = Z.ROOT; NAR = ROOT / "nar"
sys.path.insert(0, str(NAR / "17_lung/code"))
import run_lung as RL                                   # frame(), learn(), cox_row(), BASE, TC (training_table_v2), SL
TC, SL, BASE = RL.TC, RL.SL, RL.BASE
from obd import network as N
OUT = Z.HERE / "tables" / "preds_pos"; OUT.mkdir(parents=True, exist_ok=True)
KW = ["CELL_CYCLE", "MITOTIC", "G1_S", "G2_M", "E2F", "DNA_REPLICATION", "PRE_REPLICATIVE", "CHECKPOINT"]
N_RANDOM = 200
NORANDOM = os.environ.get('NORANDOM') == '1'
CUT = -1.2816


def martingale_target(tr):
    """Martingale residuals of a Cox model stratified by cohort with age, stage and chemo as covariates (Breslow)."""
    from statsmodels.duration.hazard_regression import PHReg
    X = tr[["age_n", "stage_n", "treated"]].astype(float).values
    r = PHReg(tr["time"].values.astype(float), X, status=tr["event"].values.astype(float), strata=tr["cohort"].values, ties="breslow").fit(disp=False)
    lp = X @ r.params
    res = np.zeros(len(tr))
    # baseline cumulative hazard per stratum (Breslow), step function evaluated at each subject's own time
    out = pd.Series(0.0, index=tr.index)
    for c, d in tr.groupby("cohort"):
        t = d["time"].values.astype(float); e = d["event"].values.astype(float); w = np.exp(lp[tr.index.get_indexer(d.index)])
        order = np.argsort(t); t_s, e_s, w_s = t[order], e[order], w[order]
        risk = np.cumsum(w_s[::-1])[::-1]
        H = np.cumsum(e_s / np.where(risk > 0, risk, np.inf))          # cumulative hazard increments at each sorted time
        # with ties the Breslow estimate is the same at equal times; evaluate at each subject's time (last tied index)
        pos = np.searchsorted(t_s, t, side="right") - 1
        H_t = H[pos]
        m = e - w * H_t
        out.loc[d.index] = (m - m.mean()) / (m.std() or 1.0)
    return out


def build_frames():
    P, Zm = TC.load()
    ADJ = [c for c, d in P[(P["cancer_group"] == "lung") & P["treated"].notna()].groupby("cohort")
           if d["treated"].nunique() == 2 and d["os_event"].notna().any()]
    allfr = RL.frame(P, ADJ)
    lungc = P[P["cancer_group"] == "lung"]
    cand = []
    for c, d in lungc.groupby("cohort"):
        if d["os_event"].notna().sum() + d["rfs_event"].notna().sum() < 15:
            continue
        if d["drugs"].fillna("").str.contains("ICI|ERLOTINIB|SORAFENIB").any():
            continue
        cand.append(c)
    un = RL.frame(P, cand, treated_known=False)
    un = un[~un["treated"].eq(1)]
    rec = un.groupby("cohort")["treated"].transform(lambda s: s.notna().any())
    un = un[~(rec & un["treated"].isna())]
    tests = [c for c in sorted(un["cohort"].unique()) if int(un.loc[un["cohort"] == c, "event"].sum()) >= 15]
    return P, Zm, ADJ, allfr, un, tests


def adj_cols(t):
    cols_c = [x for x in BASE if (x != "stage_n" or (t["stage_n"] != 0).mean() > 0.8)]
    if "age_n" in cols_c and (t["age_n"] != 0).mean() < 0.8:
        cols_c.remove("age_n")
    return cols_c


def main(only=None):
    t0 = time.time()
    P, Zm, ADJ, allfr, un, tests = build_frames()
    feats = list(Zm.columns)
    prox = N.precomputed_proximity()["CISPLATIN"].reindex(feats)
    short = [i for i, f in enumerate(feats) if pd.notna(prox.iloc[i]) and prox.iloc[i] <= CUT]
    prior = (-prox.fillna(0).values).astype(float)
    prol = [i for i, f in enumerate(feats) if any(k in f for k in KW)]
    short_names = [feats[i] for i in short]
    print("test cohorts", len(tests), "ADJ cohorts", ADJ, "shortlist", len(short), flush=True)
    names = (["Cox_ridge_all", "Cox_ridge_shortlist", "Proliferation_score"] + Z.MODEL_NAMES)
    if only:
        names = [n for n in names if n in only]
    scores = {n: pd.Series(np.nan, index=un.index) for n in names}
    hps, cache, rnd_rows = [], {}, []
    rng = np.random.default_rng(Z.SEED)
    for c in tests:
        te = un[un["cohort"] == c]
        tr = allfr[allfr["cohort"] != c]
        key = tuple(sorted(tr["cohort"].unique()))
        Zt = Zm.loc[te.index, feats].values
        if key not in cache:
            cache[key] = {}
            y = martingale_target(tr)
            Xtr = Zm.loc[tr.index, feats].values; grp = tr["cohort"].values
            cache[key].update(tr=tr, y=y.values, Xtr=Xtr, grp=grp)
            print(f"--- training set without {c}: {len(tr)} patients, {len(key)} cohorts ({time.time()-t0:.0f}s)", flush=True)
        K = cache[key]
        for n in names:
            if n in K.get("done", {}):
                continue
            if n == "Cox_ridge_all":
                a, g = RL.learn(K["tr"], Zm, feats); K.setdefault("fn", {})[n] = (lambda Xt, a=a: Xt @ a.values)
            elif n == "Cox_ridge_shortlist":
                a, g = RL.learn(K["tr"], Zm, short_names); K.setdefault("fn", {})[n] = (lambda Xt, a=a: Xt[:, short] @ a.values)
            elif n == "Proliferation_score":
                from scipy.stats import spearmanr
                r = spearmanr(K["Xtr"][:, prol].mean(1), K["y"])[0]
                sg = -1.0 if r < 0 else 1.0
                K.setdefault("fn", {})[n] = (lambda Xt, sg=sg: sg * Xt[:, prol].mean(1))
            else:
                f, hp, best = Z.fit_model(n, K["Xtr"], K["y"], feats, short, groups=K["grp"], prior=prior)
                K.setdefault("fn", {})[n] = f
                hps.append(dict(training_without=";".join(c_ for c_ in (set(un.cohort) | set(ADJ)) if c_ not in key) or "none", model=n, inner_cv_spearman=best, **hp))
                print(f"   {n} {hp} inner={best:.3f} ({time.time()-t0:.0f}s)", flush=True)
            K.setdefault("done", {})[n] = True
        for n in names:
            scores[n].loc[te.index] = np.asarray(K["fn"][n](Zt)).ravel()
        # random pathway sets of the size of the shortlist (ridge, inner CV), one set of 200 draws per test cohort, same draws for all
        if not NORANDOM:
            if "rand_fn" not in K:
                K["rand_fn"] = []
                for i in range(N_RANDOM):
                    ix = np.sort(rng.choice(len(feats), len(short), replace=False))
                    f, hp, _ = Z.fit_model("L2_shortlist", K["Xtr"], K["y"], feats, ix, groups=K["grp"], prior=prior)
                    K["rand_fn"].append(f)
            rs = np.column_stack([np.asarray(f(Zt)).ravel() for f in K["rand_fn"]])
            rnd_rows.append((c, pd.DataFrame(rs, index=te.index)))
    sc = pd.DataFrame(scores)
    ev = un.loc[sc.index, ["cohort", "time", "event", "endpoint", "age_n", "stage_n", "treated"]]
    ev = ev.rename(columns={"time": "T", "event": "E"})
    rand = pd.concat([d for _, d in rnd_rows]) if rnd_rows else pd.DataFrame()
    res = dict(task="lung_untreated_prognosis", kind="surv_pooled", orig="Cox_ridge_all", scores=sc, eval=ev, hps=pd.DataFrame(hps), random=rand,
               tests=tests, adj=ADJ, n_short=len(short), n_feat=len(feats), label="Lung untreated patients: prognosis score, leave-one-cohort-out")
    pickle.dump(res, open(OUT / "lung_untreated_prognosis.pkl", "wb"))
    print("saved", f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--models", default="")
    a = ap.parse_args()
    main(only=[m for m in a.models.split(",") if m] or None)
