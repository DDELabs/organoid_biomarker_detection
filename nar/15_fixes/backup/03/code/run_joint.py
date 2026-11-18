"""Organoid + patient joint model, every unit, both feature sets (network shortlist, all pathways).

Usage: python run_joint.py [UNIT ...] [--features shortlist|all|both] [--jobs 2]
Writes tables/<unit>__<features>__{eval,interaction,scores,meta}.tsv and cache/<unit>__<features>__fit.pkl
"""
import argparse
import pickle
import time

import numpy as np
import pandas as pd

import common as C
from units import UNITS


# ----------------------------------------------------------------------------------------------- data per unit
def shortlist(name, columns):
    if name.startswith("paper:"):
        z = pd.read_csv(C.REPO / "data" / "coad_blca_organoid_drugs_zscore_result_reactome.txt", sep="\t", index_col=0)
        z = z[name.split(":")[1]].reindex(columns).dropna()
        return list(z.index[z <= C.CUT]), "Kong et al. 2020 proximity table (STRING>700), z <= -1.28"
    co = pd.read_csv(C.REPO / "robust" / "results" / "survlearn" / "coefficients.tsv.gz", sep="\t")
    p = [x for x in co.loc[co["drug"] == name, "pathway"] if x in columns]
    return p, f"STRING v12 proximity shortlist of {name} (robust/results/survlearn, z <= -1.28)"


def unit_data(u):
    P, S = C.patients()
    coh = list(u["cohorts"])
    rev = {}
    for a, b in u.get("merge", {}).items():
        rev.setdefault(b, []).append(a)
    orig = sum([rev.get(c, [c]) for c in coh], [])
    d = P[P["cohort"].isin(orig)].copy()
    Z = C.zscore_within(S.loc[d.index], d["cohort"])        # within-cohort standardisation (original cohorts)
    d["unit_cohort"] = d["cohort"].replace(u.get("merge", {}))
    rows = []
    for c, cfg in u["cohorts"].items():
        x = d[d["unit_cohort"] == c].copy()
        drugs = x["drugs"].fillna("")
        tr = drugs.str.contains(cfg["treated"], regex=True)
        cmp = drugs.str.contains(cfg["comparison"], regex=True) & ~tr if cfg.get("comparison") else pd.Series(False, index=x.index)
        if cfg.get("no_stage4"):
            st = x["stage"].map(C.stage_num)
            keep = ~(st == 4)
            tr, cmp = tr & keep, cmp & keep
        x["group"] = np.where(tr, "treated", np.where(cmp, "comparison", "other"))
        x = x[x["group"] != "other"]
        t, e, ep = C.endpoint(x[x["group"] == "treated"]) if tr.sum() else (None, None, None)
        if ep is not None:
            a = ep.lower()
            x["time"], x["event"] = x[f"{a}_months"].astype(float), x[f"{a}_event"].astype(float)
        else:
            x["time"], x["event"] = np.nan, np.nan
        x["endpoint"] = ep
        rows.append(x)
    D = pd.concat(rows)
    D["has_surv"] = D["time"].notna() & D["event"].notna() & (D["time"] > 0)
    D["has_resp"] = D["responder"].notna()
    return D, Z.loc[D.index], Z


def lab_data(u, Zref):
    Xs, ys, names = [], [], []
    for name, drug, lineage in u["lab"]:
        r = C.lab_set(name, drug, lineage)
        if r is None or len(r[1]) < 4:
            continue
        sc, y, grp = r
        X = C.zscore_within(sc.reindex(columns=Zref.columns), grp).fillna(0.0)   # within-set standardisation
        Xa, ndrop = C.align_models(X.values, Zref.values)
        Xs.append(pd.DataFrame(Xa, index=[f"{name}:{i}" for i in X.index], columns=Zref.columns))
        ys.append(pd.Series(C.rank_normal(y.values), index=Xs[-1].index))
        names.append(dict(set=name, drug=drug, lineage=lineage, n=len(y), dropped_axes=ndrop))
    return pd.concat(Xs), pd.concat(ys), names


MIN_N, MIN_EV = 15, 8          # a cohort trains as a Cox block with >= MIN_N treated patients and >= MIN_EV events


def blocks_for(D, cohorts, feats, Z, min_n=None, min_ev=None):
    """Training blocks of the given cohorts (treated patients only)."""
    min_n, min_ev = min_n or MIN_N, min_ev or MIN_EV
    out, used = [], []
    for c in cohorts:
        x = D[(D["unit_cohort"] == c) & (D["group"] == "treated")]
        s = x[x["has_surv"]]
        if len(s) >= min_n and s["event"].sum() >= min_ev:
            out.append(C.PatientBlock("cox", Z.loc[s.index, feats].values, s["time"], s["event"], strata=np.repeat(c, len(s))))
            used.append((c, "survival", s.index))
            continue
        r = x[x["has_resp"]]
        if len(r) and r["responder"].sum() >= 5 and (1 - r["responder"]).sum() >= 5:
            out.append(C.PatientBlock("logit", Z.loc[r.index, feats].values, y=r["responder"].values))
            used.append((c, "response", r.index))
    return out, used


def trainable(D, min_n=None, min_ev=None):
    min_n, min_ev = min_n or MIN_N, min_ev or MIN_EV
    ok = []
    for c in D["unit_cohort"].unique():
        x = D[(D["unit_cohort"] == c) & (D["group"] == "treated")]
        s = x[x["has_surv"]]
        r = x[x["has_resp"]]
        if (len(s) >= min_n and s["event"].sum() >= min_ev) or (len(r) and r["responder"].sum() >= 5 and (1 - r["responder"]).sum() >= 5):
            ok.append(c)
    return ok


# ----------------------------------------------------------------------------------------------- fitting
def fit_all(Xo, yo, blocks, seed):
    """Patient-only and joint models with inner-CV tuning. Returns dict of (beta, info)."""
    out = {}
    (w0, l0), res0 = C.tune(None, None, blocks, ws=(0.0,), seed=seed)
    out["patient-only"] = (C.fit_joint(None, None, blocks, 0.0, l0)[0], dict(w=0.0, lam=l0))
    (wj, lj), resj = C.tune(Xo, yo, blocks, ws=C.WS, seed=seed)
    out["joint"] = (C.fit_joint(Xo, yo, blocks, wj, lj)[0], dict(w=wj, lam=lj))
    lf = max((k for k in resj if k[0] == 0.5), key=lambda k: resj[k])[1]
    out["joint (w=0.5)"] = (C.fit_joint(Xo, yo, blocks, 0.5, lf)[0], dict(w=0.5, lam=lf))
    return out


def run_unit(u, features, seed=C.SEED):
    t0 = time.time()
    D, ZD, Zall = unit_data(u)
    Xo_full, yo, labinfo = lab_data(u, Zall)
    if features == "shortlist":
        feats, how = shortlist(u["short"], list(Zall.columns))
    else:
        feats, how = list(Zall.columns), "all 673 Reactome pathways"
    Xo = Xo_full[feats].values
    yo_v = yo.values
    beta_org, alpha_org = C.organoid_only(Xo, yo_v)
    mn = u.get("min_train", (None, None))
    train_c = trainable(D, *mn)
    models = ["organoid-only", "patient-only", "joint", "joint (w=0.5)"]
    scores = pd.DataFrame(np.nan, index=D.index, columns=models)
    scores["organoid-only"] = ZD[feats].values @ beta_org
    fold_info = []
    if len(train_c) >= 2:
        design = "leave-one-cohort-out"
        for c in train_c:
            blocks, used = blocks_for(D, [x for x in train_c if x != c], feats, ZD)
            fits = fit_all(Xo, yo_v, blocks, seed)
            idx = D.index[D["unit_cohort"] == c]
            for m, (b, info) in fits.items():
                scores.loc[idx, m] = ZD.loc[idx, feats].values @ b
                fold_info.append(dict(heldout=c, model=m, train_cohorts=";".join(f"{k}:{t}" for k, t, _ in used), **info))
        # cohorts too small to train on: frozen model trained on every trainable cohort
        rest = [c for c in D["unit_cohort"].unique() if c not in train_c]
        full_fits = None
        if rest:
            blocks, used = blocks_for(D, train_c, feats, ZD)
            full_fits = fit_all(Xo, yo_v, blocks, seed)
            for c in rest:
                idx = D.index[D["unit_cohort"] == c]
                for m, (b, info) in full_fits.items():
                    scores.loc[idx, m] = ZD.loc[idx, feats].values @ b
                    fold_info.append(dict(heldout=c, model=m, train_cohorts="all trainable", **info))
    else:
        design = "5-fold cross-fitting by patient within the single cohort (5 repeats)"
        c = train_c[0] if train_c else D["unit_cohort"].iloc[0]
        acc = {m: pd.Series(0.0, index=D.index) for m in models[1:]}
        cnt = pd.Series(0.0, index=D.index)
        cntc = 0
        for rep in range(5):
            rng = np.random.default_rng(seed + rep)
            tr_idx = D.index[(D["group"] == "treated")]
            fold = pd.Series(rng.permutation(len(tr_idx)) % 5, index=tr_idx)
            for f in range(5):
                Dtr = D.drop(fold.index[fold == f])
                blocks, used = blocks_for(Dtr, [c], feats, ZD, *mn)
                fits = fit_all(Xo, yo_v, blocks, seed + rep)
                te = fold.index[fold == f]
                comp = D.index[D["group"] != "treated"]
                for m, (b, info) in fits.items():
                    acc[m].loc[te] += ZD.loc[te, feats].values @ b
                    acc[m].loc[comp] += ZD.loc[comp, feats].values @ b
                    if rep == 0:
                        fold_info.append(dict(heldout=f"{c} fold {f}", model=m, train_cohorts=c, **info))
                cnt.loc[te] += 1
                cntc += 1
        cnt.loc[D.index[D["group"] != "treated"]] = cntc
        for m in models[1:]:
            scores[m] = acc[m] / cnt
    ev, it = evaluate(D, scores, models)
    meta = dict(unit=u["id"], features=features, n_features=len(feats), feature_source=how, design=design,
                trainable_cohorts=";".join(train_c), organoid_alpha=alpha_org,
                lab=";".join(f"{x['set']}:{x['drug']}{'/' + x['lineage'] if x['lineage'] else ''} n={x['n']} dropped_axes={x['dropped_axes']}" for x in labinfo),
                n_lab=int(len(yo)), seconds=round(time.time() - t0, 1))
    return D, scores, ev, it, pd.DataFrame(fold_info), meta, dict(feats=feats, beta_org=beta_org)


# ----------------------------------------------------------------------------------------------- evaluation
def evaluate(D, scores, models):
    ev, it = [], []
    for c in D["unit_cohort"].unique():
        x = D[D["unit_cohort"] == c]
        t = x[x["group"] == "treated"]
        for m in models:
            s = scores.loc[t.index, m]
            if s.isna().all():
                continue
            ts = t[t["has_surv"]]
            if len(ts) >= 10 and ts["event"].sum() >= 5:
                ss = s.loc[ts.index]
                hr = C.cox_hr(ss, ts["time"], ts["event"])
                cc, num, den = C.concord(ss, ts["time"], ts["event"])
                lr, _ = C.logrank_median(ss, ts["time"], ts["event"])
                ev.append(dict(cohort=c, model=m, outcome=f"survival ({ts['endpoint'].iloc[0]})", n=len(ts),
                               events=int(ts["event"].sum()), C=cc, num=num, den=den, logrank_p=lr, **hr))
            tr = t[t["has_resp"]]
            if len(tr) >= 10 and tr["responder"].sum() >= 3 and (1 - tr["responder"]).sum() >= 3:
                a = C.auc_test(s.loc[tr.index], tr["responder"])
                ev.append(dict(cohort=c, model=m, outcome="response", n=len(tr), events=int(tr["responder"].sum()),
                               C=a["AUROC"], num=a["num"], den=a["den"], p=a["p"]))
            # interaction with the comparison group
            cmp = x[x["group"] == "comparison"]
            if len(cmp) >= 10:
                both = pd.concat([t, cmp])
                arm = (both["group"] == "treated").astype(float)
                sb = scores.loc[both.index, m]
                bs = both[both["has_surv"]]
                if len(bs) and bs.loc[arm.loc[bs.index] == 1, "event"].sum() >= 5 and bs.loc[arm.loc[bs.index] == 0, "event"].sum() >= 5:
                    r = C.interaction_cox(sb.loc[bs.index], arm.loc[bs.index], bs["time"], bs["event"])
                    sc_cmp = C.cox_hr(sb.loc[bs.index][arm.loc[bs.index] == 0], bs.loc[arm.loc[bs.index] == 0, "time"],
                                      bs.loc[arm.loc[bs.index] == 0, "event"])
                    it.append(dict(cohort=c, model=m, outcome=f"survival ({bs['endpoint'].iloc[0]})",
                                   n_treated=int(arm.loc[bs.index].sum()), n_comparison=int((1 - arm.loc[bs.index]).sum()),
                                   events_treated=int(bs.loc[arm.loc[bs.index] == 1, "event"].sum()),
                                   events_comparison=int(bs.loc[arm.loc[bs.index] == 0, "event"].sum()),
                                   ratio=r["int_HR"], low=r["int_low"], high=r["int_high"], p=r["int_p"],
                                   HR_in_comparison=sc_cmp["HR"], HR_in_comparison_p=sc_cmp["p"], measure="interaction HR"))
                br = both[both["has_resp"]]
                if len(br) and br["responder"].sum() >= 5 and (1 - br["responder"]).sum() >= 5 and arm.loc[br.index].nunique() == 2:
                    try:
                        r = C.interaction_logit(sb.loc[br.index], arm.loc[br.index], br["responder"])
                        it.append(dict(cohort=c, model=m, outcome="response", n_treated=int(arm.loc[br.index].sum()),
                                       n_comparison=int((1 - arm.loc[br.index]).sum()),
                                       events_treated=int(br.loc[arm.loc[br.index] == 1, "responder"].sum()),
                                       events_comparison=int(br.loc[arm.loc[br.index] == 0, "responder"].sum()),
                                       ratio=r["int_OR"], low=r["int_low"], high=r["int_high"], p=r["int_p"],
                                       HR_in_comparison=r["score_OR_in_comparison"], measure="interaction OR"))
                    except Exception as exc:  # perfect separation in tiny arms
                        print("interaction logit failed", c, m, exc)
    ev = pd.DataFrame(ev)
    # pooled over held-out cohorts (survival: Cox stratified by cohort on within-cohort z scores; C summed within cohort)
    pooled = []
    for m in models:
        for kind in ("survival", "response"):
            e = ev[(ev["model"] == m) & ev["outcome"].str.startswith(kind)]
            if not len(e):
                continue
            row = dict(cohort="POOLED", model=m, outcome=kind, n=int(e["n"].sum()), events=int(e["events"].sum()),
                       C=e["num"].sum() / e["den"].sum(), num=e["num"].sum(), den=e["den"].sum(), n_cohorts=len(e))
            if kind == "survival" and len(e):
                idx, st = [], []
                for c in e["cohort"]:
                    t = D[(D["unit_cohort"] == c) & (D["group"] == "treated") & D["has_surv"]]
                    idx += list(t.index)
                    st += [c] * len(t)
                s = scores.loc[idx, m]
                s = s.groupby(np.array(st)).transform(lambda v: (v - v.mean()) / (v.std() or 1))
                hr = C.cox_hr(s, D.loc[idx, "time"], D.loc[idx, "event"], strata=np.array(st))
                row.update(hr)
            pooled.append(row)
    ev = pd.concat([ev, pd.DataFrame(pooled)], ignore_index=True)
    return ev, pd.DataFrame(it)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("units", nargs="*")
    ap.add_argument("--features", default="both")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    us = [u for u in UNITS if not a.units or u["id"] in a.units]
    feats = ["shortlist", "all"] if a.features == "both" else [a.features]
    jobs = [(u, f) for u in us for f in feats]

    def one(u, f):
        D, sc, ev, it, fi, meta, extra = run_unit(u, f)
        tag = f"{u['id']}__{f}"
        ev.assign(unit=u["id"], features=f).to_csv(C.TABLES / f"{tag}__eval.tsv", sep="\t", index=False)
        it.assign(unit=u["id"], features=f).to_csv(C.TABLES / f"{tag}__interaction.tsv", sep="\t", index=False)
        fi.assign(unit=u["id"], features=f).to_csv(C.TABLES / f"{tag}__folds.tsv", sep="\t", index=False)
        pd.Series(meta).to_csv(C.TABLES / f"{tag}__meta.tsv", sep="\t", header=False)
        out = D[["unit_cohort", "cohort", "group", "time", "event", "endpoint", "responder", "has_surv", "has_resp"]].join(sc)
        out.to_csv(C.CACHE / f"{tag}__scores.tsv.gz", sep="\t")
        with open(C.CACHE / f"{tag}__extra.pkl", "wb") as fh:
            pickle.dump(extra, fh)
        print(tag, meta["seconds"], "s", flush=True)
        return tag

    from joblib import Parallel, delayed
    Parallel(n_jobs=a.jobs)(delayed(one)(u, f) for u, f in jobs)


if __name__ == "__main__":
    main()
