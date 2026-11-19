"""Statistics for the W3 positives benchmark: reads tables/preds_pos/*.pkl, writes
tables/positives_all_models.tsv (model x task: metric, 95% CI, p, paired difference vs the original method),
tables/positives_per_cohort.tsv, tables/positives_model_comparison.tsv.
Paired bootstrap: the same resampled patients (within cohort) for every model.  Seeds fixed.
Usage: python pos_analyze.py [task ...]
"""
import os, sys, pickle, glob, time
os.environ["OMP_NUM_THREADS"] = "1"
from pathlib import Path
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd
from statsmodels.duration.hazard_regression import PHReg
import pos_stats as S

HERE = CODE.parent; TAB = HERE / "tables"; PRED = TAB / "preds_pos"
SEED = 20261119
FAM = {"L2": "zoo (673/all pathways)", "L1": "zoo (all pathways)", "ElasticNet": "zoo (all pathways)", "SVM_linear": "zoo (all pathways)", "SVM_RBF": "zoo (all pathways)",
       "RandomForest": "zoo (all pathways)", "XGBoost": "zoo (all pathways)", "MLP": "zoo (all pathways)", "NetSmoothL2": "zoo (all pathways)", "GCN_prior": "zoo (all pathways)"}


def family(m, orig):
    if m == orig: return "original method"
    if m.startswith(("Paper_", "Cox_ridge", "Joint_original", "PatientOnly", "Corr_signature")): return "original-type method"
    if m == "Proliferation_score": return "simple baseline"
    if m.endswith("shortlist"): return "zoo (network shortlist)"
    return "zoo (all pathways)"


def std(s):
    s = np.asarray(s, float); sd = s.std()
    return (s - s.mean()) / (sd if sd > 1e-12 else 1.0)


def cox_adj(T, E, s, cov=None, strata=None, ties="efron"):
    z = std(s)
    X = z[:, None] if cov is None or cov.shape[1] == 0 else np.column_stack([z, cov])
    try:
        r = PHReg(T, X, status=E, strata=strata, ties=ties).fit(disp=False)
        return float(r.params[0]), float(r.bse[0]), float(r.pvalues[0])
    except Exception:
        return np.nan, np.nan, np.nan


def iv_pool(co, se):
    co, se = np.asarray(co, float), np.asarray(se, float); ok = np.isfinite(co) & np.isfinite(se) & (se > 0)
    if ok.sum() == 0: return np.nan, np.nan
    w = 1 / se[ok] ** 2; b = (w * co[ok]).sum() / w.sum(); return float(b), float(np.sqrt(1 / w.sum()))


def ci_hr(b, se): return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)


def pnorm(b, se):
    from scipy.stats import norm
    return float(2 * norm.sf(abs(b / se))) if np.isfinite(b) and se > 0 else np.nan


# ------------------------------------------------------------------------------------------------ task specs
def spec_surv_single(r):
    ev = r["eval"]; T, E = ev["T"].values.astype(float), ev["E"].values.astype(int)
    sc = r["scores"].copy()
    def prim(s, ix): return S.cindex(T[ix], E[ix], s[ix])
    def eff(s, ix):
        return S.cox(T[ix], E[ix], s[ix])["coef"]
    def point(s):
        c = S.cox(T, E, s); sp = S.split_stats(T, E, s)
        return dict(est=S.cindex(T, E, s), hr=c["hr"], hr_lo=c["lo"], hr_hi=c["hi"], p_cox=c["p"], lr_p=sp["lr_p"], split_hr=sp["split_hr"], split_lo=sp["split_lo"], split_hi=sp["split_hi"])
    return dict(scores=sc, strata=None, prim=prim, eff=eff, point=point, n=len(ev), events=int(E.sum()), metric="C-index",
                p_kind=("log-rank (median split)" if r["task"].startswith(("bladder", "ovarian")) else "Cox (HR per SD)"), B=500)


def spec_bin(r):
    ev = r["eval"]; y = ev["label"].values.astype(int); sign = r["sign"]
    sc = r["scores"] * sign
    def prim(s, ix): return S.auc_fast(y[ix], s[ix])
    def point(s):
        a = S.auc_stats(y, s)
        return dict(est=a["auc"], p_cox=np.nan, p_one_sided=a["p"])
    return dict(scores=sc, strata=y, prim=prim, eff=None, point=point, n=len(ev), events=int(y.sum()), metric="AUC (mutant vs wild type)",
                p_kind="one-sided Mann-Whitney", B=1000)


def spec_lung(r):
    ev = r["eval"]; ev = ev[ev["cohort"].isin(r["tests"])]
    sc = r["scores"].loc[ev.index]
    import pos_lung as PL
    coh = ev["cohort"].values; T = ev["T"].values.astype(float); E = ev["E"].values.astype(int)
    codes = pd.factorize(coh)[0]; cohs = list(pd.factorize(coh)[1])
    covs = {}
    for g, c in enumerate(cohs):
        d = ev[ev["cohort"] == c]; cols = PL.adj_cols(d)
        covs[g] = ev.loc[ev["cohort"] == c, cols].values if cols else None
    pos = {g: np.where(codes == g)[0] for g in range(len(cohs))}
    loc = {g: {int(p): i for i, p in enumerate(pos[g])} for g in pos}
    def prim(s, ix):
        n = d_ = 0.0
        cg = codes[ix]
        for g in range(len(cohs)):
            sub = ix[cg == g]
            a, b = S.cparts(T[sub], E[sub], s[sub]); n += a; d_ += b
        return n / d_ if d_ else np.nan
    def eff(s, ix, ties="breslow"):
        co, se = [], []; cg = codes[ix]
        for g in range(len(cohs)):
            sub = ix[cg == g]
            cv = None if covs[g] is None else covs[g][[loc[g][int(p)] for p in sub]]
            b, s_, _ = cox_adj(T[sub], E[sub], s[sub], cv, ties=ties); co.append(b); se.append(s_)
        return iv_pool(co, se)[0]
    def point(s):
        co, se = [], []
        cg = codes
        for g in range(len(cohs)):
            sub = pos[g]; cv = covs[g]
            b, s_, _ = cox_adj(T[sub], E[sub], s[sub], cv, ties="breslow"); co.append(b); se.append(s_)
        b, sb = iv_pool(co, se); hr, lo, hi = ci_hr(b, sb)
        return dict(est=prim(s, np.arange(len(s))), hr=hr, hr_lo=lo, hr_hi=hi, p_cox=pnorm(b, sb))
    percoh = lambda s: [dict(cohort=cohs[g], n=len(pos[g]), events=int(E[pos[g]].sum()),
                             **dict(zip(("coef", "se", "p"), cox_adj(T[pos[g]], E[pos[g]], s[pos[g]], covs[g], ties="breslow"))),
                             c_index=S.cindex(T[pos[g]], E[pos[g]], s[pos[g]])) for g in range(len(cohs))]
    return dict(scores=sc, strata=coh, prim=prim, eff=eff, point=point, n=len(ev), events=int(E.sum()), metric="C-index (pooled over cohorts)",
                p_kind="Cox HR per SD, inverse-variance pooled", B=200, percoh=percoh, cohorts=cohs, codes=codes, T=T, E=E)


def spec_ov_surv(r):
    ev = r["eval"]; ev = ev[ev["has_surv"].astype(bool)]
    sc = r["scores"].loc[ev.index].copy()
    coh = ev["cohort"].values; T = ev["T"].values.astype(float); E = ev["E"].values.astype(int)
    codes = pd.factorize(coh)[0]; cohs = list(pd.factorize(coh)[1])
    def prim(s, ix):
        n = d_ = 0.0; cg = codes[ix]
        for g in range(len(cohs)):
            sub = ix[cg == g]; a, b = S.cparts(T[sub], E[sub], s[sub]); n += a; d_ += b
        return n / d_ if d_ else np.nan
    def zs(s, ix):
        z = np.empty(len(ix)); cg = codes[ix]
        for g in range(len(cohs)):
            m = cg == g; z[m] = std(s[ix][m])
        return z
    def eff(s, ix):
        z = zs(s, ix)
        try:
            return float(PHReg(T[ix], z[:, None], status=E[ix], strata=codes[ix], ties="efron").fit(disp=False).params[0])
        except Exception:
            return np.nan
    def point(s):
        ix = np.arange(len(s)); z = zs(s, ix)
        try:
            m = PHReg(T, z[:, None], status=E, strata=codes, ties="efron").fit(disp=False); b, se, p = float(m.params[0]), float(m.bse[0]), float(m.pvalues[0])
        except Exception:
            b = se = p = np.nan
        hr, lo, hi = ci_hr(b, se)
        return dict(est=prim(s, ix), hr=hr, hr_lo=lo, hr_hi=hi, p_cox=p)
    percoh = lambda s: [dict(cohort=cohs[g], n=int((codes == g).sum()), events=int(E[codes == g].sum()), c_index=S.cindex(T[codes == g], E[codes == g], s[codes == g])) for g in range(len(cohs))]
    return dict(scores=sc, strata=coh, prim=prim, eff=eff, point=point, n=len(ev), events=int(E.sum()), metric="C-index (pooled over cohorts)",
                p_kind="Cox HR per SD (stratified by cohort)", B=300, percoh=percoh, cohorts=cohs, codes=codes, T=T, E=E)


def spec_ov_resp(r):
    ev = r["eval"]; ev = ev[ev["has_resp"].astype(bool)]
    sc = r["scores"].loc[ev.index].copy()
    coh = ev["cohort"].values; y = (1 - ev["responder"].values.astype(float)).astype(int)         # 1 = non-response, score = risk
    codes = pd.factorize(coh)[0]; cohs = list(pd.factorize(coh)[1])
    from scipy.stats import rankdata
    def parts(yy, ss):
        n1 = yy.sum(); n0 = len(yy) - n1
        if n1 < 2 or n0 < 2: return 0.0, 0.0
        rk = rankdata(ss); return rk[yy == 1].sum() - n1 * (n1 + 1) / 2, n1 * n0
    def prim(s, ix):
        n = d_ = 0.0; cg = codes[ix]
        for g in range(len(cohs)):
            sub = ix[cg == g]; a, b = parts(y[sub], s[sub]); n += a; d_ += b
        return n / d_ if d_ else np.nan
    def point(s):
        from scipy.stats import norm
        ix = np.arange(len(s)); a = prim(s, ix)
        # stratified (van Elteren-like) test: z from pooled within-cohort Mann-Whitney U
        U = V = 0.0
        for g in range(len(cohs)):
            m = codes == g; yy, ss = y[m], s[m]; n1 = yy.sum(); n0 = len(yy) - n1
            if n1 < 2 or n0 < 2: continue
            u, _ = parts(yy, ss); U += u - n1 * n0 / 2; V += n1 * n0 * (len(yy) + 1) / 12
        p = float(norm.sf(U / np.sqrt(V))) if V > 0 else np.nan
        return dict(est=a, p_cox=np.nan, p_one_sided=p)
    percoh = lambda s: [dict(cohort=cohs[g], n=int((codes == g).sum()), events=int(y[codes == g].sum()), auc=prim(s, np.where(codes == g)[0])) for g in range(len(cohs))]
    return dict(scores=sc, strata=coh, prim=prim, eff=None, point=point, n=len(ev), events=int(y.sum()), metric="AUC for non-response (pooled within cohorts)",
                p_kind="one-sided stratified Mann-Whitney", B=500, percoh=percoh, cohorts=cohs, codes=codes, y=y)


# ------------------------------------------------------------------------------------------------ main loop
def analyse(tag, label, spec, orig, rand_metric=None, rand_note=""):
    t0 = time.time()
    sc = spec["scores"]; n = len(sc)
    models = [m for m in sc.columns if sc[m].notna().all()]
    sc = sc[models]
    M = sc.values
    bt = S.paired_boot(sc, {k: v for k, v in (("prim", spec["prim"]), ("eff", spec["eff"])) if v is not None}, n, spec["strata"], B=spec["B"], seed=SEED)
    rows, pc = [], []
    for m in models:
        s = M[:, models.index(m)]
        pt = spec["point"](s)
        lo, hi = S.ci(bt["prim"][m])
        d = bt["prim"][m] - bt["prim"][orig]
        row = dict(task=tag, label=label, model=m, family=family(m, orig), n=spec["n"], events=spec["events"], metric=spec["metric"], est=pt["est"], lo=lo, hi=hi,
                   p_test=pt.get("p_one_sided", pt.get("lr_p") if spec["p_kind"].startswith("log-rank") else pt.get("p_cox")), p_kind=spec["p_kind"],
                   hr=pt.get("hr", np.nan), hr_lo=pt.get("hr_lo", np.nan), hr_hi=pt.get("hr_hi", np.nan), p_cox=pt.get("p_cox", np.nan), lr_p=pt.get("lr_p", np.nan),
                   split_hr=pt.get("split_hr", np.nan), split_lo=pt.get("split_lo", np.nan), split_hi=pt.get("split_hi", np.nan),
                   diff_vs_orig=float(np.nanmean(bt["prim"][m]) - np.nanmean(bt["prim"][orig])) if False else None)
        # paired difference on the point estimates, CI from the paired bootstrap
        d0 = pt["est"] - spec["point"](M[:, models.index(orig)])["est"]
        dl, dh = S.ci(d)
        row.update(diff_vs_orig=d0 if m != orig else np.nan, diff_lo=dl if m != orig else np.nan, diff_hi=dh if m != orig else np.nan, diff_p=S.boot_p(d) if m != orig else np.nan)
        if "eff" in bt:
            de = bt["eff"][m] - bt["eff"][orig]
            row.update(dlogHR=(np.log(pt["hr"]) - np.log(spec["point"](M[:, models.index(orig)])["hr"])) if m != orig and np.isfinite(pt.get("hr", np.nan)) else np.nan,
                       dlogHR_lo=S.ci(de)[0] if m != orig else np.nan, dlogHR_hi=S.ci(de)[1] if m != orig else np.nan, dlogHR_p=S.boot_p(de) if m != orig else np.nan)
        rows.append(row)
        if "percoh" in spec:
            for x in spec["percoh"](s):
                pc.append(dict(task=tag, model=m, **x))
    df = pd.DataFrame(rows)
    df["q_diff"] = S.bh(df["diff_p"].values)
    # random pathway sets (200 draws) of the shortlist size, ridge: median and 2.5-97.5% range of the same metric
    if rand_metric is not None and len(rand_metric):
        rm = np.asarray(rand_metric, float); rm = rm[np.isfinite(rm)]
        o_est = df.loc[df.model == orig, "est"].iloc[0]
        df = pd.concat([df, pd.DataFrame([dict(task=tag, label=label, model="RandomSets_median", family="random pathway sets", n=spec["n"], events=spec["events"], metric=spec["metric"],
                                             est=float(np.median(rm)), lo=float(np.percentile(rm, 2.5)), hi=float(np.percentile(rm, 97.5)),
                                             p_test=float((rm >= o_est).mean()), p_kind="share of random sets with metric >= original method",
                                             diff_vs_orig=float(np.median(rm) - o_est), note=f"{len(rm)} draws; lo/hi = 2.5-97.5% range across draws. {rand_note}")])], ignore_index=True)
    print(tag, f"{time.time()-t0:.0f}s", flush=True)
    return df, pd.DataFrame(pc)


def run_all(only=None):
    out, pcs, hp_out = [], [], []
    for f in sorted(glob.glob(str(PRED / "*.pkl"))):
        r = pickle.load(open(f, "rb"))
        name = Path(f).stem
        if only and not any(o in name for o in only):
            continue
        if name == "lung_untreated_prognosis":
            spec = spec_lung(r)
            rm = None
            if len(r["random"]):
                rr = r["random"].loc[spec["scores"].index]
                rm = [spec["prim"](rr.iloc[:, j].values, np.arange(len(rr))) for j in range(rr.shape[1])]
            df, pc = analyse(name, r["label"], spec, r["orig"], rm, "Random sets use 47 pathways (cisplatin shortlist size).")
            hp = r["hps"].assign(task=name)
        elif name.startswith("ovarian_joint_"):
            uid = r["unit"]
            for kind, sp in (("survival", spec_ov_surv(r)), ("response", spec_ov_resp(r))):
                tag = f"ovarian_{uid.split('_')[1].lower()}_joint_{kind}"
                rm = None
                if kind == "survival":
                    ev = r["eval"]
                    ids = sp["scores"].index
                    rr = pd.concat([r["random"][c] for c in r["random"]]).reindex(ids)
                    rm = [sp["prim"](rr.iloc[:, j].values, np.arange(len(rr))) for j in range(rr.shape[1])]
                else:
                    ids = sp["scores"].index
                    rr = pd.concat([r["random"][c] for c in r["random"]]).reindex(ids)
                    rm = [sp["prim"](rr.iloc[:, j].values, np.arange(len(rr))) for j in range(rr.shape[1])]
                lab = f"Ovarian {uid.split('_')[1].lower()} joint organoid+patient model, {kind}, leave-one-cohort-out"
                df, pc = analyse(tag, lab, sp, "Joint_original", rm, f"Random sets use {r['n_short']} pathways (network shortlist size).")
                out.append(df); pcs.append(pc)
            hp_out.append(r["hps"].assign(task=name)); continue
        else:
            spec = spec_bin(r) if r["kind"] == "bin" else spec_surv_single(r)
            rm = r["random"]["metric"].values if len(r["random"]) else None
            df, pc = analyse(name, r["label"], spec, r["orig"], rm, f"Random sets use {r['n_short']} pathways (network shortlist size).")
            hp = r["hps"].assign(task=name)
        out.append(df); pcs.append(pc); hp_out.append(hp)
    A = pd.concat(out, ignore_index=True)
    sfx = os.environ.get("POS_TAG", "" if not only else "_partial")
    A.to_csv(TAB / f"positives_all_models{sfx}.tsv", sep="\t", index=False)
    if pcs:
        pd.concat(pcs, ignore_index=True).to_csv(TAB / f"positives_per_cohort{sfx}.tsv", sep="\t", index=False)
    pd.concat(hp_out, ignore_index=True).to_csv(TAB / f"positives_chosen_hyperparameters{sfx}.tsv", sep="\t", index=False)
    return A


if __name__ == "__main__":
    run_all(sys.argv[1:] or None)
