"""Kong et al. 2020 organoid-only strategy, every step of the paper, for one (model set, drug, TCGA cohort).

Steps (python/single_pathway_prediction.py and python/multiple_pathway_prediction.py, unchanged in substance):
 1. proximal pathways: network proximity z <= cut-off (primary -1.2816);
 2. pathway scores of the lab models (ssGSEA), z-scored across models;
 3. Ridge (RidgeCV cv=3, alpha 0.1..0.9), linear SVR, LinearRegression fitted on drug response (lower = more sensitive);
 4. pathway with the largest |coefficient| = biomarker (rank 1); patients scored as coef x z-scored pathway score
    (z-scored over all patients with survival), median split: <= median = predicted responder;
 5. log-rank test on overall survival in patients given the drug, and in patients with no drug record (untreated).
Added (all reported, none used to choose anything): Cox HR with 95% CI, events, 5-year survival, leave-one-out
predicted vs observed response, multi-pathway k = 2..10, paper-style bootstrap of random pathway sets (Fig 4e/f),
whole-pipeline permutation of the lab-model response (1000 shuffles), proximity cut-off sensitivity.
"""
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import pearsonr, spearmanr, ttest_ind  # noqa: E402
from sklearn.linear_model import LinearRegression, RidgeCV  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from sklearn.svm import SVR  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402
from statsmodels.duration.survfunc import SurvfuncRight, survdiff  # noqa: E402

PRIMARY = -1.2816
CUTOFFS = [-1.0, -1.1824, -1.2816, -1.645, -2.1925]
MODELS = ("Ridge", "SVR", "LinearRegression")
OK = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "green": "#009E73", "sky": "#56B4E9",
      "purple": "#CC79A7", "grey": "#777777", "black": "#000000", "yellow": "#F0E442"}
MCOL = {"Ridge": OK["blue"], "SVR": OK["verm"], "LinearRegression": OK["green"]}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})


def short(p, n=48):
    s = p.replace("REACTOME_", "").replace("_", " ").lower()
    return s if len(s) <= n else s[: n - 1] + "…"


def fit(ml, X, y):
    if ml == "Ridge":
        return RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(X, y)
    if ml == "SVR":
        return SVR(kernel="linear").fit(X, y)
    return LinearRegression().fit(X, y)


def coefs(ml, X, y):
    y = (np.asarray(y, float) - np.mean(y)) / (np.std(y) or 1.0)   # z-scored response: same Ridge/LR ranking, usable SVR
    c = fit(ml, X, y).coef_
    return np.asarray(c[0] if ml == "SVR" else c, dtype=float)


# ------------------------------------------------------------------ survival helpers
def surv_test(score, T, E):
    """Median split (score > median = predicted non-responder). Returns dict with p, HR (non-resp vs resp), CI ..."""
    g = (score > np.median(score)).astype(int)
    out = {"n": int(len(T)), "events": int(E.sum()), "n_nonresp": int(g.sum()), "n_resp": int((1 - g).sum()),
           "events_nonresp": int(E[g == 1].sum()), "events_resp": int(E[g == 0].sum())}
    if g.sum() in (0, len(g)) or E.sum() == 0:
        out.update(p=np.nan, hr=np.nan, hr_lo=np.nan, hr_hi=np.nan, surv5_resp=np.nan, surv5_nonresp=np.nan)
        return out, g
    out["p"] = float(survdiff(T, E, g)[1])
    try:
        r = PHReg(T, g.reshape(-1, 1).astype(float), status=E).fit(disp=0)
        b, se = float(r.params[0]), float(r.bse[0])
        out.update(hr=np.exp(b), hr_lo=np.exp(b - 1.96 * se), hr_hi=np.exp(b + 1.96 * se))
    except Exception:
        out.update(hr=np.nan, hr_lo=np.nan, hr_hi=np.nan)
    for lab, m in (("resp", g == 0), ("nonresp", g == 1)):
        sf = SurvfuncRight(T[m], E[m])
        t, s = sf.surv_times, sf.surv_prob
        out[f"surv5_{lab}"] = float(s[t <= 60][-1]) if (t <= 60).any() else 1.0
    return out, g


def km_draw(ax, T, E, g, title, res):
    for lab, m, col in (("Predicted responder", g == 0, OK["blue"]), ("Predicted non-responder", g == 1, OK["verm"])):
        sf = SurvfuncRight(T[m], E[m])
        t = np.r_[0, sf.surv_times]
        s = np.r_[1, sf.surv_prob]
        se = np.r_[0, sf.surv_prob_se]
        ax.step(t, s, where="post", color=col, lw=1.6, label=f"{lab} (n={int(m.sum())}, deaths={int(E[m].sum())})")
        ax.fill_between(t, np.clip(s - 1.96 * se, 0, 1), np.clip(s + 1.96 * se, 0, 1), step="post", color=col, alpha=0.12)
        cens = T[m][E[m] == 0]
        ax.plot(cens, [s[np.searchsorted(t, c, side="right") - 1] for c in cens], "|", color=col, ms=4, alpha=0.6)
    ax.axvline(60, color=OK["grey"], ls="--", lw=0.8)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Fraction alive")
    hr = res.get("hr", np.nan)
    ax.set_title(f"{title}\nlog-rank p = {res['p']:.3g}; HR = {hr:.2f} ({res['hr_lo']:.2f}-{res['hr_hi']:.2f})")
    ax.legend(loc="lower left", frameon=False, fontsize=7)


def panel_letter(ax, s):
    ax.text(-0.14, 1.08, s, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


# ------------------------------------------------------------------ the analysis
class Analysis:
    def __init__(self, set_id, set_label, project, drug, model_drug, osc, tsc, y, surv, treated, untreated, z,
                 z_source, outdir, response_label, n_perm=1000, n_boot=2000, isogenic=None, note="", seed=0,
                 cutoffs=CUTOFFS, model_word="organoids"):
        self.__dict__.update(locals())
        del self.__dict__["self"]
        self.outdir = Path(outdir)
        (self.outdir / "figures").mkdir(parents=True, exist_ok=True)
        (self.outdir / "tables").mkdir(parents=True, exist_ok=True)
        self.tag = f"{project}_{drug}".replace(" ", "_")
        # shared pathway universe and data
        self.universe = sorted(set(osc.index) & set(tsc.index))
        self.samples = [s for s in y.index if s in osc.columns]
        self.y = y.loc[self.samples].astype(float)
        pats = [p for p in tsc.columns if p in surv.index]
        self.P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[self.universe, pats].T.values), index=pats,
                              columns=self.universe)
        self.tp = [p for p in pats if p in treated]
        self.up = [p for p in pats if p in untreated]
        # negative follow-up days occur in the biotab files; set to 0 (log-rank unaffected, Cox needs >= 0)
        self.Tt, self.Et = np.clip(surv.loc[self.tp, "T"].values.astype(float), 0, None), surv.loc[self.tp, "E"].values.astype(int)
        self.Tu, self.Eu = np.clip(surv.loc[self.up, "T"].values.astype(float), 0, None), surv.loc[self.up, "E"].values.astype(int)
        self.Xall = pd.DataFrame(StandardScaler().fit_transform(osc.loc[self.universe, self.samples].T.values),
                                 index=self.samples, columns=self.universe)
        self._pcache = {}

    # p-value of one pathway with positive coefficient (sign only flips direction, not the partition)
    def path_stats(self, j, which="t"):
        key = (j, which)
        if key not in self._pcache:
            P = self.P.loc[self.tp if which == "t" else self.up, j].values
            T, E = (self.Tt, self.Et) if which == "t" else (self.Tu, self.Eu)
            r, _ = surv_test(P, T, E) if len(T) > 3 else ({"p": np.nan, "hr": np.nan}, None)
            self._pcache[key] = (r["p"], r["hr"])
        return self._pcache[key]

    def feats(self, cutoff):
        zz = self.z.dropna()
        return [p for p in zz.index[zz <= cutoff] if p in self.universe]

    def rank1(self, feats, y=None, ml="Ridge"):
        y = self.y.values if y is None else y
        X = self.Xall[feats].values
        c = coefs(ml, X, y)
        order = np.argsort(-np.abs(c))
        return c, order

    def run(self):
        rng = np.random.default_rng(self.seed)
        F = self.feats(PRIMARY)
        self.F = F
        res = {"set": self.set_id, "set_label": self.set_label, "project": self.project, "drug": self.drug,
               "model_drug": self.model_drug, "n_models": len(self.samples), "proximity_source": self.z_source,
               "n_pathways_universe": len(self.universe), "n_proximal": len(F), "response_metric": self.response_label,
               "note": self.note}
        if len(F) < 2 or len(self.samples) < 5 or len(self.tp) < 5:
            res["skipped"] = f"too little data (proximal={len(F)}, models={len(self.samples)}, treated={len(self.tp)})"
            self.res = res
            return res
        X, y = self.Xall[F].values, self.y.values
        self.coef = {ml: coefs(ml, X, y) for ml in MODELS}
        top = {ml: F[int(np.argmax(np.abs(c)))] for ml, c in self.coef.items()}
        res["biomarker"] = top["Ridge"]
        res["same_top_all_models"] = len(set(top.values())) == 1
        res.update({f"top_{ml}": top[ml] for ml in MODELS})
        # patient tests rank 1, every model
        self.km = {}
        rows = []
        for ml in MODELS:
            c = self.coef[ml]
            i = int(np.argmax(np.abs(c)))
            for grp, pats, T, E in (("treated", self.tp, self.Tt, self.Et), ("untreated", self.up, self.Tu, self.Eu)):
                if len(pats) < 6:
                    rows.append({"model": ml, "group": grp, "n": len(pats), "p": np.nan})
                    continue
                score = c[i] * self.P.loc[pats, F[i]].values
                r, g = surv_test(score, T, E)
                r.update(model=ml, group=grp, pathway=F[i], coef=c[i])
                rows.append(r)
                self.km[(ml, grp)] = (T, E, g, r)
        self.rank1_table = pd.DataFrame(rows)
        rt = self.rank1_table.set_index(["model", "group"])
        for grp in ("treated", "untreated"):
            if ("Ridge", grp) in rt.index:
                for k in ("n", "events", "p", "hr", "hr_lo", "hr_hi", "surv5_resp", "surv5_nonresp",
                          "events_nonresp", "events_resp"):
                    res[f"{grp}_{k}"] = rt.loc[("Ridge", grp)].get(k, np.nan)
        res["direction_expected"] = bool(res.get("treated_hr", np.nan) > 1)
        # leave-one-out predicted vs observed
        loo = {}
        for ml in MODELS:
            pred = np.empty(len(y))
            for k in range(len(y)):
                tr = np.arange(len(y)) != k
                sc = StandardScaler().fit(self.osc.loc[F, np.array(self.samples)[tr]].T.values)
                Xtr = sc.transform(self.osc.loc[F, np.array(self.samples)[tr]].T.values)
                Xte = sc.transform(self.osc.loc[F, [self.samples[k]]].T.values)
                mu, sd = y[tr].mean(), y[tr].std() or 1.0
                pred[k] = fit(ml, Xtr, (y[tr] - mu) / sd).predict(Xte)[0] * sd + mu
            loo[ml] = pred
            rs, ps = spearmanr(pred, y)
            rp, pp = pearsonr(pred, y) if np.std(pred) > 0 else (np.nan, np.nan)
            res[f"loo_{ml}_spearman"], res[f"loo_{ml}_spearman_p"] = rs, ps
            res[f"loo_{ml}_pearson"], res[f"loo_{ml}_pearson_p"] = rp, pp
        self.loo = pd.DataFrame(loo, index=self.samples).assign(observed=y)
        # multi-pathway k = 2..10
        mk = []
        for ml in MODELS:
            c = self.coef[ml]
            order = np.argsort(-np.abs(c))
            for k in range(1, min(10, len(F)) + 1):
                idx = order[:k]
                for grp, pats, T, E in (("treated", self.tp, self.Tt, self.Et), ("untreated", self.up, self.Tu, self.Eu)):
                    if len(pats) < 6:
                        continue
                    score = self.P.loc[pats, [F[i] for i in idx]].values @ c[idx]
                    r, _ = surv_test(score, T, E)
                    r.update(model=ml, k=k, group=grp)
                    mk.append(r)
        self.multik = pd.DataFrame(mk)
        mt = self.multik[(self.multik.group == "treated") & (self.multik.k >= 2)]
        res["multik_min_p"] = float(mt.p.min())
        res["multik_min_setting"] = f"{mt.loc[mt.p.idxmin(), 'model']} k={int(mt.loc[mt.p.idxmin(), 'k'])}"
        res["multik_n_p_below_0.05"] = int((mt.p < 0.05).sum())
        res["multik_n_tests"] = int(len(mt))
        r2 = mt[(mt.model == "Ridge") & (mt.k == 2)]
        res["ridge_k2_p"] = float(r2.p.iloc[0]) if len(r2) else np.nan
        # permutation null: shuffle lab-model response, whole rank-1 Ridge pipeline
        pobs = res["treated_p"]
        perm = []
        for _ in range(self.n_perm):
            c, order = self.rank1(F, rng.permutation(y))
            j = F[order[0]]
            p, hr = self.path_stats(j)
            perm.append((p, hr if c[order[0]] > 0 else 1 / hr if hr else np.nan, j))
        self.perm = pd.DataFrame(perm, columns=["p", "hr", "pathway"])
        res["perm_empirical_p"] = (np.sum(self.perm.p <= pobs + 1e-12) + 1) / (self.n_perm + 1)
        res["perm_frac_p_below_0.05"] = float((self.perm.p < 0.05).mean())
        # paper-style bootstrap of random pathway sets (Fig 4e/f)
        uni = np.array(self.universe)
        iso = self.isogenic
        boots = []
        for it in range(self.n_boot):
            rf = list(rng.choice(uni, size=len(F), replace=False))
            c, order = self.rank1(rf)
            top_rank, top_feat = None, None
            for rk, oi in enumerate(order, 1):
                j = rf[oi]
                p, _ = self.path_stats(j)
                okk = p < 0.05
                if okk and iso is not None:
                    okk = iso_pass(iso, j, c[oi])
                if okk:
                    top_rank, top_feat = rk, j
                    break
            p1, hr1 = self.path_stats(rf[order[0]])
            boots.append({"iteration": it + 1, "top_rank": top_rank, "feature": top_feat, "rank1_p": p1,
                          "rank1_pathway": rf[order[0]]})
        self.boot = pd.DataFrame(boots)
        # observed rank: first proximal pathway (by |Ridge coef|) meeting the same criteria
        c = self.coef["Ridge"]
        order = np.argsort(-np.abs(c))
        obs_rank = None
        for rk, oi in enumerate(order, 1):
            p, _ = self.path_stats(F[oi])
            okk = p < 0.05 and (iso is None or iso_pass(iso, F[oi], c[oi]))
            if okk:
                obs_rank = rk
                break
        res["boot_criteria"] = "survival p<0.05" + (" and isogenic t-test p<0.05 in the predicted direction" if iso is not None else "")
        res["boot_observed_rank"] = obs_rank
        tr = self.boot.top_rank.astype(float)
        res["boot_empirical_p_rank"] = float(((tr <= obs_rank).sum()) / self.n_boot) if obs_rank else np.nan
        res["boot_empirical_p_rank1p"] = float((np.sum(self.boot.rank1_p <= pobs + 1e-12) + 1) / (self.n_boot + 1))
        # z cut-off sensitivity
        cut = []
        self.cut_km = {}
        for zc in self.cutoffs:
            f = self.feats(zc)
            row = {"cutoff": zc, "n_proximal": len(f)}
            if len(f) >= 2:
                c, order = self.rank1(f)
                j = f[order[0]]
                row["pathway"] = j
                for grp, pats, T, E in (("treated", self.tp, self.Tt, self.Et), ("untreated", self.up, self.Tu, self.Eu)):
                    if len(pats) < 6:
                        continue
                    r, g = surv_test(c[order[0]] * self.P.loc[pats, j].values, T, E)
                    for k in ("p", "hr", "hr_lo", "hr_hi", "n", "events"):
                        row[f"{grp}_{k}"] = r[k]
                    if grp == "treated":
                        self.cut_km[zc] = (T, E, g, r, j)
            cut.append(row)
        self.cut = pd.DataFrame(cut)
        for zc in (-2.1925, -1.1824):
            rr = self.cut[self.cut.cutoff == zc]
            res[f"cut{zc}_p"] = float(rr["treated_p"].iloc[0]) if len(rr) and "treated_p" in rr else np.nan
        self.res = res
        self.save()
        self.figure()
        self.figure_cutoffs()
        return res

    # ------------------------------------------------------------------ outputs
    def save(self):
        t = self.outdir / "tables"
        pd.DataFrame({"pathway": self.F, **{f"coef_{ml}": self.coef[ml] for ml in MODELS},
                      "proximity_z": [self.z[f] for f in self.F]}).sort_values("coef_Ridge", key=abs, ascending=False) \
            .to_csv(t / f"{self.tag}_coefficients.tsv", sep="\t", index=False)
        self.rank1_table.to_csv(t / f"{self.tag}_rank1_tests.tsv", sep="\t", index=False)
        self.multik.to_csv(t / f"{self.tag}_multipathway_k.tsv", sep="\t", index=False)
        self.loo.to_csv(t / f"{self.tag}_loo_predictions.tsv", sep="\t")
        self.perm.to_csv(t / f"{self.tag}_permutation_null.tsv.gz", sep="\t", index=False)
        self.boot.to_csv(t / f"{self.tag}_bootstrap_random_pathways.tsv.gz", sep="\t", index=False)
        self.cut.to_csv(t / f"{self.tag}_zcutoff_sensitivity.tsv", sep="\t", index=False)
        with open(t / f"{self.tag}_summary.json", "w") as f:
            json.dump({k: (v.item() if hasattr(v, "item") else v) for k, v in self.res.items()}, f, indent=1, default=str)

    def figure(self):
        r, F = self.res, self.F
        fig, axes = plt.subplots(4, 3, figsize=(16, 19.5))
        A = axes.ravel()
        for ax, L in zip(A, "abcdefghijkl"):
            panel_letter(ax, L)
        # a proximity
        ax = A[0]
        zz = self.z.dropna().reindex(self.universe).dropna().sort_values()
        top = zz.head(15)[::-1]
        ax.barh(range(len(top)), top.values, color=[OK["blue"] if v <= PRIMARY else OK["grey"] for v in top.values])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([short(p, 40) for p in top.index], fontsize=7.5)
        ax.axvline(PRIMARY, color=OK["verm"], ls="--", lw=1)
        ax.set_xlabel("Network proximity z-score (lower = closer)")
        ax.set_title(f"Pathways closest to {self.drug.lower()} targets\n{len(F)} of {len(zz)} pathways have z <= -1.2816")
        # b heatmap
        ax = A[1]
        cr = pd.Series(self.coef["Ridge"], index=F)
        hp = cr.abs().sort_values(ascending=False).index[:25]
        hp = cr[hp].sort_values().index
        order = self.y.sort_values().index
        M = self.Xall.loc[order, hp].T.values
        im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_yticks(range(len(hp)))
        ax.set_yticklabels([short(p, 34) for p in hp], fontsize=6.5 if len(hp) > 18 else 7.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([str(s)[:10] for s in order], rotation=90, fontsize=6.5 if len(order) > 20 else 7.5)
        ax.set_title(f"Pathway activity in {self.model_word} (z-score)\ncolumns: most sensitive -> least sensitive")
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        # c Ridge coefficients
        ax = A[2]
        cc = pd.DataFrame(self.coef, index=F)
        topc = cc["Ridge"].abs().sort_values(ascending=False).index[:10][::-1]
        yy = np.arange(len(topc))
        for k, ml in enumerate(MODELS):
            v = cc.loc[topc, ml].values
            v = v / (np.abs(cc[ml]).max() or 1)
            ax.barh(yy + (k - 1) * 0.27, v, height=0.27, color=MCOL[ml], label=ml)
        ax.set_yticks(yy)
        ax.set_yticklabels([short(p, 40) for p in topc], fontsize=7.5)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("Coefficient (each model scaled to max |coef| = 1)")
        ax.set_title("Top-10 pathways by Ridge |coefficient|\n(positive = higher activity, more resistant)")
        ax.legend(frameon=False, loc="lower right")
        # d LOO
        ax = A[3]
        lo, hi = self.y.min(), self.y.max()
        for ml in MODELS:
            ax.scatter(self.loo["observed"], self.loo[ml], s=22, color=MCOL[ml], alpha=0.8,
                       label=f"{ml}: Spearman r={r[f'loo_{ml}_spearman']:.2f} (p={r[f'loo_{ml}_spearman_p']:.2g})")
        ax.plot([lo, hi], [lo, hi], color=OK["grey"], ls=":", lw=1)
        ax.set_xlabel(f"Observed response ({self.response_label})")
        ax.set_ylabel("Predicted (leave-one-out)")
        ax.set_title(f"Leave-one-out prediction in {len(self.samples)} {self.model_word}")
        ax.legend(frameon=False, fontsize=7)
        # e,f KM
        for ax, grp in ((A[4], "treated"), (A[5], "untreated")):
            if ("Ridge", grp) in self.km:
                T, E, g, rr = self.km[("Ridge", grp)]
                km_draw(ax, T, E, g, f"{self.project} patients {'given ' + self.drug.lower() if grp == 'treated' else 'with no drug record'}"
                        f"\nbiomarker: {short(r['biomarker'], 44)}", rr)
            else:
                ax.text(0.5, 0.5, f"fewer than 6 {grp} patients", ha="center", transform=ax.transAxes)
                ax.set_axis_off()
        # g multi-k
        ax = A[6]
        mk = self.multik
        for ml in MODELS:
            for grp, ls in (("treated", "-"), ("untreated", "--")):
                d = mk[(mk.model == ml) & (mk.group == grp)].sort_values("k")
                if len(d):
                    ax.plot(d.k, -np.log10(d.p), ls, marker="o", ms=3.5, color=MCOL[ml], label=f"{ml}, {grp}")
        ax.axhline(-np.log10(0.05), color=OK["grey"], ls=":", lw=1)
        ax.set_xlabel("Number of top pathways used (k)")
        ax.set_ylabel("-log10 log-rank p")
        ax.set_title("Multi-pathway models, k = 1..10\n(dotted line: p = 0.05)")
        ax.legend(frameon=False, fontsize=6.5, ncol=2)
        # h bootstrap rank
        ax = A[7]
        tr = self.boot.top_rank
        nf = int(tr.isna().sum())
        vals = tr.dropna().astype(int)
        mx = len(F)
        cnt = vals.value_counts().reindex(range(1, mx + 1), fill_value=0)
        ax.bar(cnt.index, cnt.values, color=OK["sky"])
        if r["boot_observed_rank"]:
            ax.axvline(r["boot_observed_rank"], color=OK["verm"], lw=1.5, label=f"proximal pathways: rank {r['boot_observed_rank']}")
            ax.legend(frameon=False)
        ax.set_xlabel("Best rank meeting the criteria (random pathway sets)")
        ax.set_ylabel("Iterations")
        ax.set_title(f"Bootstrap, {self.n_boot} random sets of {len(F)} pathways\n"
                     f"criteria: {r['boot_criteria'].replace(' and isogenic', '\nand isogenic')}; none found in {nf}\nempirical p = {r['boot_empirical_p_rank']:.3g}",
                     fontsize=8.5)
        # i bootstrap p hist
        ax = A[8]
        ax.hist(-np.log10(self.boot.rank1_p.dropna()), bins=40, color=OK["sky"])
        ax.axvline(-np.log10(r["treated_p"]), color=OK["verm"], lw=1.5, label="proximal pathways (observed)")
        ax.set_xlabel("-log10 p of the rank-1 pathway")
        ax.set_ylabel("Iterations")
        ax.set_title(f"Bootstrap: rank-1 pathway of random sets\nfraction with p <= observed = {r['boot_empirical_p_rank1p']:.3g}")
        ax.legend(frameon=False)
        # j permutation
        ax = A[9]
        ax.hist(-np.log10(self.perm.p.dropna()), bins=40, color=OK["purple"], alpha=0.8)
        ax.axvline(-np.log10(r["treated_p"]), color=OK["verm"], lw=1.5, label="real response")
        ax.set_xlabel("-log10 p (rank-1 pathway)")
        ax.set_ylabel("Shuffles")
        ax.set_title(f"Permutation null: {self.n_perm} shuffles of the {self.model_word}' response\n"
                     f"empirical p = {r['perm_empirical_p']:.3g}")
        ax.legend(frameon=False)
        # k z cut-off
        ax = A[10]
        cu = self.cut.dropna(subset=["treated_p"]) if "treated_p" in self.cut else self.cut.iloc[:0]
        ax.plot(cu.cutoff, -np.log10(cu.treated_p), "o-", color=OK["verm"], label="treated")
        if "untreated_p" in cu:
            ax.plot(cu.cutoff, -np.log10(cu.untreated_p), "s--", color=OK["blue"], label="untreated")
        for _, row in cu.iterrows():
            ax.annotate(f"{int(row.n_proximal)}", (row.cutoff, -np.log10(row.treated_p)), textcoords="offset points",
                        xytext=(0, 6), ha="center", fontsize=7)
        ax.axhline(-np.log10(0.05), color=OK["grey"], ls=":", lw=1)
        ax.axvline(PRIMARY, color=OK["grey"], ls="--", lw=0.8)
        ax.set_xlabel("Proximity z cut-off (numbers = pathways kept)")
        ax.set_ylabel("-log10 log-rank p (Ridge rank 1)")
        ax.set_title("Cut-off sensitivity\n(paper -1.2816; user's -2.1925 BLCA, -1.1824 LUAD)")
        ax.legend(frameon=False)
        # l forest
        ax = A[11]
        ft = self.rank1_table.dropna(subset=["p"]).reset_index(drop=True)
        for i, row in ft.iterrows():
            col = MCOL[row.model]
            ax.errorbar(row.hr, i, xerr=[[row.hr - row.hr_lo], [row.hr_hi - row.hr]], fmt="o" if row.group == "treated" else "s",
                        color=col, mfc=col if row.group == "treated" else "white", capsize=2)
        ax.set_yticks(range(len(ft)))
        ax.set_yticklabels([f"{m}, {g} (p={p:.2g})" for m, g, p in zip(ft.model, ft.group, ft.p)], fontsize=7.5)
        ax.axvline(1, color="k", lw=0.6)
        ax.set_xscale("log")
        from matplotlib.ticker import FixedLocator, NullFormatter
        ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5, 10, 20]))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("HR, predicted non-resp. vs resp. (log scale)")
        ax.set_title("Rank-1 pathway of each model\n(HR > 1: predicted non-responders die sooner)")
        ax.invert_yaxis()
        fig.suptitle(f"{self.set_label} -> TCGA-{self.project}, {self.drug.lower()} ({self.model_word}: {self.model_drug.lower()}, "
                     f"n={len(self.samples)})", fontsize=13, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        stem = self.outdir / "figures" / f"{self.tag}_full"
        fig.savefig(f"{stem}.png", dpi=300)
        fig.savefig(f"{stem}.pdf")
        plt.close(fig)
        self.write_caption(stem)

    def write_caption(self, stem):
        r = self.res
        t = lambda k: r.get(k, np.nan)  # noqa: E731
        txt = f"""# {self.set_label} -> TCGA-{self.project}, {self.drug.lower()}

All numbers below come from `tables/{self.tag}_*.tsv` and `{self.tag}_summary.json`.
Lab models: {len(self.samples)} {self.model_word}, response to {self.model_drug.lower()} measured as {self.response_label}
(lower = more sensitive). The response is z-scored before model fitting (this does not change the Ridge or
LinearRegression pathway ranking; it lets linear SVR, whose default margin of 0.1 would otherwise swallow AUC-scale
differences, return non-zero coefficients). {self.note}

**a.** Network proximity. Each bar is a Reactome pathway; a lower z-score means the pathway sits closer to the drug's
target proteins in the protein network than random proteins would. Source: {self.z_source}. Blue bars pass the
paper's cut-off (z <= -1.2816, red dashed line). {r['n_proximal']} pathways pass and are used as features.
**b.** Activity of the most important proximal pathways (rows) in each lab model (columns, ordered from most to least
sensitive). Red = higher than average, blue = lower.
**c.** Coefficients of the three regression models (Ridge, linear SVR, LinearRegression). The pathway with the largest
Ridge coefficient is the biomarker: `{short(r['biomarker'], 80)}`. Same top pathway in all three models: {r['same_top_all_models']}.
**d.** Leave-one-out test inside the lab models: each model is predicted from the others. Spearman r: Ridge
{t('loo_Ridge_spearman'):.2f} (p = {t('loo_Ridge_spearman_p'):.3g}), SVR {t('loo_SVR_spearman'):.2f} (p = {t('loo_SVR_spearman_p'):.3g}),
LinearRegression {t('loo_LinearRegression_spearman'):.2f} (p = {t('loo_LinearRegression_spearman_p'):.3g}).
**e.** Kaplan-Meier curves for TCGA patients given {self.drug.lower()}, split at the median of the biomarker score.
n = {t('treated_n')}, deaths = {t('treated_events')}, log-rank p = {t('treated_p'):.3g}, hazard ratio (HR, non-responders vs
responders) = {t('treated_hr'):.2f} (95% CI {t('treated_hr_lo'):.2f}-{t('treated_hr_hi'):.2f}). 5-year survival: predicted
responders {t('treated_surv5_resp'):.2f}, predicted non-responders {t('treated_surv5_nonresp'):.2f}.
**f.** The same split in patients with no drug record (comparison group): n = {t('untreated_n')}, deaths =
{t('untreated_events')}, p = {t('untreated_p'):.3g}, HR = {t('untreated_hr'):.2f} ({t('untreated_hr_lo'):.2f}-{t('untreated_hr_hi'):.2f}).
**g.** Multi-pathway models (top k pathways, weighted by their coefficients). Smallest p over the {t('multik_n_tests')}
treated tests with k = 2..10: {t('multik_min_p'):.3g} ({r.get('multik_min_setting')}); {t('multik_n_p_below_0.05')} of them have p < 0.05.
Ridge k = 2: p = {t('ridge_k2_p'):.3g}. Solid lines = treated, dashed = untreated.
**h-i.** Bootstrap as in the paper's Fig. 4e/f: {self.n_boot} random sets of {r['n_proximal']} pathways replace the proximal ones.
h: best coefficient rank at which a pathway meets "{r['boot_criteria']}"; the proximal set reaches it at rank
{r['boot_observed_rank']}; empirical p (random sets at least as good) = {t('boot_empirical_p_rank'):.3g}. i: p-value of the rank-1 pathway;
fraction of random sets with p at or below the observed = {t('boot_empirical_p_rank1p'):.3g}.
**j.** Permutation null: the lab-model responses were shuffled {self.n_perm} times and the whole rank-1 Ridge procedure rerun.
Empirical p = {t('perm_empirical_p'):.3g}; {100 * t('perm_frac_p_below_0.05'):.1f}% of shuffles reach p < 0.05.
**k.** The same analysis with other proximity cut-offs (numbers above points = pathways kept). Cut-off -2.1925
(user's BLCA choice): p = {t('cut-2.1925_p'):.3g}; -1.1824 (user's LUAD choice): p = {t('cut-1.1824_p'):.3g}.
**l.** Hazard ratios (95% CI) for the rank-1 pathway of each regression model; filled = treated, open = untreated.
"""
        open(f"{stem}.caption.md", "w").write(txt)

    def figure_cutoffs(self):
        ks = [k for k in self.cutoffs if k in self.cut_km]
        if not ks:
            return
        fig, axes = plt.subplots(1, len(ks), figsize=(4.6 * len(ks), 4.6))
        axes = np.atleast_1d(axes)
        for ax, zc, L in zip(axes, ks, "abcdefg"):
            T, E, g, rr, j = self.cut_km[zc]
            km_draw(ax, T, E, g, f"z <= {zc} ({int(self.cut.loc[self.cut.cutoff == zc, 'n_proximal'].iloc[0])} pathways)\n{short(j, 40)}", rr)
            ax.legend(loc="lower left", frameon=False, fontsize=6.5)
            panel_letter(ax, L)
        fig.suptitle(f"Supplement: proximity cut-off sensitivity, {self.set_label} -> TCGA-{self.project} {self.drug.lower()}", y=1.0)
        fig.tight_layout()
        stem = self.outdir / "figures" / f"{self.tag}_supp_zcutoff_KM"
        fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{stem}.pdf", bbox_inches="tight")
        plt.close(fig)
        lines = [f"# Cut-off sensitivity, {self.set_label} -> TCGA-{self.project} {self.drug.lower()}", "",
                 "Each panel repeats the rank-1 Ridge analysis in treated patients with a different proximity z cut-off.",
                 "-1.2816 is the paper's choice; -2.1925 and -1.1824 are the cut-offs used in the user's own BLCA and LUAD notebooks.", ""]
        for zc in ks:
            T, E, g, rr, j = self.cut_km[zc]
            lines.append(f"- z <= {zc}: biomarker {short(j, 70)}; n = {rr['n']}, deaths = {rr['events']}, log-rank p = {rr['p']:.3g}, "
                         f"HR = {rr['hr']:.2f} ({rr['hr_lo']:.2f}-{rr['hr_hi']:.2f}).")
        open(f"{stem}.caption.md", "w").write("\n".join(lines) + "\n")


def iso_pass(iso, pathway, coef):
    """Paper Fig 4 criterion: two-sided t-test p<0.05 between resistant and sensitive isogenic lines, in the direction the
    organoid model implies (coef > 0: higher activity = more resistant)."""
    if pathway not in iso["nes"].index:
        return False
    a = iso["nes"].loc[pathway, iso["resistant"]].values.astype(float)
    b = iso["nes"].loc[pathway, iso["sensitive"]].values.astype(float)
    t, p = ttest_ind(a, b)
    return bool(p < 0.05 and np.sign(t) == np.sign(coef))


def wrap(s, n=90):
    return "\n".join(textwrap.wrap(s, n))
