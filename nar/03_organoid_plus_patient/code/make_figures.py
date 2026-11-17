"""Figures and summary tables from the run_joint / learning_curves / align_figs outputs.

Main: figures/main_<cancer>.png (composite per cancer, main unit, network shortlist).
Supplementary: figures/supp/unit_<unit>.png (per drug: KM per held-out cohort x model, shortlist vs all, interaction).
Cross-cancer: figures/forest_cross_cancer.png; tables/summary_main.tsv, tables/summary_all.tsv, tables/interaction_all.tsv.
"""
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import common as C
import plotstyle as PS
from units import MAIN_UNIT, UNITS

MODELS = ["organoid-only", "patient-only", "joint"]
UNIT = {u["id"]: u for u in UNITS}


def load(unit, feats):
    t = f"{unit}__{feats}"
    ev = pd.read_csv(C.TABLES / f"{t}__eval.tsv", sep="\t")
    it = pd.read_csv(C.TABLES / f"{t}__interaction.tsv", sep="\t") if (C.TABLES / f"{t}__interaction.tsv").stat().st_size > 5 else pd.DataFrame()
    sc = pd.read_csv(C.CACHE / f"{t}__scores.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(C.TABLES / f"{t}__meta.tsv", sep="\t", header=None, index_col=0)[1]
    return ev, it, sc, meta


def km_split(ax, sc, model, title, cohorts=None):
    """KM of held-out treated patients split at the within-cohort median of the score."""
    d = sc[(sc["group"] == "treated") & sc["has_surv"] & sc[model].notna()]
    if cohorts is not None:
        d = d[d["unit_cohort"].isin(cohorts)]
    if len(d) < 10:
        ax.text(0.5, 0.5, "no survival data", ha="center", transform=ax.transAxes); ax.set_axis_off(); return None
    hi = d.groupby("unit_cohort")[model].transform(lambda v: v > v.median()).astype(bool)
    for flag, col, lab in ((True, PS.CB[2], "high score (predicted benefit)"), (False, PS.CB[3], "low score")):
        x = d[hi == flag]
        t, s = C.SV.km_curve(x["time"], x["event"])
        ax.step(t, s, where="post", color=col, lw=1.5, label=f"{lab}, n={len(x)}, events={int(x['event'].sum())}")
    p = C.SV.logrank(d.loc[hi, "time"], d.loc[hi, "event"], d.loc[~hi, "time"], d.loc[~hi, "event"])
    ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("event-free probability")
    ax.set_title(f"{title}\nlog-rank p = {p:.3g}", fontsize=8)
    ax.legend(fontsize=7, loc="lower left", frameon=False)
    return p


def conc_bars(ax, ev, title):
    e = ev[ev["model"].isin(MODELS)].copy()
    e["lab"] = e["cohort"] + "\n" + e["outcome"].str.replace("survival ", "").str.replace("response", "resp.")
    cats = list(dict.fromkeys(e.loc[e["cohort"] != "POOLED", "lab"])) + list(dict.fromkeys(e.loc[e["cohort"] == "POOLED", "lab"]))
    x = np.arange(len(cats))
    for j, m in enumerate(MODELS):
        v = e[e["model"] == m].set_index("lab")["C"].reindex(cats)
        ax.bar(x + (j - 1) * 0.27, v.values, 0.27, color=PS.MODEL_COL[m], label=m)
    ax.axhline(0.5, color="k", lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels(cats, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("held-out C-index / AUROC"); ax.set_ylim(0.2, max(0.85, float(e["C"].max()) + 0.08))
    ax.set_xlim(-0.6, len(cats) - 0.4)
    ax.set_title(title, fontsize=8); ax.legend(fontsize=7, ncol=3, loc="upper left", frameon=False)


def hr_forest(ax, ev, title):
    if "HR" not in ev:
        ev = ev.assign(HR=np.nan, HR_low=np.nan, HR_high=np.nan)
    e = ev[ev["model"].isin(MODELS) & ev["outcome"].str.startswith("survival") & ev["HR"].notna()].copy()
    cats = list(dict.fromkeys(e.loc[e["cohort"] != "POOLED", "cohort"])) + (["POOLED"] if (e["cohort"] == "POOLED").any() else [])
    if not cats:
        ax.text(0.5, 0.5, "no held-out survival cohort\n(no public survival times)", ha="center", transform=ax.transAxes); ax.set_axis_off(); return
    for j, m in enumerate(MODELS):
        for i, c in enumerate(cats):
            r = e[(e["model"] == m) & (e["cohort"] == c)]
            if not len(r):
                continue
            r = r.iloc[0]
            y = i + (j - 1) * 0.25
            ax.errorbar(r["HR"], y, xerr=[[r["HR"] - r["HR_low"]], [r["HR_high"] - r["HR"]]], fmt="o" if c != "POOLED" else "D",
                        color=PS.MODEL_COL[m], ms=4, lw=1, capsize=2, label=m if i == 0 else None)
    ax.axvline(1, color="k", lw=0.8, ls="--"); ax.set_xscale("log")
    ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("HR per SD of score (<1 = higher score, longer survival)")
    ax.set_title(title, fontsize=8); ax.legend(fontsize=7, frameon=False, loc="best")


def int_forest(ax, it, title):
    if it is None or not len(it):
        ax.text(0.5, 0.5, "no comparison group\nin this cancer", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    i2 = it[it["model"].isin(MODELS)].copy()
    i2["lab"] = i2["cohort"] + " " + i2["outcome"].str.replace("survival ", "").str.replace("response", "resp.")
    cats = list(dict.fromkeys(i2["lab"]))
    for j, m in enumerate(MODELS):
        for i, c in enumerate(cats):
            r = i2[(i2["model"] == m) & (i2["lab"] == c)]
            if not len(r):
                continue
            r = r.iloc[0]
            ax.errorbar(r["ratio"], i + (j - 1) * 0.25, xerr=[[r["ratio"] - r["low"]], [r["high"] - r["ratio"]]], fmt="s",
                        color=PS.MODEL_COL[m], ms=4, lw=1, capsize=2, label=m if i == 0 else None)
    ax.axvline(1, color="k", lw=0.8, ls="--"); ax.set_xscale("log")
    ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("treatment x score interaction\n(HR, or OR for response)")
    ax.set_title(title, fontsize=8); ax.legend(fontsize=7, frameon=False, loc="best")


def lc_panel(ax, unit):
    f = C.TABLES / f"{unit}__learning_curve.tsv"
    if not f.exists():
        ax.set_axis_off(); return
    lc = pd.read_csv(f, sep="\t")
    for axis, ls, mk in (("lab models", "-", "o"), ("patients", "--", "s")):
        for m in ("organoid-only", "patient-only", "joint"):
            x = lc[(lc["axis"] == axis) & (lc["model"] == m)]
            if not len(x):
                continue
            g = x.groupby("fraction")["concordance"]
            mu, lo, hi = g.mean(), g.quantile(0.1), g.quantile(0.9)
            lab = f"{m}: vary {axis}"
            ax.plot(mu.index * 100, mu.values, ls=ls, marker=mk, color=PS.MODEL_COL[m], ms=3, lw=1.2, label=lab)
            ax.fill_between(mu.index * 100, lo.values, hi.values, color=PS.MODEL_COL[m], alpha=0.15, lw=0)
    ax.axhline(0.5, color="k", lw=0.8, ls=":")
    ax.set_xlabel("% of lab models (solid) or training patients (dashed) used")
    ax.set_ylabel("pooled held-out concordance")
    ax.set_title("Learning curves (mean, 10th-90th percentile band)", fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="best")


def align_panel(axes, cancer):
    f = C.CACHE / f"align_coords_{cancer}.pkl"
    if not f.exists():
        for a in axes:
            a.set_axis_off()
        return
    d = pd.read_pickle(f)
    lab = d["labels"]
    for a, key, ttl in ((axes[0], "raw", "Before: raw scores (PCA)"), (axes[1], "aligned", "After: standardised + aligned (PCA)")):
        Y = d["coords"][key][0]
        t = lab == "tumours"
        a.scatter(Y[t, 0], Y[t, 1], s=2, c=PS.GREY, alpha=0.3, lw=0, label="tumours", rasterized=True)
        for k, s in enumerate(pd.unique(lab[~t])):
            m = lab == s
            a.scatter(Y[m, 0], Y[m, 1], s=10, c=PS.CB[k % len(PS.CB)], edgecolor="k", lw=0.3, label=s, zorder=3)
        a.set_xticks([]); a.set_yticks([]); a.set_xlabel("PC1"); a.set_ylabel("PC2"); a.set_title(ttl, fontsize=8)
    axes[0].legend(fontsize=6, frameon=False, loc="best", markerscale=1.3)


def response_panel(ax, sc, model, title):
    """For response-only units (liver): responder rate by score half in each arm."""
    d = sc[sc["has_resp"] & sc[model].notna()].copy()
    if not len(d):
        ax.set_axis_off(); return
    d["hi"] = d.groupby("unit_cohort")[model].transform(lambda v: v > v.median())
    rows = d.groupby(["group", "hi"])["responder"].agg(["mean", "size", "sum"]).reset_index()
    x = 0
    ticks, labs = [], []
    for g, col in (("treated", PS.CB[5]), ("comparison", PS.GREY)):
        for h in (False, True):
            r = rows[(rows["group"] == g) & (rows["hi"] == h)]
            if not len(r):
                continue
            r = r.iloc[0]
            ax.bar(x, r["mean"], color=col if h else "white", edgecolor=col, hatch=None if h else "//", lw=1.2)
            ax.text(x, r["mean"] + 0.01, f"{int(r['sum'])}/{int(r['size'])}", ha="center", fontsize=7)
            ticks.append(x); labs.append(f"{'drug' if g == 'treated' else 'control'}\n{'high' if h else 'low'}")
            x += 1
        x += 0.5
    ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=7); ax.set_ylabel("responder fraction")
    ax.set_title(title, fontsize=8)


def main_fig(cancer):
    unit = MAIN_UNIT[cancer]
    u = UNIT[unit]
    ev, it, sc, meta = load(unit, "shortlist")
    fig = plt.figure(figsize=(12, 10.5))
    gs = fig.add_gridspec(3, 3, hspace=0.75, wspace=0.42)
    a0, a1 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    align_panel((a0, a1), cancer)
    a2 = fig.add_subplot(gs[0, 2]); lc_panel(a2, unit)
    a3 = fig.add_subplot(gs[1, :2]); conc_bars(a3, ev, f"Held-out patients: concordance by cohort ({meta['n_features']} network-proximal pathways)")
    a4 = fig.add_subplot(gs[1, 2]); hr_forest(a4, ev, "Held-out treated patients: HR per SD")
    a5 = fig.add_subplot(gs[2, 0])
    has_surv = sc["has_surv"].any() and (sc["group"] == "treated").any()
    if sc.loc[sc["group"] == "treated", "has_surv"].sum() >= 10:
        km_split(a5, sc, "joint", "Joint model, held-out treated (pooled)")
    else:
        response_panel(a5, sc, "joint", "Joint model: responders by score half")
    a6 = fig.add_subplot(gs[2, 1])
    if sc.loc[sc["group"] == "treated", "has_surv"].sum() >= 10:
        km_split(a6, sc, "organoid-only", "Organoid-only model, held-out treated (pooled)")
    else:
        response_panel(a6, sc, "organoid-only", "Organoid-only model: responders by score half")
    a7 = fig.add_subplot(gs[2, 2]); int_forest(a7, it, "Drug vs comparison: does the score\nmatter more with the drug?")
    PS.letters([a0, a1, a2, a3, a4, a5, a6, a7], x=-0.1)
    fig.suptitle(f"{cancer.capitalize()} - {u['drug']}: lab models + patients", fontsize=12, y=0.995)
    PS.save(fig, C.FIGS / f"main_{cancer}", caption=main_caption(cancer, u, meta))


def main_caption(cancer, u, meta):
    return f"""**{cancer.capitalize()}, {u['drug']}: a model that learns from lab models and patients together.**
Lab models: {meta['lab']}. Pathway features: {meta['n_features']} network-proximal Reactome pathways ({meta['feature_source']}).
Patients enter training only through cross-fitting ({meta['design']}): every score shown was made by a model that never saw that patient's cohort (or, with one cohort, that patient).
(a) Lab models (coloured) and tumours (grey) on the first two principal components of the raw pathway scores. (b) The same after within-cohort standardisation and alignment (model-only directions removed, then covariance matched to tumours; see alignment_{cancer}.png).
(c) Learning curves: pooled held-out concordance when only part of the lab models (solid lines) or part of the training patients (dashed lines) is used; bands = 10th-90th percentile over random subsets. 0.5 = no better than chance.
(d) Held-out concordance in every cohort (C-index for survival, AUROC for response; 0.5 = chance, higher = better) for the organoid-only model (orange), the patient-only model (light blue) and the joint model (green). POOLED counts pairs within cohorts only.
(e) Hazard ratio (HR) per standard deviation of the score in held-out treated patients, with 95% confidence interval. HR below 1 means patients with a higher score lived longer without an event.
(f) Kaplan-Meier curves (fraction still event-free over time) of held-out treated patients split at the within-cohort median of the joint score; (g) the same for the organoid-only score. For liver (no survival times public) the panels show the fraction of responders in the high and low score halves in the sorafenib and placebo arms.
(h) Interaction test where a comparison group exists: ratio below 1 means the score separates outcomes more strongly in drug-treated patients than in the comparison group (untreated, placebo or another drug), which is what a drug-specific marker should do.
"""


def supp_fig(unit):
    u = UNIT[unit]
    ev, it, sc, meta = load(unit, "shortlist")
    try:
        ev_all, it_all, sc_all, meta_all = load(unit, "all")
    except FileNotFoundError:
        ev_all = None
    surv_c = [c for c in sc["unit_cohort"].unique()
              if ((sc["unit_cohort"] == c) & (sc["group"] == "treated") & sc["has_surv"]).sum() >= 10]
    nrow = 2 + len(surv_c)
    fig, ax = plt.subplots(nrow, 3, figsize=(11.5, 3.1 * nrow))
    ax = np.atleast_2d(ax)
    conc_bars(ax[0, 0], ev, "Network shortlist: held-out concordance")
    if ev_all is not None:
        conc_bars(ax[0, 1], ev_all, "All 673 pathways: held-out concordance")
    else:
        ax[0, 1].set_axis_off()
    int_forest(ax[0, 2], it, "Interaction (network shortlist)")
    hr_forest(ax[1, 0], ev, "HR per SD (network shortlist)")
    if ev_all is not None:
        hr_forest(ax[1, 1], ev_all, "HR per SD (all pathways)")
        int_forest(ax[1, 2], it_all, "Interaction (all pathways)")
    else:
        ax[1, 1].set_axis_off(); ax[1, 2].set_axis_off()
    for i, c in enumerate(surv_c):
        for j, m in enumerate(MODELS):
            km_split(ax[2 + i, j], sc, m, f"{c}: {m}", cohorts=[c])
    PS.letters(ax.ravel()[: 6 + 3 * len(surv_c)], x=-0.12)
    fig.suptitle(f"Supplementary: {unit} ({u['cancer']}, {u['drug']})", fontsize=11)
    fig.tight_layout()
    PS.save(fig, C.FIGS / "supp" / f"unit_{unit}", caption=f"""**Supplementary figure for {unit} ({u['cancer']}, {u['drug']}).**
Lab models: {meta['lab']}. Design: {meta['design']}.
(a) Held-out concordance per cohort using the {meta['n_features']} network-proximal pathways (C-index for survival, AUROC for response; 0.5 = chance) for organoid-only (orange), patient-only (light blue) and joint (green) models. (b) The same using all 673 pathways. (c) Treatment x score interaction (ratio < 1: the score matters more in drug-treated patients than in the comparison group), network shortlist.
(d, e) Hazard ratio per SD of the score in held-out treated patients with 95% CI, shortlist and all pathways. (f) Interaction, all pathways.
Remaining rows: Kaplan-Meier curves in each held-out survival cohort, treated patients split at the cohort median score, one column per model (shortlist features). Log-rank p in each title.
""")


def summary_tables():
    rows, irows = [], []
    for u in UNITS:
        for f in ("shortlist", "all"):
            try:
                ev, it, sc, meta = load(u["id"], f)
            except FileNotFoundError:
                continue
            ev = ev.assign(unit=u["id"], cancer=u["cancer"], drug=u["drug"], features=f, n_features=meta["n_features"], design=meta["design"])
            rows.append(ev)
            if len(it):
                irows.append(it.assign(unit=u["id"], cancer=u["cancer"], drug=u["drug"], features=f))
    S = pd.concat(rows, ignore_index=True)
    S.drop(columns=["num", "den"], errors="ignore").to_csv(C.TABLES / "summary_all.tsv", sep="\t", index=False)
    I = pd.concat(irows, ignore_index=True) if irows else pd.DataFrame()
    I.to_csv(C.TABLES / "interaction_all.tsv", sep="\t", index=False)
    M = S[(S["cohort"] == "POOLED") & S["model"].isin(MODELS + ["joint (w=0.5)"])]
    M = M[["cancer", "unit", "drug", "features", "n_features", "model", "outcome", "n_cohorts", "n", "events", "C", "HR", "HR_low", "HR_high", "p", "design"]]
    M.to_csv(C.TABLES / "summary_main.tsv", sep="\t", index=False)
    return S, I


def forest(S, I):
    M = S[(S["cohort"] == "POOLED") & S["model"].isin(MODELS) & (S["features"] == "shortlist")]
    surv = M[M["outcome"] == "survival"]
    units = [u["id"] for u in UNITS if u["id"] in set(M["unit"])]
    fig, ax = plt.subplots(1, 3, figsize=(13, 0.45 * len(units) + 2.5), gridspec_kw=dict(width_ratios=[1.2, 1, 1.2]))
    for j, m in enumerate(MODELS):
        for i, un in enumerate(units):
            r = surv[(surv["model"] == m) & (surv["unit"] == un)]
            if len(r):
                r = r.iloc[0]
                ax[0].errorbar(r["HR"], i + (j - 1) * 0.25, xerr=[[r["HR"] - r["HR_low"]], [r["HR_high"] - r["HR"]]], fmt="o",
                               color=PS.MODEL_COL[m], ms=4, lw=1, capsize=2, label=m if i == 0 else None)
            cc = M[(M["model"] == m) & (M["unit"] == un)]
            for k, (_, rr) in enumerate(cc.iterrows()):
                ax[1].plot(rr["C"], i + (j - 1) * 0.25, "o" if rr["outcome"] == "survival" else "^", color=PS.MODEL_COL[m], ms=5,
                           label=(f"{m} ({rr['outcome']})" if i == 0 or not ax[1].get_legend_handles_labels()[1].count(f"{m} ({rr['outcome']})") == 0 else None))
    ylab = [f"{un}\n({UNIT[un]['drug']})" for un in units]
    for a in ax[:2]:
        a.set_yticks(range(len(units))); a.set_yticklabels(ylab, fontsize=7); a.invert_yaxis()
    ax[0].axvline(1, color="k", ls="--", lw=0.8); ax[0].set_xscale("log")
    from matplotlib.ticker import NullFormatter
    ax[0].set_xticks([0.5, 0.75, 1, 1.5, 2, 3]); ax[0].set_xticklabels(["0.5", "0.75", "1", "1.5", "2", "3"])
    ax[0].xaxis.set_minor_formatter(NullFormatter())
    ax[0].set_xlabel("pooled HR per SD in held-out treated patients\n(<1 = higher score, longer survival)")
    ax[0].set_title("Survival in held-out treated patients", fontsize=9); ax[0].legend(fontsize=7, frameon=False)
    ax[1].axvline(0.5, color="k", ls="--", lw=0.8); ax[1].set_xlabel("pooled held-out concordance\n(circle = C-index, triangle = AUROC)")
    ax[1].set_title("Concordance", fontsize=9)
    h, l = ax[1].get_legend_handles_labels()
    keep = dict(zip(l, h))
    ax[1].legend(keep.values(), keep.keys(), fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    # interaction: joint model, every cohort with a comparison
    J = I[(I["features"] == "shortlist") & (I["model"] == "joint")].copy() if len(I) else pd.DataFrame()
    if len(J):
        J["lab"] = J["unit"] + " | " + J["cohort"] + " " + J["outcome"].str.replace("survival ", "").str.replace("response", "resp.")
        for i, (_, r) in enumerate(J.iterrows()):
            ax[2].errorbar(r["ratio"], i, xerr=[[r["ratio"] - r["low"]], [r["high"] - r["ratio"]]], fmt="s", color=PS.MODEL_COL["joint"], ms=4, capsize=2)
            ax[2].text(1.02, i, f"p={r['p']:.2g}", transform=ax[2].get_yaxis_transform(), fontsize=6.5, va="center")
        ax[2].set_yticks(range(len(J))); ax[2].set_yticklabels(J["lab"], fontsize=6.5); ax[2].invert_yaxis()
        ax[2].axvline(1, color="k", ls="--", lw=0.8); ax[2].set_xscale("log")
        ax[2].set_xlabel("treatment x score interaction (joint model)\nHR (survival) or OR (response)")
        ax[2].set_title("Drug-specific test (drug vs comparison)", fontsize=9)
    PS.letters(ax, x=-0.25)
    fig.tight_layout()
    PS.save(fig, C.FIGS / "forest_cross_cancer", caption="""**Cross-cancer summary (network-proximal pathways).**
(a) Pooled hazard ratio per SD of the score in held-out treated patients (Cox model stratified by cohort; 95% CI), for the organoid-only (orange), patient-only (light blue) and joint (green) models, one row per cancer x drug unit. HR below 1 means a higher score goes with longer event-free survival.
(b) Pooled held-out concordance (pairs compared only within a cohort): circles are C-index for survival, triangles are AUROC for response. 0.5 = chance.
(c) Treatment x score interaction of the joint model in every cohort that has a comparison group (untreated, placebo, or another drug). Below 1 means the score separates outcomes more in drug-treated patients than in the comparison group. p-values printed at the right.
""")


def main():
    PS.setup()
    S, I = summary_tables()
    forest(S, I)
    for c in MAIN_UNIT:
        main_fig(c)
    for u in UNITS:
        supp_fig(u["id"])


if __name__ == "__main__":
    main()
