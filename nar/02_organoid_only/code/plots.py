"""All figures, captions and tables from the pickles written by run_all.py.

Per set x cohort x drug: a 12-panel supplementary figure and a Kaplan-Meier grid (1..10 pathways, treated and
untreated). Per set: a composite main figure and a proximity heatmap. Across sets: summary table + forest plot.
"""
import pickle
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from common import ALT_CUTOFFS, DATA, PAPER_CUTOFF, ROOT, SETS  # noqa: E402
from engine import MLS, km_curve  # noqa: E402

RES = DATA / "results"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
# Okabe-Ito colour-blind-safe palette
BLUE, ORANGE, GREEN, VERM, SKY, PURPLE, GREY, YELLOW = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
                                                        "#CC79A7", "#999999", "#F0E442")
MLC = {"Ridge": BLUE, "SVR": ORANGE, "LinearRegression": GREEN}
COHORT_LABEL = {"COAD_authors": "COAD (authors' 2019 files)", "COAD": "COAD (current GDC)", "READ": "READ (current GDC)",
                "COADREAD": "COAD+READ pooled (current GDC)", "PAAD": "PAAD (current GDC)", "LIHC": "LIHC (current GDC)"}


def short(p, n=48):
    p = str(p).replace("REACTOME_", "").replace("_", " ").lower()
    p = p[0].upper() + p[1:] if p else p
    return p if len(p) <= n else p[: n - 1] + "…"


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2g}" if p < 0.001 else f"{p:.3f}"


def hr_ok(h, lo, hi):
    """Cox fits with (near) zero deaths in one group do not converge; such HRs are reported as not estimable."""
    return bool(np.isfinite(h) and np.isfinite(lo) and np.isfinite(hi) and hi / max(lo, 1e-300) < 1e4)


def fmt_hr(h, lo, hi):
    if not hr_ok(h, lo, hi):
        return "HR not estimable (0 deaths in one group)"
    return f"HR {h:.2f} (95% CI {lo:.2f}-{hi:.2f})"


def save(fig, path, caption):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    path.with_suffix(".caption.md").write_text(caption)


def letter(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom", ha="left")


def load_all():
    out = {}
    for f in sorted(RES.glob("*.pkl")):
        with open(f, "rb") as fh:
            out[f.stem] = pickle.load(fh)
    return out


def rank1(r, who="treated", ml="Ridge"):
    e = r["ev"]
    return e[(e.model == ml) & (e.k == 1)].iloc[0]


# --------------------------------------------------------------------------- panels
def p_prox(ax, r, top=25):
    z = r["z"].sort_values()
    zz = z.iloc[:max(top, 1)]
    cols = [BLUE if v <= PAPER_CUTOFF else GREY for v in zz.values]
    ax.barh(range(len(zz)), zz.values, color=cols)
    ax.set_yticks(range(len(zz)))
    ax.set_yticklabels([short(p, 55) for p in zz.index])
    ax.invert_yaxis()
    ax.axvline(PAPER_CUTOFF, color="k", ls="--", lw=1)
    ax.set_xlabel("Network proximity z-score (more negative = closer to drug targets)")
    m = r["meta"]
    ax.set_title(f"{len(r['feats'])} of {len(z)} pathways proximal (z <= {PAPER_CUTOFF}); targets of "
                 f"{m['proximity_column'].lower()}; {min(top, len(z))} closest shown", loc="left")


def p_nes(ax, r, fig, cap=30):
    nes = r["nes"]
    coef = r["coefs"]["Ridge"].abs().sort_values(ascending=False)
    rows = list(coef.index[:cap])
    y = r["y"].sort_values()
    M = nes.loc[rows, y.index]
    M = M.sub(M.mean(1), axis=0).div(M.std(1).replace(0, 1), axis=0)
    im = ax.imshow(M.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([short(p, 45) for p in rows])
    if len(y) <= 40:
        ax.set_xticks(range(len(y)))
        ax.set_xticklabels(y.index, rotation=90)
    else:
        ax.set_xticks([])
    ax.set_xlabel(f"Lab models, ordered from most sensitive (left) to most resistant (right); n = {len(y)}")
    ax.set_title(f"Pathway activity (ssGSEA score, z-scored per row) of the proximal pathways"
                 f"{' (top ' + str(cap) + ' by |Ridge coefficient|)' if len(coef) > cap else ''}", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("z-score")
    # response strip
    ax2 = ax.inset_axes([0, 1.01, 1, 0.04])
    ax2.imshow(y.values[None, :], aspect="auto", cmap="viridis")
    ax2.set_xticks([]), ax2.set_yticks([0]), ax2.set_yticklabels([r["meta"]["metric"]])


def p_loo(ax, r, ml):
    pred, rho, p = r["loo"][ml]
    y = r["y"]
    ax.scatter(y.values, pred.loc[y.index].values, s=18, color=MLC[ml], edgecolor="k", lw=0.3)
    lo, hi = np.nanmin([y.min(), pred.min()]), np.nanmax([y.max(), pred.max()])
    ax.plot([lo, hi], [lo, hi], color=GREY, ls=":", lw=1)
    ax.set_xlabel(f"Observed {r['meta']['metric']}")
    ax.set_ylabel("Leave-one-out predicted")
    ax.set_title(f"{ml}: Spearman rho = {rho:.2f}, p = {fmt_p(p)} (n = {len(y)})", loc="left")


def p_km(ax, T, E, g, title, short_title=False):
    st = None
    for grp, col, lab in ((0, BLUE, "Predicted responder"), (1, VERM, "Predicted non-responder")):
        m = g == grp
        t, s, lo, hi = km_curve(T[m], E[m])
        ax.step(t, s, where="post", color=col, lw=1.6, label=f"{lab} (n={m.sum()}, deaths={E[m].sum()})")
        ax.fill_between(t, lo, hi, step="post", color=col, alpha=0.12, lw=0)
        cen = T[m][E[m] == 0]
        if len(cen):
            idx = np.searchsorted(t, cen, side="right") - 1
            ax.plot(cen, s[idx], "|", color=col, ms=5, mew=0.8)
    from engine import surv_stats
    st = surv_stats(T, E, g)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Fraction alive")
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.98, 0.98, f"log-rank p = {fmt_p(st['p'])}\n{fmt_hr(st['HR'], st['HR_lo'], st['HR_hi'])}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.set_title(title, loc="left")
    return st


def p_k(ax, r):
    e = r["ev"]
    for ml in MLS:
        s = e[e.model == ml]
        ax.plot(s.k, -np.log10(s.treated_p), "-o", color=MLC[ml], ms=4, label=f"{ml}, treated")
        if "untreated_p" in s:
            ax.plot(s.k, -np.log10(s.untreated_p), "--", marker="s", mfc="white", color=MLC[ml], ms=3.5,
                    label=f"{ml}, untreated")
    ax.axhline(-np.log10(0.05), color="k", ls=":", lw=1)
    ax.text(1, -np.log10(0.05), " p = 0.05", va="bottom", fontsize=8)
    ax.set_xlabel("Number of top pathways in the model (1 = single-pathway biomarker)")
    ax.set_ylabel("-log10 log-rank p")
    ax.set_xticks(range(1, int(e.k.max()) + 1))
    ax.legend(frameon=False, ncol=2, fontsize=7.5 if len(MLS) > 2 else 8)
    ax.set_title("Single- and multi-pathway models: p vs number of pathways", loc="left")


def p_null(ax, null, obs, emp, title, what):
    v = -np.log10(np.clip(null, 1e-300, 1))
    ax.hist(v, bins=40, color=SKY, edgecolor="white")
    ax.axvline(-np.log10(obs), color=VERM, lw=2, label=f"real result (p = {fmt_p(obs)})")
    ax.set_xlabel("-log10 log-rank p")
    ax.set_ylabel(f"Count ({len(null):,} {what})")
    ax.legend(frameon=False, loc="upper right")
    ax.text(0.98, 0.80, f"share as good or better:\n{fmt_p(emp)}", transform=ax.transAxes, ha="right", va="top")
    ax.set_title(title, loc="left")


def p_rank(ax, r):
    m = len(r["feats"])
    br = r["boot_rank"]
    vals = np.arange(1, m + 2)
    cnt = np.bincount(br, minlength=m + 2)[1:m + 2]
    ax.bar(vals, cnt, color=SKY, width=0.9)
    o = r["obs_rank"]
    ax.axvline(o, color=VERM, lw=2, label=f"real proximal pathways: {'rank ' + str(o) if o <= m else 'none found'}")
    ax.set_yscale("log")
    ax.set_xlabel(f"First rank meeting the criteria ({m + 1} = none found)")
    ax.set_ylabel(f"Random sets (of {len(br):,}; log scale)")
    crit = "survival p<0.05 & isogenic p<0.05" if r.get("iso_tab") is not None else "survival p<0.05"
    ax.set_title(f"Paper's Fig. 4e statistic ({crit})\nshare as good or better: {fmt_p(r['boot_rank_emp_p'])}", loc="left")
    ax.legend(frameon=False, loc="upper center")


def p_cut(ax, r):
    c = r["cut"].copy()
    x = np.arange(len(c))
    ax.bar(x - 0.2, -np.log10(c["treated_p"].fillna(1)), 0.4, color=BLUE, label="treated")
    ax.bar(x + 0.2, -np.log10(c["untreated_p"].fillna(1)), 0.4, color=GREY, label="untreated")
    ax.axhline(-np.log10(0.05), color="k", ls=":", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:g}\n({int(n)} pw)" + ("\npaper" if v == PAPER_CUTOFF else "")
                        for v, n in zip(c.cutoff, c.n_pathways)])
    ax.set_xlabel("Proximity z cut-off (number of proximal pathways)")
    ax.set_ylabel("-log10 log-rank p (rank-1, Ridge)")
    ax.legend(frameon=False)
    ax.set_title("Sensitivity to the proximity cut-off", loc="left")


def title_of(r):
    m = r["meta"]
    d = m["model_drug"] if m["model_drug"] == m["patient_drug"] else f"{m['model_drug']} (models) / {m['patient_drug']} (patients)"
    return f"{m['set_label']} | {COHORT_LABEL.get(m['cohort'], m['cohort'])} | {d}"


# --------------------------------------------------------------------------- per combo
def supp_figure(tag, r, outdir):
    m = r["meta"]
    fig = plt.figure(figsize=(17, 32))
    gs = GridSpec(6, 3, figure=fig, hspace=0.55, wspace=0.55, height_ratios=[1.2, 1.3, 1, 1, 1, 0.8], top=0.975, bottom=0.02)
    axes = {}
    axes["a"] = fig.add_subplot(gs[0, :2]); p_prox(axes["a"], r)
    axes["b"] = fig.add_subplot(gs[0, 2]); p_cut(axes["b"], r)
    axes["c"] = fig.add_subplot(gs[1, :2]); p_nes(axes["c"], r, fig)
    for L, ml, pos in (("d", "Ridge", gs[1, 2]), ("e", "SVR", gs[2, 0]), ("f", "LinearRegression", gs[2, 1])):
        axes[L] = fig.add_subplot(pos); p_loo(axes[L], r, ml)
    axes["g"] = fig.add_subplot(gs[2, 2]); p_k(axes["g"], r)
    T, E, g, used = r["km"][("treated", 1)]
    axes["h"] = fig.add_subplot(gs[3, 0])
    st = p_km(axes["h"], T, E, g, "Drug-treated patients, rank-1 pathway\n" + textwrap.fill(short(used[0], 80), 50))
    axes["i"] = fig.add_subplot(gs[3, 1])
    su = None
    if ("untreated", 1) in r["km"] and len(r["km"][("untreated", 1)][0]):
        Tu, Eu, gu, _ = r["km"][("untreated", 1)]
        su = p_km(axes["i"], Tu, Eu, gu, "Untreated patients (no drug record), same pathway\n" +
                  textwrap.fill(short(used[0], 80), 50))
    obs = rank1(r)["treated_p"]
    axes["j"] = fig.add_subplot(gs[3, 2])
    p_null(axes["j"], r["boot"], obs, r["boot_emp_p"], "Bootstrap (paper Fig. 4e/f): random pathway sets\n"
           f"of the same size ({len(r['feats'])}), Ridge, rank-1", "random sets")
    axes["k"] = fig.add_subplot(gs[4, 0])
    p_null(axes["k"], r["perm_rank1"], obs, r["perm_rank1_emp_p"], "Shuffled lab-model response, rank-1 Ridge", "shuffles")
    axes["l"] = fig.add_subplot(gs[4, 1])
    best = r["ev"]["treated_p"].min()
    p_null(axes["l"], r["perm_best"], best, r["perm_best_emp_p"], "Shuffled response, best of all 3 models x\n"
           "1-10 pathways (the 'pick the best' search)", "shuffles")
    axn = fig.add_subplot(gs[4, 2]); p_rank(axn, r)
    axl = fig.add_subplot(gs[5, :]); axl.axis("off")
    txt = (f"Lab models: {m['n_models']}\nProximal pathways: {m['n_proximal']} of {m['n_pathways_total']}\n"
           f"Treated patients: {m['n_treated']} ({m['treated_events']} deaths)\n"
           f"Untreated patients: {m['n_untreated']}\n\nRank-1 Ridge biomarker:\n{textwrap.fill(short(used[0], 90), 40)}\n\n"
           f"Treated: p = {fmt_p(st['p'])}, {fmt_hr(st['HR'], st['HR_lo'], st['HR_hi'])}\n"
           + (f"Untreated: p = {fmt_p(su['p'])}, {fmt_hr(su['HR'], su['HR_lo'], su['HR_hi'])}\n" if su else "")
           + f"Bootstrap share: {fmt_p(r['boot_emp_p'])}\nShuffle share (rank-1): {fmt_p(r['perm_rank1_emp_p'])}\n"
           f"Best of 30 tests: p = {fmt_p(best)}; shuffle share {fmt_p(r['perm_best_emp_p'])}")
    axl.text(0, 1, txt, va="top", fontsize=9, family="DejaVu Sans")
    for L, ax in axes.items():
        letter(ax, L)
    letter(axn, "m")
    letter(axl, "n")
    fig.suptitle(title_of(r), fontsize=12, y=0.995)
    cap = caption_supp(r, st, su, used[0], best)
    extra = "_SuppFig11" if m["patient_drug"] == "CETUXIMAB" else ""
    save(fig, outdir / f"Fig2a_Fig2d_Fig4e_SuppFig17_panels_{m['cohort']}_{m['model_drug']}{extra}", cap)
    return st, su


def caption_supp(r, st, su, pw, best):
    m = r["meta"]
    return f"""# Supplementary figure: {title_of(r)}

All panels follow the organoid-only strategy of Kong et al. 2020 (network proximity -> pathway scores ->
machine learning on lab models -> median split of drug-treated patients -> survival test).

- **a** Network proximity. Each bar is one Reactome pathway; the bar length is how close the pathway sits to the
  drug's protein targets in the STRING network (z-score; more negative = closer). Blue bars pass the paper's
  cut-off (z <= {PAPER_CUTOFF}); {m['n_proximal']} of {m['n_pathways_total']} pathways pass. Proximity values are the
  authors' published table for {m['proximity_column'].lower()}{' (oxaliplatin has no protein target, so the platinum drug cisplatin is used as its proxy)' if m['patient_drug'] == 'OXALIPLATIN' else ''}.
- **b** Same analysis with other cut-offs (the user's -1.182405 and -2.1925376, and -1.6449 = 5%). Bars show
  -log10 p of the rank-1 pathway in treated (blue) and untreated (grey) patients; the dotted line is p = 0.05.
- **c** Pathway activity (ssGSEA scores) of the proximal pathways in the {m['n_models']} lab models, ordered by drug
  response ({m['metric']}; lower = more sensitive). Red = higher than average for that pathway, blue = lower.
- **d-f** Leave-one-out test inside the lab models: each model is predicted by a model trained on all the others
  (Ridge, SVR, LinearRegression, all proximal pathways). Spearman rho measures agreement (1 = perfect order).
- **g** p-value (as -log10 p; higher = smaller p) for 1 to 10 top pathways for each learner, in treated (solid)
  and untreated (dashed) patients.
- **h** Kaplan-Meier curves for drug-treated patients split at the median predicted response using the rank-1
  Ridge pathway ({short(pw, 90)}). Log-rank p = {fmt_p(st['p'])}; n = {st['n']}; deaths = {st['events']}
  ({st['events_resp']} in predicted responders, {st['events_nonresp']} in predicted non-responders);
  {fmt_hr(st['HR'], st['HR_lo'], st['HR_hi'])} (non-responders vs responders; > 1 is the expected direction).
- **i** The same pathway and split rule in patients of the same cancer with no drug record ("untreated"), to see
  whether the split is drug-specific. {('Log-rank p = ' + fmt_p(su['p']) + '; n = ' + str(su['n']) + '; deaths = ' + str(su['events']) + '; ' + fmt_hr(su['HR'], su['HR_lo'], su['HR_hi']) + '.') if su else 'No untreated patients.'}
- **j** The paper's bootstrap (Fig. 4e/f): {len(r['boot']):,} random pathway sets of the same size, each run through
  Ridge -> rank-1 -> log-rank. Share of random sets with p as small or smaller than the real one = {fmt_p(r['boot_emp_p'])}.
- **k** The lab-model drug response shuffled {len(r['perm_rank1']):,} times and the whole rank-1 analysis repeated.
  Share of shuffles as good or better = {fmt_p(r['perm_rank1_emp_p'])}.
- **l** The same shuffles, but keeping the best of all 30 tests (3 learners x 1-10 pathways), which mirrors picking
  the best result. Real best p = {fmt_p(best)}; share of shuffles as good or better = {fmt_p(r['perm_best_emp_p'])}.
- **m** The paper's own bootstrap statistic (Fig. 4e): for each random set, the first rank (by |coefficient|) whose
  pathway alone gives p < 0.05 in treated patients{' and a t-test p < 0.05 between 5-FU-sensitive and -resistant isogenic HCT-8 cells' if r.get('iso_tab') is not None else ''}.
  Real proximal pathways: first such rank = {r['obs_rank'] if r['obs_rank'] <= len(r['feats']) else 'none found'};
  share of random sets reaching that rank or better = {fmt_p(r['boot_rank_emp_p'])} (paper, colon 5-FU: 0.0012).
- **n** Key numbers. The rank-1 pathway's own proximity z = {r['z'].get(pw, float('nan')):.4f} (cut-off {PAPER_CUTOFF}).
- Panels d-f stand in for the paper's Supp. Fig. 17 (a 60/10/30 split of the organoids); leave-one-out uses every
  model once as test data, which is more stable at these sample sizes.

Data: {m['data']}. Lab-model response metric: {m['metric']}.
"""


def km_grid(tag, r, outdir):
    m = r["meta"]
    ks = sorted(k for (w, k) in r["km"] if w == "treated")
    hasu = any(w == "untreated" for (w, k) in r["km"]) and len(r["km"].get(("untreated", 1), [[]])[0])
    nrow = 4 if hasu else 2
    fig, axs = plt.subplots(nrow, 5, figsize=(22, 4.6 * nrow))
    axs = np.atleast_2d(axs)
    rows = []
    for wi, who in enumerate(["treated", "untreated"] if hasu else ["treated"]):
        for j, k in enumerate(range(1, 11)):
            ax = axs[wi * 2 + j // 5, j % 5]
            if k not in ks:
                ax.axis("off")
                continue
            T, E, g, used = r["km"][(who, k)]
            st = p_km(ax, T, E, g, f"{who.capitalize()}, top {k} pathway{'s' if k > 1 else ''} (Ridge)")
            ax.legend(loc="lower left", frameon=False, fontsize=7.5)
            rows.append((who, k, st))
            letter(ax, "abcdefghijklmnopqrst"[wi * 10 + j])
    fig.suptitle(title_of(r) + " - Kaplan-Meier curves for 1-10 pathways", fontsize=12)
    fig.tight_layout()
    lines = "\n".join(f"- {who}, {k} pathway(s): p = {fmt_p(s['p'])}, n = {s['n']}, deaths = {s['events']}, "
                      f"{fmt_hr(s['HR'], s['HR_lo'], s['HR_hi'])}" for who, k, s in rows)
    cap = f"""# Kaplan-Meier grid: {title_of(r)}

Patients are split at the median of the predicted drug response from a Ridge model trained on the lab models,
using the top 1 to 10 pathways (ranked by absolute coefficient) among the proximal pathways
(paper's multi-pathway analysis). Blue = predicted responders, vermillion = predicted non-responders; shaded bands
are 95% confidence intervals; ticks are censored patients. {'Top two rows: drug-treated patients. Bottom two rows: untreated patients of the same cancer (same models, same rule).' if hasu else 'Drug-treated patients.'}
SVR and LinearRegression results for every k are in the table `all_tests_{m['cohort']}_{m['model_drug']}.tsv`.

{lines}
"""
    save(fig, outdir / f"SuppFig5_SourceData-multipathway_KMgrid_{m['cohort']}_{m['model_drug']}", cap)


def coef_figure(tag, r, outdir, cap=40):
    """Paper Fig. 2b (Ridge coefficients), Supp. Fig. 2 (LinearRegression, SVR) and Supp. Fig. 3 (agreement)."""
    m = r["meta"]
    C = pd.DataFrame(r["coefs"])
    order = C["Ridge"].abs().sort_values(ascending=False).index
    shown = list(order[:cap])
    fig = plt.figure(figsize=(20, 0.24 * len(shown) + 7))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.24 * len(shown) + 1, 4.2], hspace=0.25, wspace=0.9)
    for j, ml in enumerate(MLS):
        ax = fig.add_subplot(gs[0, j])
        v = C.loc[shown, ml]
        top = C[ml].abs().idxmax()
        ax.barh(range(len(v)), v.abs().values, color=[ORANGE if p == top else MLC[ml] for p in v.index])
        ax.set_yticks(range(len(v)))
        ax.set_yticklabels([short(p, 45) + (" (-)" if c < 0 else " (+)") for p, c in zip(v.index, v.values)])
        ax.invert_yaxis()
        ax.set_xlabel("|coefficient| (predictive ability)")
        ax.set_title(f"{ml}: rank-1 = {short(top, 40)}", loc="left")
        letter(ax, "abc"[j])
    pairs = [("Ridge", "LinearRegression"), ("Ridge", "SVR"), ("LinearRegression", "SVR")]
    for j, (a, b) in enumerate(pairs):
        ax = fig.add_subplot(gs[1, j])
        ax.scatter(C[a], C[b], s=14, color=PURPLE)
        r2 = np.corrcoef(C[a], C[b])[0, 1] ** 2 if len(C) > 2 else np.nan
        ax.set_xlabel(f"{a} coefficient"); ax.set_ylabel(f"{b} coefficient")
        ax.set_title(f"R2 = {r2:.3f}", loc="left")
        letter(ax, "def"[j])
    fig.suptitle(title_of(r) + " - organoid model coefficients", fontsize=12)
    top = {ml: C[ml].abs().idxmax() for ml in MLS}
    save(fig, outdir / f"Fig2b_SuppFig2_SuppFig3_coefficients_{m['cohort']}_{m['model_drug']}", f"""# Model coefficients: {title_of(r)}

(a-c) Absolute regression coefficients ("predictive ability" in the paper) of every proximal pathway
({len(C)}{', top ' + str(cap) + ' shown' if len(C) > cap else ''}) for Ridge (paper Fig. 2b), LinearRegression and linear SVR (paper Supp. Fig. 2),
trained on {m['n_models']} lab models ({m['metric']}). (+) = higher pathway score goes with higher (more
resistant) response; (-) = with lower (more sensitive). Orange bar = rank-1 pathway. Rank-1: Ridge
{short(top['Ridge'], 80)}; SVR {short(top['SVR'], 80)}; LinearRegression {short(top['LinearRegression'], 80)}.
(d-f) Agreement of the coefficients between learners (paper Supp. Fig. 3; R2 = squared Pearson correlation).
""")


def tables(tag, r, tdir):
    m = r["meta"]
    s = f"{m['cohort']}_{m['model_drug']}"
    tdir.mkdir(parents=True, exist_ok=True)
    r["ev"].to_csv(tdir / f"all_tests_{s}.tsv", sep="\t", index=False)
    pd.DataFrame(r["coefs"]).assign(z=r["z"].reindex(r["feats"])).sort_values("Ridge", key=abs, ascending=False) \
        .to_csv(tdir / f"coefficients_{s}.tsv", sep="\t")
    pd.DataFrame({"observed": r["y"], **{f"loo_{ml}": r["loo"][ml][0] for ml in MLS}}).to_csv(tdir / f"loo_{s}.tsv", sep="\t")
    r["cut"].to_csv(tdir / f"cutoff_sensitivity_{s}.tsv", sep="\t", index=False)
    r["z"].sort_values().to_frame("proximity_z").assign(proximal=lambda d: d.proximity_z <= PAPER_CUTOFF) \
        .to_csv(tdir / f"proximity_{s}.tsv", sep="\t")
    pd.DataFrame({"bootstrap_random_pathways_p": r["boot"]}).describe().to_csv(tdir / f"bootstrap_summary_{s}.tsv", sep="\t")
    np.savetxt(tdir / f"bootstrap_p_{s}.txt.gz", r["boot"])
    np.savetxt(tdir / f"permutation_p_{s}.txt.gz", np.c_[r["perm_rank1"], r["perm_best"]], header="rank1_ridge_p\tbest_of_30_p")


def summary_row(r):
    m = r["meta"]
    a = rank1(r)
    e = r["ev"]
    row = {"set": m["set_label"], "folder": m["folder"], "cohort": m["cohort"], "model_drug": m["model_drug"],
           "patient_drug": m["patient_drug"], "proximity_from": m["proximity_column"], "metric": m["metric"],
           "n_models": m["n_models"], "n_proximal": m["n_proximal"], "rank1_pathway": a["pathways"],
           "treated_n": a["treated_n"], "treated_events": a["treated_events"],
           "treated_events_resp": a["treated_events_resp"], "treated_events_nonresp": a["treated_events_nonresp"],
           "treated_HR": a["treated_HR"], "treated_HR_lo": a["treated_HR_lo"], "treated_HR_hi": a["treated_HR_hi"],
           "treated_p": a["treated_p"]}
    for kk in ("n", "events", "HR", "HR_lo", "HR_hi", "p"):
        row[f"untreated_{kk}"] = a.get(f"untreated_{kk}", np.nan)
    for ml in ("SVR", "LinearRegression"):
        b = rank1(r, ml=ml)
        row[f"rank1_pathway_{ml}"] = b["pathways"]
        row[f"treated_p_{ml}"] = b["treated_p"]
    bi = e["treated_p"].idxmin()
    row.update(best_p_30_tests=e.loc[bi, "treated_p"], best_setting=f"{e.loc[bi, 'model']} k={e.loc[bi, 'k']}",
               tests_p_below_0_05=int((e.treated_p < 0.05).sum()),
               bootstrap_share=r["boot_emp_p"], shuffle_share_rank1=r["perm_rank1_emp_p"],
               shuffle_share_best=r["perm_best_emp_p"],
               **{f"loo_rho_{ml}": r["loo"][ml][1] for ml in MLS}, **{f"loo_p_{ml}": r["loo"][ml][2] for ml in MLS},
               borderline_n=bool(m.get("borderline", False)))
    return row


# --------------------------------------------------------------------------- per set
def main_figure(folder, items, outdir):
    n = len(items)
    fig = plt.figure(figsize=(30, 5.0 * n))
    gs = GridSpec(n, 6, figure=fig, hspace=0.75, wspace=0.6, width_ratios=[1.5, 1, 1.15, 1.15, 1.1, 1])
    L = iter("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    lines = []
    for i, (tag, r) in enumerate(items):
        m = r["meta"]
        ax = fig.add_subplot(gs[i, 0]); p_prox(ax, r, top=12); ax.set_title(
            f"{m['model_drug']} | {COHORT_LABEL.get(m['cohort'], m['cohort'])}\n{len(r['feats'])} proximal pathways "
            f"(z <= {PAPER_CUTOFF}), 12 closest", loc="left"); letter(ax, next(L))
        ax = fig.add_subplot(gs[i, 1]); p_loo(ax, r, "Ridge"); ax.set_title(
            f"Ridge leave-one-out\nrho = {r['loo']['Ridge'][1]:.2f}, p = {fmt_p(r['loo']['Ridge'][2])}", loc="left"); letter(ax, next(L))
        T, E, g, used = r["km"][("treated", 1)]
        ax = fig.add_subplot(gs[i, 2]); st = p_km(ax, T, E, g, "Treated, rank-1:\n" + textwrap.fill(short(used[0], 70), 42)); letter(ax, next(L))
        ax = fig.add_subplot(gs[i, 3])
        su = None
        if len(r["km"].get(("untreated", 1), [[]])[0]):
            Tu, Eu, gu, _ = r["km"][("untreated", 1)]
            su = p_km(ax, Tu, Eu, gu, "Untreated, same pathway")
        letter(ax, next(L))
        ax = fig.add_subplot(gs[i, 4]); p_k(ax, r); ax.set_title("p vs number of pathways", loc="left")
        ax.legend(frameon=False, fontsize=7, ncol=1); letter(ax, next(L))
        ax = fig.add_subplot(gs[i, 5]); p_null(ax, r["perm_rank1"], rank1(r)["treated_p"], r["perm_rank1_emp_p"],
                                               "Shuffled response (1,000x)", "shuffles"); letter(ax, next(L))
        lines.append(f"- Row {i + 1}, {m['model_drug']} in {COHORT_LABEL.get(m['cohort'], m['cohort'])}: "
                     f"{m['n_models']} models, {m['n_proximal']} proximal pathways; rank-1 pathway "
                     f"'{short(used[0], 90)}'; treated n = {st['n']}, deaths = {st['events']}, p = {fmt_p(st['p'])}, "
                     f"{fmt_hr(st['HR'], st['HR_lo'], st['HR_hi'])}"
                     + (f"; untreated n = {su['n']}, deaths = {su['events']}, p = {fmt_p(su['p'])}, "
                        f"{fmt_hr(su['HR'], su['HR_lo'], su['HR_hi'])}" if su else "")
                     + f"; Ridge leave-one-out rho = {r['loo']['Ridge'][1]:.2f}; shuffle share = {fmt_p(r['perm_rank1_emp_p'])}.")
    fig.suptitle(f"{items[0][1]['meta']['set_label']}: the paper's organoid-only strategy, drug by drug", fontsize=14, y=1.0)
    cap = f"""# Main figure: {items[0][1]['meta']['set_label']}

Each row is one drug (and patient cohort). Columns, left to right: (1) the 12 Reactome pathways closest to the
drug's targets in the protein network (blue = passes the paper's cut-off z <= {PAPER_CUTOFF}); (2) how well a Ridge
model trained on the lab models predicts a left-out model (Spearman rho); (3) survival of drug-treated TCGA patients
split at the median prediction of the rank-1 pathway (blue = predicted responders, vermillion = predicted
non-responders); (4) the same split in untreated patients of the same cancer; (5) p-value for 1-10 pathways and the
three learners (solid = treated, dashed = untreated; dotted line p = 0.05); (6) where the real p-value sits among
1,000 analyses with the lab-model response shuffled. HR = hazard ratio of predicted non-responders vs responders
(> 1 is the expected direction).

{chr(10).join(lines)}
"""
    save(fig, outdir / f"Fig2d_main_composite_{folder}", cap)


def prox_heatmap(folder, items, outdir, cap_rows=60):
    zs = {}
    for tag, r in items:
        d = r["meta"]["patient_drug"]
        if d not in zs:
            zs[d] = r["z"]
    Z = pd.DataFrame(zs)
    prox = Z[(Z <= PAPER_CUTOFF).any(axis=1)]
    order = prox.min(axis=1).sort_values().index
    shown = order[:cap_rows]
    M = prox.loc[shown]
    fig, ax = plt.subplots(figsize=(3 + 1.1 * M.shape[1], 1.5 + 0.2 * len(M)))
    im = ax.imshow(M.clip(-5, 2).values, aspect="auto", cmap="PuOr", vmin=-5, vmax=5)
    for (i, j), v in np.ndenumerate(M.values):
        if np.isfinite(v) and v <= PAPER_CUTOFF:
            ax.text(j, i, "•", ha="center", va="center", fontsize=8, color="k")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([c.lower() + (" (cisplatin proxy)" if c == "OXALIPLATIN" else "") for c in M.columns], rotation=45, ha="right")
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels([short(p, 60) for p in M.index])
    cb = fig.colorbar(im, ax=ax, fraction=0.04)
    cb.set_label("proximity z (negative = close)")
    ax.set_title(f"Pathways proximal (z <= {PAPER_CUTOFF}, dot) to at least one drug; "
                 f"{len(shown)} of {len(prox)} shown", loc="left")
    fig.tight_layout()
    (outdir.parent / "tables").mkdir(parents=True, exist_ok=True)
    Z.to_csv(outdir.parent / "tables" / "proximity_all_drugs.tsv", sep="\t")
    cap = f"""# Proximity heatmap: {items[0][1]['meta']['set_label']}

Rows: Reactome pathways that pass the paper's proximity cut-off (z <= {PAPER_CUTOFF}) for at least one of the drugs
tested in this set ({len(prox)} pathways; the {len(shown)} closest are shown). Columns: drugs. Colour: network
proximity z-score between the drug's protein targets and the pathway's genes in the STRING network (authors'
published table; purple/negative = closer than random). A dot marks a proximal pathway. Oxaliplatin has no protein
target, so cisplatin (same platinum class) is used. Proximal pathways per drug: {', '.join(f'{c.lower()} {int((Z[c] <= PAPER_CUTOFF).sum())}' for c in Z.columns)}.
Full table: `tables/proximity_all_drugs.tsv`.
"""
    save(fig, outdir / "supp_proximity_heatmap", cap)


# --------------------------------------------------------------------------- cross set
def forest(S, path):
    S = S.reset_index(drop=True)
    n = len(S)
    fig = plt.figure(figsize=(16, 0.42 * n + 2))
    gs = GridSpec(1, 2, width_ratios=[1.1, 1.6], wspace=0.02)
    ax = fig.add_subplot(gs[0, 1])
    axt = fig.add_subplot(gs[0, 0], sharey=ax)
    y = np.arange(n)[::-1]
    for i, r in S.iterrows():
        yy = y[i]
        if hr_ok(r.treated_HR, r.treated_HR_lo, r.treated_HR_hi):
            ax.errorbar(r.treated_HR, yy + 0.12, xerr=[[r.treated_HR - r.treated_HR_lo], [r.treated_HR_hi - r.treated_HR]],
                        fmt="o", color=BLUE, ms=5, capsize=2, lw=1.2)
        else:
            ax.text(1, yy + 0.12, "treated HR not estimable (0 deaths in one group)", fontsize=7, color=BLUE, va="center", ha="center")
        if hr_ok(r.untreated_HR, r.untreated_HR_lo, r.untreated_HR_hi):
            ax.errorbar(r.untreated_HR, yy - 0.18, xerr=[[r.untreated_HR - r.untreated_HR_lo], [r.untreated_HR_hi - r.untreated_HR]],
                        fmt="s", mfc="white", color=GREY, ms=4, capsize=2, lw=1)
    ax.axvline(1, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 30)
    ax.set_xlabel("Hazard ratio, predicted non-responders vs responders (log scale; > 1 = expected direction)")
    ax.plot([], [], "o", color=BLUE, label="drug-treated patients")
    ax.plot([], [], "s", mfc="white", color=GREY, label="untreated patients (same biomarker)")
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 1.0 + 1.2 / n))
    plt.setp(ax.get_yticklabels(), visible=False)
    axt.set_xlim(0, 1)
    axt.axis("off")
    hdr = ["Set | cohort | drug", "n (deaths)", "p treated", "p untr.", "shuffle"]
    xs = [0, 0.56, 0.71, 0.83, 0.94]
    for x, h in zip(xs, hdr):
        axt.text(x, n - 0.2, h, fontweight="bold", fontsize=8.5)
    for i, r in S.iterrows():
        yy = y[i]
        lab = f"{r.folder.split('_')[1]} | {r.cohort.replace('_authors', ' (authors)')} | {r.model_drug.lower()}"
        if r.borderline_n:
            lab += " *"
        vals = [lab, f"{int(r.treated_n)} ({int(r.treated_events)})", fmt_p(r.treated_p), fmt_p(r.untreated_p),
                fmt_p(r.shuffle_share_rank1)]
        for x, v in zip(xs, vals):
            axt.text(x, yy, v, va="center", fontsize=8.5)
    ax.set_ylim(-0.8, n - 0.2 + 0.6)
    fig.suptitle("All lab-model sets x drugs: rank-1 pathway (Ridge), median split, overall survival", fontsize=12)
    cap = f"""# Cross-set forest plot

One row per lab-model set x patient cohort x drug ({n} rows). Each row uses the paper's single-pathway biomarker
(rank-1 pathway of a Ridge model trained on the lab models) and splits TCGA patients at the median prediction.
Filled blue circles: drug-treated patients; open grey squares: untreated patients of the same cancer with the same
biomarker and rule. Lines: 95% confidence intervals of the hazard ratio (Cox model). A hazard ratio above 1 means
predicted non-responders died sooner, which is the direction the method expects. "Shuffle" is the share of 1,000
analyses with the lab-model response shuffled that gave a p-value as small or smaller. * = 12-14 treated patients
(below the ~15 target, kept and flagged). Numbers: `tables/cross_set_summary.tsv`.
"""
    save(fig, path, cap)


def main():
    allr = load_all()
    rows = []
    by_folder = {}
    for tag, r in allr.items():
        m = r["meta"]
        d = ROOT / m["folder"]
        supp_figure(tag, r, d / "figures")
        km_grid(tag, r, d / "figures")
        coef_figure(tag, r, d / "figures")
        tables(tag, r, d / "tables")
        rows.append(summary_row(r))
        by_folder.setdefault(m["folder"], []).append((tag, r))
        print("plotted", tag, flush=True)
    order = {"COAD_authors": 0, "COAD": 1, "READ": 2, "COADREAD": 3, "PAAD": 4, "LIHC": 5}
    for folder, items in by_folder.items():
        items = sorted(items, key=lambda t: (order[t[1]["meta"]["cohort"]], t[1]["meta"]["model_drug"]))
        main_items = [t for t in items if t[1]["meta"]["cohort"] in ("COAD_authors", "COAD", "PAAD", "LIHC")]
        main_figure(folder, main_items, ROOT / folder / "figures")
        extra = [t for t in items if t[1]["meta"]["cohort"] in ("READ", "COADREAD")]
        if extra:
            fig_dir = ROOT / folder / "figures"
            main_figure(folder, extra, fig_dir / "_tmp")
            for ext in ("png", "pdf", "caption.md"):
                (fig_dir / "_tmp" / f"Fig2d_main_composite_{folder}.{ext}").rename(fig_dir / f"Fig2d_supp_composite_READ_and_pooled.{ext}")
            (fig_dir / "_tmp").rmdir()
        prox_heatmap(folder, items, ROOT / folder / "figures")
    S = pd.DataFrame(rows)
    fo = {v[0]: i for i, v in enumerate(SETS.values())}
    S["_o"] = S.folder.map(fo) * 10 + S.cohort.map(order)
    S = S.sort_values(["_o", "model_drug"]).drop(columns="_o")
    (ROOT / "tables").mkdir(exist_ok=True)
    S.to_csv(ROOT / "tables" / "cross_set_summary.tsv", sep="\t", index=False)
    forest(S, ROOT / "figures" / "cross_set_forest")
    # cut-off sensitivity across sets
    cut = pd.concat([r["cut"].assign(folder=r["meta"]["folder"], cohort=r["meta"]["cohort"], drug=r["meta"]["model_drug"])
                     for r in allr.values()])
    cut.to_csv(ROOT / "tables" / "cross_set_cutoff_sensitivity.tsv", sep="\t", index=False)
    cutoff_overview(cut, ROOT / "figures" / "cross_set_cutoff_sensitivity")


def cutoff_overview(cut, path):
    cut = cut.copy()
    cut["row"] = cut.folder.str.split("_").str[1] + " | " + cut.cohort.str.replace("_authors", " (authors)") + " | " + cut.drug.str.lower()
    P = cut.pivot_table(index="row", columns="cutoff", values="treated_p", sort=False)
    N = cut.pivot_table(index="row", columns="cutoff", values="n_pathways", sort=False)
    fig, ax = plt.subplots(figsize=(8, 0.38 * len(P) + 1.5))
    im = ax.imshow(-np.log10(P.values.astype(float)), aspect="auto", cmap="viridis", vmin=0, vmax=3)
    for (i, j), v in np.ndenumerate(P.values):
        if np.isfinite(v):
            ax.text(j, i, f"{fmt_p(v)}\n({int(N.values[i, j])} pw)", ha="center", va="center", fontsize=7,
                    color="white" if -np.log10(v) < 1.8 else "k")
    ax.set_xticks(range(P.shape[1]))
    ax.set_xticklabels([f"{c:g}" + (" (paper)" if c == PAPER_CUTOFF else "") for c in P.columns])
    ax.set_yticks(range(len(P)))
    ax.set_yticklabels(P.index)
    ax.set_xlabel("Proximity z cut-off")
    fig.colorbar(im, ax=ax, fraction=0.04).set_label("-log10 p (treated, rank-1 Ridge)")
    ax.set_title("Cut-off sensitivity: treated-patient p and number of proximal pathways", loc="left")
    fig.tight_layout()
    save(fig, path, f"""# Cut-off sensitivity across all sets

Each cell is the log-rank p of the rank-1 Ridge pathway in drug-treated patients when the network-proximity
cut-off is changed. Columns: the user's alternative cut-offs (-1.182405, -2.1925376), the 5% value (-1.6449) and
the paper's value ({PAPER_CUTOFF}). In brackets: how many pathways pass that cut-off. Brighter = smaller p.
Numbers: `tables/cross_set_cutoff_sensitivity.tsv`.
""")


if __name__ == "__main__":
    main()
