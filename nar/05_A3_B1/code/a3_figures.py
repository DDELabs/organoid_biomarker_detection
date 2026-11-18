"""A3 figures: master forest, QQ + observed/expected counts, adjusted vs unadjusted HR, interaction forest,
power vs events. Reads tables/A3_master_results.tsv."""
import numpy as np
import pandas as pd
from scipy.stats import beta, binom

from common import OK, TAB, letter, mdhr, plot_style, power_logrank, save_fig

plt = plot_style()
M = pd.read_csv(TAB / "A3_master_results.tsv", sep="\t")


def label(r):
    s = r["analysis_id"].split("__")
    setn = s[0].replace("_vandewetering", " vdW").replace("_lee2018", " Lee").replace("_cellline_substitute", " cell lines") \
        .replace("_alshihabi2024", " AlShihabi").replace("_vias2023", " Vias").replace("PDX_gao2015_", "PDX ") \
        .replace("_tiriac", " Tiriac").replace("_shi", " Shi").replace("_ji", " Ji").replace("_broutier", " Broutier") \
        .replace("_isella_pdx", " Isella PDX")
    coh = s[1].replace("COAD_authors", "COAD (2019 files)")
    return f"{setn} | {coh} | {s[2].lower()} (n={int(r['treated_n'])}, d={int(r['treated_events'])})"


M["label"] = M.apply(label, axis=1)
CLIP = (0.02, 50)


def forest(ax, d, est, lo, hi, colour, xlab):
    y = np.arange(len(d))[::-1]
    for yi, (_, r) in zip(y, d.iterrows()):
        e, l, h = r[est], r[lo], r[hi]
        c = colour(r)
        if not np.isfinite(e) or e > 1e3 or e < 1e-3:
            ax.plot([1], [yi], marker="x", color=c, ms=5)
            ax.text(CLIP[1] * 0.9, yi, "not estimable", fontsize=6.5, va="center", ha="right", color=OK["grey"])
            continue
        ax.plot([max(l, CLIP[0]), min(h, CLIP[1])], [yi, yi], color=c, lw=1.2)
        ax.plot([e], [yi], "s", color=c, ms=4)
    ax.axvline(1, color="k", lw=0.7)
    ax.set_xscale("log")
    ax.set_xlim(*CLIP)
    ax.set_yticks(y)
    ax.set_yticklabels(d["label"], fontsize=6.5)
    ax.set_xlabel(xlab)


def dircol(r):
    return OK["blue"] if str(r["direction"]).startswith("expected") else OK["verm"]


# ---------------------------------------------------------------- 1 master forest
fig, ax = plt.subplots(figsize=(8.5, 12))
forest(ax, M, "hr", "hr_lo", "hr_hi", dircol, "Hazard ratio, predicted non-responders vs responders (log scale)")
for yi, (_, r) in zip(np.arange(len(M))[::-1], M.iterrows()):
    ax.text(CLIP[1] * 1.15, yi, f"p={r['treated_p']:.2g}  q={r['BH_q']:.2f}", fontsize=6.3, va="center")
ax.plot([], [], "s", color=OK["blue"], label="expected direction (HR > 1)")
ax.plot([], [], "s", color=OK["verm"], label="opposite direction (HR < 1)")
ax.legend(loc="lower left", frameon=False, bbox_to_anchor=(0, 1.0), ncol=2)
save_fig(fig, "A3_Fig1_master_forest", """
**Master forest plot of all 57 organoid-only analyses.** Each row is one lab-model set x drug x TCGA cohort
(n = treated patients, d = deaths). The square is the hazard ratio (HR) of predicted non-responders versus predicted
responders (median split of the frozen rank-1 pathway score); the line is its 95% confidence interval (clipped to
0.02-50). Blue = HR > 1 (the expected direction: predicted non-responders die sooner); vermilion = HR < 1. "not
estimable" = no deaths in one group (colour then follows which group had more deaths). Right: log-rank p and the
Benjamini-Hochberg q-value across all 57 tests.""")
plt.close(fig)

# ---------------------------------------------------------------- 2 QQ + counts
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axs[0]
fams = [("treated_p", "treated log-rank (rank-1)", OK["blue"]), ("untreated_p", "untreated log-rank", OK["orange"]),
        ("int_p", "score x treatment interaction", OK["green"]), ("perm_p", "response-shuffle p", OK["purple"])]
n = len(M)
i = np.arange(1, n + 1)
ex = -np.log10(i / (n + 1))
ax.fill_between(ex, -np.log10(beta.ppf(0.975, i, n - i + 1)), -np.log10(beta.ppf(0.025, i, n - i + 1)),
                color="0.85", label="95% band under no effect")
for col, lab, c in fams:
    p = np.sort(M[col].dropna().values)
    e = -np.log10(np.arange(1, len(p) + 1) / (len(p) + 1))
    ax.plot(e, -np.log10(p), "o", ms=3.5, color=c, label=lab)
ax.plot([0, ex.max()], [0, ex.max()], "k--", lw=0.8)
ax.set_xlabel("Expected -log10 p (uniform, no real effect)")
ax.set_ylabel("Observed -log10 p")
ax.legend(frameon=False, fontsize=7)
letter(ax, "a")
ax = axs[1]
C = pd.read_csv(TAB / "A3_expected_by_chance.tsv", sep="\t")
xx = np.arange(len(C))
ax.bar(xx - 0.2, C["observed_p_below_0.05"], 0.4, color=OK["blue"], label="observed p < 0.05")
ax.bar(xx + 0.2, C["expected_by_chance"], 0.4, color=OK["grey"], label="expected by chance (5%)")
lo, hi = binom.ppf(0.025, 57, 0.05), binom.ppf(0.975, 57, 0.05)
ax.errorbar(xx + 0.2, C["expected_by_chance"], yerr=[C["expected_by_chance"] - lo, hi - C["expected_by_chance"]],
            fmt="none", ecolor="k", lw=0.8, capsize=2)
for k, r in C.iterrows():
    ax.text(k - 0.2, r["observed_p_below_0.05"] + 0.2, f"{int(r['observed_p_below_0.05'])}", ha="center", fontsize=7)
ax.set_xticks(xx)
ax.set_xticklabels(["treated\nlog-rank", "response\nshuffle", "random-\npathway", "untreated\nlog-rank",
                    "interaction", "multi-pathway\nk=10"], fontsize=7, rotation=35, ha="right")   # audit I23: no overlap
ax.set_ylabel("Number of the 57 analyses")
ax.set_ylim(0, 14.5)
ax.legend(frameon=False, fontsize=7, loc="upper left")
letter(ax, "b")
save_fig(fig, "A3_Fig2_QQ_expected_by_chance", """
**How many small p-values would chance alone give?** (a) QQ plot: each dot is one of the 57 analyses, sorted by
p-value; the x-axis is where it would sit if no analysis had a real effect. Dots above the grey band are smaller
than chance would give. Blue = drug-treated patients (the main test), orange = untreated patients (same split, a
check for general prognosis), green = score x treatment interaction, purple = the response-shuffle p. (b) Observed
number of analyses with p < 0.05 versus the 2.85 expected by chance (5% of 57); black bar = 95% range of a
binomial(57, 0.05). The 57 analyses share patients (for example colon and colon + rectum) so they are not fully
independent; the band is a guide, not an exact test.""")
plt.close(fig)

# ---------------------------------------------------------------- 3 adjusted vs unadjusted
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
for ax, (a, b, t, L) in zip(axs, [("hr", "adj_HR", "Median split (group) HR", "a"),
                                  ("cont_HR", "adj_cont_HR", "Continuous score, HR per SD", "b")]):
    d = M[np.isfinite(M[a]) & np.isfinite(M[b]) & (M[a] < 1e3) & (M[b] < 1e3) & (M[b] > 1e-3)]
    cc = [dircol(r) for _, r in d.iterrows()]
    sz = 8 + 60 * d["treated_events"] / d["treated_events"].max()
    ax.scatter(d[a], d[b], s=sz, c=cc, alpha=0.8, edgecolor="k", lw=0.3)
    ax.set_xscale("log"); ax.set_yscale("log")
    lim = (min(d[a].min(), d[b].min()) * 0.8, max(d[a].max(), d[b].max()) * 1.25)
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.axhline(1, color=OK["grey"], lw=0.5); ax.axvline(1, color=OK["grey"], lw=0.5)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(f"Unadjusted {t}")
    ax.set_ylabel(f"Adjusted {t}")
    rho = np.corrcoef(np.log(d[a]), np.log(d[b]))[0, 1]
    ax.set_title(f"{len(d)} estimable analyses; r(log HR) = {rho:.2f}")
    letter(ax, L)
save_fig(fig, "A3_Fig3_adjusted_vs_unadjusted", """
**Does adjusting for clinical factors change the hazard ratios?** Each dot is one analysis (dot size = number of
deaths). x = hazard ratio without adjustment; y = hazard ratio from a Cox model adjusted for age, sex, stage,
tumour purity (ABSOLUTE), MSI status (colon/rectum) and the regimen (indicators for the two most common partner
drugs, e.g. oxaliplatin for 5-FU). Missing covariate values were set to the median with a "missing" flag. Dots on
the dashed line are unchanged by adjustment. (a) median-split groups, as in the paper; (b) the continuous score per
standard deviation. Blue = HR > 1 unadjusted (expected direction), vermilion = HR < 1. Analyses with very few
deaths are unstable after adjustment (wide intervals, see the master table).""")
plt.close(fig)

# ---------------------------------------------------------------- 4 interaction forest
fig, ax = plt.subplots(figsize=(8.5, 12))
d = M.copy()
forest(ax, d, "int_HR", "int_lo", "int_hi",
       lambda r: OK["blue"] if r["int_HR"] > 1 else OK["verm"],
       "Interaction HR per SD of score: treated vs untreated (log scale)")
for yi, (_, r) in zip(np.arange(len(d))[::-1], d.iterrows()):
    ax.text(CLIP[1] * 1.15, yi, f"p={r['int_p']:.2g}", fontsize=6.3, va="center")
save_fig(fig, "A3_Fig4_interaction_forest", """
**Is the score linked to survival more strongly in drug-treated than in untreated patients?** For each analysis a
Cox model was fitted on treated plus untreated (no drug record) TCGA patients with the continuous score (per SD), the
treatment group and their product. The plotted value is the interaction hazard ratio (the score's HR in treated
patients divided by its HR in untreated patients). Above 1 (blue) = the score marks risk more in treated patients, as
a drug-specific (predictive) marker should. Below 1 (vermilion) = the opposite. Lines = 95% CI.""")
plt.close(fig)

# ---------------------------------------------------------------- 5 power vs events
fig, axs = plt.subplots(1, 2, figsize=(10, 4.3))
D = np.arange(2, 400)
ax = axs[0]
for hr, c in ((1.25, OK["grey"]), (1.5, OK["blue"]), (2.0, OK["green"]), (3.0, OK["orange"])):
    ax.plot(D, [power_logrank(x, hr) for x in D], color=c, label=f"true HR = {hr}")
ev = M["treated_events"].values
ax.scatter(ev, [power_logrank(x, 1.5) for x in ev], s=14, color=OK["blue"], edgecolor="k", lw=0.3, zorder=3,
           label="the 57 analyses (at HR 1.5)")
ax.axhline(0.8, color="k", ls=":", lw=0.8)
ax.set_xscale("log")
ax.set_xlabel("Deaths among treated patients (log scale)")
ax.set_ylabel("Power (two-sided alpha 0.05)")
ax.legend(frameon=False, fontsize=7)
letter(ax, "a")
ax = axs[1]
ax.plot(D, [mdhr(x) for x in D], color="k")
ax.scatter(ev, M["MDHR_80"], s=14, c=[dircol(r) for _, r in M.iterrows()], edgecolor="k", lw=0.3, zorder=3)
for k, (_, r) in enumerate(M[M["treated_p"] < 0.05].sort_values("treated_events").iterrows()):   # audit I23: staggered labels
    ax.annotate(r["analysis_id"].split("__")[-1].lower(), (r["treated_events"], r["MDHR_80"]), fontsize=7,
                xytext=(4, 10 if k % 2 == 0 else -12), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", lw=0.4, color="0.4"))
ax.set_xscale("log"); ax.set_yscale("log")
from matplotlib.ticker import NullFormatter
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([1.25, 1.5, 2, 3, 5, 10, 20, 50]); ax.set_yticklabels(["1.25", "1.5", "2", "3", "5", "10", "20", "50"])
ax.set_xlabel("Deaths among treated patients (log scale)")
ax.set_ylabel("Smallest HR detectable with 80% power")
letter(ax, "b")
save_fig(fig, "A3_Fig5_power_vs_events", """
**How much could each analysis detect?** (a) Power of the median-split log-rank test (Schoenfeld formula, equal
groups, two-sided 0.05) against the number of deaths, for true hazard ratios 1.25-3. Dots = the 57 analyses at a true
HR of 1.5; the dotted line marks 80% power. (b) The smallest hazard ratio each analysis could detect with 80% power.
Labelled dots are the 7 analyses with p < 0.05. Most analyses could only detect very large effects (HR >= 3).""")
plt.close(fig)
print("done")
