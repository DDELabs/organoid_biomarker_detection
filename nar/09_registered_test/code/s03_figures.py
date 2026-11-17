"""Figures for Stages 1-2 (each with a plain-language caption file)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.duration.survfunc import SurvfuncRight, survdiff  # noqa: E402

import common as C  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42})
T = pd.read_csv(C.TAB / "stage1_primary.tsv", sep="\t")
SX = pd.read_csv(C.TAB / "stage1_sensitivity.tsv", sep="\t")
NC = pd.read_csv(C.TAB / "stage1_negative_controls.tsv", sep="\t")
SP = pd.read_csv(C.TAB / "stage1_specificity_matrix.tsv", sep="\t")
PC = pd.read_csv(C.TAB / "stage1_per_cancer.tsv", sep="\t")
MA = pd.read_csv(C.TAB / "stage1_meta_analysis.tsv", sep="\t")
NU = pd.read_csv(C.TAB / "stage1_random_null_z.tsv.gz", sep="\t")
E2 = pd.read_csv(C.TAB / "stage2_controlled.tsv", sep="\t")
E1 = pd.read_csv(C.TAB / "stage2_single_arm_supportive.tsv", sep="\t")
G2 = pd.read_csv(C.TAB / "stage1_gate_G2.tsv", sep="\t")
G3 = pd.read_csv(C.TAB / "stage2_gate_G3.tsv", sep="\t")
L = C.LABEL


def letter(ax, s):
    ax.text(-0.12, 1.04, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def forest(ax, labels, est, lo, hi, colors=None, ref=1.0, xlabel="HR", markers=None):
    y = np.arange(len(labels))[::-1]
    for i, (e, a, b) in enumerate(zip(est, lo, hi)):
        c = CB if colors is None else colors[i]
        ax.plot([a, b], [y[i], y[i]], color=c, lw=1.3)
        ax.plot(e, y[i], "s" if markers is None else markers[i], color=c, ms=5)
    ax.axvline(ref, color=C.CB["grey"], ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_ylim(-0.5, len(labels) - 0.5)
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    vals = [v for v in list(lo) + list(hi) if np.isfinite(v)]
    ticks = [t for t in (0.1, 0.25, 0.5, 1, 2, 4, 10, 25) if min(vals) * 0.8 <= t <= max(vals) * 1.25]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel(xlabel)


CB = C.CB["blue"]


# ---------------------------------------------------------------- 1 flow
def fig_flow():
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.axis("off")
    n_tcga = int(T["PFI_n"].max())
    plan_ts = open(C.HERE / "ANALYSIS_PLAN.timestamp").read().split()[0]
    boxes = [
        (0.5, 0.92, f"Stage 0 (no patient outcomes)\n35 signed pathways for 7 drugs, frozen 2026-11-17T13:03:40Z\nsha256 {C.FROZEN_SHA[:16]}..."),
        (0.5, 0.76, f"Analysis plan written before outcomes were opened\nANALYSIS_PLAN.md ({plan_ts})\n7 primary tests (PFI interaction), Holm correction"),
        (0.5, 0.58, "Stage 1: TCGA, pooled per drug, stratified by cancer\n"
                    f"{len(T)} drugs, {int(T['n_projects'].sum())} drug x cancer strata, up to {n_tcga} patients per drug\n"
                    "score x exposure Cox; IPW, purity, active comparator; 1,000 random sets"),
        (0.5, 0.42, f"Gate G2 per drug: {int(G2['G2_pass'].sum())} of {len(G2)} pass"),
        (0.5, 0.25, f"Stage 2: external cohorts, one-sided in the frozen direction\n"
                    f"{int((E2['role'] == 'primary').sum())} primary controlled tests in {E2.loc[E2['role'] == 'primary', 'cohort'].nunique()} cohorts"
                    f" + {len(E1)} single-arm supportive tests"),
        (0.5, 0.09, f"Gate G3 per drug: {int(G3['G3_pass'].sum())} of {len(G3)} pass"),
    ]
    for x, y, s in boxes:
        ax.add_patch(FancyBboxPatch((x - 0.42, y - 0.06), 0.84, 0.12, boxstyle="round,pad=0.01", fc="#EAF2FA" if "Gate" not in s else "#FDF0E0",
                                    ec=C.CB["blue"] if "Gate" not in s else C.CB["orange"], lw=1.2, transform=ax.transAxes))
        ax.text(x, y, s, ha="center", va="center", fontsize=8.5, transform=ax.transAxes)
    for (x, y, _), (_, y2, _) in zip(boxes[:-1], boxes[1:]):
        ax.annotate("", xy=(x, y2 + 0.065), xytext=(x, y - 0.065), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="black", lw=1))
    C.save_fig(fig, "fig1_preregistration_flow", """
**Figure 1. How the test was pre-registered.** Each box is one step, top to bottom. The 35 pathway hypotheses were
fixed (frozen) from cell-line and network data only, before any patient outcome was used. The file's fingerprint
(sha256) was checked by the code at every run. Then the analysis plan was written, still without outcomes. Stage 1
tests each drug's score in TCGA patients, pooled across cancers. Stage 2 tests the same frozen score in outside
patient cohorts that have a comparison group. Orange boxes are the decision gates (G2, G3) and how many drugs passed.
""")


# ---------------------------------------------------------------- 2 pooled forest
def fig_pooled():
    fig, axes = plt.subplots(1, 2, figsize=(10, 9), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    labs, est, lo, hi, cols = [], [], [], [], []
    for _, r in T.iterrows():
        labs += [f"{L[r['drug']]} PFI", f"{L[r['drug']]} OS"]
        est += [r["PFI_HR_int"], r["OS_HR_int"]]
        lo += [r["PFI_lo"], r["OS_lo"]]
        hi += [r["PFI_hi"], r["OS_hi"]]
        cols += [C.CB["blue"], C.CB["sky"]]
    forest(ax, labs, est, lo, hi, cols, xlabel="Interaction HR per SD (score x drug)")
    for i, (_, r) in enumerate(T.iterrows()):
        ax.text(1.02, 1 - (2 * i + 0.5) / (2 * len(T)), f"Holm p={r['PFI_holm_p']:.2g}", transform=ax.transAxes, fontsize=8, va="center")
    ax.set_title("Primary (PFI) and secondary (OS)")
    letter(ax, "a")
    ax = axes[1]
    an = ["primary", "S1 no purity adjustment", "S2 IPW (bootstrap SE)", "S3 active comparator (other systemic therapy)",
          "S4 + proliferation score and its interaction", "S5 5-FU or capecitabine"]
    short = {"primary": "primary", "S1 no purity adjustment": "no purity", "S2 IPW (bootstrap SE)": "IPW",
             "S3 active comparator (other systemic therapy)": "active comparator",
             "S4 + proliferation score and its interaction": "+ proliferation", "S5 5-FU or capecitabine": "5-FU or capecitabine"}
    pal = [C.CB["blue"], C.CB["orange"], C.CB["green"], C.CB["red"], C.CB["purple"], C.CB["black"]]
    labs, est, lo, hi, cols = [], [], [], [], []
    for _, r in T.iterrows():
        for k, a in enumerate(an):
            if a == "primary":
                e, l_, h_ = r["PFI_HR_int"], r["PFI_lo"], r["PFI_hi"]
            else:
                s = SX[(SX["drug"] == r["drug"]) & (SX["analysis"] == a)]
                if s.empty:
                    continue
                e, l_, h_ = s.iloc[0][["HR_int", "lo", "hi"]]
            labs.append(f"{L[r['drug']]}: {short[a]}")
            est.append(e); lo.append(l_); hi.append(h_); cols.append(pal[k])
    forest(ax, labs, est, lo, hi, cols, xlabel="Interaction HR per SD, PFI")
    ax.tick_params(axis="y", labelsize=8)
    ax.set_title("Sensitivity analyses (PFI)")
    letter(ax, "b")
    fig.tight_layout()
    C.save_fig(fig, "fig2_stage1_pooled_interaction_forest", """
**Figure 2. Stage 1: does the frozen score predict a worse outcome specifically under the drug?** Each square is a
hazard ratio (HR) for the interaction between the score and receiving the drug, per 1 standard deviation (SD) of
score, with its 95% confidence interval (line). TCGA patients are pooled across cancers, with cancer as a stratum
(each cancer has its own baseline risk). The frozen hypothesis predicts HR above 1 (dashed line = no interaction).
**a**, primary endpoint PFI (progression-free interval, dark blue) and secondary OS (overall survival, light blue);
Holm-adjusted p for the 7 primary tests on the right. **b**, the same PFI interaction under the pre-specified
sensitivity analyses: without purity adjustment, inverse-probability weighting (IPW; reweights patients so treated
and untreated look alike on age, stage, sex, purity and cancer), comparator limited to patients who got some other
drug, adding a proliferation score, and (5-FU only) counting capecitabine as 5-FU.
""")


# ---------------------------------------------------------------- 3 KM
def fig_km():
    d = pd.read_csv(C.DATA / "tcga_stage1_patients.tsv.gz", sep="\t", index_col=0)
    fig, axes = plt.subplots(7, 2, figsize=(7.2, 15.5))
    rows = []
    for i, drug in enumerate(C.DRUGS):
        projs = T.loc[T["drug"] == drug, "projects"].iloc[0].split(";")
        dd = d[d["project"].isin(projs)].copy()
        med = dd.groupby("project")[f"score_{drug}"].transform("median")
        dd["high"] = (dd[f"score_{drug}"] > med).astype(int)
        st = dd[f"start_{drug}"]
        ex = dd[st.notna() & (st < dd["PFI_time"])].copy()
        ex["t"] = (ex["PFI_time"] - ex[f"start_{drug}"]) / 30.44
        co = dd[st.isna()].copy()
        co["t"] = co["PFI_time"] / 30.44
        for j, (g, name) in enumerate([(ex, "treated with drug (from drug start)"), (co, "comparator (from diagnosis)")]):
            ax = axes[i, j]
            for h, col, lab in [(0, C.CB["blue"], "low score"), (1, C.CB["orange"], "high score (predicted resistant)")]:
                gg = g[g["high"] == h]
                sf = SurvfuncRight(gg["t"].values, gg["PFI_event"].values)
                ax.step(np.r_[0, sf.surv_times], np.r_[1, sf.surv_prob], where="post", color=col, lw=1.3,
                        label=f"{lab} (n={len(gg)}, ev={int(gg['PFI_event'].sum())})")
            chi, p = survdiff(g["t"].values, g["PFI_event"].values, g["high"].values)
            ax.set_xlim(0, 60)
            ax.set_ylim(0, 1.02)
            ax.set_title(f"{L[drug]}: {name}\nlog-rank p={p:.2g}", fontsize=8.5)
            ax.legend(loc="upper right", fontsize=8, frameon=False)
            if i == 6:
                ax.set_xlabel("Months")
            if j == 0:
                ax.set_ylabel("PFI probability")
            letter(ax, "abcdefghijklmn"[2 * i + j])
            rows.append({"drug": drug, "group": name, "n": len(g), "n_high": int(g["high"].sum()), "events": int(g["PFI_event"].sum()),
                         "logrank_chi2": chi, "logrank_p": p})
    fig.tight_layout()
    pd.DataFrame(rows).to_csv(C.TAB / "stage1_km_logrank.tsv", sep="\t", index=False)
    C.save_fig(fig, "fig3_km_high_vs_low_score", """
**Figure 3. Kaplan–Meier curves of progression-free interval (PFI) by frozen score, per drug (TCGA, pooled
cancers).** Orange = score above the median of the patient's own cancer type (predicted more resistant); blue =
below. Left column: patients who received the drug before any progression, timed from the drug start. Right column:
patients of the same cancer types who never received the drug, timed from diagnosis. If the score were drug-specific,
the orange and blue curves would separate on the left (orange lower) more than on the right. These curves are
descriptive; the formal test is the interaction model in Figure 2. n = patients, ev = PFI events; curves cut at 60
months.
""")


# ---------------------------------------------------------------- 4 null
def fig_null():
    fig, axes = plt.subplots(2, 4, figsize=(9.5, 4.8))
    for k, drug in enumerate(C.DRUGS):
        ax = axes.flat[k]
        z = NU[drug].values
        r = T[T["drug"] == drug].iloc[0]
        ax.hist(z, bins=40, color=C.CB["sky"], edgecolor="white")
        ax.axvline(r["PFI_z"], color=C.CB["red"], lw=1.8)
        ax.set_title(f"{L[drug]}\nobserved z={r['PFI_z']:.2f}, null p={r['random_null_p_one_sided']:.3f}", fontsize=8.5)
        ax.set_xlabel("Interaction z (random 5-pathway sets)")
        letter(ax, "abcdefg"[k])
    axes.flat[-1].axis("off")
    fig.tight_layout()
    C.save_fig(fig, "fig4_random_signature_null", """
**Figure 4. Is the frozen score better than a random score?** For each drug, 1,000 random scores were built, each
from 5 pathways picked at random with random signs (the same size as the frozen hypothesis). Each was tested with
exactly the same Stage 1 model. The blue histogram is the spread of their interaction z-values (z = estimate divided
by its standard error; positive = the direction the hypothesis predicts). The red line is the frozen score. "null p"
is the share of random scores with a z at least as large (one-sided). A value near 0.5 means the frozen score does
no better than a typical random set.
""")


# ---------------------------------------------------------------- 5 negative controls
def fig_negctl():
    fig, axes = plt.subplots(1, 2, figsize=(11, 7), gridspec_kw={"width_ratios": [1, 1.1]})
    ax = axes[0]
    labs, est, lo, hi, cols = [], [], [], [], []
    for drug in C.DRUGS:
        r = T[T["drug"] == drug].iloc[0]
        labs.append(f"{L[drug]}: frozen score"); est.append(r["PFI_HR_int"]); lo.append(r["PFI_lo"]); hi.append(r["PFI_hi"]); cols.append(C.CB["blue"])
        for ctl, col in [("proliferation", C.CB["orange"]), ("tumour purity", C.CB["green"])]:
            s = NC[(NC["drug"] == drug) & (NC["control"] == ctl)].iloc[0]
            labs.append(f"{L[drug]}: {ctl}"); est.append(s["HR_int"]); lo.append(s["lo"]); hi.append(s["hi"]); cols.append(col)
    forest(ax, labs, est, lo, hi, cols, xlabel="Interaction HR per SD, PFI")
    ax.tick_params(axis="y", labelsize=8)
    ax.set_title("Negative controls")
    letter(ax, "a")
    ax = axes[1]
    M = SP.pivot(index="score", columns="exposure", values="z").reindex(index=C.DRUGS, columns=C.DRUGS)
    v = np.nanmax(np.abs(M.values))
    im = ax.imshow(M.values, cmap="PuOr_r", vmin=-v, vmax=v)
    ax.set_xticks(range(7)); ax.set_xticklabels([L[x] for x in C.DRUGS], rotation=45, ha="right")
    ax.set_yticks(range(7)); ax.set_yticklabels([L[x] for x in C.DRUGS])
    for i in range(7):
        for j in range(7):
            ax.text(j, i, f"{M.values[i, j]:.1f}", ha="center", va="center", fontsize=8,
                    fontweight="bold" if i == j else "normal")
    ax.set_xlabel("Drug exposure tested"); ax.set_ylabel("Frozen score used")
    plt.colorbar(im, ax=ax, fraction=0.046, label="Interaction z (+ = predicted direction)")
    ax.set_title("Specificity: every score x every drug")
    letter(ax, "b")
    fig.tight_layout()
    C.save_fig(fig, "fig5_negative_controls", """
**Figure 5. Negative controls and drug specificity (TCGA, PFI).** **a**, the frozen score (blue) next to two control
scores tested with the same model: a proliferation score (orange; how fast tumour cells divide) and tumour purity
(green; the share of cancer cells in the sample). The controls are not expected to interact with any one drug; if
they did as strongly as the frozen score, the frozen score would just be tracking these general traits. **b**, each
row is one drug's frozen score, each column one drug's exposure; the number is the interaction z-value (positive =
higher score goes with relatively worse PFI under that drug). The diagonal (bold) is the pre-registered test; a
drug-specific score would stand out on the diagonal.
""")


# ---------------------------------------------------------------- 6 per-cancer
def fig_percancer():
    fig, axes = plt.subplots(4, 2, figsize=(8.5, 11))
    for k, drug in enumerate(C.DRUGS):
        ax = axes.flat[k]
        g = PC[(PC["drug"] == drug) & PC["fitted"]]
        m = MA[MA["drug"] == drug]
        labs = [f"{r.project} (n={r.n}, ev {int(r.events_exposed)}/{int(r.events_unexposed)})" for r in g.itertuples()]
        est, lo, hi = list(g["HR_int"]), list(g["lo"]), list(g["hi"])
        cols, mk = [C.CB["blue"]] * len(g), ["s"] * len(g)
        if len(m):
            m = m.iloc[0]
            labs.append(f"Random-effects (I2={m['I2']:.0%})"); est.append(m["HR_re"]); lo.append(m["lo_re"]); hi.append(m["hi_re"])
            cols.append(C.CB["red"]); mk.append("D")
        skipped = PC[(PC["drug"] == drug) & ~PC["fitted"]]["project"].tolist()
        forest(ax, labs, est, lo, hi, cols, markers=mk, xlabel="Interaction HR per SD, PFI")
        ax.set_title(f"{L[drug]}" + (f"  (too few events: {', '.join(skipped)})" if skipped else ""), fontsize=8.5)
        letter(ax, "abcdefg"[k])
    axes.flat[-1].axis("off")
    fig.tight_layout()
    C.save_fig(fig, "fig6_per_cancer_forest_meta", """
**Figure 6. Stage 1 per cancer type, with a random-effects meta-analysis.** For each drug, the score x drug
interaction hazard ratio (HR, per 1 SD of score, PFI) is estimated separately in each TCGA cancer type (blue squares,
95% confidence intervals). "ev a/b" = PFI events while on the drug / while not on it. The red diamond combines the
cancers with a random-effects meta-analysis (this allows the true effect to differ between cancers); I2 is the share
of the spread between cancers that is beyond chance. HR above 1 = the direction the frozen hypothesis predicts.
Cancers with fewer than 5 events in either group are listed in the title and not plotted.
""")


# ---------------------------------------------------------------- 7 external forest
def fig_external():
    g = E2[E2["fitted"].eq(True)].copy()
    g["label"] = [f"{L[r.drug]} | {r.cohort} {r.endpoint} ({r.role}{', rand.' if r.randomised else ''})" for r in g.itertuples()]
    # put everything on "HR-like" scale where >1 = expected direction: OR_int is inverted
    g["e"] = np.where(g["kind"] == "tte", g["estimate"], 1 / g["estimate"])
    g["l"] = np.where(g["kind"] == "tte", g["lo"], 1 / g["hi"])
    g["h"] = np.where(g["kind"] == "tte", g["hi"], 1 / g["lo"])
    s = E1.copy()
    s["label"] = [f"{L[r.drug]} | {r.cohort} {r.endpoint} (single arm)" for r in s.itertuples()]
    s["e"] = np.where(s["kind"] == "tte", s["estimate"], 1 / s["estimate"])
    s["l"] = np.where(s["kind"] == "tte", s["lo"], 1 / s["hi"])
    s["h"] = np.where(s["kind"] == "tte", s["hi"], 1 / s["lo"])
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 10), gridspec_kw={"width_ratios": [1.3, 1]})
    cols = [C.CB["blue"] if r == "primary" else (C.CB["sky"] if r == "secondary" else C.CB["grey"]) for r in g["role"]]
    forest(axes[0], list(g["label"]), list(g["e"]), list(g["l"]), list(g["h"]), cols,
           xlabel="Effect in predicted direction (HR_int, or 1/OR_int)")
    for i, r in enumerate(g.itertuples()):
        axes[0].text(1.01, 1 - (i + 0.5) / len(g), f"p1={r.p_one_sided:.2g}, pow={r.power_obs:.2f}", transform=axes[0].transAxes,
                     fontsize=8, va="center")
    axes[0].tick_params(axis="y", labelsize=8)
    axes[0].set_title("Controlled cohorts: score x arm")
    letter(axes[0], "a")
    forest(axes[1], list(s["label"]), list(s["e"]), list(s["l"]), list(s["h"]), [C.CB["green"]] * len(s),
           xlabel="Effect in predicted direction (HR, or 1/OR)")
    for i, r in enumerate(s.itertuples()):
        axes[1].text(1.01, 1 - (i + 0.5) / len(s), f"p1={r.p_one_sided:.2g}", transform=axes[1].transAxes, fontsize=8, va="center")
    axes[1].tick_params(axis="y", labelsize=8)
    axes[1].set_title("Single-arm cohorts (supportive)")
    letter(axes[1], "b")
    fig.tight_layout()
    C.save_fig(fig, "fig7_external_cohort_forest", """
**Figure 7. Stage 2: the frozen scores in outside patient cohorts.** **a**, cohorts with a comparison arm. Each line
is the interaction between the score and being in the drug arm, with a 95% confidence interval. Survival endpoints
give a hazard ratio (HR); response endpoints give an odds ratio (OR), shown here as 1/OR so that, for both, values
above 1 mean "the direction the frozen hypothesis predicts" (a higher score is worse under the drug than in the
comparator). Dark blue = pre-specified primary test, light blue = secondary endpoint, grey = exploratory.
"rand." = randomised trial. p1 = one-sided p in the predicted direction; pow = power to detect an effect of 1.43
(HR) or 0.7 (OR) per SD at one-sided α 0.025. **b**, cohorts where everyone got the drug (no comparator): the score's
plain association with outcome. These cannot separate "predicts drug benefit" from "predicts outcome in general".
""")


# ---------------------------------------------------------------- 8 summary heatmap
def fig_summary():
    rows = []
    for _, r in T.iterrows():
        rows.append({"drug": r["drug"], "col": "TCGA pooled PFI", "z": r["PFI_z"], "p": r["PFI_p_one_sided"]})
        rows.append({"drug": r["drug"], "col": "TCGA pooled OS", "z": r["OS_z"], "p": stats.norm.sf(r["OS_z"])})
    for r in E2[E2["fitted"].eq(True) & E2["role"].isin(["primary", "exploratory"])].itertuples():
        z = r.b / r.se * (1 if r.kind == "tte" else -1)
        rows.append({"drug": r.drug, "col": f"{r.cohort} {r.endpoint}", "z": z, "p": r.p_one_sided})
    for r in E1.itertuples():
        z = r.b / r.se * (1 if r.kind == "tte" else -1)
        rows.append({"drug": r.drug, "col": f"{r.cohort} {r.endpoint} (1-arm)", "z": z, "p": r.p_one_sided})
    D = pd.DataFrame(rows)
    D.to_csv(C.TAB / "summary_drug_by_cohort.tsv", sep="\t", index=False)
    cols = list(dict.fromkeys(D["col"]))
    M = D.pivot_table(index="drug", columns="col", values="z").reindex(index=C.DRUGS, columns=cols)
    Pm = D.pivot_table(index="drug", columns="col", values="p").reindex(index=C.DRUGS, columns=cols)
    fig, ax = plt.subplots(figsize=(12, 4.6))
    v = 3
    im = ax.imshow(M.values, cmap="PuOr_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(7)); ax.set_yticklabels([L[x] for x in C.DRUGS])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M.values[i, j]):
                ax.text(j, i, f"{M.values[i, j]:.1f}" + ("*" if Pm.values[i, j] < 0.025 else ""), ha="center", va="center", fontsize=8)
    ax.set_facecolor("#F2F2F2")
    plt.colorbar(im, ax=ax, fraction=0.03, label="z in predicted direction")
    C.save_fig(fig, "fig8_summary_heatmap", """
**Figure 8. Summary: every test of every frozen score, drug by cohort.** Each cell is a z-value (estimate divided by
its standard error) oriented so that positive (brown) = the direction the frozen hypothesis predicts and negative
(purple) = the opposite direction. Grey cells = no data for that drug in that cohort. "*" = one-sided p < 0.025.
The first two columns are the TCGA pooled interaction tests (Stage 1); the next are controlled external cohorts
(Stage 2, primary or exploratory endpoint); "(1-arm)" columns are single-arm cohorts, which test only a plain
association with outcome.
""")


if __name__ == "__main__":
    for f in (fig_flow, fig_pooled, fig_km, fig_null, fig_negctl, fig_percancer, fig_external, fig_summary):
        f()
        print(f.__name__, "done", flush=True)
