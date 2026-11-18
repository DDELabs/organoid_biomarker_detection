"""B1 figures: per-cancer external forest plots, KM curves in the main cohorts, summary heatmap."""
import numpy as np
import pandas as pd
from statsmodels.duration.survfunc import SurvfuncRight

from common import DATA, OK, TAB, letter, logrank, plot_style, save_fig

plt = plot_style()
R = pd.read_csv(TAB / "B1_external_results.tsv", sep="\t")
SC = pd.read_csv(DATA / "b1_patient_scores.tsv.gz", sep="\t")
SIGLAB = {"COAD_vandewetering__COAD_authors__FLUOROURACIL": "vdW 5-FU (2019 files)",
          "COAD_vandewetering__COAD__FLUOROURACIL": "vdW 5-FU (current)",
          "PDX_gao2015_CRC__COAD__FLUOROURACIL": "Novartis PDX 5-FU",
          "COAD_vandewetering__COAD_authors__OXALIPLATIN": "vdW oxaliplatin (2019 files)",
          "COAD_vandewetering__COAD__OXALIPLATIN": "vdW oxaliplatin (current)",
          "COAD_vandewetering__COAD_authors__IRINOTECAN": "vdW irinotecan (2019 files)",
          "COAD_vandewetering__COAD__IRINOTECAN": "vdW irinotecan (current)",
          "BLCA_lee2018_paper__BLCA__CISPLATIN": "Lee cisplatin (paper NES)",
          "BLCA_lee2018_current__BLCA__CISPLATIN": "Lee cisplatin (current)",
          "OV_vias2023__OV__CARBOPLATIN": "Vias platinum (oxaliplatin screen)",
          "OV_vias2023__OV__PACLITAXEL": "Vias paclitaxel",
          "LIHC_ji__LIHC__SORAFENIB": "Ji sorafenib", "LIHC_broutier__LIHC__SORAFENIB": "Broutier sorafenib",
          "LUAD_cellline_substitute__LUAD__CISPLATIN": "LUAD cell-line cisplatin",
          "LUAD_cellline_substitute__LUAD__CARBOPLATIN": "LUAD cell-line carboplatin"}
R["siglab"] = R["analysis_id"].map(SIGLAB) + " [" + R["kind"] + "]"
R["testlab"] = R["cohort"] + " " + R["endpoint"] + ": " + R["comparison"]
CANCER = {"colorectal": ["CRC_5FU", "CRC_OXALIPLATIN", "CRC_IRINOTECAN"], "bladder": ["BLCA_CISPLATIN"],
          "ovarian": ["OV_PLATINUM", "OV_TAXANE"], "liver and lung": ["LIHC_SORAFENIB", "LUNG_PLATINUM"]}
CL = (0.1, 10)


def forest(ax, d, est, lo, hi, title):
    y = np.arange(len(d))[::-1]
    for yi, (_, r) in zip(y, d.iterrows()):
        c = OK["blue"] if r["kind"] == "rank1" else OK["orange"]
        if not np.isfinite(r[est]):
            ax.plot([1], [yi], "x", color=c)
            continue
        ax.plot([max(r[lo], CL[0]), min(r[hi], CL[1])], [yi, yi], color=c, lw=1)
        ax.plot([np.clip(r[est], *CL)], [yi], "D" if r["has_comparator"] else "s", color=c, ms=3.5)
        ax.text(CL[1] * 1.1, yi, f"{r[est]:.2f} p={r[est.replace('est', 'p') if est != 'est' else 'p']:.2g}",
                fontsize=5.8, va="center")
    ax.axvline(1, color="k", lw=0.6)
    ax.set_xscale("log"); ax.set_xlim(*CL)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} | {b}" for a, b in zip(d["testlab"], d["siglab"])], fontsize=5.8)
    ax.set_title(title, fontsize=8.5)


for i, (cancer, groups) in enumerate(CANCER.items(), 1):
    d = R[R.group.isin(groups)].reset_index(drop=True)
    h = max(3.5, 0.13 * len(d) + 1.2)
    fig, axs = plt.subplots(1, 2, figsize=(13, h), sharey=True)
    forest(axs[0], d, "est", "lo", "hi", "Treated patients: HR or OR (non-response) per SD of score")
    di = d.copy()
    forest(axs[1], di, "int_est", "int_lo", "int_hi", "Interaction with comparator arm (treated / comparator)")
    axs[1].tick_params(labelleft=False)
    letter(axs[0], "a", x=-0.02); letter(axs[1], "b", x=-0.02)
    axs[0].set_xlabel("Ratio per SD (log scale); > 1 = expected direction")
    axs[1].set_xlabel("Interaction ratio (log scale); > 1 = score matters more with the drug")
    axs[0].plot([], [], "s", color=OK["blue"], label="single rank-1 pathway")
    axs[0].plot([], [], "s", color=OK["orange"], label="multi-pathway (top-10 Ridge)")
    axs[0].legend(frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.04), ncol=2)
    save_fig(fig, f"B1_Fig{i}_forest_{cancer.replace(' ', '_')}", f"""
**External validation, {cancer}.** Each row = one frozen organoid-trained signature (right label; [rank1] = the
single top pathway, [top10] = the multi-pathway Ridge model, coefficients exactly as saved by the organoid-only
agents) applied to one external cohort and endpoint (left label). Nothing was refitted on patients. Pathway values
are within-sample rank scores (comparable across microarray and RNA-seq platforms), standardised within each
cohort; the score is standardised again (per SD). (a) Effect in drug-treated patients: hazard ratio (HR) for
survival endpoints (OS, RFS, DFS, PFS) or odds ratio (OR) of NOT responding for response endpoints. Values above 1
are the pre-specified direction (higher predicted resistance = worse outcome). (b) Interaction with the comparator
arm (only where a comparator exists; diamond markers in a): the score's effect in treated patients divided by its
effect in comparator patients; above 1 = drug-specific in the expected way. Lines = 95% CI (clipped 0.1-10). Numbers
= estimate and Wald p. Full table: tables/B1_external_results.tsv.""")
    plt.close(fig)

# ---------------------------------------------------------------- KM in the main cohorts
MAIN = [("GSE39582", "RFS", "5-FU-based vs no chemo (stage II-III)", "COAD_vandewetering__COAD_authors__FLUOROURACIL|rank1"),
        ("GSE103479", "OS", "5-FU adjuvant vs surgery only", "COAD_vandewetering__COAD_authors__FLUOROURACIL|rank1"),
        ("GSE72970", "PFS", "FOLFOX vs FOLFIRI", "COAD_vandewetering__COAD__OXALIPLATIN|rank1"),
        ("GSE48276", "OS", "neoadjuvant MVAC vs none", "BLCA_lee2018_paper__BLCA__CISPLATIN|rank1"),
        ("GSE9891", "OS", "platinum+taxane vs platinum alone", "OV_vias2023__OV__PACLITAXEL|rank1"),
        ("GSE42127", "OS", "adjuvant platinum vs observation", "LUAD_cellline_substitute__LUAD__CISPLATIN|rank1")]
fig, axs = plt.subplots(len(MAIN), 2, figsize=(9.5, 3.0 * len(MAIN)))
L = iter("abcdefghijklmnop")
for row, (coh, ep, comp, sig) in enumerate(MAIN):
    d = SC[(SC.cohort == coh) & (SC.endpoint == ep.lower()) & (SC.comparison == comp) & (SC.sig_id == sig)]
    for col, arm in enumerate((1, 0)):
        ax = axs[row, col]
        letter(ax, next(L), x=-0.15)
        a = d[d.arm1 == arm]
        if len(a) < 5:
            ax.set_axis_off(); continue
        g = (a.score > a.score.median()).astype(int).values
        for lab, m, c in (("low score (predicted sensitive)", g == 0, OK["blue"]),
                          ("high score (predicted resistant)", g == 1, OK["verm"])):
            sf = SurvfuncRight(a["T"].values[m], a["E"].values[m])
            ax.step(np.r_[0, sf.surv_times], np.r_[1, sf.surv_prob], where="post", color=c,
                    label=f"{lab}: n={m.sum()}, events={int(a['E'].values[m].sum())}")
        p = logrank(a["T"], a["E"], g)
        armname = comp.split(" vs ")[0] if arm == 1 else (comp.split(" vs ")[1] if " vs " in comp else "comparator")
        ax.set_title(f"{coh} {ep}, {armname}\n{SIGLAB[sig.split('|')[0]]}; log-rank p = {p:.2g}", fontsize=7.5)
        ax.set_ylim(0, 1.03); ax.set_xlabel("Months"); ax.set_ylabel(f"{ep} probability")
        ax.legend(frameon=False, fontsize=6, loc="lower left")
fig.tight_layout()
save_fig(fig, "B1_Fig5_KM_main_cohorts", """
**Kaplan-Meier curves in the main external cohorts.** Left column = patients who got the drug; right column = the
comparator arm of the same cohort (no chemotherapy, other regimen, or platinum without taxane). Patients are split at
the median of the frozen rank-1 organoid pathway score within each arm: vermilion = high score (predicted resistant),
blue = low score (predicted sensitive). If the signature were drug-specific, the curves would separate on the left
(vermilion lower) but not on the right. Signatures: colorectal 5-FU = the paper's BH3-only pathway (van de Wetering
organoids, authors' files); FOLFOX = current-data oxaliplatin signature; bladder = Lee 2018 cisplatin (paper NES);
ovarian = Vias 2023 paclitaxel; lung = lung cell-line cisplatin (the organoid substitute).""")
plt.close(fig)

# STORM response bars
d = SC[(SC.cohort == "GSE109211")]
fig, axs = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for ax, sig, L2 in zip(axs, ["LIHC_ji__LIHC__SORAFENIB|rank1", "LIHC_broutier__LIHC__SORAFENIB|rank1"], "ab"):
    x = d[d.sig_id == sig]
    vals, labs = [], []
    for arm, an in ((1, "sorafenib"), (0, "placebo")):
        a = x[x.arm1 == arm]
        g = a.score > a.score.median()
        for gg, gn in ((False, "low"), (True, "high")):
            vals.append(a.loc[g == gg, "responder"].mean()); labs.append(f"{an}\n{gn} score\n(n={int((g == gg).sum())})")
    ax.bar(range(4), vals, color=[OK["blue"], OK["verm"], OK["sky"], OK["orange"]])
    for k, v in enumerate(vals):
        ax.text(k, v + 0.01, f"{v:.2f}", ha="center", fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=7)
    rr = R[(R.sig_id == sig)].iloc[0]
    ax.set_title(f"{SIGLAB[sig.split('|')[0]]} [rank1]\nOR non-resp {rr['est']:.2f} (p={rr['p']:.2g}); interaction p={rr['int_p']:.2g}",
                 fontsize=7.5)
    ax.set_ylabel("Share of responders (STORM definition)")
    letter(ax, L2)
save_fig(fig, "B1_Fig6_STORM_sorafenib", """
**Liver: STORM trial (GSE109211, sorafenib vs placebo, randomised).** Share of responders (as defined in the
dataset, recurrence-based) in each arm, split at the median of the frozen rank-1 sorafenib organoid pathway score
(Ji 2023 organoids in a, Broutier 2017 in b). A drug-specific signature would show fewer responders with a high score
in the sorafenib arm only. Title: odds ratio of non-response per SD in the sorafenib arm and the score x arm
interaction p.""")
plt.close(fig)

# ---------------------------------------------------------------- heatmap signature x cohort
R["zdir"] = np.where(R["has_comparator"], R["int_z"], R["z"])
H = R.pivot_table(index="siglab", columns="testlab", values="zdir", aggfunc="first")
order_c = R.drop_duplicates("testlab").sort_values(["group", "cohort"])["testlab"]
order_r = R.drop_duplicates("siglab").sort_values(["group", "siglab"])["siglab"]
H = H.reindex(index=order_r, columns=order_c)
fig, ax = plt.subplots(figsize=(14, 7.5))
im = ax.imshow(H.values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        v = H.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(v) > 2 else "k")
ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns, rotation=75, ha="right", fontsize=6)
ax.set_yticks(range(H.shape[0])); ax.set_yticklabels(H.index, fontsize=6.5)
cb = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cb.set_label("Signed z (> 0 = expected direction)")
save_fig(fig, "B1_Fig7_heatmap_signature_by_cohort", """
**Summary heatmap: every frozen signature in every matching cohort.** Rows = signatures (organoid set, drug,
[rank1] single pathway or [top10] multi-pathway). Columns = cohort, endpoint and comparison. Colour and number = the
signed z-statistic of the primary test: the score x arm interaction where a comparator arm exists, otherwise the
effect in treated patients (log HR or log OR of non-response divided by its standard error). Red (> 0) = the
pre-specified direction; blue = the opposite. |z| > 1.96 corresponds to p < 0.05. Blank = signature not matched to
that cohort's drug.""")
plt.close(fig)
print("done")
