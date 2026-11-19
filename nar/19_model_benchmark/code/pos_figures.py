"""Figures + model-comparison table for the W3 positives benchmark (reads tables/positives_all_models.tsv)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parents[1]; TAB = HERE / "tables"; FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "pdf.fonttype": 42, "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
OK = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "green": "#009E73", "sky": "#56B4E9", "purple": "#CC79A7", "grey": "#777777"}
SHORT = {  # column labels (task id -> short label)
    "bladder_logrank_paperdata": "Bladder cisplatin\nsurvival (paper data)", "bladder_logrank_currentdata": "Bladder cisplatin\nsurvival (current data)",
    "bladder_ERCC2_paperdata": "Bladder ERCC2\n(paper data)", "bladder_ERCC2_currentdata": "Bladder ERCC2\n(current data)",
    "colon_cetuximab_BRAF": "Colon cetuximab\nvs BRAF V600E", "rectal_GSE87211_5FU": "Rectal GSE87211\n5-FU, DFS",
    "ovarian_olaparib_GSE63885": "Ovarian olaparib\nGSE63885", "ovarian_olaparib_TCGAOV": "Olaparib sig.\nTCGA-OV (context)", "ovarian_olaparib_GSE9891": "Olaparib sig.\nGSE9891 (context)",
    "lung_untreated_prognosis": "Lung untreated\nprognosis", "ovarian_platinum_joint_survival": "Ovarian platinum\njoint, survival", "ovarian_platinum_joint_response": "Ovarian platinum\njoint, response",
    "ovarian_paclitaxel_joint_survival": "Ovarian paclitaxel\njoint, survival", "ovarian_paclitaxel_joint_response": "Ovarian paclitaxel\njoint, response"}
ORDER = list(SHORT)
ORIG_LABEL = {"Paper_AA_pathway": "paper amino-acid pathway", "Paper_rank1": "paper rank-1 pathway", "Paper_top10_frozen": "frozen 5-FU top-10 signature",
              "Corr_signature": "correlation signature", "Cox_ridge_all": "Cox ridge (673)", "Joint_original": "joint model"}
MODEL_ROWS = ["L2", "L1", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest", "XGBoost", "MLP", "NetSmoothL2", "GCN_prior",
              "L2_shortlist", "L1_shortlist", "ElasticNet_shortlist", "SVM_linear_shortlist", "Proliferation_score", "RandomSets_median"]


def save(fig, stem, caption):
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight"); fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    open(FIG / f"{stem}.caption.md", "w").write(caption.strip() + "\n"); plt.close(fig)


def comparison(A):
    rows = []
    for m, d in A[A.model.isin(MODEL_ROWS)].groupby("model"):
        dd = d.dropna(subset=["diff_vs_orig"])
        better = int((dd.diff_vs_orig > 0).sum()); worse = int((dd.diff_vs_orig < 0).sum())
        sig_b = int(((dd.diff_vs_orig > 0) & (dd.diff_p < 0.05)).sum()); sig_w = int(((dd.diff_vs_orig < 0) & (dd.diff_p < 0.05)).sum())
        pos = dd[~dd.task.str.contains("TCGAOV|GSE9891")]       # without the two context cohorts
        rows.append(dict(model=m, tasks=len(dd), mean_metric=d.est.mean(), mean_diff_vs_original=dd.diff_vs_orig.mean(), tasks_better=better, tasks_worse=worse,
                         tasks_significantly_better=sig_b, tasks_significantly_worse=sig_w,
                         sign_test_p=binomtest(better, better + worse).pvalue if better + worse else np.nan,
                         mean_diff_without_context=pos.diff_vs_orig.mean()))
    C = pd.DataFrame(rows).sort_values("mean_diff_vs_original", ascending=False)
    C.to_csv(TAB / "positives_model_comparison.tsv", sep="\t", index=False)
    return C


def heatmap(A):
    tasks = [t for t in ORDER if t in set(A.task)]
    rows = ["ORIGINAL"] + MODEL_ROWS
    V = np.full((len(rows), len(tasks)), np.nan); Sg = np.full(V.shape, "", object)
    for j, t in enumerate(tasks):
        d = A[A.task == t].set_index("model")
        orig = d.index[d.family == "original method"][0]
        V[0, j] = d.loc[orig, "est"]
        for i, m in enumerate(MODEL_ROWS, 1):
            if m in d.index:
                V[i, j] = d.loc[m, "est"]
                if m != "RandomSets_median" and pd.notna(d.loc[m, "diff_p"]) and d.loc[m, "diff_p"] < 0.05:
                    Sg[i, j] = "+" if d.loc[m, "diff_vs_orig"] > 0 else "-"
    fig, ax = plt.subplots(figsize=(0.95 * len(tasks) + 2.6, 0.36 * len(rows) + 1.8))
    im = ax.imshow(V - 0.5, cmap="PuOr_r", vmin=-0.25, vmax=0.25, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(tasks)):
            if np.isfinite(V[i, j]):
                ax.text(j, i, f"{V[i,j]:.2f}{Sg[i,j]}", ha="center", va="center", fontsize=6.8, color="black")
    ax.set_xticks(range(len(tasks))); ax.set_xticklabels([SHORT[t] for t in tasks], rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(["Original method"] + [r.replace("_", " ") for r in MODEL_ROWS])
    ax.axhline(0.5, color="black", lw=1.2); ax.axhline(14.5 + 0.0, color="grey", lw=0.6)
    ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("metric minus 0.5 (0 = chance)")
    ax.set_title("a   All models on every positive result (C-index or AUC)", loc="left", fontsize=9)
    save(fig, "fig_positives_heatmap", """
# All models on every positive result (heatmap)

Each column is one positive result and each row is one model. The number in a cell is the main score on held-out data: the C-index
(for survival: 0.5 = chance, 1 = perfect ordering of who lives longer) or the AUC (for mutation and response tasks: 0.5 = chance, 1 = perfect).
Purple-brown means above chance, orange means below chance. The first row is the original method of that result (its name is in the
bar-plot of the task and in the report). A plus or minus sign next to a number means the model differs from the original method with a
paired bootstrap p below 0.05 (plus = better, minus = worse); no sign means no clear difference. 'RandomSets median' is the median of
200 random pathway sets of the same size as the network shortlist. The two olaparib columns marked context are extra cohorts where the original signature was not
positive. Source table: tables/positives_all_models.tsv.
""")


def bars(A):
    for t in ORDER:
        d = A[A.task == t].copy()
        if d.empty: continue
        d = d[~d.model.isin(["Paper_top10"])] if t not in ("rectal_GSE87211_5FU",) else d
        d = d.sort_values("est")
        orig = d.model[d.family == "original method"].iloc[0]
        fig, ax = plt.subplots(figsize=(6.4, 0.27 * len(d) + 1.6))
        col = [OK["verm"] if m == orig else (OK["grey"] if f in ("simple baseline", "random pathway sets") else (OK["orange"] if f == "original-type method" else OK["blue"]))
               for m, f in zip(d.model, d.family)]
        y = np.arange(len(d))
        ax.barh(y, d.est - 0.5, left=0.5, color=col, height=0.7)
        ok = d.lo.notna()
        ax.errorbar(d.est[ok], y[ok.values], xerr=[(d.est - d.lo)[ok], (d.hi - d.est)[ok]], fmt="none", ecolor="black", lw=0.8, capsize=1.5)
        ax.axvline(0.5, color="black", lw=0.8)
        ov = d.est[d.model == orig].iloc[0]; ax.axvline(ov, color=OK["verm"], lw=0.8, ls="--")
        ax.set_yticks(y); ax.set_yticklabels([(m + " (original)") if m == orig else m.replace("_", " ") for m in d.model])
        lab = d.metric.iloc[0]
        ax.set_xlabel(f"{lab}   (0.5 = chance; bars: 95% CI)")
        ax.set_title(f"{t.replace('_', ' ')}\n{d.label.iloc[0]}", loc="left", fontsize=8.5)
        for yi, (_, r) in zip(y, d.iterrows()):
            if pd.notna(r.diff_p) and r.model != orig:
                ax.text(ax.get_xlim()[1], yi, f" d={r.diff_vs_orig:+.3f}, p={r.diff_p:.2g}", va="center", fontsize=6.3, clip_on=False)
        best = d[~d.model.isin(["RandomSets_median"])].sort_values("est").iloc[-1]
        cap = f"""
# {SHORT[t].replace(chr(10), ' ')}: every model against the original method

{d.label.iloc[0]}. Each bar is one model. The bar shows the main score on held-out data ({lab}; 0.5 = chance), the black line is its 95% confidence interval
(bootstrap over patients). The red bar is the original method ({ORIG_LABEL.get(orig, orig)}); the orange bars are other published-style methods of the same study;
blue bars are the models added here (trained with settings chosen only inside the training data); grey bars are simple baselines and the median of 200 random pathway sets (with its 2.5-97.5% range as the line).
The text at the right gives the paired difference to the original method and its bootstrap p-value (same patients for both).
Original method: {ov:.3f}. Highest bar: {best.model.replace('_', ' ')} at {best.est:.3f}. Source: tables/positives_all_models.tsv.
"""
        save(fig, f"fig_bar_{t}", cap)


if __name__ == "__main__":
    A = pd.read_csv(TAB / "positives_all_models.tsv", sep="\t")
    C = comparison(A)
    print(C.round(3).to_string())
    heatmap(A); bars(A)
