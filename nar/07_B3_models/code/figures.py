"""All B3 figures (300 dpi PNG + PDF + plain-language caption files).

Usage: python figures.py [all|heatmap|bars|fig3|lc|align|random|main]
"""
import sys

import numpy as np
import pandas as pd

import b3common as B
from bench import MODELS
from tables import load_all, long_table, model_summary, pair_label

plt = B.plot_style()
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

FAMCOL = {"paper": B.OK["blue"], "paper comparison": B.OK["grey"], "modern": B.OK["verm"], "alignment": B.OK["green"]}
ORDER = [m for m, _ in MODELS]
FAM = dict(MODELS)
SHORT = {"Paper: single network pathway": "Paper: 1 network pathway", "Paper: top-2 network pathways": "Paper: top-2 pathways",
         "Ridge: all network pathways": "Ridge, all network pathways", "Random pathway sets (median of 100)": "Random pathway sets",
         "Fingerprint MLP: cell lines only (zero-shot)": "Fingerprint MLP, cell lines only",
         "Fingerprint MLP: fine-tuned on lab models": "Fingerprint MLP, fine-tuned",
         "Deep learning, MOLI-style (genes)": "Deep learning (MOLI-style)",
         "Paper: single network pathway + alignment": "1 network pathway + align",
         "Ridge: all network pathways + alignment": "Ridge network + align"}


def sh(m):
    return SHORT.get(m, m)


def letter(ax, s, x=-0.12, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def km(ax, T, E, risk, title):
    from obd import survival as SV
    hi = risk > np.median(risk)
    for g, lab, col in ((~hi, "predicted responders", B.OK["blue"]), (hi, "predicted non-responders", B.OK["verm"])):
        x, y = SV.km_curve(T[g], E[g])
        ax.step(x, y, where="post", color=col, lw=1.4, label=f"{lab} (n={g.sum()})")
    ms, _ = B.median_split(T, E, risk)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("months")
    ax.set_ylabel("overall survival")
    ax.set_title(title, fontsize=8.5)
    hr = ms["split_HR"]
    hrs = (f"HR = {hr:.2f} ({ms['split_HR_lo']:.2f}-{ms['split_HR_hi']:.2f})" if np.isfinite(hr) and hr < 1000
           else "HR not estimable (no deaths in one group)")
    ax.text(0.03, 0.04, f"log-rank p = {ms['lr_p']:.3g}\n{hrs}",
            transform=ax.transAxes, fontsize=7.5)
    ax.legend(fontsize=6.5, loc="upper right", frameon=False)


# --------------------------------------------------------------------------------------------------- heatmaps
def task_matrix(S, cells):
    cs = cells.pivot(index="model", columns="task", values="mean_rho") if cells is not None else None
    M = pd.DataFrame(index=[m for m in ORDER if m in set(S.model)])
    St = S.set_index("model")
    M["Lab models: held-out Spearman"] = St["loo_rho_mean"]
    M["Patients: C-index - 0.5"] = St["C_mean"] - 0.5
    M["Patients: log HR per SD"] = np.log(St["HR_geomean"])
    if cs is not None:
        mp = {"Paper: single network pathway": "Paper: single network pathway", "Paper: top-2 network pathways": "Paper: top-2 network pathways",
              "Ridge: all network pathways": "Ridge: all network pathways", "Elastic net (all pathways)": "Elastic net (pathways + fingerprint)",
              "XGBoost (all pathways)": "XGBoost (pathways + fingerprint)", "Random forest (all pathways)": "Random forest (pathways + fingerprint)",
              "Fingerprint MLP: cell lines only (zero-shot)": "Fingerprint MLP (cell lines)",
              "Fingerprint MLP: fine-tuned on lab models": "Fingerprint MLP (cell lines)", "GCN on STRING (genes)": "GCN on STRING (per drug)"}
        for task in ("cell-blind", "drug-blind"):
            col = f"Cell lines: {task} Spearman"
            M[col] = [cs.loc[mp[m], task] if m in mp and mp[m] in cs.index and task in cs.columns else np.nan for m in M.index]
        # pooled EN with fingerprint is the drug-blind version of the elastic net
        if "Elastic net (pathways + fingerprint)" in cs.index:
            M.loc["Elastic net (all pathways)", "Cell lines: drug-blind Spearman"] = cs.loc["Elastic net (pathways + fingerprint)", "drug-blind"]
    return M


def heatmap(res, L, S, cells):
    M = task_matrix(S, cells)
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    V = M.values.astype(float)
    norm = TwoSlopeNorm(vmin=-0.25, vcenter=0, vmax=0.25)
    im = ax.imshow(np.clip(V, -0.25, 0.25), cmap="PuOr", norm=norm, aspect="auto")
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            ax.text(j, i, "n/a" if np.isnan(V[i, j]) else f"{V[i, j]:+.3f}", ha="center", va="center", fontsize=7,
                    color="#555555" if np.isnan(V[i, j]) else "black")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([c.replace(": ", ":\n") for c in M.columns], fontsize=7.5)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels([sh(m) for m in M.index], fontsize=7.5)
    for t, m in zip(ax.get_yticklabels(), M.index):
        t.set_color(FAMCOL[FAM[m]])
    ax.xaxis.tick_top()
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("score (0 = no signal; purple = expected direction)", fontsize=7.5)
    ax.spines[:].set_visible(False)
    M.to_csv(B.TABLES / "benchmark_heatmap_values.tsv", sep="\t")
    B.savefig(fig, "B3_benchmark_heatmap", f"""
**Benchmark heatmap: every model (rows) on every task (columns).** Numbers are averages over the {L.pair.nunique()} lab-model set x drug pairs
that have at least 15 drug-treated TCGA patients (cell-line columns: averages over drug x screen units of GDSC and CTRPv2).
Columns: (1) Spearman correlation between predicted and measured drug response in lab models left out of training
(leave-one-out for sets of up to 40 models, else 10-fold); (2) C-index minus 0.5 in treated TCGA patients (0 = coin toss;
positive = patients predicted resistant died sooner); (3) log hazard ratio per standard deviation of the predicted resistance
(positive = expected direction); (4) cell-blind and (5) drug-blind cross-validation on cell lines (Spearman). Row label colour:
blue = the paper's network method, grey = the paper's comparison methods, red = modern models, green = with the
lab-model-to-tumour alignment of nar/03. "n/a" = the model cannot be run on that task (for example, the paper's method needs
response data for the drug, so it has no drug-blind version). Values are clipped at +/-0.25 for colour only.""")
    # full model x pair heatmaps (supplement)
    for metric, lab, stem, lim in (("C_index", "C-index in treated patients", "B3_supp_heatmap_Cindex", 0.25),
                                   ("loo_rho", "held-out Spearman in lab models", "B3_supp_heatmap_LOO", 0.8)):
        P = L.pivot_table(index="model", columns="pair", values=metric).reindex([m for m in ORDER if m in set(L.model)])
        pairs = [pair_label(b, d) for b, d in B.combos() if pair_label(b, d) in P.columns]
        P = P[pairs]
        V = P.values - (0.5 if metric == "C_index" else 0)
        fig, ax = plt.subplots(figsize=(14, 6.8))
        im = ax.imshow(np.clip(V, -lim, lim), cmap="PuOr", norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), aspect="auto")
        for i in range(V.shape[0]):
            for j in range(V.shape[1]):
                if np.isfinite(P.values[i, j]):
                    ax.text(j, i, f"{P.values[i, j]:.2f}", ha="center", va="center", fontsize=4.6)
        ax.set_xticks(range(len(pairs)))
        ax.set_xticklabels(pairs, rotation=90, fontsize=6)
        ax.set_yticks(range(P.shape[0]))
        ax.set_yticklabels([sh(m) for m in P.index], fontsize=6.5)
        for t, m in zip(ax.get_yticklabels(), P.index):
            t.set_color(FAMCOL[FAM[m]])
        fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01).set_label(lab + (" - 0.5" if metric == "C_index" else ""), fontsize=7)
        B.savefig(fig, stem, f"""
**{lab.capitalize()} for every model (rows) and every lab-model set x drug pair (columns).** Printed numbers are the raw values;
colour shows the distance from no signal ({'0.5' if metric == 'C_index' else '0'}). Purple = expected direction (resistant-predicted
patients die sooner / predictions rank lab models correctly), orange = the opposite. Row label colours as in the main heatmap.""")
        plt.close("all")


# --------------------------------------------------------------------------------------------------- bars
def bars(S, cells):
    fig, axs = plt.subplots(1, 5, figsize=(16, 6.2), sharey=True)
    idx = [m for m in ORDER if m in set(S.model)]
    St = S.set_index("model").loc[idx]
    yy = np.arange(len(idx))[::-1]
    cols = [FAMCOL[FAM[m]] for m in idx]
    spec = [("loo_rho_mean", "loo_rho_lo", "loo_rho_hi", "lab models: held-out Spearman", 0),
            ("C_mean", "C_lo", "C_hi", "patients: C-index", 0.5),
            ("HR_geomean", "HR_lo", "HR_hi", "patients: HR per SD (geometric mean)", 1)]
    for ax, (m, lo, hi, lab, ref) in zip(axs, spec):
        v = St[m].values
        ax.barh(yy, v - ref, left=ref, color=cols, alpha=0.85, height=0.7)
        ax.errorbar(v, yy, xerr=[v - St[lo].values, St[hi].values - v], fmt="none", ecolor="black", lw=0.8, capsize=2)
        ax.axvline(ref, color="black", lw=0.7)
        ax.set_title(lab, fontsize=8.5)
    if cells is not None:
        for ax, task in zip(axs[3:], ("cell-blind", "drug-blind")):
            c = cells[cells.task == task].set_index("model")
            M = task_matrix(S, cells)[f"Cell lines: {task} Spearman"].reindex(idx)
            v = M.values
            ax.barh(yy, np.nan_to_num(v), color=cols, alpha=0.85, height=0.7)
            # CI from the pooled cell-line summary
            mp = {}
            for m in idx:
                for cm in c.index:
                    if np.isfinite(M[m]) and abs(c.loc[cm, "mean_rho"] - M[m]) < 1e-12:
                        mp[m] = cm
            lo = np.array([c.loc[mp[m], "lo"] if m in mp else np.nan for m in idx])
            hi = np.array([c.loc[mp[m], "hi"] if m in mp else np.nan for m in idx])
            ok = np.isfinite(v)
            ax.errorbar(v[ok], yy[ok], xerr=[v[ok] - lo[ok], hi[ok] - v[ok]], fmt="none", ecolor="black", lw=0.8, capsize=2)
            for y_, vv in zip(yy, v):
                if not np.isfinite(vv):
                    ax.text(0.01, y_, "n/a", va="center", fontsize=6.5, color="#555555")
            ax.axvline(0, color="black", lw=0.7)
            ax.set_title(f"cell lines: {task} Spearman", fontsize=8.5)
    axs[0].set_yticks(yy)
    axs[0].set_yticklabels([sh(m) for m in idx], fontsize=7.5)
    for t, m in zip(axs[0].get_yticklabels(), idx):
        t.set_color(FAMCOL[FAM[m]])
    for i, ax in enumerate(axs):
        letter(ax, "abcde"[i], x=-0.05)
    B.savefig(fig, "B3_task_bars", f"""
**Each task separately, with 95% confidence intervals.** (a) Held-out Spearman correlation in lab models, mean over the
{int(S.n_pairs.max())} set x drug pairs (error bar = 95% CI of the mean across pairs). (b) C-index in drug-treated TCGA patients
(0.5 = coin toss). (c) Hazard ratio per standard deviation of predicted resistance, geometric mean over pairs (1 = no link; above 1 =
expected direction). (d, e) Cell-line cross-validation (GDSC + CTRPv2): cells never seen in training (d) or drugs never seen in
training (e); mean Spearman over drug x screen units with 95% CI. Bar colours: blue = paper's network method, grey = paper's comparison
methods, red = modern models, green = with alignment. Pairs share some patients (for example the paper-data and current-data
bladder sets), so the CIs across pairs are approximate.""")
    plt.close("all")


# --------------------------------------------------------------------------------------------------- Fig 3 replicas
FIG3 = [("Paper: single network pathway", "network-proximal pathway (paper)"), ("Whole transcriptome (ridge)", "b  whole transcriptome"),
        ("All pathways (ridge)", "c  all pathways"), ("First-neighbour genes (ridge)", "d  first network neighbours"),
        ("Bolis correlation genes (ridge)", "e  Bolis correlation genes"), ("Deep learning, MOLI-style (genes)", "f  deep learning (MOLI-style)")]


def fig3(res):
    for (b, d), r in res.items():
        fig, axs = plt.subplots(2, 3, figsize=(11, 6.6))
        for ax, (m, title), lt in zip(axs.ravel(), FIG3, "abcdef"):
            if m not in r["risk"]:
                ax.text(0.5, 0.5, "not run", ha="center")
                continue
            km(ax, r["T"], r["E"], np.asarray(r["risk"][m]), title.split("  ")[-1])
            letter(ax, lt)
        fig.suptitle(f"Paper Fig. 3 replica: {pair_label(b, d)} (treated TCGA patients, n = {len(r['T'])}, deaths = {int(r['E'].sum())})",
                     fontsize=9.5)
        fig.tight_layout()
        info = r["info"]
        bol = info.get("Bolis correlation genes (ridge)", {})
        pw = info.get("Paper: single network pathway", {}).get("pathways", ["-"])
        g = r["moli_grid"]
        sel = g[g.selected].iloc[0]
        B.savefig(fig, f"fig3_replica/Fig3_{b}_{d}", f"""
**Paper Fig. 3 (b-f) replica for {pair_label(b, d)}.** Each panel trains one model on the {r['n_models']} lab models and applies it to the
{len(r['T'])} TCGA patients given {d.lower()} ({int(r['E'].sum())} deaths). Patients are split at the median predicted resistance;
blue = predicted responders, red = predicted non-responders; log-rank p and hazard ratio (non-responders vs responders, 95% CI).
(a) the paper's network method: the network-proximal pathway with the largest ridge weight ({B.short(pw[0], 70)}).
(b) ridge on every gene. (c) ridge on every Reactome pathway. (d) ridge on the drug targets and their direct STRING neighbours.
(e) ridge on genes picked by the Bolis et al. correlation filter (Spearman p < 0.05 with response in all 10 random halves of
the lab models; {bol.get('n_genes', '?')} genes). (f) MOLI-style deep network (encoder + classifier, cross-entropy + triplet loss),
lab models split at the median response; best of 32 hyper-parameter settings by 3-fold CV AUC ({sel['average_AUC']:.2f}; Supp Table 4).""")
        plt.close("all")


# --------------------------------------------------------------------------------------------------- learning curves
def learning_curves():
    f1, f2 = B.TABLES / "learning_curve_lab_models.tsv", B.TABLES / "learning_curve_cell_lines.tsv"
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    cmap = plt.get_cmap("tab10")
    if f1.exists():
        D = pd.read_csv(f1, sep="\t")
        # per pair: mean over folds and draws, then mean over pairs (x = fraction of the available training models)
        g = D.groupby(["model", "bundle", "drug", "frac"])[["rho_heldout", "C_index"]].mean().reset_index()
        for i, m in enumerate([m for m in ORDER if m in set(g.model)]):
            s = g[g.model == m].groupby("frac")[["rho_heldout", "C_index"]].agg(["mean", "sem"])
            for ax, col in ((axs[0], "rho_heldout"), (axs[1], "C_index")):
                ax.errorbar(s.index, s[(col, "mean")], yerr=1.96 * s[(col, "sem")], marker="o", ms=3, lw=1.2, capsize=2,
                            color=cmap(i % 10), label=sh(m))
        axs[0].set_ylabel("Spearman on held-out lab models")
        axs[1].set_ylabel("C-index in treated patients")
        axs[1].axhline(0.5, color="black", lw=0.6)
        axs[0].axhline(0, color="black", lw=0.6)
        for ax in axs[:2]:
            ax.set_xlabel("fraction of training lab models used")
        axs[0].legend(fontsize=6.5, frameon=False)
        npairs = g.groupby(["bundle", "drug"]).ngroups
    if f2.exists():
        D2 = pd.read_csv(f2, sep="\t")
        for i, m in enumerate(D2.model.unique()):
            s = D2[D2.model == m].groupby("n_cells")["mean_rho"].agg(["mean", "min", "max"])
            axs[2].plot(s.index, s["mean"], marker="o", ms=3, label=m, color=[B.OK["verm"], B.OK["purple"], B.OK["sky"]][i % 3])
            axs[2].fill_between(s.index, s["min"], s["max"], alpha=0.2, color=[B.OK["verm"], B.OK["purple"], B.OK["sky"]][i % 3])
        axs[2].set_xscale("log")
        axs[2].set_xlabel("number of training cell lines")
        axs[2].set_ylabel("cell-blind Spearman (mean over drugs)")
        axs[2].legend(fontsize=6.5, frameon=False)
    for i, ax in enumerate(axs):
        letter(ax, "abc"[i])
    B.savefig(fig, "B3_learning_curves", f"""
**Learning curves: how performance changes with the amount of training data.** (a, b) The larger lab-model sets
({npairs if f1.exists() else '?'} set x drug pairs, e.g. 192 colorectal PDX for cetuximab, 61 liver organoids for sorafenib, PDX and
lung cell lines): 20% of the models are held out (5 times), and each model is trained on 25%, 50%, 75% or 100% of the rest
(one random draw per outer fold). (a) Spearman correlation on the held-out lab models; (b) C-index in the drug-treated TCGA patients.
Points = mean over pairs, bars = 95% CI. (c) Cell lines: cell-blind Spearman of the pooled models (with drug fingerprint) when
trained on 5% to 100% of the training cell lines (one cell-blind fold; band = range over 3 draws).""")
    plt.close("all")


# --------------------------------------------------------------------------------------------------- alignment + random sets
def alignment_fig(L):
    pairs = [("Paper: single network pathway", "Paper: single network pathway + alignment"),
             ("Ridge: all network pathways", "Ridge: all network pathways + alignment"),
             ("Elastic net (all pathways)", "Elastic net + alignment"), ("Random forest (all pathways)", "Random forest + alignment"),
             ("XGBoost (all pathways)", "XGBoost + alignment")]
    fig, axs = plt.subplots(2, 5, figsize=(15, 6))
    for j, (a, b) in enumerate(pairs):
        A, Bm = L[L.model == a].set_index("pair"), L[L.model == b].set_index("pair")
        k = A.index.intersection(Bm.index)
        for i, (met, lab, ref) in enumerate((("C_index", "C-index", 0.5), ("loo_rho", "held-out Spearman", 0))):
            ax = axs[i, j]
            ax.scatter(A.loc[k, met], Bm.loc[k, met], s=12, color=B.OK["green"], alpha=0.8)
            lim = [np.nanmin([A.loc[k, met].min(), Bm.loc[k, met].min()]) - 0.03, np.nanmax([A.loc[k, met].max(), Bm.loc[k, met].max()]) + 0.03]
            ax.plot(lim, lim, color="black", lw=0.6)
            ax.axhline(ref, color="grey", lw=0.5, ls=":")
            ax.axvline(ref, color="grey", lw=0.5, ls=":")
            ax.set_xlabel(f"{lab}, no alignment", fontsize=7.5)
            ax.set_ylabel(f"{lab}, with alignment", fontsize=7.5)
            if i == 0:
                ax.set_title(sh(a), fontsize=8)
            letter(ax, "abcdefghij"[i * 5 + j], x=-0.25)
    fig.tight_layout()
    B.savefig(fig, "B3_alignment_with_without", """
**With versus without the lab-model-to-tumour alignment (B3 item h).** Each dot is one set x drug pair. x = model trained on the
plain z-scored lab-model pathway scores; y = the same model trained after the alignment of nar/03 (remove lab-model-specific
principal axes, then match the tumour covariance with CORAL). Top row: C-index in drug-treated TCGA patients. Bottom row: Spearman on
held-out lab models. Dots above the diagonal = alignment helped. Numbers: tables/alignment_with_vs_without.tsv.""")
    plt.close("all")


def random_fig(res):
    keys = [k for k in res if k[0] in ("COAD_paper", "BLCA_paper")]
    keys = [k for k in keys if k[1] in ("FLUOROURACIL", "CISPLATIN")] + [k for k in res if k[0] not in ("COAD_paper", "BLCA_paper")][:0]
    fig, axs = plt.subplots(len(keys), 3, figsize=(11, 3.2 * len(keys)), squeeze=False)
    for i, k in enumerate(keys):
        r = res[k]
        rr = r["random_sets"]
        main = {x["model"]: x for x in r["rows"]}["Paper: single network pathway"]
        for j, (col, lab) in enumerate((("C_index", "C-index (treated patients)"), ("lr_p", "log-rank p (median split)"),
                                        ("loo_rho", "held-out Spearman (lab models)"))):
            ax = axs[i, j]
            v = rr[col].dropna()
            if col == "lr_p":
                v = -np.log10(v)
            ax.hist(v, bins=20, color=B.OK["grey"], alpha=0.8)
            obs = main[col] if col != "lr_p" else -np.log10(main[col])
            ax.axvline(obs, color=B.OK["blue"], lw=2)
            share = (rr[col] >= main[col]).mean() if col != "lr_p" else (rr[col] <= main[col]).mean()
            ax.set_xlabel(("-log10 " if col == "lr_p" else "") + lab)
            ax.set_title(f"{pair_label(*k)}\nrandom sets as good or better: {share:.0%}", fontsize=8)
            letter(ax, "abcdef"[i * 3 + j])
    fig.tight_layout()
    B.savefig(fig, "B3_random_pathway_sets", """
**Paper network pathway (blue line) versus 100 random pathway sets of the same size (grey histogram).** For each random set the
same rank-1 ridge procedure is run. Columns: C-index in drug-treated TCGA patients, -log10 log-rank p of the median split, and
Spearman on held-out lab models. The title gives the share of random sets that do at least as well as the network set.
All pairs: tables/random_sets_vs_paper_pathway.tsv.""")
    plt.close("all")


def main_fig(res, L, S, cells):
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1], hspace=0.45, wspace=0.55)
    # a: compact heatmap
    ax = fig.add_subplot(gs[0, :2])
    M = task_matrix(S, cells)
    V = M.values.astype(float)
    ax.imshow(np.clip(V, -0.25, 0.25), cmap="PuOr", norm=TwoSlopeNorm(vmin=-0.25, vcenter=0, vmax=0.25), aspect="auto")
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            ax.text(j, i, "n/a" if np.isnan(V[i, j]) else f"{V[i, j]:+.2f}", ha="center", va="center", fontsize=6.3)
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([c.split(": ")[1] if ": " in c else c for c in M.columns], fontsize=6.8, rotation=25, ha="left")
    ax.xaxis.tick_top()
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels([sh(m) for m in M.index], fontsize=6.8)
    for t, m in zip(ax.get_yticklabels(), M.index):
        t.set_color(FAMCOL[FAM[m]])
    letter(ax, "a", x=-0.45)
    # b: C-index bars
    idx = list(M.index)
    St = S.set_index("model").loc[idx]
    yy = np.arange(len(idx))[::-1]
    ax = fig.add_subplot(gs[0, 2])
    v = St["C_mean"].values
    ax.barh(yy, v - 0.5, left=0.5, color=[FAMCOL[FAM[m]] for m in idx], height=0.7)
    ax.errorbar(v, yy, xerr=[v - St["C_lo"], St["C_hi"] - v], fmt="none", ecolor="black", lw=0.7, capsize=1.5)
    ax.axvline(0.5, color="black", lw=0.6)
    ax.set_yticks([])
    ax.set_xlabel("C-index, treated patients")
    letter(ax, "b")
    ax = fig.add_subplot(gs[0, 3])
    v = St["loo_rho_mean"].values
    ax.barh(yy, v, color=[FAMCOL[FAM[m]] for m in idx], height=0.7)
    ax.errorbar(v, yy, xerr=[v - St["loo_rho_lo"], St["loo_rho_hi"] - v], fmt="none", ecolor="black", lw=0.7, capsize=1.5)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks([])
    ax.set_xlabel("held-out Spearman, lab models")
    letter(ax, "c")
    # d-g: KM headline pairs
    k1, k2 = ("COAD_paper", "FLUOROURACIL"), ("BLCA_paper", "CISPLATIN")
    best = S[S.family == "modern"].sort_values("C_mean").model.iloc[-1] if (S.family == "modern").any() else None
    panels = [(k1, "Paper: single network pathway"), (k1, best), (k2, "Paper: single network pathway"), (k2, best)]
    for j, ((k, m), lt) in enumerate(zip(panels, "defg")):
        ax = fig.add_subplot(gs[1, j])
        if k in res and m in res[k]["risk"]:
            r = res[k]
            km(ax, r["T"], r["E"], np.asarray(r["risk"][m]), f"{pair_label(*k)}\n{sh(m)}")
        letter(ax, lt, x=-0.3, y=1.12)
    B.savefig(fig, "B3_main_figure", f"""
**Main figure: the paper's network method next to its comparison methods and modern models, on identical splits.**
(a) Summary heatmap (rows = models, columns = tasks; numbers as in the full benchmark heatmap). Lab-model task = held-out Spearman;
patient tasks = C-index minus 0.5 and log HR per SD in drug-treated TCGA patients; cell-line tasks = cell-blind and drug-blind CV.
Purple = expected direction. (b) Mean C-index with 95% CI over the {L.pair.nunique()} set x drug pairs. (c) Mean held-out Spearman
in lab models with 95% CI. (d-g) Kaplan-Meier curves for the paper's two headline pairs (colorectal organoids -> 5-FU-treated TCGA-COAD
patients; bladder organoids -> cisplatin-treated TCGA-BLCA), paper's single network pathway versus the modern model with the
highest mean C-index across pairs ({sh(best) if best else '-'}). Row colours: blue = paper, grey = paper's comparison methods,
red = modern, green = with alignment.""")
    plt.close("all")


def main(what="all"):
    res = load_all()
    L = long_table(res)
    S = model_summary(L)
    fc = B.TABLES / "cellline_cv_summary.tsv"
    cells = pd.read_csv(fc, sep="\t") if fc.exists() else None
    (B.FIGS / "fig3_replica").mkdir(exist_ok=True)
    if what in ("all", "heatmap"):
        heatmap(res, L, S, cells)
    if what in ("all", "bars"):
        bars(S, cells)
    if what in ("all", "fig3"):
        fig3(res)
    if what in ("all", "lc"):
        learning_curves()
    if what in ("all", "align"):
        alignment_fig(L)
    if what in ("all", "random"):
        random_fig(res)
    if what in ("all", "main"):
        main_fig(res, L, S, cells)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
