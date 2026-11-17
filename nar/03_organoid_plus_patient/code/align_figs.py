"""Alignment of lab models (organoids, PDX, cell lines) with patient tumours, per cancer.

Three stages: (1) raw rank scores; (2) within-cohort / within-set standardisation; (3) stage 2 plus
the Celligner-like step (contrastive removal of model-only axes, then CORAL to the tumour covariance).
Figures: figures/alignment_<cancer>.png/.pdf (PCA and UMAP, before/after); table tables/alignment_metrics.tsv.
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors

import common as C
import plotstyle as PS
from run_joint import unit_data
from units import UNITS

CANCERS = ["colorectal", "bladder", "pancreas", "liver", "ovarian", "lung", "sarcoma"]


def cancer_sets(cancer):
    us = [u for u in UNITS if u["cancer"] == cancer]
    sets = []
    for u in us:
        for name, drug, lin in u["lab"]:
            if (name, lin) not in sets:
                sets.append((name, lin))
        for name, lin in u.get("align_extra", []):
            if (name, lin) not in sets:
                sets.append((name, lin))
    # tumours: every patient of every unit of this cancer (union)
    P, S = C.patients()
    idx, cohort = [], []
    for u in us:
        D, ZD, Zall = unit_data(u)
        for i in Zall.index:
            if i not in idx:
                idx.append(i)
    d = P.loc[idx]
    return sets, d


def stages(cancer):
    sets, d = cancer_sets(cancer)
    P, S = C.patients()
    raw_p = S.loc[d.index].fillna(S.loc[d.index].mean()).fillna(0)
    z_p = C.zscore_within(S.loc[d.index], d["cohort"]).fillna(0)
    cols = list(S.columns)
    raw_m, z_m, a_m, lab = [], [], [], []
    rng = np.random.default_rng(C.SEED)
    for name, lin in sets:
        sc = C.lab_expression_scores(name, lin).reindex(columns=cols)
        if len(sc) > 300:                                   # keep the figure readable: random 300 models
            sc = sc.iloc[np.sort(rng.choice(len(sc), 300, replace=False))]
        raw_m.append(sc.fillna(raw_p.mean()))
        zm = C.zscore_within(sc, pd.Series("all", index=sc.index)).fillna(0)
        z_m.append(zm)
        a, nd = C.align_models(zm.values, z_p.values)
        a_m.append(pd.DataFrame(a, index=zm.index, columns=cols))
        lab += [PS.SETNAME.get(name, name) + (f" ({lin})" if lin else "")] * len(sc)
    out = {}
    for key, Xp, Xm in (("raw", raw_p, raw_m), ("standardised", z_p, z_m), ("aligned", z_p, a_m)):
        X = np.vstack([Xp.values] + [m.values for m in Xm])
        out[key] = X
    labels = np.array(["tumours"] * len(d) + lab)
    cohorts = np.r_[d["cohort"].values, np.array(lab)]
    return out, labels, cohorts, len(d)


def metrics(X, labels, n_p, seed=C.SEED):
    y = (labels != "tumours").astype(int)
    Z = PCA(20, random_state=seed).fit_transform(X)
    k = min(5, int(y.sum()))
    pr = cross_val_predict(LogisticRegression(max_iter=3000, C=1.0, class_weight="balanced"), Z, y,
                           cv=StratifiedKFold(k, shuffle=True, random_state=seed), method="decision_function")
    auc = roc_auc_score(y, pr)
    nn = NearestNeighbors(n_neighbors=16).fit(Z)
    _, ind = nn.kneighbors(Z[y == 1])
    frac = (y[ind[:, 1:]] == 0).mean()
    return auc, frac, 1 - y.mean()


def main():
    import matplotlib.pyplot as plt
    import umap
    PS.setup()
    rows = []
    for cancer in CANCERS:
        st, labels, cohorts, n_p = stages(cancer)
        fig, ax = plt.subplots(2, 3, figsize=(11, 7.2))
        coords = {}
        for j, key in enumerate(("raw", "standardised", "aligned")):
            X = st[key]
            auc, frac, expct = metrics(X, labels, n_p)
            rows.append(dict(cancer=cancer, stage=key, n_tumours=n_p, n_models=int((labels != "tumours").sum()),
                             model_vs_tumour_AUROC=auc, tumour_share_of_model_neighbours=frac, tumour_share_expected=expct))
            pc = PCA(2, random_state=C.SEED).fit(X)
            Y = pc.transform(X)
            U = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=C.SEED, n_jobs=1).fit_transform(PCA(30, random_state=C.SEED).fit_transform(X))
            coords[key] = (Y, U, pc.explained_variance_ratio_)
            for i, (Q, nm) in enumerate(((Y, "PC"), (U, "UMAP"))):
                a = ax[i, j]
                t = labels == "tumours"
                a.scatter(Q[t, 0], Q[t, 1], s=3, c=PS.GREY, alpha=0.35, lw=0, label=f"tumours (n={t.sum()})", rasterized=True)
                for k, s in enumerate(pd.unique(labels[~t])):
                    m = labels == s
                    a.scatter(Q[m, 0], Q[m, 1], s=14, c=PS.CB[k % len(PS.CB)], edgecolor="k", lw=0.3, label=f"{s} (n={m.sum()})", zorder=3)
                if nm == "PC":
                    ev = pc.explained_variance_ratio_
                    a.set_xlabel(f"PC1 ({ev[0]:.0%})"); a.set_ylabel(f"PC2 ({ev[1]:.0%})")
                else:
                    a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
                a.set_xticks([]); a.set_yticks([])
                title = {"raw": "1. Raw pathway scores", "standardised": "2. Standardised within cohort / set",
                         "aligned": "3. Standardised + aligned"}[key]
                a.set_title(f"{title}\nmodel-vs-tumour AUROC {auc:.2f}" if i == 0 else title, fontsize=9)
        ax[0, 0].legend(fontsize=7, loc="best", frameon=True, markerscale=1.5)
        PS.letters(ax.ravel())
        fig.suptitle(f"{cancer.capitalize()}: lab models and tumours on a common scale", fontsize=11)
        fig.tight_layout()
        PS.save(fig, C.FIGS / f"alignment_{cancer}",
                caption=(f"**Alignment of lab models with patient tumours ({cancer}).** Each dot is one sample. Grey dots are "
                         f"patient tumours from all cohorts used for {cancer}; coloured dots are lab models (organoids, PDX or cell lines). "
                         "Top row (a-c): principal component analysis (PCA, the two directions of largest spread). Bottom row (d-f): UMAP "
                         "(a map that keeps close neighbours close). Column 1: raw pathway scores - lab models sit apart from tumours and "
                         "cohorts measured on different platforms also sit apart. Column 2: every pathway standardised inside each cohort "
                         "and inside each lab set (mean 0, SD 1). Column 3: standardised and then aligned - model-only directions are "
                         "removed (directions with much more spread in lab models than in tumours, for example missing immune and stromal "
                         "cells) and the lab-model covariance is matched to the tumour covariance (CORAL). The title number "
                         "'model-vs-tumour AUROC' is how well a classifier can tell lab models from tumours (1 = perfectly, 0.5 = cannot); "
                         "lower means better mixing."))
        pd.to_pickle(dict(coords=coords, labels=labels, cohorts=cohorts), C.CACHE / f"align_coords_{cancer}.pkl")
        print(cancer, "done", flush=True)
    pd.DataFrame(rows).to_csv(C.TABLES / "alignment_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
