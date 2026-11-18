"""Collect bench pickles into tables: long benchmark table, per-model summary over pairs, Supp Table 4 (DL grid),
Bolis gene lists, random-set comparison, alignment comparison."""
import pickle

import numpy as np
import pandas as pd

import b3common as B
from bench import MODELS, OUT

LABEL = {"COAD_paper": "COAD organoids (paper data)", "BLCA_paper": "BLCA organoids (paper data)",
         "COADREAD_vdW": "CRC organoids", "COAD_isella": "CRC PDX (Isella)", "PAAD_tiriac": "PAAD org. (Tiriac)",
         "PAAD_shi": "PAAD org. (Shi)", "LIHC_ji": "LIHC org. (Ji)", "BLCA_current": "BLCA organoids",
         "OV_vias": "OV organoids", "SARC_alshihabi": "SARC organoids", "LUAD_cells": "LUAD cell lines",
         "PDX_BRCA": "BRCA PDX", "PDX_CRC": "CRC PDX (Gao)", "PDX_NSCLC": "NSCLC PDX", "PDX_PDAC": "PDAC PDX",
         "PDX_CM": "melanoma PDX"}


def pair_label(b, d):
    return f"{LABEL.get(b, b)} | {d.lower()}"


def load_all():
    res = {}
    for b, d in B.combos():
        f = OUT / f"{b}__{d}.pkl"
        if f.exists():
            res[(b, d)] = pickle.load(open(f, "rb"))
    return res


def long_table(res):
    rows = []
    for (b, d), r in res.items():
        for row in r["rows"]:
            rows.append(dict(bundle=b, drug=d, pair=pair_label(b, d), n_models=r["n_models"], cv=r["cv"],
                             n_proximal=r["n_proximal"], **row))
    L = pd.DataFrame(rows)
    L.to_csv(B.TABLES / "benchmark_long.tsv", sep="\t", index=False)
    return L


def model_summary(L):
    order = [m for m, _ in MODELS]
    fam = dict(MODELS)
    rows = []
    for m in order:
        s = L[L.model == m]
        if s.empty:
            continue
        def mci(v):
            v = v.dropna()
            if len(v) < 2:
                return np.nan, np.nan, np.nan
            return v.mean(), v.mean() - 1.96 * v.std(ddof=1) / np.sqrt(len(v)), v.mean() + 1.96 * v.std(ddof=1) / np.sqrt(len(v))
        r, rl, rh = mci(s["loo_rho"])
        c, cl, ch = mci(s["C_index"])
        lhr = np.log(s["HR"])
        h, hl, hh = mci(lhr)
        rows.append(dict(model=m, family=fam[m], n_pairs=int(s["C_index"].notna().sum()),
                         loo_rho_mean=r, loo_rho_lo=rl, loo_rho_hi=rh,
                         C_mean=c, C_lo=cl, C_hi=ch, HR_geomean=np.exp(h), HR_lo=np.exp(hl), HR_hi=np.exp(hh),
                         share_C_above_half=(s["C_index"] > 0.5).mean(),
                         n_lr_p05_expected=int(((s["lr_p"] < 0.05) & (s["split_HR"] > 1)).sum()),
                         n_lr_p05_opposite=int(((s["lr_p"] < 0.05) & (s["split_HR"] < 1)).sum()),
                         n_int_p05=int((s["int_p"] < 0.05).sum()),
                         median_loo_rho=s["loo_rho"].median(), median_C=s["C_index"].median()))
    S = pd.DataFrame(rows)
    S.to_csv(B.TABLES / "benchmark_model_summary.tsv", sep="\t", index=False)
    return S


def supp_table4(res):
    rows = []
    for (b, d), r in res.items():
        g = r["moli_grid"].copy()
        g.insert(0, "pair", pair_label(b, d))
        rows.append(g)
    G = pd.concat(rows)
    G.to_csv(B.TABLES / "SuppTable4_DL_grid_all_pairs.tsv", sep="\t", index=False)
    sel = G[G.selected].drop(columns="selected")
    sel.to_csv(B.TABLES / "SuppTable4_DL_selected.tsv", sep="\t", index=False)
    cols = ["batch_size", "nodes", "lr", "dropout", "weight_decay", "gamma", "epochs", "n_folds", "margin", "average_AUC"]
    md = ["# Supplementary Table 4. Deep-learning (MOLI-style) hyper-parameters and cross-validated AUC", "",
          "Layout follows the paper's source data (MOESM3 sheet 'Figure 3f'): one row per hyper-parameter setting;",
          "average_AUC = mean AUC over 3 stratified folds of the lab models (label = response above the median).",
          "The setting with the highest average AUC (bold) was used for the patient prediction.", ""]
    for key in [k for k in res if k[0] in ("COAD_paper", "BLCA_paper") and k[1] in ("FLUOROURACIL", "CISPLATIN")]:
        g = res[key]["moli_grid"]
        md += [f"## {pair_label(*key)} (all 32 settings)", "", "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for _, x in g.iterrows():
            cells = [f"{x[c]:.3f}" if c == "average_AUC" else str(x[c]) for c in cols]
            if x["selected"]:
                cells = [f"**{c}**" for c in cells]
            md.append("| " + " | ".join(cells) + " |")
        md.append("")
    md += ["## Selected setting for every pair", "", "| pair | " + " | ".join(cols) + " |", "|" + "---|" * (len(cols) + 1)]
    for _, x in sel.iterrows():
        md.append(f"| {x['pair']} | " + " | ".join(f"{x[c]:.3f}" if c == "average_AUC" else str(x[c]) for c in cols) + " |")
    open(B.TABLES / "SuppTable4_DL_hyperparameters.md", "w").write("\n".join(md) + "\n")


def extras(res):
    rows, bol = [], []
    for (b, d), r in res.items():
        rr = r["random_sets"]
        main = {x["model"]: x for x in r["rows"]}.get("Paper: single network pathway", {})
        if main:
            rows.append(dict(pair=pair_label(b, d), paper_C=main.get("C_index"), paper_lr_p=main.get("lr_p"),
                             paper_loo_rho=main.get("loo_rho"), random_C_median=rr["C_index"].median(),
                             share_random_C_ge_paper=(rr["C_index"] >= main.get("C_index", np.nan)).mean(),
                             share_random_lr_p_le_paper=(rr["lr_p"] <= main.get("lr_p", np.nan)).mean(),
                             share_random_loo_ge_paper=(rr["loo_rho"] >= main.get("loo_rho", np.nan)).mean()))
        gi = r["info"].get("Bolis correlation genes (ridge)", {})
        bol.append(dict(pair=pair_label(b, d), n_genes=gi.get("n_genes"), genes=",".join(gi.get("genes", [])[:60])))
    pd.DataFrame(rows).to_csv(B.TABLES / "random_sets_vs_paper_pathway.tsv", sep="\t", index=False)
    pd.DataFrame(bol).to_csv(B.TABLES / "bolis_selected_genes.tsv", sep="\t", index=False)


def alignment_table(L):
    pairs = [("Paper: single network pathway", "Paper: single network pathway + alignment"),
             ("Ridge: all network pathways", "Ridge: all network pathways + alignment"),
             ("Elastic net (all pathways)", "Elastic net + alignment"),
             ("Random forest (all pathways)", "Random forest + alignment"), ("XGBoost (all pathways)", "XGBoost + alignment")]
    rows = []
    for a, b in pairs:
        A = L[L.model == a].set_index("pair")
        Bm = L[L.model == b].set_index("pair")
        j = A.index.intersection(Bm.index)
        dC = (Bm.loc[j, "C_index"] - A.loc[j, "C_index"]).dropna()
        dR = (Bm.loc[j, "loo_rho"] - A.loc[j, "loo_rho"]).dropna()
        from scipy.stats import wilcoxon
        rows.append(dict(model=a, n_pairs=len(dC), mean_C_without=A.loc[j, "C_index"].mean(), mean_C_with=Bm.loc[j, "C_index"].mean(),
                         mean_delta_C=dC.mean(), wilcoxon_p_C=wilcoxon(dC).pvalue if len(dC) > 5 else np.nan,
                         mean_loo_without=A.loc[j, "loo_rho"].mean(), mean_loo_with=Bm.loc[j, "loo_rho"].mean(),
                         mean_delta_loo=dR.mean(), wilcoxon_p_loo=wilcoxon(dR).pvalue if len(dR) > 5 else np.nan))
    pd.DataFrame(rows).to_csv(B.TABLES / "alignment_with_vs_without.tsv", sep="\t", index=False)


def main():
    res = load_all()
    L = long_table(res)
    model_summary(L)
    supp_table4(res)
    extras(res)
    alignment_table(L)
    print(len(res), "pairs")
    return res, L


if __name__ == "__main__":
    main()
