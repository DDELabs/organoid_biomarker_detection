"""G3 (+ G5a organoid part): BH3-only pathway score vs BH3-mimetic response in the 19 van de Wetering 2015 organoids.

Scores: (1) 'current' = raw ssGSEA ES recomputed here on the curated GSE64392 expression (MSigDB v6.1 Reactome,
sets with >= 5 genes in the organoid gene universe); (2) 'authors' = the authors' committed organoid NES
(python/results/COAD/organoid/...normalized.es.txt).
Drugs: navitoclax (BCL-2/BCL-XL/BCL-W inhibitor), obatoclax (pan-BCL-2 family inhibitor, also off-target effects),
5-FU and oxaliplatin (contrast). Response = AUC (primary; lower = more sensitive) and ln IC50.
Nulls: every other Reactome pathway; 2,000 random gene sets with the same number of genes as the BH3-only set.
"""
import sys

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats

from gcommon import (BH3, FIG, OI, REPO, SEED, TAB, VDW, panel, savefig, spearman_ci, ssgsea_fast, style, fmt_p,
                     version_rows)

sys.path.insert(0, str(REPO / "robust"))
from obd import cohorts as C  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

DRUGS = ["NAVITOCLAX", "OBATOCLAX", "FLUOROURACIL", "OXALIPLATIN"]
LAB = {"NAVITOCLAX": "navitoclax", "OBATOCLAX": "obatoclax", "FLUOROURACIL": "5-FU", "OXALIPLATIN": "oxaliplatin"}
N_RAND = 2000


def load():
    e = C.to_canonical(pd.read_csv(VDW / "expression.tsv.gz", sep="\t", index_col=0))
    gs = gene_sets(("REACTOME",))
    gs = {k: [g for g in v if g in set(e.index)] for k, v in gs.items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    cur = ssgsea_fast(e, gs)
    a = pd.read_csv(REPO / "python/results/COAD/organoid/reactome_ssgsea_result/gseapy.samples.normalized.es.txt",
                    sep="\t", comment="#", index_col=0)
    a.columns = [c.upper().replace("P19A", "P19TA").replace("P19B", "P19TB").replace("P24A", "P24TA")
                 .replace("P24B", "P24TB") for c in a.columns]
    r = pd.read_csv(VDW / "response.tsv", sep="\t")
    auc = r.pivot_table(index="sample", columns="drug", values="response", aggfunc="median")
    ic = r.pivot_table(index="sample", columns="drug", values="ln_ic50", aggfunc="median")
    return e, gs, cur, a, auc, ic


def main():
    e, gs, cur, auth, auc, ic = load()
    S = {"current": cur, "authors": auth}
    M = {"AUC": auc, "ln IC50": ic}
    rows = []
    for sname, sc in S.items():
        for mname, R in M.items():
            samp = [s for s in R.index if s in sc.columns]
            for d in DRUGS:
                y = R.loc[samp, d]
                st = spearman_ci(sc.loc[BH3, samp], y)
                # pathway null: every pathway's rho with this drug
                allr = np.array([stats.spearmanr(sc.loc[p, samp], y)[0] for p in sc.index])
                st.update(score=sname, metric=mname, drug=d, n_pathways=len(allr),
                          share_pathways_rho_le=float(np.mean(allr <= st["rho"])),
                          share_pathways_absrho_ge=float(np.mean(np.abs(allr) >= abs(st["rho"]))))
                rows.append(st)
    T = pd.DataFrame(rows)
    # random gene-set null (current scores only; the authors' NES cannot be recomputed for random sets)
    rng = np.random.default_rng(SEED)
    k = len(gs[BH3])
    pool = np.array(sorted(set(e.index)))
    rsets = {f"rand{i}": list(rng.choice(pool, k, replace=False)) for i in range(N_RAND)}
    rsc = ssgsea_fast(e, rsets)
    samp = [s for s in auc.index if s in cur.columns]
    null = {}
    for mname, R in M.items():
        for d in DRUGS:
            y = R.loc[samp, d].values
            null[(mname, d)] = np.array([stats.spearmanr(rsc.loc[s, samp].values, y)[0] for s in rsc.index])
    T["share_random_sets_rho_le"] = [float(np.mean(null[(r.metric, r.drug)] <= r.rho)) if r.score == "current"
                                     else np.nan for r in T.itertuples()]
    T["share_random_sets_absrho_ge"] = [float(np.mean(np.abs(null[(r.metric, r.drug)]) >= abs(r.rho)))
                                        if r.score == "current" else np.nan for r in T.itertuples()]
    T["n_random_sets"] = N_RAND
    T["genes_in_set"] = k
    T = T[["score", "metric", "drug", "n", "rho", "lo", "hi", "boot_lo", "boot_hi", "p", "n_pathways",
           "share_pathways_rho_le", "share_pathways_absrho_ge", "share_random_sets_rho_le",
           "share_random_sets_absrho_ge", "n_random_sets", "genes_in_set"]]
    T.to_csv(TAB / "G3_bh3_vs_drugs_organoids.tsv", sep="\t", index=False)

    # BH3-only vs every drug (current, AUC) -> where do the BH3 mimetics rank?
    alld = []
    for d in auc.columns:
        y = auc.loc[samp, d]
        if y.notna().sum() >= 10:
            r_, p_ = stats.spearmanr(cur.loc[BH3, samp][y.notna()], y.dropna())
            alld.append(dict(drug=d, n=int(y.notna().sum()), rho=r_, p=p_))
    A = pd.DataFrame(alld).sort_values("rho").reset_index(drop=True)
    A["rank_most_negative"] = np.arange(1, len(A) + 1)
    A.to_csv(TAB / "G3_bh3_vs_all_drugs_organoids_AUC.tsv", sep="\t", index=False)

    # single genes (context for G7): BCL2 family expression vs navitoclax
    g_rows = []
    for g in ["BCL2", "BCL2L1", "MCL1", "BCL2L2", "BAX", "BAK1", "BCL2L11", "PMAIP1", "BBC3", "BAD", "BID", "BMF",
              "E2F1", "TP53"]:
        if g in e.index:
            for d in ("NAVITOCLAX", "OBATOCLAX"):
                st = spearman_ci(e.loc[g, samp], auc.loc[samp, d])
                g_rows.append(dict(gene=g, drug=d, **{k_: st[k_] for k_ in ("n", "rho", "lo", "hi", "p")}))
    pd.DataFrame(g_rows).to_csv(TAB / "G3_single_genes_vs_bh3_mimetics_organoids.tsv", sep="\t", index=False)
    pd.DataFrame(version_rows([VDW / "expression.tsv.gz", VDW / "response.tsv",
                               REPO / "python/results/COAD/organoid/reactome_ssgsea_result/gseapy.samples.normalized.es.txt",
                               REPO / "data/msigdb.v6.1.symbols.gmt.txt"])).to_csv(TAB / "G3_data_versions.tsv", sep="\t", index=False)
    rr = pd.read_csv(VDW / "response.tsv", sep="\t")
    rng_rows = []
    for d in DRUGS:
        q = rr[rr.drug == d]
        rng_rows.append(dict(drug=d, n=len(q), max_conc_uM=q.max_conc_uM.iloc[0], auc_min=q.response.min(),
                             auc_median=q.response.median(), auc_max=q.response.max(), auc_sd=q.response.std(),
                             ln_ic50_min=q.ln_ic50.min(), ln_ic50_max=q.ln_ic50.max(),
                             n_ic50_above_max_conc=int((q.ln_ic50 > np.log(q.max_conc_uM)).sum())))
    pd.DataFrame(rng_rows).to_csv(TAB / "G3_response_range_organoids.tsv", sep="\t", index=False)
    figure(T, A, cur, auc, null, samp)
    print(T[(T.metric == "AUC")][["score", "drug", "n", "rho", "lo", "hi", "p", "share_pathways_rho_le",
                                  "share_random_sets_rho_le"]].to_string())
    print(A.head(8).to_string())


def figure(T, A, cur, auc, null, samp):
    plt = style()
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.2))
    cols = {"NAVITOCLAX": OI["blue"], "OBATOCLAX": OI["sky"], "FLUOROURACIL": OI["orange"], "OXALIPLATIN": OI["red"]}
    x = cur.loc[BH3, samp]
    for ax, d, L in zip(axs[0], ["NAVITOCLAX", "OBATOCLAX", "FLUOROURACIL"], "abc"):
        y = auc.loc[samp, d]
        ax.scatter(x, y, s=28, color=cols[d], edgecolor="k", lw=0.4)
        b = np.polyfit(x, y, 1)
        xx = np.linspace(x.min(), x.max(), 10)
        ax.plot(xx, np.polyval(b, xx), color="k", lw=0.8, ls="--")
        r = T[(T.score == "current") & (T.metric == "AUC") & (T.drug == d)].iloc[0]
        ax.set_title(f"{LAB[d]}: rho = {r.rho:.2f} (95% CI {r.lo:.2f} to {r.hi:.2f})\np = {fmt_p(r.p)}, n = {r.n}")
        ax.set_xlabel("BH3-only pathway score (ssGSEA)")
        ax.set_ylabel(f"{LAB[d]} AUC (lower = more sensitive)")
        panel(ax, L)
    # d forest
    ax = axs[1, 0]
    sub = T.copy()
    sub["lab"] = sub.drug.map(LAB) + ", " + sub.score + ", " + sub.metric
    sub = sub.iloc[::-1].reset_index(drop=True)
    for i, r in sub.iterrows():
        ax.errorbar(r.rho, i, xerr=[[r.rho - r.lo], [r.hi - r.rho]], fmt="o", color=cols[r.drug], ms=4,
                    mfc="white" if r.score == "authors" else cols[r.drug], capsize=2)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub.lab, fontsize=6.5)
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("Spearman rho with 95% CI")
    ax.set_title("All scores x metrics\n(filled = current ssGSEA, open = authors' NES)")
    panel(ax, "d", -0.55)
    # e null
    ax = axs[1, 1]
    nl = null[("AUC", "NAVITOCLAX")]
    ax.hist(nl, bins=40, color=OI["grey"], alpha=0.8)
    r = T[(T.score == "current") & (T.metric == "AUC") & (T.drug == "NAVITOCLAX")].iloc[0]
    ax.axvline(r.rho, color=OI["blue"], lw=2)
    ax.set_xlabel("Spearman rho with navitoclax AUC")
    ax.set_ylabel("random gene sets")
    ax.set_title(f"{len(nl)} random {int(r.genes_in_set)}-gene sets\nshare with rho <= BH3-only: "
                 f"{r.share_random_sets_rho_le:.3f}; all pathways: {r.share_pathways_rho_le:.3f}")
    panel(ax, "e")
    # f all drugs
    ax = axs[1, 2]
    c = [cols.get(d, OI["grey"]) for d in A.drug]
    ax.bar(A.rank_most_negative, A.rho, color=c, width=0.9)
    for d in ["NAVITOCLAX", "OBATOCLAX", "FLUOROURACIL", "OXALIPLATIN"]:
        rr = A[A.drug == d].iloc[0]
        ax.annotate(f"{LAB[d]} ({int(rr.rank_most_negative)})", (rr.rank_most_negative, rr.rho),
                    xytext=(rr.rank_most_negative + 4, rr.rho + (-0.12 if rr.rho < 0 else 0.12)), fontsize=7,
                    arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.set_xlabel(f"drugs ranked by rho (1 = most negative), {len(A)} drugs")
    ax.set_ylabel("rho, BH3-only score vs drug AUC")
    ax.set_title("BH3-only score vs every drug in the screen")
    panel(ax, "f")
    fig.tight_layout()
    rn = T[(T.score == "current") & (T.metric == "AUC")].set_index("drug")
    savefig(fig, "G3_bh3_vs_bh3_mimetics_organoids", f"""
# BH3-only pathway score and drug response in the 19 colorectal organoids (van de Wetering 2015)

**What is shown.** Each dot is one organoid. The x-axis is the score of the Reactome pathway "Activation of BH3-only
proteins" (ssGSEA, recomputed here from the published expression). The y-axis is the drug response as area under
the dose-response curve (AUC): lower AUC means the organoid is more sensitive. In cell lines, a higher BH3-only score
went with more sensitivity to navitoclax (a negative rho).

- **a** navitoclax (blocks BCL-2, BCL-XL and BCL-W): rho = {rn.loc['NAVITOCLAX','rho']:.2f} (95% CI {rn.loc['NAVITOCLAX','lo']:.2f} to {rn.loc['NAVITOCLAX','hi']:.2f}), p = {fmt_p(rn.loc['NAVITOCLAX','p'])}.
- **b** obatoclax (pan-BCL-2-family blocker with other off-target effects): rho = {rn.loc['OBATOCLAX','rho']:.2f}, p = {fmt_p(rn.loc['OBATOCLAX','p'])}.
- **c** 5-FU, for contrast: rho = {rn.loc['FLUOROURACIL','rho']:.2f}, p = {fmt_p(rn.loc['FLUOROURACIL','p'])}.
- **d** All drug x score x metric combinations with 95% confidence intervals (Fisher z). Filled dots use our
  recomputed score, open dots the authors' own organoid scores (NES). "ln IC50" is the log of the half-maximal
  inhibitory concentration (lower = more sensitive).
- **e** Random-pathway check: {int(rn.loc['NAVITOCLAX','n_random_sets'])} random gene sets of the same size ({int(rn.loc['NAVITOCLAX','genes_in_set'])} genes); grey =
  their correlations with navitoclax, blue line = the BH3-only set. The share of random sets with a correlation at
  least as negative is {rn.loc['NAVITOCLAX','share_random_sets_rho_le']:.3f}; among all {int(rn.loc['NAVITOCLAX','n_pathways'])} Reactome pathways it is {rn.loc['NAVITOCLAX','share_pathways_rho_le']:.3f}.
- **f** The BH3-only score correlated with every drug in the screen; bars are sorted from most negative to most
  positive. Numbers in brackets are ranks.

Note on range: navitoclax was tested up to 2 uM only; organoid AUCs lie between 0.94 and 1.00, so almost no
organoid was killed in the tested range (see tables/G3_response_range_organoids.tsv). This leaves little
variation for any marker to explain.

With 19 organoids, a 95% interval spans roughly +/-0.45 around the estimate, so only large correlations can be told
apart from zero.
""")


if __name__ == "__main__":
    main()
