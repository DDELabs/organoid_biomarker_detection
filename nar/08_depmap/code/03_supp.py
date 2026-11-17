"""Supplementary analyses S1-S8 for the cell-line validation of the 5-FU BH3-only biomarker.

python -I nar/08_depmap/code/03_supp.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import stats  # noqa: E402

from common import BH3, CACHE, OI, SEED, TAB, load, panel, savefig, style  # noqa: E402
from data import SCREENS, crispr, gdsc_matrix, mutations, screen  # noqa: E402
from stats import partial_spearman, spearman  # noqa: E402

from obd import scoring  # noqa: E402
from obd.cohorts import to_canonical  # noqa: E402
from obd.reference import canonical_symbol, gene_sets  # noqa: E402

plt = style()
FU = "FLUOROURACIL"
COL = {"GDSC1": OI["blue"], "GDSC2": OI["orange"], "CTRP": OI["green"], "PRISM": OI["purple"]}
EXPR = {"GDSC1": "expr_cmp", "GDSC2": "expr_cmp", "CTRP": "expr_ctrp", "PRISM": "expr_prism"}
BH3_GENES = ["BCL2L11", "BBC3", "PMAIP1", "BID", "BAD", "BMF", "HRK", "BIK"]


def reactome_current():
    out = {}
    for line in open(CACHE / "reactome_gmt" / "ReactomePathways.gmt"):
        p = line.rstrip("\n").split("\t")
        if p[0] in ("Activation of BH3-only proteins",
                    "BH3-only proteins associate with and inactivate anti-apoptotic BCL-2 members"):
            out[p[0]] = sorted({canonical_symbol(g) or g for g in p[2:]})
    return out


def bh3_sets():
    s = {"MSigDB v6.1 Reactome (used by paper)": gene_sets(("REACTOME",))[BH3]}
    for k, v in reactome_current().items():
        s[f"Reactome current: {k[:40]}"] = v
    s["8 core BH3-only genes (no pathway context)"] = [canonical_symbol(g) or g for g in BH3_GENES]
    return s


# ------------------------------------------------------------------ S1 matching flow
def s1():
    rows = []
    m = load("models")
    K = crispr()
    mut = mutations()
    for n in SCREENS:
        S = screen(n)
        screened = S.resp[FU].dropna().index
        withexpr = screened.intersection(S.scores.columns)
        col = [i for i in withexpr if i in S.colon]
        colk = [i for i in col if isinstance(S.sidm.get(i), str) and S.sidm[i] in K.columns]
        colm = [i for i in col if isinstance(S.sidm.get(i), str) and S.sidm[i] in mut.index]
        rows.append({"screen": n, "all_lines_in_screen": len(S.resp), "lines_with_5FU": len(screened),
                     "5FU_and_RNAseq": len(withexpr), "colorectal_5FU_RNAseq": len(col),
                     "colorectal_plus_CRISPR": len(colk), "colorectal_plus_mutation_data": len(colm),
                     "bladder_lines_any_drug": len(S.bladder)})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S1_data_matching_flow.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(9, 4))
    cols = list(t.columns[1:-1])
    for j, n in enumerate(SCREENS):
        v = t[t.screen == n][cols].values[0]
        ax.plot(range(len(cols)), v, "-o", color=COL[n], label=n)
        for i, x in enumerate(v):
            ax.text(i, x * 1.08, str(int(x)), ha="center", fontsize=7, color=COL[n])
    ax.set_yscale("log")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=20, ha="right")
    ax.set_ylabel("Number of cell lines (log scale)")
    ax.legend(frameon=False)
    panel(ax, "a")
    savefig(fig, "S1_data_matching_flow", """
**Supplementary Figure S1. How many cell lines survive each matching step.** From all lines in each screen, to
lines with a 5-FU measurement, to lines that also have RNA-seq (Cell Model Passports for GDSC; CCLE 21Q4 for CTRPv2
and PRISM), then colorectal lines only, and finally colorectal lines that also have CRISPR (Project Score) or
mutation data (Cell Model Passports). GDSC lines are matched by Sanger model ID; CTRPv2/PRISM lines by DepMap ID and
then to Sanger model IDs through the DepMap annotation. Numbers: `tables/S1_data_matching_flow.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S2 GDSC1 vs GDSC2
def s2():
    g1, g2 = gdsc_matrix("GDSC1")[FU], gdsc_matrix("GDSC2")[FU]
    S = screen("GDSC2")
    c = screen("CTRP")
    ctrp = c.resp[FU].dropna()
    ctrp.index = [c.sidm.get(i) for i in ctrp.index]
    ctrp = ctrp[[isinstance(i, str) for i in ctrp.index]]
    ctrp = ctrp.groupby(level=0).median()
    rows = []
    for a, b, na, nb in ((g1, g2, "GDSC1", "GDSC2"), (g2, ctrp, "GDSC2", "CTRP")):
        idx = a.dropna().index.intersection(b.dropna().index)
        r, p, n = spearman(a[idx], b[idx])
        col = [i for i in idx if i in S.colon or S.lineage.get(i) == "Colorectal Carcinoma"]
        rc, pc, nc = spearman(a[col], b[col])
        rows.append({"pair": f"{na} vs {nb} (5-FU)", "n_all": n, "rho_all": r, "p_all": p, "n_colorectal": nc,
                     "rho_colorectal": rc, "p_colorectal": pc})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S2_gdsc1_vs_gdsc2_agreement.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, (a, b, na, nb), L in zip(axes, ((g1, g2, "GDSC1 ln IC50", "GDSC2 ln IC50"), (g2, ctrp, "GDSC2 ln IC50", "CTRPv2 AUC")), "ab"):
        idx = a.dropna().index.intersection(b.dropna().index)
        col = [i for i in idx if S.lineage.get(i) == "Colorectal Carcinoma"]
        ax.scatter(a[idx], b[idx], s=4, color=OI["grey"], alpha=0.5, label="all lines")
        ax.scatter(a[col], b[col], s=14, color=OI["red"], edgecolor="black", lw=0.3, label="colorectal")
        r = t.iloc[0 if L == "a" else 1]
        ax.set_title(f"5-FU: {na.split()[0]} vs {nb.split()[0]}\nall rho={r.rho_all:.2f} (n={r.n_all}); "
                     f"CRC rho={r.rho_colorectal:.2f} (n={r.n_colorectal})", fontsize=8.5)
        ax.set_xlabel(na)
        ax.set_ylabel(nb)
        panel(ax, L)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    savefig(fig, "S2_gdsc1_vs_gdsc2", """
**Supplementary Figure S2. Do the screens agree on 5-FU?** (a) 5-FU ln IC50 of the same cell lines in GDSC1 and
GDSC2 (different assay generations at the Sanger Institute). (b) GDSC2 ln IC50 against CTRPv2 AUC (Broad).
Grey = all shared lines, red = colorectal lines. Spearman rho is given for both. If two screens disagree on which
lines are resistant, no expression marker can agree with both. Numbers: `tables/S2_gdsc1_vs_gdsc2_agreement.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S3 colon vs pan-cancer
def s3():
    rows = []
    for n in SCREENS:
        S = screen(n)
        a = S.common(FU)
        lin = S.lineage.loc[a]
        r, p, k = partial_spearman(S.scores.loc[BH3, a], S.resp.loc[a, FU], lin)
        rows.append({"screen": n, "group": "pan-cancer (lineage-adjusted)", "n": k, "rho": r, "p": p})
        r, p, k = spearman(S.scores.loc[BH3, a], S.resp.loc[a, FU])
        rows.append({"screen": n, "group": "pan-cancer (unadjusted)", "n": k, "rho": r, "p": p})
        for L, g in lin.groupby(lin):
            if len(g) >= 15:
                r, p, k = spearman(S.scores.loc[BH3, g.index], S.resp.loc[g.index, FU])
                rows.append({"screen": n, "group": str(L), "n": k, "rho": r, "p": p})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S3_colon_vs_pancancer.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 4, figsize=(13, 6), sharex=True)
    for ax, n, L in zip(axes, SCREENS, "abcd"):
        g = t[t.screen == n].reset_index(drop=True)
        colon = g.group.isin(["Colorectal Carcinoma", "colorectal"])
        cols = [OI["red"] if c else (OI["black"] if "pan-cancer" in grp else OI["grey"]) for c, grp in zip(colon, g.group)]
        ax.scatter(g.rho, range(len(g)), c=cols, s=18)
        for i, r in g.iterrows():
            se = 1 / np.sqrt(max(r.n - 3, 1))
            ax.plot([np.tanh(np.arctanh(r.rho) - 1.96 * se), np.tanh(np.arctanh(r.rho) + 1.96 * se)], [i, i], color=cols[i], lw=0.8)
        ax.set_yticks(range(len(g)))
        ax.set_yticklabels([f"{x} (n={k})" for x, k in zip(g.group, g.n)], fontsize=6)
        ax.axvline(0, color="black", lw=0.5)
        ax.invert_yaxis()
        ax.set_title(n)
        ax.set_xlabel("rho, BH3-only score vs 5-FU")
        panel(ax, L)
    fig.tight_layout()
    savefig(fig, "S3_colon_vs_pancancer", """
**Supplementary Figure S3. Colorectal lines compared with other cancer types.** Spearman rho between the BH3-only
activation score and 5-FU response in each cancer type with at least 15 lines (grey), colorectal lines (red), and
all lines together with and without adjustment for cancer type (black). Lines show 95% confidence intervals.
Negative = BH3-high lines more sensitive (the organoid direction). Numbers: `tables/S3_colon_vs_pancancer.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S4 permutation / random controls
def s4():
    rng = np.random.default_rng(SEED)
    rows, dists = [], {}
    gsall = gene_sets(("REACTOME",))
    for n in SCREENS:
        S = screen(n)
        c = S.common(FU, S.colon)
        y = S.resp.loc[c, FU].values
        x = S.scores.loc[BH3, c].values
        r = stats.spearmanr(x, y)[0]
        perm = np.array([stats.spearmanr(x, rng.permutation(y))[0] for _ in range(10000)])
        allp = np.array([stats.spearmanr(S.scores.loc[p, c].values, y)[0] for p in S.scores.index])
        E = load(EXPR[n])[c]
        genes = [g for g in gsall[BH3] if g in E.index]
        pool = np.array(E.index)
        pct = (E.rank(axis=0) / (len(E) + 1)).values
        idx = {g: i for i, g in enumerate(E.index)}
        rs = []
        for _ in range(1000):
            pick = rng.choice(len(pool), len(genes), replace=False)
            rs.append(stats.spearmanr(pct[pick].mean(axis=0), y)[0])
        rs = np.array(rs)
        dists[n] = (r, perm, allp, rs)
        rows.append({"screen": n, "n": len(c), "rho_BH3": r,
                     "label_permutation_p_two_sided": float((1 + (np.abs(perm) >= abs(r)).sum()) / 10001),
                     "label_permutation_p_one_sided_negative": float((1 + (perm <= r).sum()) / 10001),
                     "reactome_pathways": len(allp), "BH3_percentile_among_all_pathways": float((allp <= r).mean()),
                     "random_gene_sets_size": len(genes),
                     "random_gene_set_p_one_sided_negative": float((1 + (rs <= r).sum()) / 1001),
                     "random_gene_set_q025": float(np.quantile(rs, 0.025)), "random_gene_set_q975": float(np.quantile(rs, 0.975))})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S4_permutation_random_controls.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(3, 4, figsize=(12, 7.5))
    for j, n in enumerate(SCREENS):
        r, perm, allp, rs = dists[n]
        for i, (d, lab) in enumerate(((perm, "shuffled 5-FU labels (10,000)"), (allp, "all Reactome pathways"),
                                     (rs, f"random gene sets of size {t.random_gene_sets_size[j]}"))):
            ax = axes[i, j]
            ax.hist(d, bins=40, color=OI["grey"])
            ax.axvline(r, color=OI["red"], lw=1.5)
            ax.set_title(f"{n}: {lab}", fontsize=8)
            ax.set_xlabel("rho with 5-FU (colorectal)")
            if j == 0:
                panel(ax, "abc"[i])
    fig.tight_layout()
    savefig(fig, "S4_permutation_and_random_controls", """
**Supplementary Figure S4. Is the BH3-only result different from chance?** Red line = observed Spearman rho of the
BH3-only score with 5-FU response in colorectal lines. Grey histograms = what chance gives: (a) the same score with
the 5-FU values shuffled 10,000 times; (b) every Reactome pathway scored the same way; (c) 1,000 random gene sets of
the same size as the BH3-only set. Numbers (p-values and percentiles): `tables/S4_permutation_random_controls.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S5 mutation adjustment
def s5():
    mut = mutations()
    rows = []
    for n in SCREENS:
        S = screen(n)
        c = [i for i in S.common(FU, S.colon) if isinstance(S.sidm.get(i), str) and S.sidm[i] in mut.index]
        d = pd.DataFrame({"y": S.resp.loc[c, FU].values, "BH3": S.scores.loc[BH3, c].values}, index=c)
        d = d.join(mut.loc[[S.sidm[i] for i in c], ["TP53", "MSI", "KRAS", "BRAF"]].set_axis(c))
        d = d.dropna()
        d["y"] = (d["y"] - d["y"].mean()) / d["y"].std()
        d["BH3"] = (d["BH3"] - d["BH3"].mean()) / d["BH3"].std()
        for model, cov in (("unadjusted", ["BH3"]), ("adjusted TP53+MSI+KRAS+BRAF", ["BH3", "TP53", "MSI", "KRAS", "BRAF"])):
            cv = [x for x in cov if d[x].nunique() > 1]
            f = sm.OLS(d["y"], sm.add_constant(d[cv])).fit()
            ci = f.conf_int().loc["BH3"]
            rows.append({"screen": n, "model": model, "n": len(d), "beta_BH3_per_SD": f.params["BH3"], "ci_low": ci[0],
                         "ci_high": ci[1], "p": f.pvalues["BH3"],
                         **{f"beta_{k}": f.params.get(k, np.nan) for k in ("TP53", "MSI", "KRAS", "BRAF")},
                         **{f"p_{k}": f.pvalues.get(k, np.nan) for k in ("TP53", "MSI", "KRAS", "BRAF")},
                         **{f"n_mut_{k}": int(d[k].sum()) for k in ("TP53", "MSI", "KRAS", "BRAF")}})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S5_mutation_adjustment.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for i, r in t.iterrows():
        ax.errorbar(r.beta_BH3_per_SD, i, xerr=[[r.beta_BH3_per_SD - r.ci_low], [r.ci_high - r.beta_BH3_per_SD]],
                    fmt="o" if r.model == "unadjusted" else "s", color=COL[r.screen])
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{r.screen}, {r.model} (n={r.n})" for r in t.itertuples()], fontsize=7)
    ax.axvline(0, color="black", lw=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Change in 5-FU response (SD) per 1 SD of BH3-only score\n(negative = BH3-high lines more sensitive)")
    panel(ax, "a")
    savefig(fig, "S5_mutation_adjustment", """
**Supplementary Figure S5. Adjusting for colorectal driver mutations.** Linear regression in colorectal lines of
standardised 5-FU response on the standardised BH3-only score, alone (circles) and together with TP53, KRAS and BRAF
mutation status and microsatellite instability (MSI) (squares). Bars are 95% confidence intervals. TP53 matters
because p53 switches on the BH3-only genes PUMA (BBC3) and NOXA (PMAIP1) after DNA damage. Mutation and MSI calls:
Cell Model Passports. Numbers including the mutation coefficients: `tables/S5_mutation_adjustment.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S6 single genes
def s6():
    rows = []
    members = sorted(set(gene_sets(("REACTOME",))[BH3]) | set(BH3_GENES))
    for n in SCREENS:
        S = screen(n)
        E = load(EXPR[n])
        c = [i for i in S.common(FU, S.colon) if i in E.columns]
        for g in members:
            if g not in E.index:
                continue
            r, p, k = spearman(E.loc[g, c], S.resp.loc[c, FU])
            rows.append({"screen": n, "gene": g, "core_BH3_only": g in BH3_GENES, "n": k, "rho": r, "p": p})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S6_single_genes.tsv", sep="\t", index=False)
    m = t.pivot(index="gene", columns="screen", values="rho")[list(SCREENS)]
    m = m.loc[[g for g in BH3_GENES if g in m.index] + [g for g in m.index if g not in BH3_GENES]]
    fig, ax = plt.subplots(figsize=(5, 7))
    im = ax.imshow(m.values, cmap="PuOr_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels([f"{g}{' *' if g in BH3_GENES else ''}" for g in m.index], fontsize=7)
    ax.set_xticks(range(4))
    ax.set_xticklabels(m.columns)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if np.isfinite(m.values[i, j]):
                ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="Spearman rho, expression vs 5-FU response")
    panel(ax, "a")
    savefig(fig, "S6_single_bh3_genes", """
**Supplementary Figure S6. Each gene on its own.** Spearman rho between the expression (log2 TPM+1) of each gene and
5-FU response in colorectal lines. Starred genes are the 8 core BH3-only genes (BIM = BCL2L11, PUMA = BBC3,
NOXA = PMAIP1, BID, BAD, BMF, HRK, BIK); the others are members of the Reactome "Activation of BH3-only proteins"
pathway (kinases, 14-3-3 proteins, TP53 family, E2F1). Negative = higher expression with more sensitivity.
Numbers with n and p: `tables/S6_single_genes.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S7 isogenic pairs
def platform_map(gpl, col):
    f = CACHE / f"gpl_{gpl}" / f"{gpl}.txt"
    lines = open(f, errors="replace").read().split("\n")
    s = lines.index("!platform_table_begin") + 1
    e = lines.index("!platform_table_end")
    import io
    t = pd.read_csv(io.StringIO("\n".join(lines[s:e])), sep="\t", low_memory=False, dtype=str)
    m = t.set_index("ID")[col].dropna().astype(str).str.split(" /// ").str[0]
    return m


def isogenic_sets():
    from obd.geo_matrix import read_expr, read_meta
    out = {}
    # GSE81005: HCT-8 parental vs 5-FU resistant, 3 replicates at 0, 24, 48 h of 5-FU
    p = CACHE / "geo_GSE81005" / "GSE81005_series_matrix.txt.gz"
    x, meta = read_expr(p), read_meta(p)
    sym = platform_map("GPL15207", "Gene Symbol")
    x = x.groupby(sym.reindex(x.index).values).median()
    for tp in ("0h", "24h", "48h"):
        cols = [g for g in meta.index if f"-{tp}-" in meta.loc[g, "title"]]
        grp = ["resistant" if "5Fu" in meta.loc[g, "title"] else "parental" for g in cols]
        out[f"GSE81005 HCT-8 ({tp} 5-FU)"] = (x[cols], pd.Series(grp, index=cols))
    # GSE196900: HCT116 and SW480 parental vs 5-FU resistant (RNA-seq counts)
    c = pd.read_csv(CACHE / "geo_GSE196900" / "GSE196900_Raw_gene_counts_matrix.txt.gz", sep="\t")
    c = c.set_index("Gene Symbol").iloc[:, 5:].apply(pd.to_numeric, errors="coerce").fillna(0)
    c = c.groupby(level=0).sum()
    cpm = np.log2(c / c.sum() * 1e6 + 1)
    for line in ("HCT116", "SW480"):
        cols = [k for k in cpm.columns if k.startswith(line)]
        out[f"GSE196900 {line}"] = (cpm[cols], pd.Series(["resistant" if "R" in k[len(line):] else "parental" for k in cols], index=cols))
    # GSE56322: one parental vs one 5-FU-resistant colon cancer line (Agilent)
    p = CACHE / "geo_GSE56322" / "GSE56322_series_matrix.txt.gz"
    x, meta = read_expr(p), read_meta(p)
    sym = platform_map("GPL14550", "GENE_SYMBOL")
    x = x.groupby(sym.reindex(x.index).values).median()
    out["GSE56322 colon line (1 vs 1)"] = (x, pd.Series(["resistant" if "resistant" in t else "parental" for t in meta["title"]], index=meta.index))
    return out


def s7():
    sets = bh3_sets()
    rows = []
    for name, (x, grp) in isogenic_sets().items():
        x = to_canonical(x.apply(pd.to_numeric, errors="coerce").dropna(how="all"))
        sc = scoring.rank_score(x, {k: [g for g in v if g in x.index] for k, v in sets.items()}, min_size=5)
        for setname in sc.index:
            v = sc.loc[setname]
            a, b = v[grp[grp == "resistant"].index], v[grp[grp == "parental"].index]
            p = stats.ttest_ind(a, b, equal_var=False).pvalue if min(len(a), len(b)) >= 2 else np.nan
            sd = v.std() if len(v) > 2 else np.nan
            rows.append({"dataset": name, "gene_set": setname, "n_resistant": len(a), "n_parental": len(b),
                         "mean_resistant": a.mean(), "mean_parental": b.mean(), "difference_res_minus_par": a.mean() - b.mean(),
                         "welch_p": p, "expected_sign": "negative (resistant has lower BH3-only activation)"})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S7_isogenic_pairs.tsv", sep="\t", index=False)
    g = t[t.gene_set.str.startswith("MSigDB")].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.barh(range(len(g)), g.difference_res_minus_par, color=[OI["blue"] if d < 0 else OI["red"] for d in g.difference_res_minus_par])
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels([f"{d} (p={p:.2g})" if np.isfinite(p) else f"{d} (no replicates)" for d, p in zip(g.dataset, g.welch_p)], fontsize=7)
    ax.axvline(0, color="black", lw=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("BH3-only score: resistant minus parental\n(blue = lower in resistant, the organoid direction)")
    panel(ax, "a")
    savefig(fig, "S7_isogenic_pairs", """
**Supplementary Figure S7. Matched 5-FU-resistant and parental colorectal lines (GEO).** For each pair we compare
the BH3-only activation score (same rank score, MSigDB v6.1 gene set) between the resistant derivative and its
parent. The organoid model expects resistant cells to have LOWER BH3-only activation (blue bars). GSE81005: HCT-8
and HCT-8/5-FU, 3 replicates each, measured before (0 h) and after 24 h and 48 h of 5-FU (Affymetrix PrimeView).
GSE196900: HCT116 and SW480 with their 5-FU-resistant sublines, 2 replicates each (RNA-seq counts, log2 CPM).
GSE56322: one parental and one resistant line, no replicates (Agilent). p = Welch t-test. Other gene-set versions:
`tables/S7_isogenic_pairs.tsv`.
""")
    plt.close(fig)
    return t


# ------------------------------------------------------------------ S8 metric + gene-set version
def s8():
    rows = []
    sets = bh3_sets()
    for n in SCREENS:
        S = screen(n)
        E = load(EXPR[n])
        c0 = S.common(FU, S.colon)
        a0 = S.common(FU)
        sc = scoring.rank_score(E[list(dict.fromkeys(a0))], {k: [g for g in v if g in E.index] for k, v in sets.items()}, min_size=5)
        metrics = [("ln_ic50", None), ("auc", "auc")] if n.startswith("GDSC") else [("AUC", None)]
        for mname, alt in metrics:
            resp = gdsc_matrix(n, "auc")[FU] if alt else S.resp[FU]
            for setname in sc.index:
                c = [i for i in c0 if i in resp.dropna().index]
                a = [i for i in a0 if i in resp.dropna().index]
                r, p, k = spearman(sc.loc[setname, c], resp[c])
                r2, p2, k2 = partial_spearman(sc.loc[setname, a], resp[a], S.lineage.loc[a])
                rows.append({"screen": n, "metric": mname, "gene_set": setname, "colon_n": k, "colon_rho": r, "colon_p": p,
                             "pancancer_n": k2, "pancancer_partial_rho": r2, "pancancer_p": p2})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "S8_metric_and_geneset_robustness.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    lab = t.screen + " " + t.metric
    ys = {k: i for i, k in enumerate(dict.fromkeys(lab))}
    mk = dict(zip(sets, "osD^"))
    for ax, col, L in zip(axes, ("colon_rho", "pancancer_partial_rho"), "ab"):
        for s_ in sets:
            g = t[t.gene_set == s_]
            ax.scatter(g[col], [ys[x] for x in (g.screen + " " + g.metric)], marker=mk[s_], s=26, label=s_,
                       edgecolor="black", lw=0.3)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(list(ys.values()))
        ax.set_yticklabels(list(ys.keys()))
        ax.set_xlabel("rho with 5-FU response" + (" (colorectal)" if L == "a" else " (pan-cancer, lineage-adjusted)"))
        panel(ax, L)
    axes[0].invert_yaxis()
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0, -0.45))
    fig.tight_layout()
    savefig(fig, "S8_metric_and_geneset_robustness", """
**Supplementary Figure S8. Does the answer depend on the response measure or the gene-set version?** Spearman rho
between BH3-only scores and 5-FU response for: GDSC ln IC50 and GDSC AUC (CTRPv2/PRISM report AUC only); and four
versions of the gene set: the MSigDB v6.1 Reactome set used by the paper (16 genes), the current Reactome
"Activation of BH3-only proteins" (R-HSA-114452, 30 genes), the current Reactome "BH3-only proteins associate with
and inactivate anti-apoptotic BCL-2 members" (R-HSA-111453), and the 8 core BH3-only genes alone.
(a) colorectal lines; (b) all lines, adjusted for cancer type. Numbers: `tables/S8_metric_and_geneset_robustness.tsv`.
""")
    plt.close(fig)
    return t


if __name__ == "__main__":
    for f in (s1, s2, s3, s4, s5, s6, s7, s8):
        t = f()
        print(f.__name__)
        print(t.round(3).to_string()[:6000])
