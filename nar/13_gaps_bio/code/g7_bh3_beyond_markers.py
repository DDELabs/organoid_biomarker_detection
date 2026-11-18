"""G7: does the BH3-only pathway score add to simple markers of BH3-mimetic sensitivity in cell lines?

Screens (nar/08_depmap derived data, reused): GDSC1, GDSC2 (ln IC50), CTRPv2 (AUC); lower = more sensitive.
Drugs: navitoclax and ABT-737 (BCL-2/BCL-XL/BCL-W), plus venetoclax (BCL-2 only), WEHI-539 (BCL-XL only),
AZD5991 (MCL1 only) for contrast.
Covariates: cancer type (lineage dummies), rank-normalised expression of BCL2, BCL2L1 (BCL-XL), MCL1, and a
proliferation score (HALLMARK_E2F_TARGETS, rank score as the Reactome scores).
All variables are rank-normalised (inverse normal) before OLS, so coefficients are Spearman-like.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from gcommon import BH3, FIG, OI, REPO, SEED, TAB, fmt_p, panel, savefig, style

sys.path.insert(0, str(REPO / "nar/08_depmap/code"))
from data import screen  # noqa: E402
from common import lineage_dummies, rank_z, load  # noqa: E402  (08_depmap common)

sys.path.insert(0, str(REPO / "robust"))
from obd import scoring  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

DRUGS = {"GDSC1": ["NAVITOCLAX", "VENETOCLAX", "WEHI-539", "OBATOCLAX"],
         "GDSC2": ["NAVITOCLAX", "ABT-737", "VENETOCLAX", "WEHI-539", "AZD5991", "OBATOCLAX"],
         "CTRP": ["NAVITOCLAX", "ABT-737", "VENETOCLAX", "OBATOCLAX"]}
GENES = ["BCL2", "BCL2L1", "MCL1"]
CORE = ["BAD", "BBC3", "BCL2L11", "BID", "BMF", "PMAIP1"]      # the BH3-only genes themselves
N_RAND = 1000


def design(S, expr, extra_scores, drug):
    s = S.common(drug)
    s = [x for x in s if x in expr.columns]
    df = pd.DataFrame({"y": S.resp.loc[s, drug], "BH3": S.scores.loc[BH3, s], "lin": S.lineage.reindex(s)})
    for g in GENES:
        df[g] = expr.loc[g, s].values
    for k, v in extra_scores.items():
        df[k] = v.reindex(s).values
    df = df.dropna()
    for c in df.columns:
        if c != "lin":
            df[c] = rank_z(df[c].values)
    return df


def ols(df, xs):
    L = lineage_dummies(df["lin"], 5)
    X = pd.concat([df[xs], L], axis=1) if xs else L
    X = sm.add_constant(X.astype(float))
    return sm.OLS(df["y"].values, X).fit()


def partial_r(df, x, covs):
    L = lineage_dummies(df["lin"], 5)
    C = pd.concat([df[covs], L], axis=1) if covs else L
    C = sm.add_constant(C.astype(float))
    rx = sm.OLS(df[x].values, C).fit().resid
    ry = sm.OLS(df["y"].values, C).fit().resid
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = len(df) - C.shape[1] - 1
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    z, se = np.arctanh(r), 1 / np.sqrt(max(dof - 1, 1))
    return r, float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se)), float(2 * stats.t.sf(abs(t), dof))


def main():
    gs = gene_sets(("HALLMARK", "REACTOME"))
    rows, rows_b, nulls = [], [], {}
    rng = np.random.default_rng(SEED)
    for sn in DRUGS:
        S = screen(sn)
        expr = load("expr_cmp" if sn.startswith("GDSC") else "expr_ctrp")
        ex = expr[[c for c in S.scores.columns if c in expr.columns]]
        extra = scoring.rank_score(ex, {"E2F": gs["HALLMARK_E2F_TARGETS"], "BH3core": CORE,
                                        "BH3_noE2F": [g for g in gs[BH3] if g not in ("E2F1", "TFDP1")]}, min_size=3)
        extra = {k: extra.loc[k] for k in extra.index}
        k = len([g for g in gs[BH3] if g in ex.index])
        pool = np.array(ex.index)
        rs = scoring.rank_score(ex, {f"r{i}": list(rng.choice(pool, k, replace=False)) for i in range(N_RAND)}, min_size=3)
        for d in DRUGS[sn]:
            df = design(S, ex, extra, d)
            n = len(df)
            for label, covs in (("lineage", []), ("+ BCL2L1", ["BCL2L1"]), ("+ BCL2, BCL2L1, MCL1", GENES),
                                ("+ proliferation (E2F)", ["E2F"]), ("+ all three genes + E2F", GENES + ["E2F"])):
                r, lo, hi, p = partial_r(df, "BH3", covs)
                rows.append(dict(screen=sn, drug=d, n=n, adjusted_for=label, var="BH3-only pathway", r=r, lo=lo, hi=hi, p=p))
            for v in GENES + ["E2F", "BH3core", "BH3_noE2F"]:
                r, lo, hi, p = partial_r(df, v, [])
                rows.append(dict(screen=sn, drug=d, n=n, adjusted_for="lineage", var=v, r=r, lo=lo, hi=hi, p=p))
            m0, m1 = ols(df, GENES + ["E2F"]), ols(df, GENES + ["E2F", "BH3"])
            f, pf, _ = m1.compare_f_test(m0)
            b = m1.params["BH3"]
            ci = m1.conf_int().loc["BH3"]
            rows_b.append(dict(screen=sn, drug=d, n=n, R2_markers=m0.rsquared, R2_markers_plus_BH3=m1.rsquared,
                               dR2=m1.rsquared - m0.rsquared, F=f, p_F=pf, beta_BH3=b, beta_lo=ci[0], beta_hi=ci[1],
                               **{f"beta_{g}": m1.params[g] for g in GENES + ["E2F"]},
                               **{f"p_{g}": m1.pvalues[g] for g in GENES + ["E2F"]}))
            # null: all 672 pathways and random gene sets, lineage-adjusted partial Spearman
            if d in ("NAVITOCLAX", "ABT-737"):
                s = df.index
                L = sm.add_constant(lineage_dummies(df["lin"], 5).astype(float))
                H = L.values @ np.linalg.pinv(L.values)
                ry = df["y"].values - H @ df["y"].values

                def pr(vals):
                    x = rank_z(vals)
                    rx = x - H @ x
                    return float(np.corrcoef(rx, ry)[0, 1])
                allp = np.array([pr(S.scores.loc[p, s].values) for p in S.scores.index])
                rnd = np.array([pr(rs.loc[q, s].values) for q in rs.index])
                obs = pr(S.scores.loc[BH3, s].values)
                nulls[(sn, d)] = (obs, allp, rnd)
                rows_b[-1].update(obs_partial=obs, share_pathways_le=float(np.mean(allp <= obs)),
                                  rank_pathways=int(np.sum(allp < obs) + 1), n_pathways=len(allp),
                                  share_random_le=float(np.mean(rnd <= obs)), n_random=len(rnd), set_size=k)
            print(sn, d, n, flush=True)
    A, B = pd.DataFrame(rows), pd.DataFrame(rows_b)
    A.to_csv(TAB / "G7_partial_correlations.tsv", sep="\t", index=False)
    B.to_csv(TAB / "G7_regression_added_value.tsv", sep="\t", index=False)
    figure(A, B, nulls)


def g5_cell_lines():
    """G5a: BH3-only score vs oxaliplatin (and 5-FU for contrast) in colorectal lines and pan-cancer (lineage-adjusted)."""
    sys.path.insert(0, str(REPO / "nar/08_depmap/code"))
    from stats import partial_spearman, spearman
    rows = []
    for sn in ("GDSC1", "GDSC2", "CTRP", "PRISM"):
        S = screen(sn)
        for d in ("OXALIPLATIN", "FLUOROURACIL"):
            if d not in S.resp:
                rows.append(dict(screen=sn, drug=d, scope="colorectal", n=0))
                continue
            s = S.common(d)
            crc = [x for x in s if x in set(S.colon)]
            r, p, n = spearman(S.scores.loc[BH3, crc], S.resp.loc[crc, d])
            z, se = np.arctanh(r), 1.06 / np.sqrt(max(n - 3, 1))
            rows.append(dict(screen=sn, drug=d, scope="colorectal", n=n, rho=r, lo=np.tanh(z - 1.96 * se),
                             hi=np.tanh(z + 1.96 * se), p=p, metric=S.metric))
            r, p, n = partial_spearman(S.scores.loc[BH3, s], S.resp.loc[s, d], S.lineage.reindex(s))
            z, se = np.arctanh(r), 1 / np.sqrt(max(n - 3, 1))
            rows.append(dict(screen=sn, drug=d, scope="pan-cancer, lineage-adjusted", n=n, rho=r, lo=np.tanh(z - 1.96 * se),
                             hi=np.tanh(z + 1.96 * se), p=p, metric=S.metric))
    pd.DataFrame(rows).to_csv(TAB / "G5_cell_lines_bh3_vs_oxaliplatin.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())


def figure(A, B, nulls):
    plt = style()
    fig = plt.figure(figsize=(12.5, 8.5))
    gsp = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gsp[0, :2])
    sub = A[(A["var"] == "BH3-only pathway") & A.drug.isin(["NAVITOCLAX", "ABT-737"])].reset_index(drop=True)
    labs = ["lineage", "+ BCL2L1", "+ BCL2, BCL2L1, MCL1", "+ proliferation (E2F)", "+ all three genes + E2F"]
    cols = [OI["grey"], OI["sky"], OI["blue"], OI["orange"], OI["red"]]
    groups = list(dict.fromkeys(zip(sub.screen, sub.drug)))
    for gi, (sn, d) in enumerate(groups):
        for li, (lab, c) in enumerate(zip(labs, cols)):
            r = sub[(sub.screen == sn) & (sub.drug == d) & (sub.adjusted_for == lab)].iloc[0]
            x = gi + (li - 2) * 0.15
            ax.errorbar(x, r.r, yerr=[[r.r - r.lo], [r.hi - r.r]], fmt="o", color=c, ms=4, capsize=2,
                        label=lab if gi == 0 else None)
    ax.axhline(0, color="grey", lw=0.7)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{s}\n{d.lower()}\n(n = {int(sub[(sub.screen == s) & (sub.drug == d)].n.iloc[0])})"
                        for s, d in groups])
    ax.set_ylabel("partial Spearman r, BH3-only score vs response\n(negative = higher score, more sensitive)")
    ax.set_title("BH3-only score after adjusting for simple markers (95% CI)")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="lower left")
    panel(ax, "a", -0.07)
    ax = fig.add_subplot(gsp[0, 2])
    vs = ["BH3-only pathway", "BCL2", "BCL2L1", "MCL1", "E2F", "BH3core", "BH3_noE2F"]
    vl = ["BH3-only pathway", "BCL2", "BCL2L1 (BCL-XL)", "MCL1", "E2F proliferation", "6 core BH3-only genes",
          "BH3-only set without E2F1/TFDP1"]
    q = A[(A.drug == "NAVITOCLAX") & (A.adjusted_for == "lineage")]
    for si, (sn, c) in enumerate(zip(["GDSC1", "GDSC2", "CTRP"], [OI["blue"], OI["green"], OI["purple"]])):
        for vi, v in enumerate(vs):
            r = q[(q.screen == sn) & (q["var"] == v)].iloc[0]
            ax.errorbar(r.r, vi + (si - 1) * 0.22, xerr=[[r.r - r.lo], [r.hi - r.r]], fmt="o", color=c, ms=3.5,
                        capsize=1.5, label=sn if vi == 0 else None)
    ax.set_yticks(range(len(vs)))
    ax.set_yticklabels(vl, fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("partial r with navitoclax (lineage-adjusted)")
    ax.set_title("Single markers vs navitoclax")
    ax.legend(frameon=False, fontsize=7)
    panel(ax, "b", -0.6)
    for j, key in enumerate([("GDSC2", "NAVITOCLAX"), ("CTRP", "NAVITOCLAX"), ("GDSC2", "ABT-737")]):
        ax = fig.add_subplot(gsp[1, j])
        obs, allp, rnd = nulls[key]
        ax.hist(rnd, bins=40, color=OI["grey"], alpha=0.7, density=True, label="random gene sets")
        ax.hist(allp, bins=40, color=OI["sky"], alpha=0.5, density=True, label="all Reactome pathways")
        ax.axvline(obs, color=OI["red"], lw=2, label="BH3-only")
        b = B[(B.screen == key[0]) & (B.drug == key[1])].iloc[0]
        ax.set_title(f"{key[0]} {key[1].lower()}: rank {int(b.rank_pathways)}/{int(b.n_pathways)} pathways;\n"
                     f"random-set share <= obs: {b.share_random_le:.3f}; added R2 = {b.dR2:.3f} (p = {fmt_p(b.p_F)})",
                     fontsize=8)
        ax.set_xlabel("lineage-adjusted partial r")
        if j == 0:
            ax.legend(frameon=False, fontsize=7)
        panel(ax, "cde"[j])
    fig.tight_layout()
    bn = B[B.drug.isin(["NAVITOCLAX", "ABT-737"])]
    lines = "\n".join(f"- {r.screen} {r.drug.lower()} (n = {r.n}): added R2 = {r.dR2:.3f}, F-test p = {fmt_p(r.p_F)}; "
                      f"BH3-only beta {r.beta_BH3:.3f} ({r.beta_lo:.3f} to {r.beta_hi:.3f}); rank {int(r.rank_pathways)} of "
                      f"{int(r.n_pathways)} pathways; random-set share {r.share_random_le:.3f}" for r in bn.itertuples())
    savefig(fig, "G7_bh3_beyond_simple_markers", f"""
# Does the BH3-only pathway tell us more than BCL2-family genes and proliferation?

Cell-line screens: GDSC1, GDSC2 (response = ln IC50) and CTRPv2 (response = AUC); lower = more sensitive. Navitoclax
and ABT-737 block BCL-2, BCL-XL and BCL-W. All values were rank-normalised; "partial r" is the correlation that is
left after removing the listed covariates (always including cancer type).

- **a** BH3-only score vs response, adjusted for cancer type only (grey), and then also for BCL2L1 (the gene for
  BCL-XL), for the three main anti-apoptotic genes, for a proliferation score (Hallmark E2F targets), and for all of
  them together (red). Bars = 95% CI.
- **b** Single markers against navitoclax (cancer type removed): each gene alone, the proliferation score, the six
  core BH3-only genes alone, and the BH3-only set without E2F1 and TFDP1.
- **c-e** How special is the BH3-only set? Grey = 1,000 random gene sets of the same size, blue = all 672 Reactome
  pathways, red line = BH3-only. Panel titles give the rank among pathways, the share of random sets that do as
  well, and the extra variance explained (R2) by adding the BH3-only score to a model with BCL2, BCL2L1, MCL1, E2F
  and cancer type.

{lines}
""")


if __name__ == "__main__":
    if "g5" in sys.argv:
        g5_cell_lines()
    else:
        main()
