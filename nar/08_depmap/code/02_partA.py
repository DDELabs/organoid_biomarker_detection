"""Part A: DepMap / cell-line validation of the colon 5-FU biomarker (BH3-only activation).

Writes tables to nar/08_depmap/tables/A_*.tsv and main figures fig1-fig4 (panels a-h).
python -I nar/08_depmap/code/02_partA.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from common import BH3, CUR, OI, SEED, TAB, lineage_dummies, panel, savefig, style  # noqa: E402
from data import SCREENS, crispr, frozen_signatures, organoid_37, screen  # noqa: E402
from stats import bh, partial_spearman, random_signature_null, sig_score, spearman  # noqa: E402

plt = style()
FU = "FLUOROURACIL"
COL = {"GDSC1": OI["blue"], "GDSC2": OI["orange"], "CTRP": OI["green"], "PRISM": OI["purple"]}
N_NULL = 1000


def fisher_ci(r, n):
    if not np.isfinite(r) or n < 5:
        return np.nan, np.nan
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


# --------------------------------------------------------------------------- (b)
def part_b():
    rows = []
    for n in SCREENS:
        S = screen(n)
        c = S.common(FU, S.colon)
        r, p, k = spearman(S.scores.loc[BH3, c], S.resp.loc[c, FU])
        lo, hi = fisher_ci(r, k)
        rows.append({"screen": n, "metric": S.metric, "lines": "colorectal", "n": k, "spearman_rho": r, "ci_low": lo,
                     "ci_high": hi, "p": p})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "A_b_bh3_vs_5fu_colon.tsv", sep="\t", index=False)
    return t


# --------------------------------------------------------------------------- (c)
def part_c():
    sigs = frozen_signatures()
    rows, nulls, preds = [], {}, {}
    for n in SCREENS:
        S = screen(n)
        c = S.common(FU, S.colon)
        y = S.resp.loc[c, FU]
        sc = S.scores[c]
        for name, w in sigs.items():
            s = sig_score(sc, w)
            k = int(sum(p in sc.index for p in w.index))
            r, p, m = spearman(s, y)
            null = random_signature_null(sc, k, y, n=N_NULL, seed=SEED + k)
            nulls[(n, name)] = null
            preds[(n, name)] = pd.DataFrame({"predicted_resistance": s, "observed": y})
            rows.append({"screen": n, "signature": name, "k_pathways": k, "n": m, "spearman_rho": r, "p": p,
                         "null_q025": float(np.quantile(null, 0.025)), "null_q975": float(np.quantile(null, 0.975)),
                         "random_null_p_one_sided": float((1 + (null >= r).sum()) / (1 + len(null)))})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "A_c_frozen_model_colon_lines.tsv", sep="\t", index=False)
    pd.concat({f"{a}|{b}": v for (a, b), v in preds.items()}).to_csv(TAB / "A_c_predictions.tsv", sep="\t")
    return t, nulls, preds


# --------------------------------------------------------------------------- (d)
def part_d():
    w = organoid_37()
    rows, per = [], []
    for n in SCREENS:
        S = screen(n)
        c = S.common(FU, S.colon)
        y = S.resp.loc[c, FU]
        agree = tested = 0
        for pw, coef in w.items():
            if pw not in S.scores.index:
                continue
            r, p, k = spearman(S.scores.loc[pw, c], y)
            ok = np.sign(r) == np.sign(coef)
            per.append({"screen": n, "pathway": pw, "organoid_coef": coef, "cell_line_rho": r, "p": p, "n": k,
                        "sign_agrees": bool(ok)})
            agree += int(ok)
            tested += 1
        bt = stats.binomtest(agree, tested, 0.5, alternative="greater")
        rows.append({"screen": n, "pathways_tested": tested, "sign_agree": agree, "fraction": agree / tested,
                     "binomial_p_one_sided": bt.pvalue})
    t = pd.DataFrame(rows)
    allp = pd.DataFrame(per)
    tot_a, tot_t = int(t.sign_agree.sum()), int(t.pathways_tested.sum())
    t = pd.concat([t, pd.DataFrame([{"screen": "all four (not independent)", "pathways_tested": tot_t, "sign_agree": tot_a,
                                     "fraction": tot_a / tot_t,
                                     "binomial_p_one_sided": stats.binomtest(tot_a, tot_t, 0.5, alternative="greater").pvalue}])])
    t.to_csv(TAB / "A_d_sign_agreement_37.tsv", sep="\t", index=False)
    allp.to_csv(TAB / "A_d_sign_agreement_per_pathway.tsv", sep="\t", index=False)
    return t, allp


# --------------------------------------------------------------------------- (e)
def part_e():
    rows = []
    for n in SCREENS:
        S = screen(n)
        for d in S.resp.columns:
            c = S.common(d, S.colon)
            a = S.common(d)
            if len(a) < 30:
                continue
            r1, p1, k1 = spearman(S.scores.loc[BH3, c], S.resp.loc[c, d]) if len(c) >= 10 else (np.nan, np.nan, len(c))
            r2, p2, k2 = partial_spearman(S.scores.loc[BH3, a], S.resp.loc[a, d], S.lineage.loc[a])
            rows.append({"screen": n, "drug": d, "colon_rho": r1, "colon_p": p1, "colon_n": k1,
                         "pancancer_partial_rho": r2, "pancancer_p": p2, "pancancer_n": k2})
    t = pd.DataFrame(rows)
    out = []
    for n, g in t.groupby("screen"):
        g = g.copy()
        g["pancancer_q"] = bh(g["pancancer_p"].values)
        for col in ("colon_rho", "pancancer_partial_rho"):
            # rank 1 = most negative rho = BH3-only activation most linked to SENSITIVITY
            g[col.replace("rho", "rank")] = g[col].rank(method="min")
        out.append(g)
    t = pd.concat(out)
    t.to_csv(TAB / "A_e_bh3_vs_all_drugs.tsv", sep="\t", index=False)
    summ = []
    for n, g in t.groupby("screen"):
        f = g[g.drug == FU]
        if f.empty:
            continue
        f = f.iloc[0]
        summ.append({"screen": n, "drugs": len(g), "5FU_colon_rho": f.colon_rho,
                     "5FU_colon_rank_of": f"{int(f.colon_rank)}/{int(g.colon_rho.notna().sum())}",
                     "5FU_colon_percentile": float(f.colon_rank / g.colon_rho.notna().sum()),
                     "5FU_pancancer_rho": f.pancancer_partial_rho,
                     "5FU_pancancer_rank_of": f"{int(f.pancancer_partial_rank)}/{len(g)}",
                     "5FU_pancancer_percentile": float(f.pancancer_partial_rank / len(g)),
                     "drugs_pancancer_q<0.05": int((g.pancancer_q < 0.05).sum()),
                     "median_pancancer_rho_all_drugs": float(g.pancancer_partial_rho.median())})
    s = pd.DataFrame(summ)
    s.to_csv(TAB / "A_e_5fu_rank_summary.tsv", sep="\t", index=False)
    return t, s


# --------------------------------------------------------------------------- (f)
def bladder_model():
    """Organoid-trained bladder cisplatin model (Lee 2018 PDOs, curated), robust stability learner."""
    from obd import models as M
    from obd import network as N
    from obd.cohorts import to_canonical
    from obd.reference import gene_sets
    from obd import scoring
    e = pd.read_csv(CUR / "bladder_lee2018" / "expression.tsv.gz", sep="\t", index_col=0)
    e = to_canonical(e)
    r = pd.read_csv(CUR / "bladder_lee2018" / "response.tsv", sep="\t")
    y = r[r.drug == "CISPLATIN"].set_index("sample")["response"]
    y = y[y.index.isin(e.columns)]
    gs = {k: [g for g in v if g in e.index] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    sc = scoring.rank_score(e[y.index], gs)
    feats = [p for p in N.proximal([N.precomputed_proximity()["CISPLATIN"]]) if p in sc.index]
    X = sc.loc[feats].T.values
    w, table = M.robust_signature(X, y.values, feats, top_k=7, n_boot=200, seed=0)
    table.to_csv(TAB / "A_f_bladder_organoid_model_table.tsv", sep="\t")
    w.rename("weight").to_csv(TAB / "A_f_bladder_organoid_signature.tsv", sep="\t")
    return w, table, len(y), len(feats)


def part_f():
    w, table, n_org, n_feat = bladder_model()
    rows, preds, nulls = [], {}, {}
    for n in SCREENS:
        S = screen(n)
        drug = "CISPLATIN" if "CISPLATIN" in S.resp else "CARBOPLATIN"
        c = S.common(drug, S.bladder)
        if len(c) < 8:
            continue
        y = S.resp.loc[c, drug]
        sc = S.scores[c]
        s = sig_score(sc, w)
        k = int(sum(p in sc.index for p in w.index))
        r, p, m = spearman(s, y)
        null = random_signature_null(sc, k, y, n=N_NULL, seed=SEED + 7)
        nulls[n] = null
        preds[n] = pd.DataFrame({"predicted_resistance": s, "observed": y})
        rows.append({"screen": n, "drug": drug, "organoids_trained": n_org, "proximal_pathways": n_feat,
                     "signature_pathways": k, "n_bladder_lines": m, "spearman_rho": r, "p": p,
                     "null_q025": float(np.quantile(null, 0.025)), "null_q975": float(np.quantile(null, 0.975)),
                     "random_null_p_one_sided": float((1 + (null >= r).sum()) / (1 + len(null)))})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "A_f_bladder_cisplatin.tsv", sep="\t", index=False)
    return t, preds, nulls, w


# --------------------------------------------------------------------------- (g)
MIMETICS = ["NAVITOCLAX", "OBATOCLAX", "VENETOCLAX", "ABT-737"]


def part_g():
    rows = []
    for n in SCREENS:
        S = screen(n)
        for d in MIMETICS:
            if d not in S.resp:
                continue
            c, a = S.common(d, S.colon), S.common(d)
            r1, p1, k1 = spearman(S.scores.loc[BH3, c], S.resp.loc[c, d])
            r2, p2, k2 = partial_spearman(S.scores.loc[BH3, a], S.resp.loc[a, d], S.lineage.loc[a])
            for scope, r, p, k in (("colorectal", r1, p1, k1), ("pan-cancer, lineage-adjusted", r2, p2, k2)):
                lo, hi = fisher_ci(r, k)
                rows.append({"screen": n, "drug": d, "scope": scope, "n": k, "rho": r, "ci_low": lo, "ci_high": hi, "p": p})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "A_g_bh3_vs_bh3_mimetics.tsv", sep="\t", index=False)
    return t


# --------------------------------------------------------------------------- (h)
GENES = ["MCL1", "BCL2L1", "BCL2", "BCL2L11", "BBC3", "PMAIP1", "BID", "BAD", "BMF", "HRK", "BIK"]


def part_h():
    K = crispr()
    rows = []
    for n in SCREENS:
        S = screen(n)
        a = [i for i in S.resp[FU].dropna().index if isinstance(S.sidm.get(i), str) and S.sidm[i] in K.columns]
        c = [i for i in a if i in S.colon]
        for g in GENES:
            if g not in K.index:
                rows.append({"screen": n, "gene": g, "note": "not in Project Score"})
                continue
            dep = pd.Series(K.loc[g, [S.sidm[i] for i in a]].values, index=a)
            r1, p1, k1 = spearman(dep.loc[c], S.resp.loc[c, FU])
            r2, p2, k2 = partial_spearman(dep, S.resp.loc[a, FU], S.lineage.loc[a])
            rows.append({"screen": n, "gene": g, "colon_rho": r1, "colon_p": p1, "colon_n": k1,
                         "pancancer_partial_rho": r2, "pancancer_p": p2, "pancancer_n": k2,
                         "median_LFC_colon": float(dep.loc[c].median()) if c else np.nan})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "A_h_crispr_vs_5fu.tsv", sep="\t", index=False)
    return t


# --------------------------------------------------------------------------- figures
def fig1(tb):
    fig = plt.figure(figsize=(11, 6.6))
    gsp = fig.add_gridspec(2, 4, height_ratios=[1.05, 1], hspace=0.45, wspace=0.38)
    ax = fig.add_subplot(gsp[0, :])
    ax.axis("off")
    panel(ax, "a")
    boxes = [
        (0.02, "19 colorectal organoids\n(van de Wetering 2015)\n5-FU IC50 + expression"),
        (0.215, "Network step: 37 Reactome\npathways close to TYMS\n(STRING proximity z <= -1.28)"),
        (0.41, "Organoid model (frozen)\nBH3-only activation (-)\nhyaluronan (+), CO2 hydration (+)"),
        (0.605, "Cell lines: GDSC1, GDSC2,\nCTRPv2, PRISM + Cell Model\nPassports RNA-seq, CRISPR"),
        (0.80, "Tests: (b) BH3 score vs 5-FU\n(c) model unchanged (d) signs\n(e-h) specificity, bladder,\nBH3 mimetics, CRISPR"),
    ]
    cols = [OI["sky"], OI["yellow"], OI["orange"], OI["green"], OI["purple"]]
    for (x, txt), col in zip(boxes, cols):
        ax.add_patch(plt.Rectangle((x, 0.12), 0.18, 0.76, transform=ax.transAxes, facecolor=col, alpha=0.35,
                                   edgecolor="black", lw=0.8))
        ax.text(x + 0.09, 0.5, txt, transform=ax.transAxes, ha="center", va="center", fontsize=8)
    for x in (0.2, 0.395, 0.59, 0.785):
        ax.annotate("", xy=(x + 0.015, 0.5), xytext=(x, 0.5), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.set_title("Workflow: an organoid-trained 5-FU biomarker tested, unchanged, in cancer cell lines", fontsize=10)
    first = True
    for j, n in enumerate(SCREENS):
        S = screen(n)
        c = S.common(FU, S.colon)
        ax = fig.add_subplot(gsp[1, j])
        if first:
            panel(ax, "b")
            first = False
        x, y = S.scores.loc[BH3, c], S.resp.loc[c, FU]
        ax.scatter(x, y, s=16, color=COL[n], edgecolor="black", lw=0.3)
        b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 10)
        ax.plot(xs, np.polyval(b, xs), color="black", lw=0.8)
        row = tb[tb.screen == n].iloc[0]
        ax.set_title(f"{n} (n={row.n})\nrho={row.spearman_rho:.2f}, p={row.p:.2g}", fontsize=9)
        ax.set_xlabel("BH3-only activation score")
        ax.set_ylabel("5-FU " + ("ln IC50" if S.metric == "ln_ic50" else "AUC") + " (higher = resistant)")
    savefig(fig, "fig1_workflow_bh3_vs_5fu", """
**Figure 1. Testing the organoid 5-FU biomarker in colorectal cancer cell lines.**
(a) Workflow. The 5-FU model was trained on 19 colorectal organoids only. Its strongest pathway is "Activation of
BH3-only proteins" (BH3-only proteins are pro-death switches such as BIM, PUMA and NOXA). In organoids, a higher
score went with a LOWER IC50, that is, more sensitivity to 5-FU. We then score the same pathway in cell lines from
four drug screens without changing anything.
(b) Each dot is one colorectal cell line. x = BH3-only activation score from RNA-seq (single-sample rank score).
y = 5-FU response (GDSC: natural log of IC50; CTRPv2 and PRISM: area under the dose-response curve). Higher y means
more resistant. The organoid model predicts a downward slope (negative rho). Spearman rho and p are shown; the
line is a least-squares fit for the eye only. Numbers: `tables/A_b_bh3_vs_5fu_colon.tsv`.
""")
    plt.close(fig)


def fig2(tc, nulls, preds, td, allp):
    fig, axes = plt.subplots(2, 4, figsize=(11, 6.4))
    for j, n in enumerate(SCREENS):
        ax = axes[0, j]
        if j == 0:
            panel(ax, "c")
        d = preds[(n, "robust_v1")]
        ax.scatter(d.predicted_resistance, d.observed, s=16, color=COL[n], edgecolor="black", lw=0.3)
        row = tc[(tc.screen == n) & (tc.signature == "robust_v1")].iloc[0]
        ax.set_title(f"{n}: rho={row.spearman_rho:.2f}\nrandom 95%: {row.null_q025:.2f} to {row.null_q975:.2f}", fontsize=8.5)
        ax.set_xlabel("Predicted resistance\n(organoid model, unchanged)")
        ax.set_ylabel("Observed 5-FU response")
    # null band panel: all signatures x screens
    ax = axes[1, 0]
    panel(ax, "c")
    sigs = ["paper_top7", "robust_v1", "v1_rank", "nit_soft", "nit_deconf"]
    for j, n in enumerate(SCREENS):
        for i, s in enumerate(sigs):
            x = i + j * (len(sigs) + 1)
            row = tc[(tc.screen == n) & (tc.signature == s)].iloc[0]
            ax.plot([x, x], [row.null_q025, row.null_q975], color=OI["grey"], lw=4, alpha=0.6, solid_capstyle="butt")
            ax.scatter([x], [row.spearman_rho], color=COL[n], edgecolor="black", lw=0.3, s=18, zorder=3)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks([2 + j * 6 for j in range(4)])
    ax.set_xticklabels(SCREENS)
    ax.set_ylabel("Spearman rho, predicted vs observed")
    ax.set_title("5 frozen signatures (dots)\nvs random-pathway 95% band (grey)", fontsize=8.5)
    # (d)
    ax = axes[1, 1]
    panel(ax, "d")
    g = allp[allp.screen == "GDSC2"]
    ax.scatter(g.organoid_coef, g.cell_line_rho, c=[OI["blue"] if a else OI["red"] for a in g.sign_agrees], s=16,
               edgecolor="black", lw=0.3)
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Organoid weight (+ = resistance)")
    ax.set_ylabel("Colon-line rho with 5-FU\n(+ = resistance), GDSC2")
    ax.set_title("Blue = same sign, red = opposite", fontsize=8.5)
    ax = axes[1, 2]
    t = td[td.screen.isin(SCREENS)]
    ax.bar(range(len(t)), t.fraction, color=[COL[s] for s in t.screen], edgecolor="black", lw=0.5)
    ax.axhline(0.5, color="black", ls="--", lw=0.8)
    for i, r in enumerate(t.itertuples()):
        ax.text(i, r.fraction + 0.02, f"{r.sign_agree}/{r.pathways_tested}\np={r.binomial_p_one_sided:.2g}", ha="center",
                fontsize=7.5)
    ax.set_xticks(range(len(t)))
    ax.set_xticklabels(t.screen)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of pathways with same sign")
    ax.set_title("Sign agreement, 37 network-proximal\npathways (binomial test vs 0.5)", fontsize=8.5)
    axes[1, 3].axis("off")
    axes[1, 3].text(0, 0.5, "Reading guide\n\nc: a dot above the grey band\nmeans the organoid model\nbeats random pathway sets.\n\n"
                    "d: if organoid and cell-line\nbiology agree, more than half\nof the pathways should share\nthe sign (dashed line = chance).",
                    fontsize=8, va="center")
    fig.tight_layout()
    savefig(fig, "fig2_frozen_model_and_sign_agreement", """
**Figure 2. The organoid model applied unchanged to colorectal cell lines.**
(c, top row) Each dot is a colorectal cell line. x = predicted 5-FU resistance from the frozen organoid signature
"robust_v1" (BH3-only activation with a negative weight, hyaluronan uptake and CO2 hydration with positive weights);
y = observed 5-FU response (higher = more resistant). The title gives Spearman rho and the 2.5-97.5% range of rho for
1,000 random signatures with the same number of Reactome pathways and random signs ("random-pathway null band").
(c, bottom left) Same test for all five frozen signatures from the earlier external-validation work. Grey bar =
random 95% band; coloured dot = the real signature. A dot above the band would mean better than random.
(d) Sign agreement for the 37 pathways the network step kept for 5-FU. Left: organoid weight (from the robust
stability model) against the cell-line Spearman rho in GDSC2 colorectal lines. Right: fraction of pathways whose
sign agrees, per screen, with a one-sided binomial test against 50%. Numbers: `tables/A_c_*.tsv`, `tables/A_d_*.tsv`.
""")
    plt.close(fig)


def fig3(te, ts, tf, fpreds, fnulls):
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.6))
    for j, n in enumerate(["GDSC1", "GDSC2", "PRISM"]):
        ax = axes[0, j]
        if j == 0:
            panel(ax, "e")
        g = te[te.screen == n].sort_values("pancancer_partial_rho").reset_index(drop=True)
        ax.scatter(range(len(g)), g.pancancer_partial_rho, s=4, color=OI["grey"])
        hit = g[g.drug == FU]
        if not hit.empty:
            ax.scatter(hit.index, hit.pancancer_partial_rho, s=40, color=OI["red"], edgecolor="black", zorder=3, label="5-FU")
            ax.annotate(f"5-FU rank {int(hit.index[0]) + 1}/{len(g)}", (hit.index[0], hit.pancancer_partial_rho.iloc[0]),
                        xytext=(10, 14), textcoords="offset points", fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.5))
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("Drugs, sorted by association")
        ax.set_ylabel("BH3 score vs drug response\n(lineage-adjusted rho)")
        ax.set_title(f"{n}: {len(g)} drugs (negative = BH3-high lines\nmore sensitive)", fontsize=8.5)
    ax = axes[1, 0]   # audit I23: no second letter "e" here; the caption calls this "the bottom-left panel of (e)"
    g = te[te.screen == "GDSC2"].dropna(subset=["colon_rho"]).sort_values("colon_rho").reset_index(drop=True)
    ax.scatter(range(len(g)), g.colon_rho, s=4, color=OI["grey"])
    hit = g[g.drug == FU]
    ax.scatter(hit.index, hit.colon_rho, s=40, color=OI["red"], edgecolor="black", zorder=3)
    ax.annotate(f"5-FU rank {int(hit.index[0]) + 1}/{len(g)}", (hit.index[0], hit.colon_rho.iloc[0]), xytext=(10, 14),
                textcoords="offset points", fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Drugs, sorted")
    ax.set_ylabel("rho in colorectal lines")
    ax.set_title("GDSC2 colorectal lines only", fontsize=8.5)
    for j, n in enumerate(list(fpreds)[:2]):
        ax = axes[1, j + 1]
        if j == 0:
            panel(ax, "f")
        d = fpreds[n]
        row = tf[tf.screen == n].iloc[0]
        ax.scatter(d.predicted_resistance, d.observed, s=18, color=COL[n], edgecolor="black", lw=0.3)
        ax.set_title(f"Bladder lines, {n} {row.drug.lower()} (n={row.n_bladder_lines})\nrho={row.spearman_rho:.2f}, "
                     f"random 95%: {row.null_q025:.2f} to {row.null_q975:.2f}", fontsize=8.5)
        ax.set_xlabel("Predicted resistance (bladder organoid model)")
        ax.set_ylabel("Observed response (higher = resistant)")
    fig.tight_layout()
    savefig(fig, "fig3_drug_specificity_and_bladder", """
**Figure 3. Is the BH3-only signal specific to 5-FU, and does the bladder organoid model carry over?**
(e) For every drug in a screen we correlate the BH3-only activation score with drug response across all cell lines,
after removing differences between cancer types (lineage-adjusted partial Spearman rho). Drugs are sorted from the
most negative rho (BH3-high lines more sensitive) to the most positive. 5-FU is the red dot; its rank is printed.
The bottom-left panel (part of e, no separate letter) repeats this in GDSC2 colorectal lines only. Numbers: `tables/A_e_*.tsv`.
(f) A cisplatin model was trained on 9 bladder cancer organoids (Lee 2018) with the same robust stability learner
and the cisplatin-proximal Reactome pathways, then applied unchanged to bladder cancer cell lines. Each dot is a line.
The title gives Spearman rho and the random-pathway 95% band (1,000 random signatures of equal size).
Numbers: `tables/A_f_*.tsv`.
""")
    plt.close(fig)


def fig4(tg, th):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    panel(ax, "g")
    t = tg.reset_index(drop=True)
    labels = []
    for i, r in t.iterrows():
        mk = "o" if r.scope == "colorectal" else "s"
        ax.errorbar(r.rho, i, xerr=[[r.rho - r.ci_low], [r.ci_high - r.rho]], fmt=mk, color=COL[r.screen], ms=4, lw=0.8)
        labels.append(f"{r.drug.title()} {r.screen} ({'CRC' if r.scope == 'colorectal' else 'pan'}, n={r.n})")
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman rho: BH3-only score vs BH3-mimetic response\n(negative = BH3-high lines more sensitive)")
    ax.set_title("BH3-only score vs BH3 mimetics (circle = colorectal, square = pan-cancer adjusted)", fontsize=8.5)
    ax = axes[1]
    panel(ax, "h")
    m = th.dropna(subset=["colon_rho"]).pivot(index="gene", columns="screen", values="colon_rho").reindex(GENES)
    p = th.dropna(subset=["colon_rho"]).pivot(index="gene", columns="screen", values="pancancer_partial_rho").reindex(GENES)
    both = pd.concat({"CRC": m, "pan": p}, axis=1)
    both.columns = [f"{b}\n{a}" for a, b in both.columns]
    im = ax.imshow(both.values, cmap="PuOr_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(range(both.shape[1]))
    ax.set_xticklabels(both.columns, fontsize=6.5)
    ax.set_yticks(range(len(both)))
    ax.set_yticklabels(both.index)
    for i in range(both.shape[0]):
        for j in range(both.shape[1]):
            v = both.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="rho: CRISPR knockout effect vs 5-FU response")
    ax.set_title("CRISPR dependency (Project Score) vs 5-FU response", fontsize=8.5)
    fig.tight_layout()
    savefig(fig, "fig4_bh3_mimetics_and_crispr", """
**Figure 4. Mechanism checks: BH3 mimetics and CRISPR knockouts.**
(g) If a high BH3-only activation score means the cell is "primed" for apoptosis (close to the death threshold),
such lines should also respond better to BH3 mimetics, drugs that block BCL-2 family survival proteins
(navitoclax/ABT-263, obatoclax, venetoclax, ABT-737). Each row is one drug in one screen; circles = colorectal lines,
squares = all lines with cancer type removed. Bars are 95% confidence intervals (Fisher z). Negative rho = BH3-high
lines more sensitive. Numbers: `tables/A_g_bh3_vs_bh3_mimetics.tsv`.
(h) CRISPR knockout effect (Project Score, Sanger + Broad combined, batch-corrected log fold change; more negative =
the line needs the gene more) of the anti-apoptotic genes MCL1, BCL2L1 (BCL-XL), BCL2 and of BH3-only genes, against
5-FU response. Columns: screen and scope (CRC = colorectal lines, pan = all lines adjusted for cancer type). A positive
rho means lines that depend more on the gene (more negative knockout effect) are more SENSITIVE to 5-FU.
Numbers (with n per cell): `tables/A_h_crispr_vs_5fu.tsv`.
""")
    plt.close(fig)


def main():
    tb = part_b()
    print(tb.round(3).to_string())
    fig1(tb)
    tc, nulls, preds = part_c()
    print(tc.round(3).to_string())
    td, allp = part_d()
    print(td.round(3).to_string())
    fig2(tc, nulls, preds, td, allp)
    te, ts = part_e()
    print(ts.round(3).to_string())
    tf, fpreds, fnulls, w = part_f()
    print(w)
    print(tf.round(3).to_string())
    fig3(te, ts, tf, fpreds, fnulls)
    tg = part_g()
    print(tg.round(3).to_string())
    th = part_h()
    print(th.round(3).to_string())
    fig4(tg, th)


if __name__ == "__main__":
    main()
