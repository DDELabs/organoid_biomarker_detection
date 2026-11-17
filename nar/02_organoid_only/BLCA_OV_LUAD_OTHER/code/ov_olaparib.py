"""User's ovarian attempt: GSE229119, olaparib time course in 3 HGSOC organoids and 4 cell lines.

This is NOT a drug-response panel (every sample got the same drug; there is no 'sensitive vs resistant' measure),
so the paper's regression step cannot run. What can be done:
 1. score Reactome pathways (ssGSEA) in every sample; change = olaparib minus the same model's DMSO control;
 2. describe how pathways change over time (day 3, day 7, after recovery) in the organoids, highlighting the pathways
    closest to olaparib's targets (PARP1/2/3) in the paper's network;
 3. build an 'olaparib-response pathway signature' from the organoids (pathways that move the same way in all 3 at day 7)
    and ask whether it relates to platinum/PARP-relevant outcomes in TCGA-OV: overall survival and platinum-free interval
    in platinum-treated patients, and BRCA1/2 somatic mutation status.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr, ttest_1samp  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import common as K  # noqa: E402
import engine as E  # noqa: E402

OUT = K.PARENT / "OV_GSE229119_olaparib"
FIG, TAB = OUT / "figures", OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)


def main():
    raw = pd.read_csv(K.RAW / "geo" / "GSE229119_log2TPM.gene.expr.table.txt.gz", sep="\t", index_col=0)
    expr = K.C.to_canonical(raw)
    texpr = K.tcga_expr("OV")
    osc, tsc = K.paired_scores("OV_GSE229119", expr, texpr)
    # ssGSEA raw ES are on an arbitrary per-pathway scale: express every pathway in units of its SD across the 28 samples
    osc = osc.div(osc.std(axis=1).replace(0, np.nan), axis=0).dropna()
    models = sorted({c.split("_")[0] for c in osc.columns})
    times = {"D3": "olaparib day 3", "D7": "olaparib day 7", "R": "after recovery"}
    delta = {}
    for m in models:
        cols = [c for c in osc.columns if c.startswith(m + "_")]
        ctrl = [c for c in cols if "_NC_" in c][0]
        for c in cols:
            if c == ctrl:
                continue
            tp = c.split("_")[-1]
            key = "D3" if tp == "D3" else "D7" if tp == "D7" else "R"
            if m == "OVCA432" and tp == "R3":  # OVCA432 has no day 7: day3 then recover 3 / 7 days
                key = "R3"
            delta[(m, key)] = osc[c] - osc[ctrl]
    D = pd.DataFrame(delta)
    D.columns = [f"{m}|{k}" for m, k in D.columns]
    D.to_csv(TAB / "pathway_change_vs_DMSO.tsv", sep="\t")
    orgs = [m for m in models if m.startswith("ORG")]
    z, src = K.proximity("OLAPARIB")
    prox = [p for p in z.index[z <= -1.2816] if p in D.index]
    # signature: day-7 change consistent in all 3 organoids
    d7 = D[[f"{o}|D7" for o in orgs]]
    sd = d7.values.std()
    rows = []
    for p in D.index:
        v = d7.loc[p].values
        rows.append(dict(pathway=p, mean_change_d7=v.mean(), same_sign_all3=bool(np.all(v > 0) or np.all(v < 0)),
                         t_p=ttest_1samp(v, 0).pvalue, proximity_z=z.get(p, np.nan),
                         mean_change_d3=D.loc[p, [f"{o}|D3" for o in orgs]].mean(), mean_change_recovery=D.loc[p, [f"{o}|R" for o in orgs]].mean()))
    S = pd.DataFrame(rows).set_index("pathway")
    S["effect_vs_noise"] = S.mean_change_d7 / sd
    S.sort_values("mean_change_d7").to_csv(TAB / "organoid_olaparib_pathway_changes.tsv", sep="\t")
    cons = S[S.same_sign_all3]
    up = list(cons.sort_values("mean_change_d7", ascending=False).index[:15])
    dn = list(cons.sort_values("mean_change_d7").index[:15])
    up = [p for p in up if S.loc[p, "mean_change_d7"] > 0]
    dn = [p for p in dn if S.loc[p, "mean_change_d7"] < 0]
    json.dump({"up": up, "down": dn}, open(TAB / "olaparib_signature_pathways.json", "w"), indent=1)

    # ---- TCGA-OV association
    cdr = K.cdr()
    surv = K.survival_current("OV")
    pats = [p for p in tsc.columns if p in surv.index]
    P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[:, pats].T), index=pats, columns=tsc.index)
    score = P[up].mean(axis=1) - P[dn].mean(axis=1)
    prox_score = P[[p for p in prox if p in P.columns]].mean(axis=1) if prox else None
    plat = K.treated("OV", "CARBOPLATIN") | K.treated("OV", "CISPLATIN")
    tp = [p for p in pats if p in plat]
    T, Ev = np.clip(surv.loc[tp, "T"].values, 0, None), surv.loc[tp, "E"].values.astype(int)
    r_os, g_os = E.surv_test(score.loc[tp].values, T, Ev)
    c2 = cdr.reindex(tp)
    ok = c2["PFI.time"].notna() & c2["PFI"].notna()
    Tp, Ep = (c2.loc[ok, "PFI.time"].values / (365 / 12)).astype(float), c2.loc[ok, "PFI"].values.astype(int)
    r_pfi, g_pfi = E.surv_test(score.loc[c2.index[ok]].values, np.clip(Tp, 0, None), Ep)
    # platinum-resistant (progression/recurrence within 6 months of platinum: PFI.time < 183 d with event) vs sensitive
    res6 = c2.loc[ok & ((c2["PFI"] == 1) & (c2["PFI.time"] < 183) | (c2["PFI.time"] >= 183))]
    lab6 = ((res6["PFI"] == 1) & (res6["PFI.time"] < 183)).astype(int)
    a, b = score.loc[lab6.index[lab6 == 1]], score.loc[lab6.index[lab6 == 0]]
    p6 = mannwhitneyu(a, b).pvalue
    # BRCA1/2 somatic mutations (open data; germline status is controlled access)
    m = pd.read_csv(K.RAW / "TCGA-OV.somaticmutation_wxs.tsv.gz", sep="\t", usecols=["sample", "gene", "effect"])
    m = m[m["sample"].str[13:15] == "01"]
    m["patient"] = m["sample"].str[:12]
    prof = set(m.patient) & set(score.index)
    silent = {"synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "upstream_gene_variant",
              "downstream_gene_variant", "intergenic_region", "non_coding_transcript_exon_variant"}
    brca = set(m.loc[m.gene.isin(["BRCA1", "BRCA2"]) & ~m.effect.isin(silent), "patient"]) & prof
    sb, sw = score.loc[sorted(brca)], score.loc[sorted(prof - brca)]
    pb = mannwhitneyu(sb, sw).pvalue
    res = {"n_organoids": len(orgs), "n_models": len(models), "proximal_pathways_olaparib": len(prox),
           "signature_up": len(up), "signature_down": len(dn), "noise_sd_d7": sd,
           "os_platinum": r_os, "pfi_platinum": r_pfi, "platinum_resistant_n": int(len(a)), "platinum_sensitive_n": int(len(b)),
           "score_median_resistant": float(a.median()), "score_median_sensitive": float(b.median()), "p_resistant_vs_sensitive": p6,
           "brca_mut_n": len(sb), "brca_wt_n": len(sw), "score_median_brca": float(sb.median()), "score_median_wt": float(sw.median()),
           "p_brca": pb, "proximity_source": src}
    if prox_score is not None:
        r_px, _ = E.surv_test(prox_score.loc[tp].values, T, Ev)
        res["os_platinum_proximal_mean"] = r_px
    json.dump(res, open(TAB / "tcga_ov_association.json", "w"), indent=1, default=float)

    # ---- figure
    fig = plt.figure(figsize=(17, 15))
    gs = fig.add_gridspec(3, 3)
    ax = fig.add_subplot(gs[0, :2])
    show = up[:12] + dn[:12]
    cols = [c for o in models for c in D.columns if c.startswith(o + "|")]
    im = ax.imshow(D.loc[show, cols].values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels([E.short(p, 45) for p in show], fontsize=7)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("|", "\n") for c in cols], fontsize=7, rotation=90)
    ax.set_title("Change in pathway score after olaparib (minus the model's own DMSO control; units = SD of that pathway across all 28 samples)\nrows: the signature (12 up, 12 down at day 7 in all 3 organoids)")
    plt.colorbar(im, ax=ax, fraction=0.02)
    E.panel_letter(ax, "a")
    ax = fig.add_subplot(gs[0, 2])
    zz = z.sort_values().head(12)[::-1]
    ax.barh(range(len(zz)), [S.loc[p, "mean_change_d7"] if p in S.index else np.nan for p in zz.index],
            color=[E.OK["verm"] if (p in S.index and S.loc[p, "same_sign_all3"]) else E.OK["grey"] for p in zz.index])
    ax.set_yticks(range(len(zz)))
    ax.set_yticklabels([f"{E.short(p, 32)} (z={zz[p]:.1f})" for p in zz.index], fontsize=7)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Mean change at day 7 in 3 organoids")
    ax.set_title(f"Pathways closest to PARP1/2/3\n({len(prox)} with z <= -1.2816); red = same direction in all 3")
    E.panel_letter(ax, "b")
    ax = fig.add_subplot(gs[1, 0])
    for p, col in zip(up[:3] + dn[:3], [E.OK["verm"]] * 3 + [E.OK["blue"]] * 3):
        vals = [D.loc[p, [f"{o}|{k}" for o in orgs]].mean() for k in ("D3", "D7", "R")]
        ax.plot([0, 1, 2], [0] + vals[:2], "o-", color=col, alpha=0.8)
        ax.plot([2, 3], [vals[1], vals[2]], "o--", color=col, alpha=0.8)
        ax.text(3.05, vals[2], E.short(p, 28), fontsize=6.5, va="center")
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["DMSO", "day 3", "day 7", "recovered"])
    ax.set_xlim(-0.2, 4.6)
    ax.set_ylabel("Mean change in 3 organoids")
    ax.set_title("Time course of the 3 strongest up and down pathways")
    E.panel_letter(ax, "c")
    ax = fig.add_subplot(gs[1, 1])
    E.km_draw(ax, T, Ev, g_os, "TCGA-OV, platinum-treated: overall survival\nsplit by olaparib-response signature score", r_os)
    E.panel_letter(ax, "d")
    ax = fig.add_subplot(gs[1, 2])
    E.km_draw(ax, np.clip(Tp, 0, None), Ep, g_pfi, "Same patients: progression-free interval\n(high score = 'non-responder' label)", r_pfi)
    ax.set_ylabel("Fraction without progression")
    ax.set_xlabel("Progression-free interval (months)")
    E.panel_letter(ax, "e")
    ax = fig.add_subplot(gs[2, 0])
    ax.boxplot([b.values, a.values], widths=0.5, showfliers=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"platinum-sensitive\n(n={len(b)})", f"platinum-resistant\n(n={len(a)})"])
    ax.set_ylabel("Olaparib-response signature score")
    ax.set_title(f"Progression within 6 months vs later\nMann-Whitney p = {p6:.3g}")
    E.panel_letter(ax, "f")
    ax = fig.add_subplot(gs[2, 1])
    ax.boxplot([sw.values, sb.values], widths=0.5, showfliers=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"BRCA1/2 WT\n(n={len(sw)})", f"BRCA1/2 somatic mut\n(n={len(sb)})"])
    ax.set_ylabel("Olaparib-response signature score")
    ax.set_title(f"BRCA1/2 somatic mutation (PARP-inhibitor relevant)\nMann-Whitney p = {pb:.3g}")
    E.panel_letter(ax, "g")
    ax = fig.add_subplot(gs[2, 2])
    rho = D[[c for c in D.columns if c.endswith("|D7")]].corr(method="spearman")
    im = ax.imshow(rho.values, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(rho)))
    ax.set_xticklabels([c.split("|")[0] for c in rho.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(len(rho)))
    ax.set_yticklabels([c.split("|")[0] for c in rho.columns], fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.04)
    ax.set_title("How similar is the day-7 response between models?\n(Spearman r of pathway changes)")
    E.panel_letter(ax, "h")
    fig.tight_layout()
    stem = FIG / "OV_GSE229119_olaparib_timecourse"
    fig.savefig(f"{stem}.png", dpi=300)
    fig.savefig(f"{stem}.pdf")
    plt.close(fig)
    rho.to_csv(TAB / "day7_change_similarity_spearman.tsv", sep="\t")
    cap = f"""# Olaparib time course in ovarian organoids (GSE229119) and what it says about TCGA-OV

What the data are: 3 high-grade serous ovarian cancer organoids (2414, 2445, 3008) and 4 cell lines (Hey, MCF7 [breast],
OVCA432, SKOV3), each treated with olaparib (a PARP inhibitor) and sampled at day 3, day 7 and after the drug was removed.
There is one DMSO control per model. Every model got the drug, and no model has a 'how strongly did it respond' number,
so the paper's regression step (which needs a response value per model) cannot be run on these data.

**a.** Change in Reactome pathway score (ssGSEA) after olaparib, minus the same model's DMSO control, in units of that
pathway's standard deviation across all 28 samples. Rows are the
{len(up)} pathways that went up and {len(dn)} that went down at day 7 in all 3 organoids (the 'olaparib-response signature').
**b.** The pathways closest to olaparib's targets PARP1/2/3 in the paper's network ({len(prox)} pass z <= -1.2816) and how much
they changed at day 7.
**c.** Mean time course in the 3 organoids for the 3 strongest up and down pathways.
**d.** TCGA-OV patients given carboplatin or cisplatin, split at the median of the signature score (score = mean z of the up
pathways minus mean z of the down pathways). n = {r_os['n']}, deaths = {r_os['events']}, log-rank p = {r_os['p']:.3g}, HR (high vs
low score) = {r_os['hr']:.2f} ({r_os['hr_lo']:.2f}-{r_os['hr_hi']:.2f}).
**e.** Same split, progression-free interval (TCGA-CDR): n = {r_pfi['n']}, events = {r_pfi['events']}, p = {r_pfi['p']:.3g},
HR = {r_pfi['hr']:.2f} ({r_pfi['hr_lo']:.2f}-{r_pfi['hr_hi']:.2f}).
**f.** Score in patients who progressed within 6 months (platinum-resistant proxy, n = {len(a)}) vs later (n = {len(b)}):
medians {a.median():.3f} vs {b.median():.3f}, Mann-Whitney p = {p6:.3g}.
**g.** Score by BRCA1/2 somatic mutation (open data; germline BRCA status is controlled access): mutant n = {len(sb)}, wild-type
n = {len(sw)}, medians {sb.median():.3f} vs {sw.median():.3f}, p = {pb:.3g}.
**h.** Similarity of the day-7 change between models.

What this can answer: which pathways olaparib moves in ovarian organoids, how fast, and whether they move back after the
drug is removed; and whether tumours that already 'look like' an olaparib-treated organoid differ in platinum outcome.
What it cannot answer: which organoid or patient responds to olaparib (no response values, and no TCGA-OV patient got
olaparib), so no organoid-trained biomarker and no olaparib-treated survival test is possible.
"""
    open(f"{stem}.caption.md", "w").write(cap)
    print(json.dumps(res, default=float, indent=1)[:2000])


if __name__ == "__main__":
    main()
