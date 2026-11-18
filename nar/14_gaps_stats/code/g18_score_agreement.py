"""G18: agreement between ssGSEA (organoid-only agents, A3 survey) and the within-sample rank score (training table,
used by B1, the registered test and the joint models) in TCGA.

(a) per pathway x TCGA cohort Spearman correlation, for every pathway in the organoid-only TCGA ssGSEA caches and the
    GI result pickles; (b) for the frozen rank-1 pathway of each of the 57 analyses; (c) re-run all 57 treated tests
    (and the 7 nominal hits) with the rank score in place of ssGSEA: same pathway, same sign of the Ridge coefficient,
    same median split rule, same patients and survival.
"""
import pickle

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from g_common import A3, DATA, NAR, OK, REPO, TAB, C, load_patients, logrank_fast, master, save_fig, short_aid

TT = REPO / "robust" / "results" / "training_table"
CACHE = NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "_data" / "cache"


def rank_scores():
    pt = pd.read_csv(TT / "patients.tsv.gz", sep="\t", usecols=["sample_id", "cohort", "patient"])
    pt = pt[pt.cohort.str.startswith("TCGA-")]
    R = pd.read_csv(TT / "pathway_scores_tcga.tsv.gz", sep="\t", index_col=0)
    R = R.loc[R.index.intersection(pt.sample_id)]
    m = pt.set_index("sample_id").loc[R.index]
    R.index = pd.MultiIndex.from_arrays([m.cohort.str.replace("TCGA-", ""), m.patient])
    R = R[~R.index.duplicated()]
    return R


def main():
    R = rank_scores()
    rows = []
    # (a) caches
    for f in sorted(CACHE.glob("ssgsea_*_tcga.tsv.gz")):
        tag = f.name.replace("ssgsea_", "").replace("_tcga.tsv.gz", "")
        proj = tag.split("_")[0] if not tag.startswith("PDX") else tag.split("_")[-1]
        if proj not in R.index.get_level_values(0):
            continue
        S = pd.read_csv(f, sep="\t", index_col=0)
        S.columns = [c[:12] for c in S.columns]
        S = S.loc[:, ~S.columns.duplicated()]
        Rp = R.loc[proj]
        pats = S.columns.intersection(Rp.index)
        for pw in S.index.intersection(Rp.columns):
            r = spearmanr(S.loc[pw, pats].values, Rp.loc[pats, pw].values)[0]
            rows.append({"source": tag, "project": proj, "pathway": pw, "n": len(pats), "spearman": r})
    # GI pickles
    for f in sorted((NAR / "02_organoid_only" / "_data" / "results").glob("*.pkl")):
        d = pickle.load(open(f, "rb"))
        P = d["Pall"]
        proj = d["meta"]["cohort"]
        parts = {"COAD_authors": ["COAD"], "COADREAD": ["COAD", "READ"]}.get(proj, [proj])
        Rp = pd.concat([R.loc[p] for p in parts if p in R.index.get_level_values(0)])
        Rp = Rp[~Rp.index.duplicated()]
        pats = P.index.intersection(Rp.index)
        for pw in P.columns.intersection(Rp.columns):
            r = spearmanr(P.loc[pats, pw].values, Rp.loc[pats, pw].values)[0]
            rows.append({"source": "GI:" + f.stem, "project": proj, "pathway": pw, "n": len(pats), "spearman": r})
    A = pd.DataFrame(rows)
    # one row per project x pathway (median over sources that score the same pathway in the same cohort)
    Au = A.groupby(["project", "pathway"]).agg(spearman=("spearman", "median"), n=("n", "max"), n_sources=("source", "nunique")).reset_index()
    Au.to_csv(TAB / "G18_pathway_correlation_ssgsea_vs_rank.tsv", sep="\t", index=False)
    summ = Au.groupby("project").spearman.describe(percentiles=[.1, .25, .5, .75, .9])
    allr = Au.spearman.describe(percentiles=[.1, .25, .5, .75, .9])
    summ.loc["ALL"] = allr
    summ["share_rho_above_0.7"] = list(Au.groupby("project").spearman.apply(lambda x: (x > .7).mean())) + [(Au.spearman > .7).mean()]
    summ["share_rho_above_0.5"] = list(Au.groupby("project").spearman.apply(lambda x: (x > .5).mean())) + [(Au.spearman > .5).mean()]
    summ.to_csv(TAB / "G18_correlation_summary_by_cohort.tsv", sep="\t")
    print(summ.to_string())

    # (b)+(c)
    p = load_patients()
    M = master()
    sig = pd.concat([pd.read_csv(A3 / "_data" / f"signatures_{k}.tsv", sep="\t") for k in ("gi", "other")])
    sig1 = sig[sig["rank"] == 1].set_index("analysis_id")
    out = []
    for _, m in M.iterrows():
        aid = m.analysis_id
        pw, coef = sig1.loc[aid, "pathway"], sig1.loc[aid, "coef"]
        d = p[p.analysis_id == aid]
        parts = {"COAD_authors": ["COAD"], "COADREAD": ["COAD", "READ"]}.get(m.project, [m.project])
        Rp = pd.concat([R.loc[q] for q in parts])
        Rp = Rp[~Rp.index.duplicated()]
        r = {"analysis_id": aid, "pathway": pw, "ssgsea_treated_p": m.treated_p, "direction_ssgsea": m.direction}
        if pw not in Rp.columns:
            r["note"] = "pathway not in the 673-pathway rank-score table"
            out.append(r)
            continue
        x_ss = d.score1 / coef
        x_rk = Rp[pw].reindex(d.patient.values).values
        ok = np.isfinite(x_rk)
        r["n_with_rank_score"] = int(ok.sum())
        r["n_total"] = len(d)
        r["spearman_ssgsea_vs_rank_patients"] = spearmanr(x_ss[ok], x_rk[ok])[0]
        t = d[d.group == "treated"].copy()
        t["xr"] = np.sign(coef) * Rp[pw].reindex(t.patient.values).values
        t = t[np.isfinite(t.xr)]
        t["split_rank"] = (t.xr > t.xr.median()).astype(int)
        r["treated_n_rank"] = len(t)
        r["split_agreement"] = float((t.split == t.split_rank).mean())
        r["rank_treated_p"] = logrank_fast(t["T"], t.E, t.split_rank)
        r["ssgsea_treated_p_same_patients"] = logrank_fast(t["T"], t.E, t.split)
        z = (t.xr - t.xr.mean()) / t.xr.std()
        c = C.cox(t["T"], t.E, z.values, ["s"])
        if c is not None:
            r.update(rank_cont_HR=c.loc["s", "HR"], rank_cont_lo=c.loc["s", "lo"], rank_cont_hi=c.loc["s", "hi"], rank_cont_p=c.loc["s", "p"])
        r["ssgsea_cont_HR"] = m.cont_HR
        e1 = t.loc[t.split_rank == 1, "E"].sum() / max((t.split_rank == 1).sum(), 1)
        e0 = t.loc[t.split_rank == 0, "E"].sum() / max((t.split_rank == 0).sum(), 1)
        r["direction_rank"] = "expected" if e1 > e0 else "opposite"
        u = d[d.group == "untreated"].copy()
        if len(u):
            u["xr"] = np.sign(coef) * Rp[pw].reindex(u.patient.values).values
            u = u[np.isfinite(u.xr)]
            r["rank_untreated_p"] = logrank_fast(u["T"], u.E, (u.xr > u.xr.median()).astype(int))
        out.append(r)
    O = pd.DataFrame(out)
    O["nominal_hit_ssgsea"] = O.ssgsea_treated_p < 0.05
    O.to_csv(TAB / "G18_rerun_57_with_rank_score.tsv", sep="\t", index=False)
    cnt = pd.DataFrame([{"family": "treated p<0.05", "ssgsea": int((O.ssgsea_treated_p < .05).sum()),
                         "rank_score": int((O.rank_treated_p < .05).sum()), "n_rank": int(O.rank_treated_p.notna().sum()),
                         "both": int(((O.ssgsea_treated_p < .05) & (O.rank_treated_p < .05)).sum())},
                        {"family": "untreated p<0.05", "ssgsea": int((M.untreated_p < .05).sum()),
                         "rank_score": int((O.rank_untreated_p < .05).sum()), "n_rank": int(O.rank_untreated_p.notna().sum()), "both": np.nan}])
    cnt.to_csv(TAB / "G18_hit_counts_ssgsea_vs_rank.tsv", sep="\t", index=False)
    print(cnt)
    print(O[O.nominal_hit_ssgsea][["analysis_id", "spearman_ssgsea_vs_rank_patients", "split_agreement", "ssgsea_treated_p",
                                   "rank_treated_p", "rank_cont_HR", "direction_ssgsea", "direction_rank"]].to_string())
    print(O[["spearman_ssgsea_vs_rank_patients", "split_agreement"]].describe())
    fig_(Au, O)


def fig_(Au, O):
    plt = C.plot_style()
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    projs = sorted(Au.project.unique())
    ax[0].boxplot([Au[Au.project == q].spearman.dropna() for q in projs], vert=True, showfliers=False)
    for i, q in enumerate(projs, 1):
        v = Au[Au.project == q].spearman.dropna()
        ax[0].scatter(np.random.default_rng(1).normal(i, .07, len(v)), v, s=3, alpha=.4, color=OK["blue"])
    ax[0].set_xticks(range(1, len(projs) + 1))
    ax[0].set_xticklabels(projs, rotation=60, fontsize=7)
    ax[0].set_ylabel("Spearman correlation, ssGSEA vs rank score\n(one dot = one pathway in one TCGA cohort)")
    ax[0].axhline(0, color="k", lw=.6)
    C.letter(ax[0], "a")
    o = O.dropna(subset=["rank_treated_p"])
    hit = o.ssgsea_treated_p < .05
    ax[1].scatter(-np.log10(o.ssgsea_treated_p[~hit]), -np.log10(o.rank_treated_p[~hit]), s=14, color=OK["grey"], label="other analyses")
    ax[1].scatter(-np.log10(o.ssgsea_treated_p[hit]), -np.log10(o.rank_treated_p[hit]), s=24, color=OK["verm"], label="7 ssGSEA hits (p<0.05)")
    lim = max(3.2, -np.log10(o[["ssgsea_treated_p", "rank_treated_p"]].min().min()) + .2)
    ax[1].plot([0, lim], [0, lim], "k--", lw=.6)
    ax[1].axhline(-np.log10(.05), color=OK["orange"], lw=.6, ls=":")
    ax[1].axvline(-np.log10(.05), color=OK["orange"], lw=.6, ls=":")
    ax[1].set_xlabel("-log10 p, ssGSEA score (A3 survey)")
    ax[1].set_ylabel("-log10 p, rank score (same pathway)")
    ax[1].legend(frameon=False, fontsize=7)
    C.letter(ax[1], "b")
    ax[2].scatter(o.spearman_ssgsea_vs_rank_patients, o.split_agreement, s=14, color=OK["blue"])
    ax[2].set_xlabel("Spearman, ssGSEA vs rank score\n(frozen rank-1 pathway, analysis patients)")
    ax[2].set_ylabel("share of treated patients in the same\nmedian-split group")
    C.letter(ax[2], "c")
    fig.tight_layout()
    med = Au.spearman.median()
    cap = f"""**Figure G18. Do the two pathway scores agree in TCGA?**
Two ways of scoring a pathway in a tumour were used in this project: ssGSEA (a gene-set enrichment score; used by the
organoid-only survey of 57 analyses) and a within-sample rank score (the average rank of the pathway's genes inside one
sample; used by the external validation, the registered test and the joint models).
(a) Spearman correlation between the two scores for each pathway in each TCGA cohort (median over all pathways
{med:.2f}). (b) Treated-patient log-rank p of each of the 57 analyses with ssGSEA (x) and with the rank score of the same
pathway (y); red = the 7 analyses with ssGSEA p < 0.05; dotted lines = p 0.05. (c) For each analysis, the correlation
of the two scores for its frozen pathway and the share of treated patients who fall in the same half of the median split.
"""
    save_fig(fig, "G18_score_agreement", cap)


if __name__ == "__main__":
    main()
