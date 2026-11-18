"""G10: evidence-of-absence statistics for the 57-analysis survey.

Effect per analysis = signed log HR per SD of the frozen organoid score in treated TCGA patients (Cox, from
nar/05_A3_B1/tables/A3_master_results.tsv). Positive = expected direction (higher predicted resistance -> shorter
survival). Random-effects meta-analysis: REML tau^2, Hartung-Knapp-Sidik-Jonkman (HKSJ) CI, 95% prediction interval.
Also DerSimonian-Laird, per drug class, one-analysis-per-independent-patient-cluster sensitivity, the median-split HR,
the score x treatment interaction HR, and equivalence (TOST-style) bounds.
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, t as tdist

from g_common import OK, TAB, C, drug_class, master, save_fig, short_aid


def re_meta(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float)
    k = len(y)
    w = 1 / se ** 2
    mu_fe = np.sum(w * y) / w.sum()
    Q = np.sum(w * (y - mu_fe) ** 2)
    c = w.sum() - np.sum(w ** 2) / w.sum()
    tau_dl = max(0.0, (Q - (k - 1)) / c) if k > 1 else 0.0
    tau = tau_dl
    for _ in range(500):  # REML fixed-point
        wr = 1 / (se ** 2 + tau)
        mu = np.sum(wr * y) / wr.sum()
        new = max(0.0, np.sum(wr ** 2 * ((y - mu) ** 2 - se ** 2)) / np.sum(wr ** 2) + 1 / wr.sum())
        if abs(new - tau) < 1e-10:
            tau = new
            break
        tau = new
    wr = 1 / (se ** 2 + tau)
    mu = np.sum(wr * y) / wr.sum()
    se_mu = np.sqrt(1 / wr.sum())
    # HKSJ
    q_hk = np.sum(wr * (y - mu) ** 2) / (k - 1) if k > 1 else 1.0
    se_hk = se_mu * np.sqrt(max(q_hk, 1.0))
    tq = tdist.ppf(0.975, k - 1) if k > 1 else np.nan
    t90 = tdist.ppf(0.95, k - 1) if k > 1 else np.nan
    pi_t = tdist.ppf(0.975, k - 2) if k > 2 else np.nan
    pi_half = pi_t * np.sqrt(tau + se_hk ** 2) if k > 2 else np.nan
    wd = 1 / (se ** 2 + tau_dl)
    mu_dl = np.sum(wd * y) / wd.sum()
    se_dl = np.sqrt(1 / wd.sum())
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return {"k": k, "HR": np.exp(mu), "lo_HKSJ": np.exp(mu - tq * se_hk), "hi_HKSJ": np.exp(mu + tq * se_hk),
            "p_HKSJ": 2 * tdist.sf(abs(mu / se_hk), k - 1) if k > 1 else np.nan,
            "lo90": np.exp(mu - t90 * se_hk), "hi90": np.exp(mu + t90 * se_hk),
            "lo_z": np.exp(mu - 1.96 * se_mu), "hi_z": np.exp(mu + 1.96 * se_mu),
            "tau2_REML": tau, "tau_REML": np.sqrt(tau), "I2": I2, "Q": Q, "Q_p": chi2.sf(Q, k - 1) if k > 1 else np.nan,
            "PI_lo": np.exp(mu - pi_half) if k > 2 else np.nan, "PI_hi": np.exp(mu + pi_half) if k > 2 else np.nan,
            "HR_FE": np.exp(mu_fe), "HR_DL": np.exp(mu_dl), "lo_DL": np.exp(mu_dl - 1.96 * se_dl),
            "hi_DL": np.exp(mu_dl + 1.96 * se_dl), "tau2_DL": tau_dl}


def effects(M, col, lo, hi):
    ok = M[col].between(1e-3, 1e3) & M[lo].gt(0) & np.isfinite(M[hi])
    d = M[ok].copy()
    d["y"] = np.log(d[col])
    d["se"] = (np.log(d[hi]) - np.log(d[lo])) / (2 * 1.96)
    return d[d.se.between(1e-6, 5)]


def main():
    M = master()
    M["drug_class"] = M.drug.map(drug_class)
    cl = pd.read_csv(TAB / "G9_per_analysis_permutation.tsv", sep="\t").set_index("analysis_id")["cluster"]
    M["cluster"] = M.analysis_id.map(cl)
    rows = []
    for label, col, lo, hi in (("continuous score per SD (primary)", "cont_HR", "cont_lo", "cont_hi"),
                               ("median split", "hr", "hr_lo", "hr_hi"),
                               ("adjusted continuous score per SD", "adj_cont_HR", "adj_cont_lo", "adj_cont_hi"),
                               ("score x treatment interaction per SD", "int_HR", "int_lo", "int_hi")):
        d = effects(M, col, lo, hi)
        rows.append({"effect": label, "subset": "all analyses", **re_meta(d.y, d.se)})
        if label.startswith("continuous score per SD (primary)") or label.startswith("score x"):
            for k, g in d.groupby("drug_class"):
                if len(g) >= 2:
                    rows.append({"effect": label, "subset": f"drug class: {k}", **re_meta(g.y, g.se)})
            # one analysis per independent cluster: most events
            one = d.sort_values("treated_events", ascending=False).drop_duplicates("cluster")
            rows.append({"effect": label, "subset": f"one analysis per patient cluster ({len(one)})", **re_meta(one.y, one.se)})
            # leave out the paper-era files
            rows.append({"effect": label, "subset": "excluding the two paper-era file sets",
                         **re_meta(*d[~d.project.eq("COAD_authors") & ~d.analysis_id.str.startswith("BLCA_lee2018_paper")][["y", "se"]].T.values)})
            for k, g in d.groupby(d.project.replace({"COAD_authors": "COAD"})):
                if len(g) >= 3:
                    rows.append({"effect": label, "subset": f"TCGA cohort: {k}", **re_meta(g.y, g.se)})
    R = pd.DataFrame(rows)
    for b in (1.1, 1.25, 1.5, 2.0):
        R[f"rules_out_HR_ge_{b}"] = R["hi90"] < b      # TOST-style: one-sided 95% upper bound below b
        R[f"rules_out_HR_le_{1/b:.2f}"] = R["lo90"] > 1 / b
    R.to_csv(TAB / "G10_random_effects_meta.tsv", sep="\t", index=False)

    # per-analysis equivalence
    d = effects(M, "cont_HR", "cont_lo", "cont_hi")
    d["hi90"] = np.exp(d.y + 1.645 * d.se)
    d["lo90"] = np.exp(d.y - 1.645 * d.se)
    for b in (1.25, 1.5, 2.0):
        d[f"excludes_HR_ge_{b}"] = d.hi90 < b
        d[f"excludes_both_HR_ge_{b}_and_le_{1/b:.2f}"] = (d.hi90 < b) & (d.lo90 > 1 / b)
    d["smallest_HR_ruled_out_(one-sided 95%)"] = d.hi90
    keep = ["analysis_id", "drug", "drug_class", "project", "treated_n", "treated_events", "cont_HR", "cont_lo", "cont_hi",
            "lo90", "hi90"] + [c for c in d.columns if c.startswith("excludes")]
    d[keep].to_csv(TAB / "G10_per_analysis_equivalence.tsv", sep="\t", index=False)
    cnt = pd.DataFrame([{"bound": b, "n_analyses": len(d), "n_upper_bound_below": int((d.hi90 < b).sum()),
                         "n_two_sided_equivalent": int(((d.hi90 < b) & (d.lo90 > 1 / b)).sum())} for b in (1.25, 1.5, 2.0, 3.0)])
    cnt.to_csv(TAB / "G10_equivalence_counts.tsv", sep="\t", index=False)
    print(R[["effect", "subset", "k", "HR", "lo_HKSJ", "hi_HKSJ", "p_HKSJ", "hi90", "tau2_REML", "I2", "PI_lo", "PI_hi"]].to_string())
    print(cnt)
    figure(M, d, R)


def figure(M, d, R):
    plt = C.plot_style()
    fig, axs = plt.subplots(1, 2, figsize=(11, 10.5), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    d = d.sort_values(["drug_class", "cont_HR"]).reset_index(drop=True)
    cols = dict(zip(sorted(d.drug_class.unique()), [OK["blue"], OK["orange"], OK["green"], OK["purple"], OK["verm"], OK["sky"], OK["grey"]]))
    for i, r in d.iterrows():
        ax.plot([r.cont_lo, r.cont_hi], [i, i], color=cols[r.drug_class], lw=1)
        ax.plot(r.cont_HR, i, "s", color=cols[r.drug_class], ms=3 + 4 * np.sqrt(r.treated_events / d.treated_events.max()))
    main = R[(R.effect.str.startswith("continuous score per SD (primary)")) & (R.subset == "all analyses")].iloc[0]
    y0 = len(d) + 1
    ax.plot([main.PI_lo, main.PI_hi], [y0, y0], color=OK["grey"], lw=4, alpha=.5, label="95% prediction interval")
    ax.plot([main.lo_HKSJ, main.hi_HKSJ], [y0, y0], color=OK["black"], lw=2)
    ax.plot(main.HR, y0, "D", color=OK["black"], ms=7, label="pooled (random effects)")
    ax.set_yticks(list(range(len(d))) + [y0])
    ax.set_yticklabels([short_aid(a) for a in d.analysis_id] + ["Pooled, 57 analyses"], fontsize=6.5)
    ax.axvline(1, color="k", lw=.7, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("HR per SD of organoid score, treated patients (log scale)\n>1 = expected direction")
    ax.set_xlim(0.03, 15)
    ax.set_xticks([0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10]); ax.set_xticklabels(["0.05", "0.1", "0.25", "0.5", "1", "2", "4", "10"]); ax.minorticks_off()
    for k, c in cols.items():
        ax.plot([], [], "s", color=c, label=k)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=3, frameon=False)
    C.letter(ax, "a", x=-0.45)
    ax = axs[1]
    sub = R[R.effect.str.startswith("continuous score per SD (primary)")].reset_index(drop=True)
    sub2 = R[R.effect.ne("continuous score per SD (primary)") & R.subset.eq("all analyses")]
    S = pd.concat([sub, sub2], ignore_index=True)
    lab = [(s if e.startswith("continuous score per SD (primary)") else e + ": all") for e, s in zip(S.effect, S.subset)]
    for i, r in S.iterrows():
        ax.plot([r.PI_lo, r.PI_hi], [i, i], color=OK["grey"], lw=3.5, alpha=.4)
        ax.plot([r.lo_HKSJ, r.hi_HKSJ], [i, i], color=OK["blue"], lw=1.6)
        ax.plot(r.HR, i, "D", color=OK["blue"], ms=5)
    ax.set_yticks(range(len(S)))
    ax.set_yticklabels([f"{l} (k={int(k)})" for l, k in zip(lab, S.k)], fontsize=7)
    ax.invert_yaxis()
    ax.axvline(1, color="k", lw=.7, ls="--")
    for b in (1.25, 1.5):
        ax.axvline(b, color=OK["verm"], lw=.6, ls=":")
        ax.axvline(1 / b, color=OK["verm"], lw=.6, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(0.2, 5)
    ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"]); ax.minorticks_off()
    ax.set_xlabel("Pooled HR (blue: 95% CI, HKSJ; grey: 95% prediction interval)\nred dotted lines: HR 1.25 and 1.5 (and inverse)")
    C.letter(ax, "b", x=-0.9)
    fig.tight_layout()
    cap = f"""**Figure G10. Pooled effect of the organoid scores across the 57 analyses.**
(a) Each row is one lab-model set x drug x TCGA cohort analysis. The square is the hazard ratio (HR) per standard
deviation (SD) of the frozen organoid score in treated patients; the line is its 95% confidence interval. HR > 1 means
patients predicted to be resistant died sooner (the expected direction). Bigger squares = more deaths. Colours = drug
class. Bottom row: random-effects pooled HR {main.HR:.2f} (95% CI {main.lo_HKSJ:.2f}-{main.hi_HKSJ:.2f}, Hartung-Knapp
method) and the 95% prediction interval {main.PI_lo:.2f}-{main.PI_hi:.2f} (grey; where the true effect of a new,
similar analysis would fall).
(b) Pooled HR for subsets (drug class, TCGA cohort, one analysis per group of overlapping patients) and for other
effect measures (median split, covariate-adjusted, score x treatment interaction). Red dotted lines mark HR 1.25 and
1.5 and their inverses, used as equivalence bounds. "Random effects" = the true effect is allowed to differ between
analyses; the amount of difference is tau^2.
"""
    save_fig(fig, "G10_meta_forest", cap)


if __name__ == "__main__":
    main()
