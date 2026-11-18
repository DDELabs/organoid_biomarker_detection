"""G17: the 57-analysis survey repeated with progression-free interval (PFI, TCGA-CDR, Liu et al. 2018) instead of OS.

Same patients, same frozen scores and median splits as A3 (nar/05_A3_B1/_data/a3_patients_*.tsv.gz); only the
endpoint changes. PFI.time in days -> months (/30.44). Patients without a PFI record are dropped.
"""
import numpy as np
import pandas as pd
from scipy.stats import binom

from g_common import NAR, OK, TAB, C, load_patients, logrank_fast, master, save_fig, short_aid


def main():
    cdr = pd.read_csv(NAR / "02_organoid_only" / "_data" / "raw" / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv",
                      sep="\t", usecols=["_PATIENT", "PFI", "PFI.time", "OS", "OS.time"]).drop_duplicates("_PATIENT").set_index("_PATIENT")
    p = load_patients()
    M = master()
    rows = []
    for _, m in M.iterrows():
        d = p[p.analysis_id == m.analysis_id].copy()
        d["Tp"] = pd.to_numeric(cdr["PFI.time"].reindex(d.patient.values), errors="coerce").values / 30.44
        d["Ep"] = pd.to_numeric(cdr["PFI"].reindex(d.patient.values), errors="coerce").values
        r = {"analysis_id": m.analysis_id, "drug": m.drug, "project": m.project, "os_treated_n": m.treated_n,
             "os_treated_events": m.treated_events, "os_treated_p": m.treated_p, "os_cont_HR": m.cont_HR,
             "os_direction": m.direction, "os_untreated_p": m.untreated_p}
        for grp in ("treated", "untreated"):
            g = d[(d.group == grp) & d.Tp.notna() & d.Ep.notna()]
            if len(g) < 6:
                continue
            g = g[g.Tp >= 0]
            r[f"pfi_{grp}_n"] = len(g)
            r[f"pfi_{grp}_events"] = int(g.Ep.sum())
            r[f"pfi_{grp}_p"] = logrank_fast(g.Tp, g.Ep.astype(int), g.split)
            if grp == "treated":
                c = C.cox(g.Tp, g.Ep, g.split.values, ["s"])
                if c is not None:
                    r.update(pfi_HR=c.loc["s", "HR"], pfi_HR_lo=c.loc["s", "lo"], pfi_HR_hi=c.loc["s", "hi"])
                z = (g.score1 - g.score1.mean()) / g.score1.std()
                c = C.cox(g.Tp, g.Ep, z.values, ["s"])
                if c is not None:
                    r.update(pfi_cont_HR=c.loc["s", "HR"], pfi_cont_lo=c.loc["s", "lo"], pfi_cont_hi=c.loc["s", "hi"], pfi_cont_p=c.loc["s", "p"])
                e1 = g.loc[g.split == 1, "Ep"].mean()
                e0 = g.loc[g.split == 0, "Ep"].mean()
                hr = r.get("pfi_HR", np.nan)
                r["pfi_direction"] = ("expected" if hr > 1 else "opposite") if np.isfinite(hr) and 1e-3 < hr < 1e3 else ("expected" if e1 > e0 else "opposite")
                r["power_HR1.5_pfi"] = C.power_logrank(int(g.Ep.sum()), 1.5)
                r["MDHR_80_pfi"] = C.mdhr(int(g.Ep.sum()))
        rows.append(r)
    P = pd.DataFrame(rows)
    P["pfi_BH_q"] = C.bh(P.pfi_treated_p)
    P["pfi_Holm"] = C.holm(P.pfi_treated_p)
    P.to_csv(TAB / "G17_PFI_survey_57.tsv", sep="\t", index=False)
    out = []
    for fam, col in (("OS treated", "os_treated_p"), ("PFI treated", "pfi_treated_p"), ("OS untreated", "os_untreated_p"),
                     ("PFI untreated", "pfi_untreated_p")):
        v = P[col].dropna()
        o = int((v < .05).sum())
        out.append({"family": fam, "n_tests": len(v), "p_below_0.05": o, "expected": .05 * len(v),
                    "binomial_p": binom.sf(o - 1, len(v), .05), "BH_q_below_0.05": int((C.bh(v) < .05).sum()),
                    "min_p": v.min()})
    o = pd.DataFrame(out)
    o.loc[1, "expected_direction_among_hits"] = int(((P.pfi_treated_p < .05) & (P.pfi_direction == "expected")).sum())
    o.loc[0, "expected_direction_among_hits"] = int(((P.os_treated_p < .05) & P.os_direction.str.startswith("expected")).sum())
    ev = pd.DataFrame([{"median_events_OS": P.os_treated_events.median(), "median_events_PFI": P.pfi_treated_events.median(),
                        "median_ratio_PFI_over_OS": (P.pfi_treated_events / P.os_treated_events).median(),
                        "median_MDHR_OS": M.MDHR_80.median(), "median_MDHR_PFI": P.MDHR_80_pfi.median(),
                        "n_power80_HR1.5_PFI": int((P["power_HR1.5_pfi"] >= .8).sum()),
                        "n_same_side_p0.05_OS_and_PFI": int(((P.os_treated_p < .05) & (P.pfi_treated_p < .05)).sum()),
                        "spearman_logp_OS_PFI": P[["os_treated_p", "pfi_treated_p"]].corr("spearman").iloc[0, 1],
                        "spearman_cont_logHR_OS_PFI": np.log(P[["os_cont_HR", "pfi_cont_HR"]]).corr("spearman").iloc[0, 1],
                        "n_same_direction": int((P.os_direction.str.startswith("expected") == P.pfi_direction.eq("expected")).sum())}])
    o.to_csv(TAB / "G17_PFI_vs_OS_counts.tsv", sep="\t", index=False)
    ev.to_csv(TAB / "G17_PFI_events_power.tsv", sep="\t", index=False)
    print(o.to_string()); print(ev.T)
    print(P[(P.pfi_treated_p < .05) | (P.os_treated_p < .05)][["analysis_id", "os_treated_events", "os_treated_p", "pfi_treated_events",
                                                              "pfi_treated_p", "pfi_HR", "pfi_direction", "pfi_BH_q"]].to_string())
    plt = C.plot_style()
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    ax[0].scatter(P.os_treated_events, P.pfi_treated_events, s=14, color=OK["blue"])
    mx = P.pfi_treated_events.max() * 1.05
    ax[0].plot([0, mx], [0, mx], "k--", lw=.6)
    ax[0].set_xlabel("deaths (OS events), treated"); ax[0].set_ylabel("PFI events, treated")
    C.letter(ax[0], "a")
    h = (P.os_treated_p < .05) | (P.pfi_treated_p < .05)
    ax[1].scatter(-np.log10(P.os_treated_p), -np.log10(P.pfi_treated_p), s=14, color=np.where(h, OK["verm"], OK["grey"]))
    for _, r in P[h].iterrows():
        ax[1].annotate(short_aid(r.analysis_id), (-np.log10(r.os_treated_p), -np.log10(r.pfi_treated_p)), fontsize=5.5)
    ax[1].axhline(-np.log10(.05), color=OK["orange"], ls=":", lw=.7); ax[1].axvline(-np.log10(.05), color=OK["orange"], ls=":", lw=.7)
    ax[1].set_xlabel("-log10 p, overall survival"); ax[1].set_ylabel("-log10 p, progression-free interval")
    C.letter(ax[1], "b")
    ax[2].scatter(P.os_cont_HR, P.pfi_cont_HR, s=14, color=OK["green"])
    ax[2].set_xscale("log"); ax[2].set_yscale("log")
    ax[2].axhline(1, color="k", lw=.5); ax[2].axvline(1, color="k", lw=.5)
    ax[2].set_xlabel("HR per SD of score, OS"); ax[2].set_ylabel("HR per SD of score, PFI")
    C.letter(ax[2], "c")
    fig.tight_layout()
    r0, r1 = o.iloc[0], o.iloc[1]
    cap = f"""**Figure G17. The 57-analysis survey with progression-free interval (PFI) instead of overall survival (OS).**
PFI = time until the cancer grows, comes back, or the patient dies (TCGA Clinical Data Resource). Same patients, same
frozen organoid scores and median splits as the OS survey; only the outcome changes.
(a) Events per analysis: PFI gives more events than OS (median {ev.median_events_PFI[0]:.0f} vs {ev.median_events_OS[0]:.0f}).
(b) Treated-patient log-rank p with OS (x) and PFI (y); dotted lines = p 0.05; red/labelled = p < 0.05 with either
endpoint. OS: {int(r0['p_below_0.05'])} of 57 below 0.05; PFI: {int(r1['p_below_0.05'])} of {int(r1.n_tests)}
(about {r1.expected:.1f} expected by chance). (c) Hazard ratio per SD of the continuous score with OS and PFI.
"""
    save_fig(fig, "G17_PFI_vs_OS", cap)


if __name__ == "__main__":
    main()
