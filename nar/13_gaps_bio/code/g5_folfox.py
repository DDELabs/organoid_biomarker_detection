"""G5: the FOLFOX question. Is the paper's 5-FU marker (BH3-only pathway) a 5-FU, an oxaliplatin or a regimen marker?

(a) organoids: BH3-only vs oxaliplatin and 5-FU response (from G3 table);
(b) cell lines: colorectal and pan-cancer (from G5 cell-line table written by g7_bh3_beyond_markers.py g5);
(c) TCGA colon 5-FU patients: frozen BH3-only split within FOLFOX vs no-oxaliplatin patients, with a Cox
    score x oxaliplatin interaction (authors' 2019 files and current data);
(d) GSE39582 (Marisa 2013, stage II-III): frozen BH3-only score in 5FU, FUFOL, FOLFOX and untreated patients,
    RFS and OS, with score x regimen interactions.
Frozen score = the authors' Ridge sign x BH3-only score (higher = predicted more 5-FU-resistant).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import survdiff

from gcommon import BH3, FIG, OI, REPO, TAB, fmt_p, panel, savefig, style, version_rows
import g2_decompose as G2  # noqa: E402
from engine import Combo, split  # noqa: E402  (02 engine, path added by g2_decompose)

TT = REPO / "robust/results/training_table"


def cox(T, E, X, names):
    m = PHReg(T, X, status=E, ties="efron").fit(disp=False)
    ci = m.conf_int()
    return {f"{n}_HR": np.exp(m.params[i]) for i, n in enumerate(names)} | \
        {f"{n}_lo": np.exp(ci[i, 0]) for i, n in enumerate(names)} | \
        {f"{n}_hi": np.exp(ci[i, 1]) for i, n in enumerate(names)} | \
        {f"{n}_p": m.pvalues[i] for i, n in enumerate(names)}


def subset_stats(T, E, s, label, extra):
    out = dict(subset=label, n=len(T), events=int(E.sum()), **extra)
    if len(T) < 8 or E.sum() < 2:
        return out
    z = (s - s.mean()) / s.std()
    try:
        out.update({k.replace("x_", "HRsd_"): v for k, v in cox(T, E, z.reshape(-1, 1), ["x"]).items()})
    except Exception:
        pass
    g = split(s)
    out["events_resp"], out["events_nonresp"] = int(E[g == 0].sum()), int(E[g == 1].sum())
    try:
        out["p_logrank"] = float(survdiff(T, E, g)[1]) if 0 < E[g == 0].sum() + E[g == 1].sum() else np.nan
    except Exception:
        out["p_logrank"] = np.nan
    return out


def interaction(T, E, s, grp):
    z = (s - s.mean()) / s.std()
    X = np.column_stack([z, grp, z * grp])
    try:
        r = cox(T, E, X, ["score", "group", "inter"])
        return {k: r[k] for k in ("inter_HR", "inter_lo", "inter_hi", "inter_p")}
    except Exception:
        return {}


def tcga(I, sign):
    rows = []
    tt = I["tt"]
    for setup, F in (("authors' 2019 files", G2.AUTH), ("current data", G2.CURR)):
        osc, tsc, sur = I["org"][F["org"]], I["pat"][F["pat"]], I["sur"][F["sur"]]
        treated = set(I["treated"][F["parse"]])
        c = Combo(osc, I["resp"][F["resp"]], tsc, sur, treated, set(), I["z"])
        s = sign * c.P.loc[c.tp, BH3].values
        if setup.startswith("authors"):
            oxa = set(G2.O.authors_inputs("OXALIPLATIN")[4])
        else:
            oxa = set(tt.loc[tt.drug == "OXALIPLATIN", "patient"])
        ox = np.array([p in oxa for p in c.tp]).astype(float)
        for lab, sel in (("all 5-FU patients", np.ones(len(ox), bool)), ("5-FU + oxaliplatin (FOLFOX-type)", ox == 1),
                         ("5-FU without oxaliplatin", ox == 0)):
            # split at the median of ALL 5-FU patients (as the paper), then look inside each regimen
            g_all = split(s)
            r = subset_stats(c.T[sel], c.E[sel], s[sel], lab, dict(cohort="TCGA colon", setup=setup, endpoint="OS"))
            r["deaths_pred_resp_globalsplit"] = int(c.E[sel & (g_all == 0)].sum())
            r["deaths_pred_nonresp_globalsplit"] = int(c.E[sel & (g_all == 1)].sum())
            if lab.startswith("all"):
                r.update(interaction(c.T, c.E, s, ox))
            rows.append(r)
    return rows


def gse39582(sign):
    P = pd.read_csv(TT / "patients.tsv.gz", sep="\t")
    P = P[(P.cohort == "GSE39582")].set_index("sample_id")
    S = pd.read_csv(TT / "pathway_scores_geo.tsv.gz", sep="\t", index_col=0, usecols=["sample_id", BH3])
    P = P.join(S, how="inner")
    P = P[P.stage.isin([2, 3])]
    P["s"] = sign * (P[BH3] - P[BH3].mean()) / P[BH3].std()
    rows = []
    groups = {"5FU alone": P.arm == "5FU", "FUFOL (5-FU + leucovorin)": P.arm == "FUFOL",
              "5FU or FUFOL (no oxaliplatin)": P.arm.isin(["5FU", "FUFOL"]), "FOLFOX": P.arm == "FOLFOX",
              "no chemotherapy": P.arm == "no chemotherapy"}
    for ep, tcol, ecol in (("RFS", "rfs_months", "rfs_event"), ("OS", "os_months", "os_event")):
        for lab, sel in groups.items():
            d = P[sel & P[tcol].notna() & P[ecol].notna()]
            rows.append(subset_stats(d[tcol].values, d[ecol].values.astype(int), d.s.values, lab,
                                     dict(cohort="GSE39582 stage II-III", setup="frozen BH3-only (rank score)", endpoint=ep)))
        # interactions: FOLFOX vs no-oxaliplatin 5-FU; 5-FU-based vs untreated
        for lab, a, b in (("FOLFOX vs 5FU/FUFOL (interaction)", ["FOLFOX"], ["5FU", "FUFOL"]),
                          ("5FU/FUFOL vs no chemotherapy (interaction)", ["5FU", "FUFOL"], ["no chemotherapy"])):
            d = P[P.arm.isin(a + b) & P[tcol].notna() & P[ecol].notna()]
            r = dict(subset=lab, n=len(d), events=int(d[ecol].sum()), cohort="GSE39582 stage II-III",
                     setup="frozen BH3-only (rank score)", endpoint=ep)
            r.update(interaction(d[tcol].values, d[ecol].values.astype(int), d.s.values, d.arm.isin(a).astype(float).values))
            rows.append(r)
    return rows


def main():
    I = G2.inputs()
    T0 = pd.read_csv(TAB / "G2_decomposition_colon_5FU.tsv", sep="\t")
    sign = float(T0.frozen_sign.iloc[0])
    rows = tcga(I, sign) + gse39582(sign)
    R = pd.DataFrame(rows)
    R.insert(0, "frozen_sign", sign)
    R.to_csv(TAB / "G5_folfox_patients.tsv", sep="\t", index=False)
    pd.DataFrame(version_rows([TT / "patients.tsv.gz", TT / "pathway_scores_geo.tsv.gz"])).to_csv(
        TAB / "G5_data_versions.tsv", sep="\t", index=False)
    figure(R)


def figure(R):
    O = pd.read_csv(TAB / "G3_bh3_vs_drugs_organoids.tsv", sep="\t")
    CL = pd.read_csv(TAB / "G5_cell_lines_bh3_vs_oxaliplatin.tsv", sep="\t")
    plt = style()
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9))
    col = {"OXALIPLATIN": OI["red"], "FLUOROURACIL": OI["orange"]}
    # a organoids + b cell lines in one forest each
    ax = axs[0, 0]
    o = O[O.drug.isin(["OXALIPLATIN", "FLUOROURACIL"])].reset_index(drop=True)
    for i, r in o.iterrows():
        ax.errorbar(r.rho, i, xerr=[[r.rho - r.lo], [r.hi - r.rho]], fmt="o", color=col[r.drug], capsize=2,
                    mfc="white" if r.score == "authors" else col[r.drug])
    ax.set_yticks(range(len(o)))
    ax.set_yticklabels([f"{'oxaliplatin' if r.drug == 'OXALIPLATIN' else '5-FU'}, {r.score} score, {r.metric} (p = {fmt_p(r.p)})"
                        for r in o.itertuples()], fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("Spearman rho, BH3-only score vs response (negative = higher score, more sensitive)")
    ax.set_title("19 colorectal organoids")
    panel(ax, "a", -0.55)
    ax = axs[0, 1]
    c = CL.dropna(subset=["rho"]).reset_index(drop=True)
    for i, r in c.iterrows():
        ax.errorbar(r.rho, i, xerr=[[r.rho - r.lo], [r.hi - r.rho]], fmt="o", color=col[r.drug], capsize=2,
                    mfc="white" if r.scope.startswith("pan") else col[r.drug])
    ax.set_yticks(range(len(c)))
    ax.set_yticklabels([f"{r.screen} {'oxaliplatin' if r.drug == 'OXALIPLATIN' else '5-FU'}, {r.scope} (n = {int(r.n)}, p = {fmt_p(r.p)})"
                        for r in c.itertuples()], fontsize=6.5)
    ax.invert_yaxis()
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("Spearman rho (filled = colorectal lines, open = all lines, cancer type removed)")
    ax.set_title("Cancer cell lines")
    panel(ax, "b", -0.75)
    for ax, coh, L, title in ((axs[1, 0], "TCGA colon", "c", "TCGA colon, 5-FU-treated patients (overall survival)"),
                              (axs[1, 1], "GSE39582 stage II-III", "d", "GSE39582 stage II-III (frozen BH3-only score)")):
        d = R[(R.cohort == coh) & R.get("HRsd_HR", pd.Series(np.nan, index=R.index)).notna()].reset_index(drop=True)
        for i, r in d.iterrows():
            ax.errorbar(r.HRsd_HR, i, xerr=[[r.HRsd_HR - r.HRsd_lo], [r.HRsd_hi - r.HRsd_HR]], fmt="s", capsize=2,
                        color=OI["red"] if "FOLFOX" in r.subset or "oxaliplatin (" in r.subset else
                        (OI["grey"] if "no chemo" in r.subset else OI["blue"]), ms=4)
        ax.set_xscale("log")
        ax.axvline(1, color="grey", lw=0.7)
        ax.set_yticks(range(len(d)))
        lab = [f"{r.setup if coh.startswith('TCGA') else r.endpoint}: {r.subset} (n = {r.n}, events = {r.events}; "
               f"HR p = {fmt_p(r.HRsd_p)})" for r in d.itertuples()]
        ax.set_yticklabels(lab, fontsize=6.3)
        ax.invert_yaxis()
        ax.set_xlabel("HR per SD of predicted 5-FU resistance (frozen BH3-only), 95% CI")
        it = R[(R.cohort == coh) & R.get("inter_p", pd.Series(np.nan, index=R.index)).notna()]
        ax.set_title(title + "\n" + "; ".join(f"{(r.setup[:7] if coh.startswith('TCGA') else r.endpoint)} interaction "
                                              f"{'oxaliplatin' if coh.startswith('TCGA') else r.subset.split(' (')[0]} p = {fmt_p(r.inter_p)}"
                                              for r in it.itertuples()), fontsize=7.5)
        panel(ax, L, -1.05)
    fig.tight_layout()
    lines = "\n".join(f"- {r.cohort}, {r.setup}, {r.endpoint}, {r.subset}: n = {r.n}, events = {r.events}"
                      + (f", HR per SD {r.HRsd_HR:.2f} ({r.HRsd_lo:.2f}-{r.HRsd_hi:.2f}), p = {fmt_p(r.HRsd_p)}"
                         if np.isfinite(r.get('HRsd_HR', np.nan)) else "")
                      + (f", interaction HR {r.inter_HR:.2f} ({r.inter_lo:.2f}-{r.inter_hi:.2f}), p = {fmt_p(r.inter_p)}"
                         if np.isfinite(r.get('inter_HR', np.nan)) else "") for _, r in R.iterrows())
    savefig(fig, "G5_folfox_question", f"""
# Is the paper's 5-FU marker really a 5-FU marker, or an oxaliplatin / FOLFOX marker?

Most TCGA "5-FU" patients also had oxaliplatin (FOLFOX). The paper's marker is the "Activation of BH3-only proteins"
pathway. Here we look at it against oxaliplatin in lab models and within regimen groups in patients.

- **a** 19 colorectal organoids: correlation of the BH3-only score with oxaliplatin (red) and 5-FU (orange) response.
  Filled = our recomputed score, open = authors' score. Negative = higher score goes with more sensitivity.
- **b** Cancer cell lines (GDSC1, GDSC2, CTRPv2, PRISM): colorectal lines only (filled) and all lines with cancer type
  removed (open).
- **c** TCGA colon 5-FU-treated patients: hazard ratio per standard deviation of the frozen BH3-only prediction
  (higher = predicted more resistant; HR > 1 = the paper's direction), in all patients and separately in patients
  with and without oxaliplatin. The title gives the score x oxaliplatin interaction p (does the marker work
  differently with oxaliplatin?).
- **d** GSE39582 (French colon cancer cohort, stage II-III): the same frozen score in each adjuvant regimen and in
  untreated patients, for relapse-free (RFS) and overall survival (OS), with interaction tests.

{lines}
""")


if __name__ == "__main__":
    main()
