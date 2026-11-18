"""G12: baseline characteristics of treated vs 'no drug record' (untreated) TCGA patients in the A3 survey.

Clinical covariates exactly as in nar/05_A3_B1/code/a3_analyze.py (clinical(): TCGA-CDR age/sex/stage, ABSOLUTE
purity, MSI from Xena COADREAD + 2019 COAD biotab). One block per TCGA cohort x drug (treated set of that drug) plus a
cohort-level row (treated with any survey drug).
"""
import importlib.util

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from g_common import A3, OK, TAB, C, load_patients, save_fig

spec = importlib.util.spec_from_file_location("a3an", A3 / "code" / "a3_analyze.py")
AN = importlib.util.module_from_spec(spec)
spec.loader.exec_module(AN)


def desc(cl, pats, T=None):
    c = cl.reindex(pats)
    st = c.stage
    r = {"n": len(pats), "age_median": c.age.median(), "age_q1": c.age.quantile(.25), "age_q3": c.age.quantile(.75),
         "pct_male": 100 * c.male.mean(), "stage_known_n": int(st.notna().sum()),
         "pct_stage_III_IV": 100 * (st >= 3).sum() / max(st.notna().sum(), 1),
         "pct_stage_IV": 100 * (st >= 4).sum() / max(st.notna().sum(), 1),
         "purity_known_n": int(c.purity.notna().sum()), "purity_median": c.purity.median(),
         "msi_known_n": int(c.msi.notna().sum()), "pct_MSI_H": 100 * c.msi.mean() if c.msi.notna().any() else np.nan}
    return r, c


def compare(ct, cu):
    out = {}
    for v in ("age", "purity"):
        a, b = ct[v].dropna(), cu[v].dropna()
        out[f"p_{v}_MannWhitney"] = mannwhitneyu(a, b).pvalue if len(a) > 2 and len(b) > 2 else np.nan
    for v, f in (("male", lambda s: s == 1), ("stage", lambda s: s >= 3), ("msi", lambda s: s == 1)):
        a, b = ct[v].dropna(), cu[v].dropna()
        if len(a) > 2 and len(b) > 2 and (a.nunique() > 1 or b.nunique() > 1):
            out[f"p_{'stageIII_IV' if v == 'stage' else v}_Fisher"] = fisher_exact([[f(a).sum(), (~f(a)).sum()], [f(b).sum(), (~f(b)).sum()]])[1]
    return out


def main():
    cl = AN.clinical()
    p = load_patients()
    rows = []
    blocks = [(proj, "any survey drug", p[p.proj == proj]) for proj in sorted(p.proj.unique())]
    for (proj, drug), d in p.assign(drug=p.analysis_id.str.split("__").str[2]).groupby(["proj", "drug"]):
        blocks.append((proj, drug.lower(), d))
    for proj, drug, d in blocks:
        tp = sorted(d[d.group == "treated"].patient.unique())
        up = sorted(d[d.group == "untreated"].patient.unique())
        rt, ct = desc(cl, tp)
        ru, cu = desc(cl, up)
        row = {"cohort": proj, "treated_with": drug}
        row.update({f"treated_{k}": v for k, v in rt.items()})
        row.update({f"untreated_{k}": v for k, v in ru.items()})
        row.update(compare(ct, cu))
        rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(TAB / "G12_baseline_treated_vs_untreated.tsv", sep="\t", index=False)
    # readable version
    def f(r, g):
        s = f"{r[g+'_age_median']:.0f} ({r[g+'_age_q1']:.0f}-{r[g+'_age_q3']:.0f})"
        return s
    R = pd.DataFrame({"cohort": T.cohort, "treated with": T.treated_with,
                      "n treated / untreated": T.treated_n.astype(int).astype(str) + " / " + T.untreated_n.astype(int).astype(str),
                      "age treated": [f(r, "treated") for _, r in T.iterrows()],
                      "age untreated": [f(r, "untreated") for _, r in T.iterrows()],
                      "p age": T.p_age_MannWhitney.map(C.fp),
                      "% stage III-IV treated / untreated": T.treated_pct_stage_III_IV.round(0).astype("Int64").astype(str) + " / " + T.untreated_pct_stage_III_IV.round(0).astype("Int64").astype(str),
                      "p stage": T.get("p_stageIII_IV_Fisher").map(C.fp),
                      "purity treated / untreated": T.treated_purity_median.round(2).astype(str) + " / " + T.untreated_purity_median.round(2).astype(str),
                      "p purity": T.p_purity_MannWhitney.map(C.fp),
                      "% MSI-H treated / untreated": [f"{a:.0f} / {b:.0f}" if np.isfinite(a) and r in ("COAD", "READ", "COADREAD", "COAD_authors") else "" for a, b, r in zip(T.treated_pct_MSI_H, T.untreated_pct_MSI_H, T.cohort)],
                      "p MSI": [C.fp(x) if r in ("COAD", "READ", "COADREAD", "COAD_authors") else "" for x, r in zip(T.get("p_msi_Fisher"), T.cohort)]})
    R.to_csv(TAB / "G12_baseline_readable.tsv", sep="\t", index=False)
    coh = T[T.treated_with == "any survey drug"]
    print(R[R["treated with"] == "any survey drug"].to_string())
    print(T[T.treated_with == "any survey drug"][["cohort", "treated_pct_male", "untreated_pct_male", "p_male_Fisher"]].to_string())
    plt = C.plot_style()
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    x = np.arange(len(coh))
    w = .38
    for a, (vt, vu, lab) in zip(ax, (("treated_age_median", "untreated_age_median", "median age (years)"),
                                     ("treated_pct_stage_III_IV", "untreated_pct_stage_III_IV", "% stage III-IV (of known)"),
                                     ("treated_purity_median", "untreated_purity_median", "median tumour purity (ABSOLUTE)"))):
        a.bar(x - w / 2, coh[vt], w, color=OK["blue"], label="treated with a survey drug")
        a.bar(x + w / 2, coh[vu], w, color=OK["orange"], label="no drug record")
        a.set_xticks(x); a.set_xticklabels(coh.cohort, rotation=60, fontsize=7)
        a.set_ylabel(lab)
    ax[0].legend(frameon=False, fontsize=7, loc="lower left")
    for a, l in zip(ax, "abc"):
        C.letter(a, l)
    fig.tight_layout()
    cap = """**Figure G12. Treated patients and patients with no drug record are different groups.**
For each TCGA cohort used in the 57-analysis survey, blue = patients treated with at least one of the survey drugs,
orange = patients with no drug record (the "untreated" comparison group). (a) median age at diagnosis; (b) share of
patients with stage III or IV disease, among those with a known stage; (c) median tumour purity (share of cancer cells
in the sample, ABSOLUTE method). The untreated group is not a randomised control: in most cohorts it is older or
younger, and has more or less advanced disease, than the treated group. Exact numbers, MSI (colon/rectum) and p-values
are in tables/G12_baseline_treated_vs_untreated.tsv. COAD_authors = the 2019 colon files used by the original paper.
"""
    save_fig(fig, "G12_baseline_treated_vs_untreated", cap)


if __name__ == "__main__":
    main()
