"""G8 part 2: BRAF V600E colon tumours and the cetuximab organoid model -- random-pathway nulls and adjustment for
MSI, tumour side, stage and tumour mutation burden (TMB). Also draws the G8 figure (reads the ERCC2 results).

Score as nar/02_organoid_only/code/extras.fig5a: Ridge on the 19 van de Wetering organoids (cetuximab ln IC50),
cetuximab-proximal pathways, rank-1 pathway x coefficient = predicted resistance; V600E expected HIGHER.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd

from gcommon import FIG, NAR, O2, OI, TAB, fmt_p, panel, savefig, style, version_rows
from g8_nulls import adjust, all_pathway_null, mwu, pipeline_null, rz, tmb

sys.path.insert(0, str(O2 / "code"))
import common as O  # noqa: E402
import run_all as A  # noqa: E402
from engine import Combo  # noqa: E402

RIGHT = {"Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon"}
LEFT = {"Splenic Flexure", "Descending Colon", "Sigmoid Colon", "Rectosigmoid Junction", "Rectum"}


def covariates(T):
    c = pd.read_csv(NAR / "05_A3_B1/_data/raw/COADREAD_clinicalMatrix", sep="\t", index_col=0)
    c = c[c.index.str[13:15] == "01"]
    c.index = c.index.str[:12]
    c = c[~c.index.duplicated()]
    msi = c["CDE_ID_3226963"].map({"MSI-H": 1.0, "MSI-L": 0.0, "MSS": 0.0})
    side = c["anatomic_neoplasm_subdivision"].map(lambda s: 1.0 if s in RIGHT else (0.0 if s in LEFT else np.nan))
    st = c["pathologic_stage"].astype(str).str.replace("Stage ", "")
    stage = pd.Series(np.where(st.str.startswith(("III", "IV")), 1.0, np.where(st.str.startswith("I"), 0.0, np.nan)), index=c.index)
    cov = pd.DataFrame({"MSI_H": msi, "right_side": side, "stage_III_IV": stage})
    cov = cov.join(pd.DataFrame({"logTMB_rz": rz(np.log1p(T.values))}, index=T.index), how="inner")
    return cov


def main():
    cd = O.cdr()
    import json
    braf = json.loads((O2 / "_data/braf_v600e_COAD.json").read_text())
    mut, prof_gdc = set(braf["v600e"]), set(braf["profiled"])
    T, _ = tmb(NAR / "10_completeness/_data/raw/xena_mut/TCGA-COAD.somaticmutation_wxs.tsv.gz")
    cov = covariates(T)
    z, _ = O.proximity("CETUXIMAB")
    rows, nulls = [], {}
    for version in ("authors", "current"):
        if version == "authors":
            osc, y, tsc = O.authors_inputs("CETUXIMAB")[:3]
        else:
            osc, tsc, _, _, _ = A.cohort_inputs("vandewetering", "COAD", cd)
            y = O.set_response("vandewetering")["CETUXIMAB"]
        surall = pd.DataFrame({"months": 1.0, "event": 0}, index=list(tsc.columns))
        c = Combo(osc, y, tsc, surall, set(), set(), z)
        cf = c.coefs()["Ridge"]
        i = int(np.argmax(np.abs(cf)))
        pw = c.feats[i]
        Pz = c.P
        s = cf[i] * Pz[pw]
        prof = prof_gdc
        s = s[s.index.isin(prof)]
        p = mwu(s, mut, "greater")
        sh_all, n_all, ps_all = all_pathway_null(Pz, prof, mut, "greater", p)
        sh_pipe, ps_pipe = pipeline_null(c.X, c.y, Pz, list(c.common), len(c.feats), prof, mut, "greater", p, n=2000, seed=2)
        nulls[version] = (p, ps_all, ps_pipe)
        r = dict(version=version, model="rank-1 pathway", pathway=pw.replace("REACTOME_", ""), coef=cf[i],
                 n_proximal=len(c.feats), n_mut=int(s.index.isin(mut).sum()), n_wt=int((~s.index.isin(mut)).sum()),
                 p_one_sided=p, share_all_pathways_both_signs=sh_all, n_null_tests=n_all,
                 share_random_proximal_sized_sets=sh_pipe, n_random_sets=2000)
        for lab, cols in (("MSI", ["MSI_H"]), ("side", ["right_side"]), ("TMB", ["logTMB_rz"]),
                          ("all", ["MSI_H", "right_side", "stage_III_IV", "logTMB_rz"])):
            r.update({f"adj{lab}_{k}": v for k, v in adjust(s, mut, cov[cols], 1).items()})
        # within MSS tumours only
        mss = [x for x in s.index if cov["MSI_H"].get(x) == 0.0]
        r["MSS_only_n_mut"] = int(pd.Index(mss).isin(mut).sum())
        r["MSS_only_p_one_sided"] = mwu(s.loc[mss], mut, "greater")
        rows.append(r)
        print(version, pw, p, sh_all, sh_pipe, r["adjall_ols_p_one_sided"], flush=True)
    R = pd.DataFrame(rows)
    R.insert(0, "test", "BRAF V600E (colon, cetuximab)")
    R.to_csv(TAB / "G8_BRAF_null_and_adjustment.tsv", sep="\t", index=False)
    pd.to_pickle(nulls, HERE.parent / "_data" / "g8_braf_nulls.pkl")
    pd.DataFrame(version_rows([O2 / "_data/braf_v600e_COAD.json", NAR / "05_A3_B1/_data/raw/COADREAD_clinicalMatrix",
                               NAR / "10_completeness/_data/raw/xena_mut/TCGA-COAD.somaticmutation_wxs.tsv.gz"])).to_csv(
        TAB / "G8_BRAF_data_versions.tsv", sep="\t", index=False)
    figure()


def figure():
    E = pd.read_csv(TAB / "G8_ERCC2_null_and_adjustment.tsv", sep="\t")
    B = pd.read_csv(TAB / "G8_BRAF_null_and_adjustment.tsv", sep="\t")
    ne = pd.read_pickle(HERE.parent / "_data" / "g8_ercc2_nulls.pkl")
    nb = pd.read_pickle(HERE.parent / "_data" / "g8_braf_nulls.pkl")
    plt = style()
    fig, axs = plt.subplots(2, 3, figsize=(13, 8))
    panels = [(("paper", "paper biomarker (amino acid synthesis)"), ne, "ERCC2, paper data, amino-acid pathway"),
              (("current", "paper biomarker (amino acid synthesis)"), ne, "ERCC2, current data, amino-acid pathway"),
              ("authors", nb, "BRAF V600E, authors' data"), ("current", nb, "BRAF V600E, current data")]
    for ax, (key, nd, title), L in zip([axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]], panels, "abde"):
        p, pa, pp = nd[key]
        bins = np.linspace(0, 6, 40)
        ax.hist(np.clip(-np.log10(pa[np.isfinite(pa)]), 0, 6), bins=bins, density=True, color=OI["sky"], alpha=0.6,
                label="every pathway, both directions")
        ax.hist(np.clip(-np.log10(pp[np.isfinite(pp)]), 0, 6), bins=bins, density=True, color=OI["grey"], alpha=0.6,
                label="random proximal-sized sets -> Ridge rank 1")
        ax.axvline(-np.log10(p), color=OI["red"], lw=2, label=f"observed p = {fmt_p(p)}")
        sa = np.nanmean(pa <= p + 1e-15)
        sp = np.nanmean(pp <= p + 1e-15)
        ax.set_title(f"{title}\nshare as small: pathways {sa:.3f}; random sets {sp:.3f}", fontsize=8.5)
        ax.set_xlabel("-log10 one-sided Mann-Whitney p")
        ax.set_ylabel("density")
        ax.legend(frameon=False, fontsize=6.5)
        panel(ax, L)
    # adjusted effects
    for ax, D, labs, title, L in ((axs[0, 2], E, [("raw", None), ("TMB", "adjTMB"), ("TMB + stage", "adjTMBstage")],
                                   "ERCC2: score difference, mutant vs wild type", "c"),
                                  (axs[1, 2], B, [("raw", None), ("MSI", "adjMSI"), ("side", "adjside"), ("TMB", "adjTMB"),
                                                  ("MSI+side+stage+TMB", "adjall")],
                                   "BRAF V600E: score difference, mutant vs wild type", "f")):
        y = 0
        ticks = []
        for r in D.itertuples():
            if "rank-1" in r.model and "ERCC2" in r.test:
                continue
            for lab, pre in labs:
                if pre is None:
                    ticks.append((y, f"{r.version}: unadjusted p = {fmt_p(r.p_one_sided)}"))
                    y += 1
                    continue
                b, lo, hi, p1 = (getattr(r, f"{pre}_ols_beta_mut"), getattr(r, f"{pre}_ols_lo"), getattr(r, f"{pre}_ols_hi"),
                                 getattr(r, f"{pre}_ols_p_one_sided"))
                ax.errorbar(b, y, xerr=[[b - lo], [hi - b]], fmt="o", color=OI["blue"] if r.version != "current" else OI["orange"],
                            capsize=2, ms=4)
                ticks.append((y, f"{r.version}, + {lab}: one-sided p = {fmt_p(p1)}"))
                y += 1
        ax.set_yticks([t[0] for t in ticks])
        ax.set_yticklabels([t[1] for t in ticks], fontsize=6.5)
        ax.invert_yaxis()
        ax.axvline(0, color="grey", lw=0.7)
        ax.set_xlabel("difference in rank-normalised score (OLS, 95% CI)")
        ax.set_title(title, fontsize=8.5)
        panel(ax, L, -0.9)
    fig.tight_layout()
    e = E[E.model.str.startswith("paper")]
    lines = "\n".join(f"- ERCC2, {r.version} data: n = {r.n_mut} mutant vs {r.n_wt} wild type, p = {fmt_p(r.p_one_sided)}; "
                      f"share of pathways (both directions) as small {r.share_all_pathways_both_signs:.3f}; random proximal-sized "
                      f"sets {r.share_random_proximal_sized_sets:.3f}; TMB-adjusted p = {fmt_p(r.adjTMB_ols_p_one_sided)}"
                      for r in e.itertuples())
    lines += "\n" + "\n".join(f"- BRAF V600E, {r.version} data ({r.pathway.lower().replace('_', ' ')}): n = {r.n_mut} vs {r.n_wt}, p = "
                              f"{fmt_p(r.p_one_sided)}; pathways {r.share_all_pathways_both_signs:.3f}; random sets "
                              f"{r.share_random_proximal_sized_sets:.3f}; adjusted for MSI, side, stage and TMB p = "
                              f"{fmt_p(r.adjall_ols_p_one_sided)}; MSS tumours only p = {fmt_p(r.MSS_only_p_one_sided)}"
                              for r in B.itertuples())
    savefig(fig, "G8_mutation_concordance_null_and_adjustment", f"""
# Are the mutation checks special to the organoid pathways, and do they survive adjustment?

The paper checked its models against known mutations: bladder tumours with ERCC2 mutations should look more
cisplatin-sensitive, and colon tumours with BRAF V600E should look more cetuximab-resistant. Both reproduced. Here we
ask two questions. (1) Would any pathway do as well? (2) Is the result explained by tumour type features that travel
with the mutation (high tumour mutation burden = TMB for ERCC2; MSI-high status, right-sided colon and TMB for BRAF)?

- **a, b, d, e** Red line = observed p (higher = stronger). Blue = the same test with every Reactome pathway in both
  directions; grey = 2,000 runs of the whole pipeline with random pathway sets of the proximal size (Ridge on the
  organoids, rank-1 pathway, test). The titles give the share of each null that is at least as strong.
- **c, f** Difference in the (rank-normalised) predicted-resistance score between mutant and wild-type tumours, first
  without and then with adjustment (linear regression). For ERCC2 a negative difference is the expected direction;
  for BRAF V600E a positive difference.

{lines}
""")


if __name__ == "__main__":
    main()
