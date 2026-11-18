"""G2: colon 5-FU, 2019 (authors' inputs, p = 0.014) -> current data (p = 0.866): change one input at a time.

Factors (authors' value -> current value):
  organoid scores : authors' organoid NES -> ssGSEA recomputed on curated GSE64392 expression (own gene universe)
  patient scores  : authors' TCGA NES (2019 HTSeq FPKM-UQ) -> ssGSEA on current GDC STAR FPKM-UQ (protein-coding)
  gene universe   : each side scored on its own genes -> both sides on the shared organoid x TCGA COAD x READ genes
                    (this is the nar/02_organoid_only pipeline, cached scores reused)
  drug parsing    : authors' DrugBank-synonym split on 'AND' -> robust/obd two-pass normaliser (same 2019 file; the
                    current GDC biotab drug file is byte-identical, checked by cmp)
  survival        : authors' 2019 clinical_patient file -> TCGA-CDR OS (Liu 2018); also CDR capped at the 2019 follow-up
  pathway version : MSigDB v6.1 Reactome (2017) -> current Reactome release (name-matched pathways)
For every setting: (i) re-selected rank-1 Ridge pathway (the paper's method), (ii) the paper's BH3-only pathway
frozen (sign from the authors' Ridge fit, median split in treated patients), with log-rank p, HR (95% CI) and a
continuous Cox HR per SD.
"""
import pickle
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd

from gcommon import BH3, DATA, DEPMAP, FIG, O2, OI, REPO, TAB, VDW, fmt_p, panel, savefig, ssgsea_fast, style, version_rows

sys.path.insert(0, str(O2 / "code"))
import common as O  # noqa: E402  (nar/02_organoid_only/code/common.py)
from engine import Combo, surv_stats, split  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

DRUG = "FLUOROURACIL"
CACHE = DATA / "g2_scores.pkl"


def norm_name(s):
    return "REACTOME_" + re.sub(r"[^A-Z0-9]+", "_", s.upper()).strip("_")


def reactome_current():
    out = {}
    for line in open(DEPMAP / "reactome_gmt" / "ReactomePathways.gmt"):
        p = line.rstrip("\n").split("\t")
        out[norm_name(p[0])] = sorted(set(O.C.canonical_symbol(g) for g in p[2:]) - {None})
    return out


def build_scores():
    if CACHE.exists():
        return pickle.load(open(CACHE, "rb"))
    v61 = O.gene_sets(("REACTOME",))
    cur = reactome_current()
    keep = set(v61) | {BH3}
    curm = {k: v for k, v in cur.items() if k in keep}           # name-matched to v6.1 (proximity z exists)
    oe = O.set_expression("vandewetering")
    pc = set(O.C.canonical_symbol(g) for g in pd.read_csv(REPO / "data/2017_07_31_biomart_protein_coding_genes.txt",
                                                          sep="\t")["Gene name"].dropna()) - {None}
    te = O.tcga_expr("COAD")
    te = te.loc[[g for g in te.index if g in pc]]
    S = {}
    for tag, gs in (("v61", v61), ("cur", curm)):
        for side, ex in (("org", oe), ("pat", te)):
            g = {k: [x for x in v if x in set(ex.index)] for k, v in gs.items()}
            g = {k: v for k, v in g.items() if len(v) >= 5}
            S[(side, tag)] = ssgsea_fast(ex, g)
            print(side, tag, S[(side, tag)].shape, flush=True)
    S["n_genes"] = {"org": len(oe), "pat": len(te), "name_matched_current": len(curm)}
    pickle.dump(S, open(CACHE, "wb"))
    return S


def inputs():
    S = build_scores()
    org_a, resp_a, pat_a, sur_a, treat_a, any_a = O.authors_inputs(DRUG)
    sh = O.scores("vandewetering", ["COAD", "READ"])
    tt = O.treatment_table("COAD")
    cd = O.cdr()
    ok = cd["OS.time"].notna() & cd["OS"].notna()
    sur_c = pd.DataFrame({"months": cd.loc[ok, "OS.time"] / (365 / 12), "event": cd.loc[ok, "OS"].astype(int)})
    # CDR capped at each patient's 2019 follow-up (authors' file): a death after that date counts as censored
    j = sur_c.join(sur_a, rsuffix="_a", how="inner")
    cap = pd.DataFrame({"months": np.minimum(j.months, j.months_a),
                        "event": ((j.event == 1) & (j.months <= j.months_a + 1e-9)).astype(int)})
    resp_c = O.set_response("vandewetering")[DRUG]
    z = O.proximity(DRUG)[0]
    return dict(
        org={"authors": org_a, "own": S[("org", "v61")], "shared": sh["models"], "own_cur": S[("org", "cur")]},
        pat={"authors": pat_a, "star": S[("pat", "v61")], "shared": sh["COAD"], "star_cur": S[("pat", "cur")]},
        resp={"authors": resp_a, "current": resp_c},
        treated={"authors": treat_a, "current": set(tt.loc[tt.drug == DRUG, "patient"])},
        anyrec={"authors": any_a, "current": tt.attrs["any_record"]},
        sur={"authors": sur_a, "cdr": sur_c, "cdr_capped": cap},
        tt=tt, z=z, n_genes=S["n_genes"])


AUTH = dict(org="authors", pat="authors", parse="authors", sur="authors", resp="authors", patients="all")
CURR = dict(org="shared", pat="shared", parse="current", sur="cdr", resp="current", patients="all")


def run(I, f, frozen_sign, subset=None):
    osc, tsc = I["org"][f["org"]], I["pat"][f["pat"]]
    sur = I["sur"][f["sur"]]
    treated = set(I["treated"][f["parse"]])
    if f["patients"] == "overlap":                      # patients present in both authors' NES and current STAR
        both = set(I["pat"]["authors"].columns) & set(I["pat"]["star"].columns)
        tsc = tsc[[c for c in tsc.columns if c in both]]
    if subset is not None:
        treated &= subset
    untreated = set(tsc.columns) - set(I["anyrec"][f["parse"]])
    c = Combo(osc, I["resp"][f["resp"]], tsc, sur, treated, untreated, I["z"])
    ev, cf = c.evaluate(untreated=False)
    r1 = ev[(ev.model == "Ridge") & (ev.k == 1)].iloc[0]
    out = dict(n_proximal=len(c.feats), bh3_proximal=BH3 in c.feats, rank1=r1.pathways.replace("REACTOME_", ""),
               n=int(r1.treated_n), deaths=int(r1.treated_events), p_rank1=r1.treated_p, HR_rank1=r1.treated_HR,
               HR_rank1_lo=r1.treated_HR_lo, HR_rank1_hi=r1.treated_HR_hi,
               deaths_resp_rank1=int(r1.treated_events_resp), deaths_nonresp_rank1=int(r1.treated_events_nonresp))
    if BH3 in c.feats:
        out["bh3_rank_in_ridge"] = int(np.where(np.argsort(-np.abs(cf["Ridge"]), kind="stable") ==
                                                c.feats.index(BH3))[0][0] + 1)
        out["bh3_sign_refit"] = float(np.sign(cf["Ridge"][c.feats.index(BH3)]))
    if BH3 in c.P.columns:
        s = frozen_sign * c.P.loc[c.tp, BH3].values
        g = split(s)
        st = surv_stats(c.T, c.E, g)
        out.update(p_frozen=st["p"], HR_frozen=st["HR"], HR_frozen_lo=st["HR_lo"], HR_frozen_hi=st["HR_hi"],
                   deaths_resp_frozen=st["events_resp"], deaths_nonresp_frozen=st["events_nonresp"])
        try:
            zs = (s - s.mean()) / s.std()
            m = PHReg(c.T, zs.reshape(-1, 1), status=c.E, ties="efron").fit(disp=False)
            b, se = float(m.params[0]), float(m.bse[0])
            out.update(HRsd_frozen=np.exp(b), HRsd_lo=np.exp(b - 1.96 * se), HRsd_hi=np.exp(b + 1.96 * se),
                       p_cox_frozen=float(m.pvalues[0]))
        except Exception:
            pass
    return out


def main():
    I = inputs()
    base = run(I, AUTH, 1.0)
    # sign of the BH3-only coefficient in the authors' Ridge fit; frozen biomarker = sign * score
    sign = base.get("bh3_sign_refit", 1.0)
    base = run(I, AUTH, sign)
    rows = []

    def add(group, step, label, f, subset=None):
        r = run(I, f, sign, subset)
        r.update(group=group, step=step, label=label, **{f"f_{k}": v for k, v in f.items()})
        rows.append(r)
        print(group, step, label, r["n"], r["deaths"], fmt_p(r["p_rank1"]), fmt_p(r.get("p_frozen", np.nan)), flush=True)

    # A. cumulative chain authors -> current
    chain = [("0", "authors' setup (2019)", dict(AUTH)),
             ("1", "+ only patients with current STAR data", dict(AUTH, patients="overlap")),
             ("2", "+ organoid scores recomputed", dict(AUTH, patients="overlap", org="own")),
             ("3", "+ patient scores from current STAR FPKM-UQ", dict(AUTH, patients="overlap", org="own", pat="star")),
             ("4", "+ shared gene universe (02 pipeline)", dict(AUTH, patients="all", org="shared", pat="shared")),
             ("5", "+ organoid response from curated table", dict(AUTH, org="shared", pat="shared", resp="current")),
             ("6", "+ current drug-record parser", dict(AUTH, org="shared", pat="shared", resp="current", parse="current")),
             ("7", "+ TCGA-CDR survival = current setup", dict(CURR))]
    for s, l, f in chain:
        add("A_chain", s, l, f)
    # B. one factor at a time from the authors' setup
    one = [("org", "organoid scores recomputed", dict(AUTH, org="own")),
           ("pat", "patient scores from current STAR", dict(AUTH, pat="star")),
           ("univ", "shared gene universe, both sides", dict(AUTH, org="shared", pat="shared")),
           ("resp", "organoid response from curated table", dict(AUTH, resp="current")),
           ("parse", "current drug-record parser", dict(AUTH, parse="current")),
           ("sur", "TCGA-CDR survival", dict(AUTH, sur="cdr")),
           ("surcap", "TCGA-CDR survival capped at 2019 follow-up", dict(AUTH, sur="cdr_capped")),
           ("overlap", "only patients with current STAR data", dict(AUTH, patients="overlap"))]
    for s, l, f in one:
        add("B_one_from_authors", s, l, f)
    # C. one factor at a time put back into the current setup
    back = [("org", "authors' organoid NES", dict(CURR, org="authors")),
            ("pat", "authors' patient NES", dict(CURR, pat="authors")),
            ("orgpat", "authors' organoid + patient NES", dict(CURR, org="authors", pat="authors")),
            ("resp", "authors' organoid response file", dict(CURR, resp="authors")),
            ("parse", "authors' drug-record parser", dict(CURR, parse="authors")),
            ("sur", "authors' 2019 survival", dict(CURR, sur="authors")),
            ("surcap", "TCGA-CDR capped at 2019 follow-up", dict(CURR, sur="cdr_capped")),
            ("version", "current Reactome release (both sides)", dict(CURR, org="own_cur", pat="star_cur"))]
    for s, l, f in back:
        add("C_one_into_current", s, l, f)
    # D. frozen BH3-only on current data, by regimen (FOLFOX vs no oxaliplatin) -- also used by G5
    tt = I["tt"]
    oxa = set(tt.loc[tt.drug == "OXALIPLATIN", "patient"])
    iri = set(tt.loc[tt.drug == "IRINOTECAN", "patient"])
    allp = set(I["pat"]["shared"].columns) | set(I["pat"]["authors"].columns)
    for setup, F in (("authors", AUTH), ("current", CURR), ("current_capped", dict(CURR, sur="cdr_capped"))):
        for s, l, sub in (("all", "all 5-FU patients", None), ("folfox", "5-FU + oxaliplatin", oxa),
                          ("nooxa", "5-FU without oxaliplatin", allp - oxa),
                          ("alone", "5-FU, no oxaliplatin or irinotecan", allp - oxa - iri)):
            if setup == "authors" and sub is not None:     # authors' parser for the regimen split on authors' files
                at = I["treated"]
                syn_oxa = set(O.authors_inputs("OXALIPLATIN")[4])
                syn_iri = set(O.authors_inputs("IRINOTECAN")[4])
                sub = {"folfox": syn_oxa, "nooxa": allp - syn_oxa, "alone": allp - syn_oxa - syn_iri}[s]
            add(f"D_regimen_{setup}", s, l, F, sub)
    T = pd.DataFrame(rows)
    T.insert(0, "frozen_sign", sign)
    T.to_csv(TAB / "G2_decomposition_colon_5FU.tsv", sep="\t", index=False)
    pd.Series(I["n_genes"]).to_csv(TAB / "G2_gene_universe_sizes.tsv", sep="\t")
    pd.DataFrame(version_rows([
        REPO / "python/results/COAD/TCGA/reactome_ssgsea_result/gseapy.samples.normalized.es.txt",
        REPO / "python/results/COAD/organoid/reactome_ssgsea_result/gseapy.samples.normalized.es.txt",
        REPO / "data/TCGA_COAD/clinical_drug/dae92461-d098-40a3-997e-1a2d3e5b8dfe/nationwidechildrens.org_clinical_drug_coad.txt",
        O2 / "_data/raw/nationwidechildrens.org_clinical_drug_coad.txt",
        REPO / "data/TCGA_COAD/clinical_patient/4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt",
        O2 / "_data/raw/Survival_SupplementalTable_S1_20171025_xena_sp.tsv", O2 / "_data/raw/TCGA-COAD.star_fpkm-uq.tsv.gz",
        DEPMAP / "reactome_gmt/ReactomePathways.gmt", REPO / "data/msigdb.v6.1.symbols.gmt.txt"])).to_csv(
        TAB / "G2_data_versions.tsv", sep="\t", index=False)
    figure(T)


def figure(T):
    plt = style()
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5), gridspec_kw=dict(height_ratios=[1, 1]))
    lp = lambda p: -np.log10(p) if np.isfinite(p) and p > 0 else np.nan  # noqa: E731
    for col, grp, title, L in ((0, "A_chain", "Cumulative: authors' setup -> current setup", "a"),
                               (1, "C_one_into_current", "One authors' input put back into the current setup", "c")):
        d = T[T.group == grp].reset_index(drop=True)
        ax = axs[0, col]
        x = np.arange(len(d))
        ax.plot(x, [lp(p) for p in d.p_rank1], "-o", color=OI["orange"], label="re-selected rank-1 pathway")
        ax.plot(x, [lp(p) for p in d.p_frozen], "-s", color=OI["blue"], label="BH3-only pathway frozen")
        ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.8)
        ax.text(len(d) - 0.5, -np.log10(0.05) + 0.03, "p = 0.05", fontsize=7, ha="right", color="grey")
        for i, r in d.iterrows():
            ax.annotate(f"{r.n}/{r.deaths}", (i, lp(r.p_frozen)), xytext=(0, -11), textcoords="offset points",
                        fontsize=6.5, ha="center", color=OI["blue"])
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r.step}. {r.label}" if grp == "A_chain" else r.label for r in d.itertuples()],
                           rotation=35, ha="right", fontsize=7)
        ax.set_ylabel("-log10 log-rank p (higher = stronger split)")
        ax.set_title(title + "\n(numbers = treated patients / deaths)")
        ax.legend(frameon=False, loc="upper right")
        panel(ax, L)
        ax = axs[1, col]
        for i, r in d.iterrows():
            if np.isfinite(r.get("HRsd_frozen", np.nan)):
                ax.errorbar(r.HRsd_frozen, i, xerr=[[r.HRsd_frozen - r.HRsd_lo], [r.HRsd_hi - r.HRsd_frozen]],
                            fmt="s", color=OI["blue"], capsize=2, ms=4)
                ax.text(r.HRsd_hi * 1.05, i, f"{r.HRsd_frozen:.2f} ({r.HRsd_lo:.2f}-{r.HRsd_hi:.2f})", va="center",
                        fontsize=6.5)
        ax.set_xscale("log")
        ax.axvline(1, color="grey", lw=0.7)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([r.label for r in d.itertuples()], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Cox HR per SD of frozen BH3-only prediction (>1 = predicted non-responders die sooner)")
        ax.set_title("Frozen BH3-only pathway as a continuous score")
        panel(ax, "b" if col == 0 else "d", -0.45)
    fig.tight_layout()
    a = T[T.group == "A_chain"].reset_index(drop=True)
    lines = "\n".join(f"- Step {r.step} ({r.label}): n = {r.n}, deaths = {r.deaths}; re-selected rank-1 = "
                      f"{r.rank1.lower().replace('_', ' ')}, p = {fmt_p(r.p_rank1)}; frozen BH3-only p = "
                      f"{fmt_p(r.p_frozen)}, HR per SD {r.HRsd_frozen:.2f} ({r.HRsd_lo:.2f}-{r.HRsd_hi:.2f})"
                      for r in a.itertuples())
    savefig(fig, "G2_decomposition_colon_5FU", f"""
# What changed between the 2019 colon 5-FU result and today's data

The paper split 5-FU-treated TCGA colon patients by an organoid-trained model and got log-rank p = 0.014. On today's
data the same pipeline gives p = 0.866. Here we change one input at a time to see which change moves the result.
Orange = the paper's method run fresh at each step (the model re-selects its top pathway). Blue = the paper's own
pathway, "Activation of BH3-only proteins", kept fixed (frozen) with the direction learned in 2019.

- **a** Cumulative path from the authors' 2019 inputs (left) to the current inputs (right). Each step adds one
  current input on top of the previous ones. Height = -log10 p (the dashed line is p = 0.05).
- **b** Same steps; the frozen BH3-only score used as a continuous predictor (Cox hazard ratio per standard deviation,
  95% CI). HR > 1 means patients predicted to respond poorly die sooner, the paper's direction.
- **c, d** The reverse view: start from the current setup and put back one of the authors' inputs at a time.

Steps of panel a:
{lines}
""")


if __name__ == "__main__":
    main()
