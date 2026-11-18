"""G14: paper Supp. Fig. 15 (models trained on cancer cell lines instead of organoids, applied to TCGA patients) for
pancreas gemcitabine, liver sorafenib, ovarian platinum and paclitaxel, and sarcoma gemcitabine.

Same pipeline as nar/02_organoid_only (engine.Combo): drug-proximal Reactome pathways (paper's STRING proximity,
z <= -1.2816), ssGSEA on the shared gene universe (cell lines x TCGA project), Ridge on cell-line response, the
top three pathways each used alone, median split of treated patients, log-rank + Cox HR (overall survival, TCGA-CDR).
Cell lines: GDSC2 (release 8.5, ln IC50) because GDSC2 screened sorafenib and paclitaxel in ~950 lines (GDSC1 only
~400); expression = Cell Model Passports RNA-seq (nar/_cache/depmap/derived/expr_cmp.pkl).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd

from gcommon import DATA, DEPMAP, FIG, O2, OI, REPO, TAB, fmt_p, panel, savefig, ssgsea_fast, style, version_rows

sys.path.insert(0, str(O2 / "code"))
import common as O  # noqa: E402
from engine import Combo, km_curve, surv_stats  # noqa: E402

BRAW = O2 / "BLCA_OV_LUAD_OTHER/_data/raw"
SARC_LINES = {"Rhabdomyosarcoma", "Other Sarcomas", "Ewing's Sarcoma", "Osteosarcoma"}
SPECS = [  # project, patient drug(s), cell-line drug, lineage set, label
    ("PAAD", ["GEMCITABINE"], "GEMCITABINE", {"Pancreatic Carcinoma"}, "pancreas | gemcitabine"),
    ("LIHC", ["SORAFENIB"], "SORAFENIB", {"Hepatocellular Carcinoma"}, "liver | sorafenib"),
    ("OV", ["CARBOPLATIN", "CISPLATIN"], "CISPLATIN", {"Ovarian Carcinoma"}, "ovary | platinum (carboplatin or cisplatin)"),
    ("OV", ["PACLITAXEL"], "PACLITAXEL", {"Ovarian Carcinoma"}, "ovary | paclitaxel"),
    ("SARC", ["GEMCITABINE"], "GEMCITABINE", SARC_LINES, "sarcoma | gemcitabine"),
]


def raw_dir(project):
    return O.RAW if (O.RAW / f"nationwidechildrens.org_clinical_drug_{project.lower()}.txt").exists() else BRAW


def treated_sets(project):
    O.RAW = raw_dir(project)
    O.C.EXTERNAL = O.RAW
    O._TREAT.pop(project, None)
    tt = O.treatment_table(project)
    O.RAW = O2 / "_data/raw"
    return tt


def tcga(project):
    O.C.EXTERNAL = O2 / "_data/raw" if (O2 / f"_data/raw/TCGA-{project}.star_fpkm-uq.tsv.gz").exists() else BRAW
    return O.C.tcga_star_fpkm_uq(project)


def main():
    lines_expr = pd.read_pickle(DEPMAP / "derived/expr_cmp.pkl")
    models = pd.read_pickle(DEPMAP / "derived/models.pkl")
    g = pd.read_pickle(DEPMAP / "derived/gdsc_long.pkl")
    g = g[(g.dataset == "GDSC2") & g.ln_ic50.notna() & g.sanger_model_id.notna()]
    resp = g.groupby(["sanger_model_id", "drug"]).ln_ic50.median().unstack()
    cd = O.cdr()
    ok = cd["OS.time"].notna() & cd["OS"].notna()
    sur = pd.DataFrame({"months": cd.loc[ok, "OS.time"] / (365 / 12), "event": cd.loc[ok, "OS"].astype(int)})
    gs_all = O.gene_sets(("REACTOME",))
    rows, kms = [], []
    cache = {}
    for project, pdrugs, ldrug, lineage, label in SPECS:
        if project not in cache:
            te = tcga(project)
            genes = lines_expr.index.intersection(te.index)
            gs = {k: [x for x in v if x in set(genes)] for k, v in gs_all.items()}
            gs = {k: v for k, v in gs.items() if len(v) >= 5}
            cache[project] = (ssgsea_fast(lines_expr.loc[genes], gs), ssgsea_fast(te.loc[genes], gs), len(genes))
            print(project, cache[project][0].shape, cache[project][1].shape, flush=True)
        lsc, tsc, ng = cache[project]
        tt = treated_sets(project)
        treated = set(tt.loc[tt.drug.isin(pdrugs), "patient"])
        untreated = set(tsc.columns) - tt.attrs["any_record"]
        z, zcol = O.proximity(ldrug)
        y_all = resp[ldrug].dropna()
        y_all = y_all[y_all.index.isin(lsc.columns)]
        for tr, ids in (("all GDSC2 lines", list(y_all.index)),
                        ("lineage-matched lines", [i for i in y_all.index if models.cancer_type.get(i) in lineage])):
            y = y_all.loc[ids]
            c = Combo(lsc, y, tsc, sur, treated, untreated, z)
            cf = c.coefs()["Ridge"]
            loo = c.loo()["Ridge"] if len(c.samples) < 150 else (None, np.nan, np.nan)   # LOO only for the smaller lineage sets
            for k in (1, 2, 3):
                o = np.argsort(-np.abs(cf), kind="stable")[k - 1]
                pw = c.feats[o]
                pred = cf[o] * c.P.loc[c.tp, pw].values
                grp = (pred > np.median(pred)).astype(int)
                st = surv_stats(c.T, c.E, grp)
                predu = cf[o] * c.P.loc[c.up, pw].values
                stu = surv_stats(c.Tu, c.Eu, (predu > np.median(predu)).astype(int)) if len(c.up) > 10 else {}
                rows.append(dict(setting=label, project=project, patient_drugs="+".join(pdrugs), line_drug=ldrug,
                                 proximity_column=zcol, training=tr, n_lines=len(c.samples), n_proximal=len(c.feats),
                                 genes=ng, rank=k, pathway=pw, loo_rho_lines=loo[1], loo_p_lines=loo[2],
                                 **st, **{f"untreated_{a}": b for a, b in stu.items() if a in ("n", "events", "p", "HR")}))
                if k == 1:
                    kms.append((label, tr, len(c.samples), pw, c.T, c.E, grp, st))
            print(label, tr, len(c.samples), rows[-3]["pathway"], rows[-3]["n"], rows[-3]["events"], rows[-3]["p"], flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(TAB / "G14_SuppFig15_cell_line_models.tsv", sep="\t", index=False)
    pd.DataFrame(version_rows([DEPMAP / "derived/expr_cmp.pkl", DEPMAP / "derived/gdsc_long.pkl", DEPMAP / "derived/models.pkl"]
                              + [p for p in (O2 / "_data/raw").glob("TCGA-*.star_fpkm-uq.tsv.gz") if p.name[5:9] in ("PAAD", "LIHC")]
                              + [BRAW / "TCGA-OV.star_fpkm-uq.tsv.gz", BRAW / "TCGA-SARC.star_fpkm-uq.tsv.gz"])).to_csv(
        TAB / "G14_data_versions.tsv", sep="\t", index=False)
    figure(T, kms)


def figure(T, kms):
    plt = style()
    n = len(kms)
    fig = plt.figure(figsize=(15, 3.2 * (n // 4 + 2)))
    gsp = fig.add_gridspec(n // 4 + 2, 4)
    ax = fig.add_subplot(gsp[0, :])
    d = T.reset_index(drop=True)
    x = np.arange(len(d))
    hr = d.HR.values
    ok = np.isfinite(hr)
    ax.errorbar(x[ok], hr[ok], yerr=[hr[ok] - d.HR_lo.values[ok], d.HR_hi.values[ok] - hr[ok]], fmt="o", capsize=2,
                color=OI["blue"], ms=4)
    for i, r in d.iterrows():
        ax.text(i, (r.HR_hi if np.isfinite(r.HR_hi) else 1) * 1.15, f"p={fmt_p(r.p)}", rotation=90, fontsize=6, ha="center", va="bottom")
    ax.set_yscale("log")
    ax.axhline(1, color="grey", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.setting.split(' |')[0]} {r.line_drug.lower()[:5]}, {'all' if r.training.startswith('all') else 'lineage'}, rank {r.rank}"
                        for r in d.itertuples()], rotation=60, ha="right", fontsize=6.5)
    ax.set_ylabel("HR, predicted non-responders vs responders")
    ax.set_title("Cell-line-trained models applied to TCGA treated patients (top 3 pathways, each alone)")
    panel(ax, "a", -0.04)
    L = iter("bcdefghijklmnop")
    for i, (label, tr, nl, pw, Tt, Et, grp, st) in enumerate(kms):
        ax = fig.add_subplot(gsp[1 + i // 4, i % 4])
        for gv, c, nm in ((0, OI["blue"], "predicted responders"), (1, OI["orange"], "predicted non-responders")):
            t, s, lo, hi = km_curve(Tt[grp == gv], Et[grp == gv])
            ax.step(t, s, where="post", color=c, label=f"{nm} (n={int((grp == gv).sum())}, deaths={int(Et[grp == gv].sum())})")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("months")
        ax.set_ylabel("overall survival")
        ax.set_title(f"{label}\n{tr} (n={nl}); rank 1: {pw.replace('REACTOME_', '').lower()[:38]}\n"
                     f"log-rank p = {fmt_p(st['p'])}", fontsize=7)
        ax.legend(frameon=False, fontsize=6)
        panel(ax, next(L))
    fig.tight_layout()
    r1 = T[T["rank"] == 1]
    lines = "\n".join(f"- {r.setting}, {r.training} ({r.n_lines} lines): rank-1 {r.pathway.replace('REACTOME_', '').lower()}; "
                      f"treated n = {r.n}, deaths = {r.events}, p = {fmt_p(r.p)}"
                      + (f", HR {r.HR:.2f} ({r.HR_lo:.2f}-{r.HR_hi:.2f})" if np.isfinite(r.HR) else "")
                      + f"; inside the lines leave-one-out rho = {r.loo_rho_lines:.2f}" for r in r1.itertuples())
    savefig(fig, "G14_SuppFig15_cell_line_models", f"""
# Supp. Fig. 15 for the other cancers: models trained on cancer cell lines instead of organoids

The paper compared its organoid-trained model with the same pipeline trained on cancer cell lines (GDSC). We repeat
this for pancreas (gemcitabine), liver (sorafenib), ovary (platinum and paclitaxel) and sarcoma (gemcitabine).
Training data: GDSC2 cell lines (response = log IC50, lower = more sensitive), either all lines with RNA-seq or only
lines from the matching cancer type. Ovarian platinum: the model is trained on cisplatin (no carboplatin screen in
GDSC) and applied to patients who received carboplatin or cisplatin. Each of the top three pathways splits treated
TCGA patients at the median prediction.

- **a** Hazard ratio (predicted non-responders vs predicted responders, 95% CI) and log-rank p for every setting and
  each of the top 3 pathways. HR > 1 = the expected direction.
- **b onward** Kaplan-Meier curves for the rank-1 pathway of each setting.

{lines}
""")


if __name__ == "__main__":
    main()
