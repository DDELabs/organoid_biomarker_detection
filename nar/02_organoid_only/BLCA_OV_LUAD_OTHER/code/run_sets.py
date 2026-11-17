"""Organoid-only strategy (Kong et al. 2020) for the ovarian, sarcoma, lung (cell-line substitute) and Novartis PDX sets.

Usage: python run_sets.py OV|SARC|LUAD|PDX
Every (set, drug) runs the full engine (engine.Analysis): proximity, heatmap, 3 regressions + leave-one-out, KM rank 1
in treated and untreated, multi-pathway k, bootstrap, permutation null, cut-off sensitivity.
"""
import sys

import numpy as np
import pandas as pd

import common as K
import engine as E

PLAT = "No carboplatin/cisplatin screen in these organoids: oxaliplatin response (the only platinum screened) and the paper's cisplatin proximity are used as the platinum proxy."


def ov_sets():
    expr, r = K.vias_ovarian()
    texpr = K.tcga_expr("OV")
    osc, tsc = K.paired_scores("OV_vias2023", expr, texpr)
    resp = lambda d: r[r.drug == d].groupby("sample")["response"].median()  # noqa: E731
    base = dict(set_id="OV_vias2023", set_label="Ovarian HGSOC organoids (Vias 2023)", project="OV", osc=osc, tsc=tsc,
                surv=K.survival_current("OV"), outdir=K.PARENT / "OV_vias2023",
                response_label="AUC of 5-parameter viability curve", model_word="organoids")
    note2 = "14 organoids from 11 patients (3 patients give 2 organoids each)."
    for d, od, zd, note in [("CARBOPLATIN", "OXALIPLATIN", "CISPLATIN", PLAT), ("PACLITAXEL", "PACLITAXEL", "PACLITAXEL", ""),
                            ("CISPLATIN", "OXALIPLATIN", "CISPLATIN", PLAT), ("DOXORUBICIN", "DOXORUBICIN", "DOXORUBICIN", ""),
                            ("GEMCITABINE", "GEMCITABINE", "GEMCITABINE", ""),
                            ("DOCETAXEL", "PACLITAXEL", "DOCETAXEL", "No docetaxel screen: paclitaxel (same taxane class) response used.")]:
        yield dict(base, drug=d, model_drug=od, y=resp(od), zdrug=zd, note=(note + " " + note2).strip())


def sarc_sets():
    expr, r = K.load_preclinical("sarcoma_alshihabi2024")
    texpr = K.tcga_expr("SARC")
    osc, tsc = K.paired_scores("SARC_alshihabi2024", expr, texpr)
    resp = lambda d: r[r.drug == d].groupby("sample")["response"].median()  # noqa: E731
    base = dict(set_id="SARC_alshihabi2024", set_label="Sarcoma organoids (Al Shihabi 2024)", project="SARC", osc=osc, tsc=tsc,
                surv=K.survival_current("SARC"), outdir=K.PARENT / "SARC_alshihabi2024",
                response_label="% viability at a fixed dose", model_word="organoids")
    note = "Organoids span several sarcoma subtypes; TCGA-SARC is also mixed. Doxorubicin was not screened in these organoids."
    for d in ("GEMCITABINE", "DOCETAXEL"):
        yield dict(base, drug=d, model_drug=d, y=resp(d), zdrug=d, note=note)


def luad_sets():
    """GSE165611 drug data are not public (see report); substitute: GDSC1/2 lung adenocarcinoma cell lines (ln IC50,
    TCGA_DESC = LUAD) with CCLE expression (DepMap, via COSMIC ID), and CTRPv2 AUC for carboplatin."""
    s = pd.read_csv(K.PRECLIN / "ctrpv2_ccle" / "samples.tsv", sep="\t")
    expr, _ = K.load_preclinical("ctrpv2_ccle")
    cos2ach = {int(c): a for c, a in zip(s.COSMICID, s["sample"]) if pd.notna(c)}
    g = pd.concat([pd.read_csv(K.RAW / f"GDSC{i}_fitted.tsv.gz", sep="\t").assign(ds=f"GDSC{i}") for i in (2, 1)])
    g = g[g.TCGA_DESC == "LUAD"]
    g["ach"] = g.COSMIC_ID.map(cos2ach)
    g = g.dropna(subset=["ach"])
    luad = sorted(set(g.ach) & set(expr.columns))
    texpr = K.tcga_expr("LUAD")
    osc, tsc = K.paired_scores("LUAD_gdsc_ccle", expr[luad], texpr)

    def gd(name):
        d = g[g.DRUG_NAME.str.upper() == name]
        d = d[d.ds == "GDSC2"] if (d.ds == "GDSC2").sum() >= 15 else d[d.ds == "GDSC1"]
        return d.groupby("ach")["LN_IC50"].median(), d.ds.iloc[0] if len(d) else ""
    c = pd.read_csv(K.PRECLIN / "ctrpv2_ccle" / "response.tsv", sep="\t")
    base = dict(set_id="LUAD_cellline_substitute", set_label="Lung adenocarcinoma cell lines (GDSC/CTRP; substitute for GSE165611)",
                project="LUAD", osc=osc, tsc=tsc, surv=K.survival_current("LUAD"), outdir=K.PARENT / "LUAD_GSE165611_substitute",
                model_word="cell lines")
    for d, zd in [("CISPLATIN", "CISPLATIN"), ("PEMETREXED", "PEMETREXED"), ("PACLITAXEL", "PACLITAXEL"),
                  ("VINORELBINE", "VINORELBINE"), ("GEMCITABINE", "GEMCITABINE"), ("DOCETAXEL", "DOCETAXEL"),
                  ("ETOPOSIDE", "ETOPOSIDE"), ("ERLOTINIB", "ERLOTINIB")]:
        y, ds = gd(d)
        yield dict(base, drug=d, model_drug=d, y=y, zdrug=zd, response_label=f"ln IC50 ({ds})",
                   note=f"Lab models are {len(y)} GDSC lung adenocarcinoma cell lines ({ds}), not patient-derived cells.")
    cc = c[(c.drug == "CARBOPLATIN") & c["sample"].isin(luad)].groupby("sample")["response"].median()
    yield dict(base, drug="CARBOPLATIN", model_drug="CARBOPLATIN", y=cc, zdrug="CISPLATIN", response_label="AUC (CTRPv2)",
               note="Carboplatin response from CTRPv2 on the same GDSC/CCLE LUAD lines; cisplatin proximity used (carboplatin has no protein targets in DrugBank).")


def pdx_sets():
    expr, r = K.load_preclinical("pdx_novartis_gao2015")
    r = r[r.treatment_type == "single"]
    plan = [("BRCA", "BRCA", "PACLITAXEL", "PACLITAXEL"), ("BRCA", "BRCA", "TAMOXIFEN", "TAMOXIFEN"),
            ("BRCA", "BRCA", "TRASTUZUMAB", "TRASTUZUMAB"), ("PDAC", "PAAD", "GEMCITABINE", "GEMCITABINE"),
            ("CRC", "COAD", "FLUOROURACIL", "FLUOROURACIL"), ("CRC", "COAD", "CETUXIMAB", "CETUXIMAB"),
            ("NSCLC", "LUAD", "PACLITAXEL", "PACLITAXEL"), ("NSCLC", "LUAD", "ERLOTINIB", "ERLOTINIB"),
            ("CM", "SKCM", "DACARBAZINE", "DACARBAZINE")]
    cache = {}
    for tt, proj, d, zd in plan:
        rr = r[(r.tumour_type == tt) & (r.drug == d)]
        y = rr.groupby("sample")["response"].median()
        models = [m for m in y.index if m in expr.columns]
        if proj not in cache:
            mods = sorted(set(r.loc[r.tumour_type == tt, "sample"]) & set(expr.columns))
            texpr = K.tcga_expr(proj) if proj != "SKCM" else K.C.tcga_star_fpkm_uq("SKCM", ("01", "06"))
            cache[proj] = K.paired_scores(f"PDX_{tt}_{proj}", expr[mods], texpr)
            del texpr
        osc, tsc = cache[proj]
        yield dict(set_id=f"PDX_gao2015_{tt}", set_label=f"Novartis PDX ({tt}, Gao 2015)", project=proj, osc=osc, tsc=tsc,
                   surv=K.survival_current(proj), outdir=K.PARENT / "PDX_novartis_gao2015", drug=d, model_drug=d,
                   y=y.loc[models], zdrug=zd, response_label="BestAvgResponse (% tumour volume change)",
                   model_word="PDX models", note=f"One mouse per PDX model ({len(models)} {tt} models).")


SETS = {"OV": ov_sets, "SARC": sarc_sets, "LUAD": luad_sets, "PDX": pdx_sets}


def main(name):
    rows = []
    import json
    for cfg in SETS[name]():
        done = cfg["outdir"] / "tables" / f"{cfg['project']}_{cfg['drug']}_summary.json"
        if done.exists():  # resume after an interruption
            rows.append(json.load(open(done)))
            continue
        z, src = K.proximity(cfg.pop("zdrug"))
        anyd = K.any_drug(cfg["project"])
        surv = cfg["surv"]
        a = E.Analysis(z=z, z_source=src, treated=K.treated(cfg["project"], cfg["drug"]),
                       untreated=set(surv.index) - anyd, n_perm=1000, n_boot=2000, seed=0, **cfg)
        res = a.run()
        rows.append(res)
        print(name, cfg["drug"], {k: res.get(k) for k in ("n_models", "n_proximal", "biomarker", "treated_n",
                                                          "treated_p", "treated_hr", "untreated_p", "perm_empirical_p", "skipped")}, flush=True)
        out = cfg["outdir"] / "tables"
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out / f"{name}_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
