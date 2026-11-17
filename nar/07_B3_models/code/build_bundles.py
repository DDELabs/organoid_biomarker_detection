"""Build one cached data bundle per (lab-model set, patient cohort) pair, reusing the loaders and caches of the
earlier agents (nar/02_organoid_only/code for GI sets, nar/02_organoid_only/BLCA_OV_LUAD_OTHER/code for the rest).

A bundle holds: lab-model pathway ssGSEA (osc), patient ssGSEA (tsc), lab-model and patient gene expression on
shared genes, survival (T months, E), untreated patients, and per drug: response y, treated patients, proximity z.
Only drugs with >= 15 treated patients (with expression + survival) are kept.

Usage: python build_bundles.py GI | OTHER
"""
import importlib
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
B3 = HERE.parent
OUT = B3 / "_data" / "raw" / "bundles"
OUT.mkdir(parents=True, exist_ok=True)
NAR = B3.parent
MIN_TREATED = 15


def save(name, osc, tsc, gm, gp, surv, untreated, drugs, meta):
    pats = [p for p in tsc.columns if p in surv.index]
    keep = {}
    for d, v in drugs.items():
        tp = [p for p in pats if p in v["treated"]]
        ys = v["y"].dropna()
        ys = ys[[s for s in ys.index if s in osc.columns]]
        if len(tp) < MIN_TREATED or v["z"] is None or len(ys) < 5:
            print(name, d, "skipped", len(tp), len(ys), flush=True)
            continue
        keep[d] = dict(y=ys.astype(float), treated=tp, z=v["z"], zsrc=v.get("zsrc", ""), model_drug=v.get("model_drug", d),
                       note=v.get("note", ""))
    if not keep:
        print(name, "no drugs kept")
        return
    genes = gm.index.intersection(gp.index)
    gm, gp = gm.loc[genes], gp.loc[genes]
    gm = gm.loc[:, ~gm.columns.duplicated()]
    gp = gp.loc[:, ~gp.columns.duplicated()]
    ok = (gm.var(axis=1) > 0) & (gp.var(axis=1) > 0)
    b = dict(name=name, osc=osc.astype(np.float32), tsc=tsc[pats].astype(np.float32),
             gm=gm[ok].astype(np.float32), gp=gp.loc[ok, [p for p in pats if p in gp.columns]].astype(np.float32),
             surv=surv.loc[pats, ["T", "E"]].astype(float), untreated=[p for p in pats if p in untreated],
             drugs=keep, meta=meta)
    pickle.dump(b, open(OUT / f"{name}.pkl", "wb"), protocol=4)
    print("saved", name, {d: (len(v["y"]), len(v["treated"])) for d, v in keep.items()}, "genes", int(ok.sum()),
          "patients", len(pats), flush=True)


def gi():
    sys.path.insert(0, str(NAR / "02_organoid_only" / "code"))
    K = importlib.import_module("common")
    RA = importlib.import_module("run_all")
    cd = K.cdr()
    # --- paper data: authors' committed NES (COAD) -------------------------------------------------------------
    oexpr = K.set_expression("vandewetering")
    tex = K.tcga_expr("COAD")
    drugs = {}
    for d in ("FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN"):
        org, resp, pat, sur, treated, anydrug = K.authors_inputs(d)
        z, zc = K.proximity(d)
        drugs[d] = dict(y=resp, treated=treated, z=z, zsrc=f"paper table ({zc})")
    org.columns = [c.upper() for c in org.columns]
    oe = oexpr.copy()
    oe.columns = [c.upper() for c in oe.columns]
    tex.columns = [c[:12] for c in tex.columns]
    surv = sur.rename(columns={"months": "T", "event": "E"})
    save("COAD_paper", org, pat, oe, tex, surv, set(sur.index) - anydrug, drugs,
         dict(label="Colorectal organoids (van de Wetering 2015) - paper source data", project="COAD", models="organoids",
              metric="ln IC50", data="authors' committed organoid/TCGA NES, 2019 biotab drug + survival"))
    # --- current data -------------------------------------------------------------------------------------------
    plan = {("vandewetering", "COADREAD"): ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CETUXIMAB"],
            ("tiriac", "PAAD"): ["GEMCITABINE", "FLUOROURACIL", "OXALIPLATIN"],
            ("shi", "PAAD"): ["GEMCITABINE", "FLUOROURACIL", "OXALIPLATIN"],
            ("ji", "LIHC"): ["SORAFENIB", "GEMCITABINE", "DOXORUBICIN", "FLUOROURACIL"],
            ("isella", "COAD"): ["CETUXIMAB"]}
    names = {"vandewetering": "COADREAD_vdW", "tiriac": "PAAD_tiriac", "shi": "PAAD_shi", "ji": "LIHC_ji", "isella": "COAD_isella"}
    for (key, cohort), dl in plan.items():
        osc, tsc, sur, tt, untreated = RA.cohort_inputs(key, cohort, cd)
        resp = K.set_response(key)
        parts = ["COAD", "READ"] if cohort == "COADREAD" else [cohort]
        gp = pd.concat([K.tcga_expr(p) for p in parts], axis=1)
        gp.columns = [c[:12] for c in gp.columns]
        tsc = tsc.copy()
        tsc.columns = [c[:12] for c in tsc.columns]
        tsc = tsc.loc[:, ~tsc.columns.duplicated()]
        gm = K.set_expression(key)
        drugs = {}
        for d in dl:
            pdrug = K.DRUG_MAP.get(d, d)
            tr = set(tt.loc[tt["drug"] == pdrug, "patient"])
            if d == "CETUXIMAB" and cohort == "COADREAD":       # the earlier agents used COAD only for cetuximab
                tr = {p for p in tr if p in set(K.treatment_table("COAD")["patient"])}
            z, zc = K.proximity(pdrug)
            if d not in resp.columns:
                continue
            drugs[d] = dict(y=resp[d], treated=tr, z=z, zsrc=f"paper table ({zc})")
        surv = sur.rename(columns={"months": "T", "event": "E"})
        save(names[key], osc, tsc, gm, gp, surv, untreated, drugs,
             dict(label=K.SETS[key][4], project=cohort, models="PDX models" if key == "isella" else "organoids",
                  metric="ln IC50" if key == "vandewetering" else ("tumour volume change" if key == "isella" else "AUC"),
                  data="current GDC STAR FPKM-UQ + TCGA-CDR OS + GDC biotab drugs"))
        del gp


def other():
    sys.path.insert(0, str(NAR / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "code"))
    K = importlib.import_module("common")
    RB = importlib.import_module("run_blca")
    RS = importlib.import_module("run_sets")
    # --- bladder ------------------------------------------------------------------------------------------------
    oexpr, _ = K.load_preclinical("bladder_lee2018")
    texpr = K.tcga_expr("BLCA")
    texpr.columns = [c[:12] for c in texpr.columns]
    anyd = K.any_drug("BLCA")
    for v in ("paper", "current"):
        osc, tsc, surv = RB.paper_version() if v == "paper" else RB.current_version()
        tsc = tsc.copy()
        tsc.columns = [c[:12] for c in tsc.columns]
        tsc = tsc.loc[:, ~tsc.columns.duplicated()]
        drugs = {}
        for d, (od, pdrug, note) in RB.DRUGS.items():
            if v == "paper" and d != "CISPLATIN":
                continue
            z, src = K.proximity(pdrug)
            drugs[d] = dict(y=RB.response(od), treated=K.treated("BLCA", d), z=z, zsrc=src, model_drug=od, note=note)
        save(f"BLCA_{v}", osc, tsc, oexpr, texpr, surv, set(surv.index) - anyd, drugs,
             dict(label=f"Bladder organoids (Lee 2018) - {'paper source data' if v == 'paper' else 'current GDC data'}",
                  project="BLCA", models="organoids", metric="AUC (CoderData refit)",
                  data="MOESM3 NES + biotab survival" if v == "paper" else "current GDC"))
    del texpr
    # --- ovary, sarcoma, lung cell lines, PDX ------------------------------------------------------------------
    gexpr = {}
    expr, _ = K.vias_ovarian()
    gexpr["OV_vias2023"] = expr
    gexpr["SARC_alshihabi2024"] = K.load_preclinical("sarcoma_alshihabi2024")[0]
    gexpr["LUAD_cellline_substitute"] = K.load_preclinical("ctrpv2_ccle")[0]
    pdx = K.load_preclinical("pdx_novartis_gao2015")[0]
    for name in ("OV", "SARC", "LUAD", "PDX"):
        groups = {}
        for cfg in RS.SETS[name]():
            groups.setdefault(cfg["set_id"], []).append(cfg)
        for sid, cfgs in groups.items():
            c0 = cfgs[0]
            proj = c0["project"]
            tex = K.tcga_expr(proj) if proj != "SKCM" else K.C.tcga_star_fpkm_uq("SKCM", ("01", "06"))
            tex.columns = [c[:12] for c in tex.columns]
            gm = gexpr.get(sid, pdx)
            gm = gm[[c for c in c0["osc"].columns if c in gm.columns]]
            tsc = c0["tsc"].copy()
            tsc.columns = [c[:12] for c in tsc.columns]
            tsc = tsc.loc[:, ~tsc.columns.duplicated()]
            surv = c0["surv"]
            anyd = K.any_drug(proj)
            drugs = {}
            for cfg in cfgs:
                z, src = K.proximity(cfg["zdrug"])
                drugs[cfg["drug"]] = dict(y=cfg["y"], treated=K.treated(proj, cfg["drug"]), z=z, zsrc=src,
                                          model_drug=cfg["model_drug"], note=cfg.get("note", ""))
            short = {"OV_vias2023": "OV_vias", "SARC_alshihabi2024": "SARC_alshihabi", "LUAD_cellline_substitute": "LUAD_cells"}
            save(short.get(sid, sid.replace("PDX_gao2015", "PDX")), c0["osc"], tsc, gm, tex, surv, set(surv.index) - anyd, drugs,
                 dict(label=c0["set_label"], project=proj, models=c0.get("model_word", "organoids"),
                      metric=c0.get("response_label", ""), data="current GDC"))
            del tex


if __name__ == "__main__":
    {"GI": gi, "OTHER": other}[sys.argv[1]]()
