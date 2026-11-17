"""Compute every set x drug x cohort analysis and store results (pickles + tables).

Usage: python run_all.py [worker_index n_workers]      (2 workers max; rows are split round-robin)
"""
import pickle
import sys

import numpy as np
import pandas as pd

from common import (ALT_CUTOFFS, DATA, DRUG_MAP, PAPER_CUTOFF, SETS, authors_inputs, cdr, proximity, scores,
                    set_response, treatment_table)
from engine import Combo, MLS

OUT = DATA / "results"
MIN_TREATED = 12          # ~15 requested; 12-14 kept and flagged as borderline
import os
N_BOOT, N_PERM = int(os.environ.get("N_BOOT", 10000)), int(os.environ.get("N_PERM", 1000))


def cohort_inputs(key, cohort, cd):
    """Patient scores, survival, treated table and untreated set for a cohort (COAD, READ, PAAD, LIHC, COADREAD)."""
    projects = SETS[key][2]
    sc = scores(key, projects)
    parts = ["COAD", "READ"] if cohort == "COADREAD" else [cohort]
    tsc = pd.concat([sc[p] for p in parts], axis=1)
    tsc = tsc.loc[:, ~tsc.columns.duplicated()]
    ok = cd["OS.time"].notna() & cd["OS"].notna()
    sur = pd.DataFrame({"months": cd.loc[ok, "OS.time"] / (365 / 12), "event": cd.loc[ok, "OS"].astype(int)})
    tt = pd.concat([treatment_table(p) for p in parts])
    anyrec = set().union(*[treatment_table(p).attrs["any_record"] for p in parts])
    untreated = set(tsc.columns) - anyrec
    return sc["models"], tsc, sur, tt, untreated


def specs():
    out = []
    # A. authors' committed files, exact reproduction
    for d in ("FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CETUXIMAB"):
        out.append(("vandewetering", "COAD_authors", d, d))
    plan = {
        "vandewetering": {"COAD": ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CETUXIMAB"],
                          "READ": ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CETUXIMAB"],
                          "COADREAD": ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN"]},
        "tiriac": {"PAAD": ["GEMCITABINE", "FLUOROURACIL", "OXALIPLATIN", "SN-38", "PACLITAXEL"]},
        "shi": {"PAAD": ["GEMCITABINE", "FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "PACLITAXEL", "MITOMYCIN"]},
        "ji": {"LIHC": ["SORAFENIB", "GEMCITABINE", "DOXORUBICIN", "FLUOROURACIL", "LENVATINIB", "REGORAFENIB"]},
        "broutier": {"LIHC": ["SORAFENIB", "GEMCITABINE", "DOXORUBICIN", "CISPLATIN", "FLUOROURACIL"]},
        "isella": {"COAD": ["CETUXIMAB"]},
    }
    for key, coh in plan.items():
        for cohort, drugs in coh.items():
            for d in drugs:
                out.append((key, cohort, d, DRUG_MAP.get(d, d)))
    return out


def build(spec, cd):
    key, cohort, mdrug, pdrug = spec
    if cohort == "COAD_authors":
        org, resp, pat, sur, treated, anydrug = authors_inputs(pdrug)
        untreated = set(sur.index) - anydrug
        y = resp
        osc, tsc = org, pat
        note = "authors' committed NES + 2019 biotab files"
    else:
        osc, tsc, sur, tt, untreated = cohort_inputs(key, cohort, cd)
        treated = set(tt.loc[tt["drug"] == pdrug, "patient"])
        y = set_response(key)[mdrug] if mdrug in set_response(key).columns else None
        note = "current GDC STAR FPKM-UQ + TCGA-CDR OS + GDC biotab drug file"
    z, zcol = proximity(pdrug)
    meta = dict(set=key, folder=SETS[key][0], set_label=SETS[key][4], cohort=cohort, model_drug=mdrug,
                patient_drug=pdrug, proximity_column=zcol, data=note,
                metric="ln IC50" if (key == "vandewetering") else
                ("tumour volume change (3 wk)" if key == "isella" else "AUC"))
    if y is None or z is None:
        return None, meta
    return Combo(osc, y, tsc, sur, treated, untreated, z, PAPER_CUTOFF, meta), meta


def run(spec, cd):
    c, meta = build(spec, cd)
    tag = f"{meta['folder']}__{meta['cohort']}__{meta['model_drug']}"
    if c is None:
        meta["skip"] = "no drug response or no proximity"
        return tag, meta
    meta.update(n_models=len(c.samples), n_proximal=len(c.feats), n_treated=len(c.tp), n_untreated=len(c.up),
                treated_events=int(c.E.sum()), n_pathways_total=len(c.common))
    if len(c.tp) < MIN_TREATED or len(c.feats) < 1:
        meta["skip"] = f"{len(c.tp)} treated patients with survival and expression (< {MIN_TREATED})"
        return tag, meta
    meta["borderline"] = len(c.tp) < 15
    ev, cf = c.evaluate()
    loo = c.loo()
    # isogenic HCT-8 5-FU-sensitive/-resistant cells exist only for colorectal 5-FU (paper Fig. 4b/4e)
    iso_p, iso_tab = None, None
    if meta["patient_drug"] == "FLUOROURACIL" and meta["cohort"] in ("COAD_authors", "COAD", "READ", "COADREAD"):
        import iso
        iso_tab = iso.iso_table(iso.authors_nes() if meta["cohort"] == "COAD_authors" else iso.current_nes())[0]
        iso_p = iso_tab["p"]
    boot, boot_rank = c.bootstrap(N_BOOT, seed=0, iso_p=iso_p)
    pp, pn = c.single_pathway_p(1), c.single_pathway_p(-1)
    co = cf["Ridge"]
    signed = pd.Series([pp[f] if co[i] > 0 else pn[f] for i, f in enumerate(c.feats)], index=c.feats)
    obs_rank = c.top_rank(co, c.feats, signed, iso_p)
    r1, best = c.permutation(N_PERM, seed=1)
    cut = []
    for cu in ALT_CUTOFFS:
        f = c.set_cutoff(cu)
        if len(f) == 0:
            cut.append({"cutoff": cu, "n_pathways": 0})
            continue
        e, _ = c.evaluate(f)
        r = e[(e.model == "Ridge") & (e.k == 1)].iloc[0]
        cut.append({"cutoff": cu, "n_pathways": len(f), "rank1_pathway": r["pathways"],
                    "treated_p": r["treated_p"], "treated_HR": r["treated_HR"], "treated_HR_lo": r["treated_HR_lo"],
                    "treated_HR_hi": r["treated_HR_hi"], "untreated_p": r.get("untreated_p", np.nan),
                    "best_p_27_tests": e["treated_p"].min()})
    c.set_cutoff(PAPER_CUTOFF)
    # KM groups for rank 1..10 (Ridge), treated and untreated
    km = {}
    for who, Tm, Em in (("treated", c.T, c.E), ("untreated", c.Tu, c.Eu)):
        for k in range(1, min(10, len(c.feats)) + 1):
            g, used, _ = c.groups("Ridge", k, who=who)
            km[(who, k)] = (Tm, Em, g, used)
    obs = ev.loc[(ev.model == "Ridge") & (ev.k == 1), "treated_p"].iloc[0]
    res = dict(meta=meta, ev=ev, coefs={m: pd.Series(cf[m], index=c.feats) for m in MLS},
               loo={m: (pd.Series(v[0], index=c.samples), v[1], v[2]) for m, v in loo.items()},
               y=pd.Series(c.y, index=c.samples), boot=boot, perm_rank1=r1, perm_best=best,
               cut=pd.DataFrame(cut), km=km, nes=c.osc.loc[c.feats], z=c.z_all.reindex(c.common).dropna(),
               feats=c.feats, tp=c.tp, up=c.up, Pall=c.P[c.feats], iso_tab=iso_tab,
               single_p=pd.DataFrame({"p_pos": pp, "p_neg": pn, "z": c.z_all.reindex(c.common)}),
               boot_rank=boot_rank, obs_rank=obs_rank,
               boot_rank_emp_p=float((boot_rank <= obs_rank).mean()) if obs_rank <= len(c.feats) else np.nan,
               boot_emp_p=float((boot <= obs + 1e-12).mean()),
               perm_rank1_emp_p=float((np.sum(r1 <= obs + 1e-12) + 1) / (N_PERM + 1)),
               perm_best_emp_p=float((np.sum(best <= ev["treated_p"].min() + 1e-12) + 1) / (N_PERM + 1)))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"{tag}.pkl", "wb") as f:
        pickle.dump(res, f)
    return tag, meta


def main():
    w, nw = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
    cd = cdr()
    rows = []
    for i, s in enumerate(specs()):
        if i % nw != w:
            continue
        key, cohort, mdrug, _ = s
        done = OUT / f"{SETS[key][0]}__{cohort}__{mdrug}.pkl"
        if done.exists():                       # resume: reuse finished results
            with open(done, "rb") as fh:
                tag, meta = done.stem, pickle.load(fh)["meta"]
        else:
            tag, meta = run(s, cd)
        rows.append(meta)
        print(tag, meta.get("skip", "done"), meta.get("n_treated"), flush=True)
    pd.DataFrame(rows).to_csv(OUT / f"run_log_{w}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
