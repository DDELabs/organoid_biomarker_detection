"""Top pathway candidates for every cancer x drug pair (relaxed, integrated evidence).

Local evidence (that cancer only): response notes, survival interaction, independent cohorts of that cancer.
Backing evidence (drug level): pan-cancer patients, cell lines, organoids (from run_integrated.py).
A source supports a pathway at |z| >= 1.64 in the same direction. Ranking: number of supporting sources,
then mean |z| of the local sources. Needs robust/results/integrated/evidence.tsv.gz.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_integrated as RI  # noqa: E402
import run_survlearn as RS  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.curated import load_curated_trial, trial_catalog  # noqa: E402

OUT = RI.OUT
T = 1.64
MAP = [("breast", ["BRCA"]), ("colorectal", ["COAD", "READ"]), ("ovarian", ["OV"]), ("bladder", ["BLCA"]),
       ("oesophag", ["ESCA"]), ("esophag", ["ESCA"]), ("nsclc", ["LUAD", "LUSC"]), ("lung", ["LUAD", "LUSC"]),
       ("gastric", ["STAD"]), ("stomach", ["STAD"]), ("glio", ["GBM", "LGG"]), ("melanoma", ["SKCM"]),
       ("hepato", ["LIHC"]), ("liver", ["LIHC"]), ("renal", ["KIRC"]), ("kidney", ["KIRC"]), ("pancrea", ["PAAD"]),
       ("head", ["HNSC"]), ("cervi", ["CESC"]), ("endometr", ["UCEC"])]
NAMES = {"BRCA": "breast", "COAD": "colon", "READ": "rectal", "OV": "ovarian", "BLCA": "bladder", "ESCA": "oesophageal",
         "LUAD": "lung adeno", "LUSC": "lung squamous", "STAD": "stomach", "GBM": "glioblastoma", "LGG": "lower-grade glioma",
         "SKCM": "melanoma", "LIHC": "liver", "KIRC": "kidney", "PAAD": "pancreatic", "HNSC": "head & neck",
         "CESC": "cervical", "UCEC": "endometrial", "UCS": "uterine carcinosarcoma", "SARC": "sarcoma", "TGCT": "testicular",
         "MESO": "mesothelioma", "PRAD": "prostate", "DLBC": "lymphoma", "THYM": "thymoma", "ACC": "adrenal", "CHOL": "bile duct"}


def projects_of(cancer):
    c = str(cancer).lower()
    for k, v in MAP:
        if k in c:
            return v
    return []


def main():
    Z, cdr, long = RS.load()
    paths = list(Z.columns)
    fold = SL.folds_by_patient(cdr.index, 5, seed=0)
    cdr["g"] = SL.general_risk(Z, cdr, "time", "event", fold)
    cdr["g"] /= cdr["g"].std()
    E = pd.read_csv(OUT / "evidence.tsv.gz", sep="\t")
    cat = trial_catalog()
    cat = cat[cat["usable"].astype(str).str.lower().eq("true")]
    trials = {}
    for cid, cancer in zip(cat["cohort_id"], cat["cancer"]):
        try:
            _, clin = load_curated_trial(cid, canonical=False)
            if "responder" in clin and "drugs" in clin:
                trials[cid] = (RI.cohort_scores(cid).reindex(columns=paths).fillna(0), clin, projects_of(cancer))
        except Exception:  # noqa: BLE001
            pass
    n = long.drop_duplicates(["patient", "drug"]).groupby(["drug", "project"]).size()
    rows, pairs = [], []
    for (d, c), npat in n[n >= 15].items():
        Ed = E[E["drug"] == d].set_index("pathway").reindex(paths)
        if Ed["drug"].isna().all():
            continue
        ev = pd.DataFrame(index=paths)
        info = {"cancer": c, "drug": d, "patients_given_drug": int(npat)}
        rec = long[(long["drug"] == d) & (long["project"] == c)]
        lab = rec[rec["responder"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")
        if len(lab) >= 15 and lab["responder"].value_counts().min() >= 4:
            ev["local_notes"] = RI.corr_z(Z.loc[lab.index].values, lab["responder"].astype(float).values)
            info["notes_n"] = len(lab)
        first = set(rec.loc[rec["before_event"], "patient"])
        D = cdr[cdr["project"] == c].copy()
        D["treated"] = D.index.isin(first).astype(float)
        D["other_systemic"] = D.index.isin(set(long.loc[long["before_event"] & (long["drug"] != d), "patient"])).astype(float)
        et, eu = D.loc[D.treated == 1, "event"].sum(), D.loc[D.treated == 0, "event"].sum()
        if D["treated"].sum() >= 15 and et >= 8 and eu >= 8:
            try:
                bt, it = RI.score_test(D[D.treated == 1], Z.loc[D.index[D.treated == 1]].values)
                bu, iu = RI.score_test(D[D.treated == 0], Z.loc[D.index[D.treated == 0]].values)
                ev["local_survival"] = -(bt - bu) / np.sqrt(1 / it + 1 / iu)
                info["survival_events"] = f"{int(et)}/{int(eu)}"
            except Exception:  # noqa: BLE001
                pass
        zs, ns = [], []
        for cid, (sc, clin, projs) in trials.items():
            if c not in projs:
                continue
            got = clin["drugs"].fillna("").astype(str).str.upper().str.split(";").map(lambda x: d in x)
            sel = [s for s in clin.index[got & clin["responder"].notna()] if s in sc.index]
            y = clin.loc[sel, "responder"].astype(float)
            if len(sel) >= 15 and y.value_counts().min() >= 4:
                zs.append(RI.corr_z(sc.loc[sel].values, y.values))
                ns.append(len(sel))
        if zs:
            ev["local_external"] = RI.stouffer(zs, ns)
            info["external_n"] = int(sum(ns))
        local = [x for x in ("local_notes", "local_survival", "local_external") if x in ev]
        if not local:
            continue
        ev["pan_cancer"] = Ed[["notes", "survival", "external"]].mean(axis=1, skipna=True)
        ev["cells"], ev["organoids"], ev["proximal"] = Ed["cells"], Ed["organoids"], Ed["network_z"] <= -1.2816
        src = local + ["pan_cancer", "cells", "organoids"]
        V = ev[src].values.astype(float)
        best = None
        for sgn in (1, -1):
            cnt = np.nansum(sgn * V >= T, axis=1)
            loc = np.nansum(sgn * ev[local].values >= T, axis=1)
            strength = np.nanmean(sgn * ev[local].values, axis=1)
            cur = pd.DataFrame({"direction": sgn, "support": cnt, "local_support": loc, "local_strength": strength}, index=paths)
            best = cur if best is None else best.where((best["support"] > cur["support"]) | ((best["support"] == cur["support"]) & (best["local_strength"] >= cur["local_strength"])), cur)
        ev = ev.join(best)
        top = ev[(ev["local_support"] >= 1) & (ev["support"] >= 2)].sort_values(["support", "local_strength"], ascending=False).head(5)
        pairs.append({**info, "local_sources": len(local), "candidates": len(top), "best_support": int(top["support"].max()) if len(top) else 0})
        for rank, (p, r) in enumerate(top.iterrows(), 1):
            rows.append({**info, "rank": rank, "pathway": p, **r.to_dict()})
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "top_by_cancer_drug.tsv", sep="\t", index=False)
    P = pd.DataFrame(pairs)
    P.to_csv(OUT / "cancer_drug_pairs.tsv", sep="\t", index=False)
    f = lambda v: "" if pd.isna(v) else f"{v:+.1f}"  # noqa: E731
    with open(OUT / "TOP_BY_CANCER_DRUG.md", "w") as fh:
        fh.write("# Top pathway candidates per cancer and drug (relaxed, integrated evidence)\n\n"
                 "Numbers are z-scores; + means higher pathway activity goes with better response or more benefit. "
                 "A source supports a pathway at 1.64 or more in the row's direction. These are candidates, not validated biomarkers.\n\n"
                 "| Cancer | Drug | Patients | # | Pathway | Direction | Sources agreeing | Notes | Survival | Other cohorts | Pan-cancer | Cell lines | Organoids | Network |\n"
                 "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in R.sort_values(["cancer", "drug", "rank"]).itertuples():
            fh.write(f"| {NAMES.get(r.cancer, r.cancer)} | {r.drug.title()} | {r.patients_given_drug} | {r.rank} | "
                     f"{r.pathway.replace('REACTOME_', '').replace('_', ' ').lower()} | {'higher = better' if r.direction > 0 else 'higher = worse'} | {int(r.support)} | "
                     f"{f(getattr(r, 'local_notes', np.nan))} | {f(getattr(r, 'local_survival', np.nan))} | {f(getattr(r, 'local_external', np.nan))} | "
                     f"{f(r.pan_cancer)} | {f(r.cells)} | {f(r.organoids)} | {'yes' if r.proximal else ''} |\n")
    print(f"pairs with >= 15 patients and local data: {len(P)}; with >= 1 candidate: {(P.candidates > 0).sum()}; rows: {len(R)}")
    print(P.best_support.value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
