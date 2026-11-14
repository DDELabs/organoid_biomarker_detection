"""The paper's headline single-pathway method on every cancer x drug pair that has patients and a lab model.

For each TCGA cancer x drug with >= 15 treated patients (>= 5 deaths) and each available lab model of that tissue
(patient-derived organoids / PDX where they exist, otherwise cancer cell lines of the matching tissue):
  proximal pathways (paper's table, else STRING v12) -> ridge on lab-model drug response -> top pathway ->
  median split of treated patients -> log-rank on overall survival; plus untreated control, direction, shuffle test.
Pathway scores: ssGSEA (gseapy), as in the paper.
Usage: nohup python robust/run_single_pathway_all.py [--perm 500] &
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_paper_single_pathway as SP  # noqa: E402
import run_survlearn as RS  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd.curated import PRECLIN, load_curated_preclinical  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "single_pathway_all"
ALIAS = {"5-FLUOROURACIL": "FLUOROURACIL", "SN-38": "IRINOTECAN"}
ORGANOID = {"colorectal_vandewetering2015": ["COAD", "READ"], "bladder_lee2018": ["BLCA"], "liver_ji2023": ["LIHC"],
            "liver_broutier2017": ["LIHC"], "pancreas_shi2022": ["PAAD"], "pancreas_tiriac2018": ["PAAD"],
            "sarcoma_alshihabi2024": ["SARC"], "colorectal_pdx_isella2017": ["COAD", "READ"]}
LINEAGE = {"lung": ["LUAD", "LUSC"], "colorectal": ["COAD", "READ"], "central_nervous_system": ["GBM", "LGG"], "skin": ["SKCM"],
           "breast": ["BRCA"], "ovary": ["OV"], "pancreas": ["PAAD"], "gastric": ["STAD"], "upper_aerodigestive": ["HNSC"],
           "uterus": ["UCEC", "UCS"], "esophagus": ["ESCA"], "urinary_tract": ["BLCA"], "kidney": ["KIRC", "KIRP"],
           "liver": ["LIHC"], "soft_tissue": ["SARC"], "thyroid": ["THCA"], "prostate": ["PRAD"], "bile_duct": ["CHOL"]}
CELLS = ("ctrpv2_ccle", "prism_repurposing")


def model_scores(name, gs_all):
    def compute():
        expr, _ = load_curated_preclinical(name)
        gs = {k: [g for g in v if g in expr.index] for k, v in gs_all.items()}
        return scoring.ssgsea(expr, {k: v for k, v in gs.items() if len(v) >= 5})
    return scoring.cached(f"sp_model_{name}_ssgsea", compute)


def tcga_scores(proj, gs_all):
    def compute():
        expr = C.tcga_star_fpkm_uq(proj, ("01", "06") if proj == "SKCM" else ("01",))
        gs = {k: [g for g in v if g in expr.index] for k, v in gs_all.items()}
        sc = scoring.ssgsea(expr, {k: v for k, v in gs.items() if len(v) >= 5})
        f = C.EXTERNAL / f"TCGA-{proj}.star_fpkm-uq.tsv.gz" if hasattr(C, "EXTERNAL") else None
        return sc
    return scoring.cached(f"sp_tcga_{proj}_ssgsea", compute)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=500)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    gs_all = gene_sets(("REACTOME",))
    paper = N.precomputed_proximity()
    long, cdr = R.load_response_labels(), R.load_cdr()
    from run_atlas_relaxed import SUPPORTIVE, THERAPY
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    # lab-model responses
    models = {}
    for name in list(ORGANOID) + list(CELLS):
        r = pd.read_csv(PRECLIN / name / "response.tsv", sep="\t", low_memory=False)
        r["drug"] = r["drug"].str.upper().map(lambda x: ALIAS.get(x, x))
        models[name] = r
    # eligible patient pairs
    n = long.drop_duplicates(["patient", "drug"]).groupby(["project", "drug"])["patient"].apply(set)
    pairs = []
    for (proj, drug), pats in n.items():
        ok = [p for p in pats if p in cdr.index and pd.notna(cdr.loc[p, "OS.time"]) and cdr.loc[p, "OS.time"] > 0]
        if len(ok) >= 15 and cdr.loc[ok, "OS"].sum() >= 5:
            z = paper[drug].dropna() if drug in paper.columns else RS.proximity(drug)
            if z is None:
                continue
            srcs = []
            for name, projs in ORGANOID.items():
                if proj in projs:
                    k = models[name].loc[models[name]["drug"] == drug, "sample"].nunique()
                    if k >= 8:
                        srcs.append((name, None, "organoid/PDX"))
            for lin, projs in LINEAGE.items():
                if proj in projs:
                    for name in CELLS:
                        m = models[name]
                        k = m.loc[(m["drug"] == drug) & (m["lineage"] == lin), "sample"].nunique()
                        if k >= 8:
                            srcs.append((name, lin, "cell lines"))
            for s in srcs:
                pairs.append((proj, drug, *s, "paper table" if drug in paper.columns else "STRING v12"))
    print(f"{len(pairs)} runs over {len({(p[0], p[1]) for p in pairs})} cancer-drug pairs, {len({p[0] for p in pairs})} cancers", flush=True)
    rng = np.random.default_rng(0)
    rows = []
    for proj, drug, name, lin, kind, psrc in sorted(pairs):
        try:
            tsc, osc = tcga_scores(proj, gs_all), model_scores(name, gs_all)
            r = models[name]
            r = r[(r["drug"] == drug) & ((r["lineage"] == lin) if lin else True)].groupby("sample")["response"].median().dropna()
            samples = [s for s in r.index if s in osc.columns]
            z = paper[drug].dropna() if drug in paper.columns else RS.proximity(drug)
            feats = [p for p in z.index[z <= -1.2816] if p in osc.index and p in tsc.index]
            if len(samples) < 8 or len(feats) < 3:
                continue
            pats = [p for p in tsc.columns if p in cdr.index and pd.notna(cdr.loc[p, "OS.time"]) and pd.notna(cdr.loc[p, "OS"]) and cdr.loc[p, "OS.time"] > 0]
            P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[feats, pats].T.values), index=pats, columns=feats)
            lp = long[long["project"] == proj]
            treated, anydrug = set(lp.loc[lp["drug"] == drug, "patient"]), set(R.load_response_labels(projects=[proj])["patient"])
            tp, up = [p for p in pats if p in treated], [p for p in pats if p not in anydrug]
            te = lambda ps: (cdr.loc[ps, "OS.time"].values / (365 / 12), cdr.loc[ps, "OS"].values.astype(int))  # noqa: E731
            X = StandardScaler().fit_transform(osc.loc[feats, samples].T.values)
            res = SP.analyse(f"{proj} {drug}", X, r.loc[samples].values, feats, P.loc[tp].values, *te(tp), P.loc[up].values, *te(up), a.perm, rng)
            res.update(cancer=proj, drug=drug, model=name + (f" [{lin}]" if lin else ""), model_type=kind, proximity=psrc)
            rows.append(res)
            print(proj, drug, name, kind, "p=%.3f" % res["p_treated"], "right" if res["direction_correct"] else "wrong", "shuffle=%.3f" % res["shuffled_as_good"], flush=True)
            pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False)
        except Exception as exc:  # noqa: BLE001
            print(proj, drug, name, "FAILED", str(exc)[:100], flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
