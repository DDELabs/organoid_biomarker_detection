"""Pathway biomarkers in the cancers TCGA cannot cover: liver, kidney, glioma/GBM, melanoma.

Uses the public treated cohorts in data/curated/trials (Reactome rank scores, z-scored in cohort).

  response cohorts   per-pathway logistic score-test z (+ = higher in responders), BH q, 200-bootstrap
                     top-10 stability; top 5 response + 5 resistance pathways
                     GSE104580 liver TACE | EMTAB3267 kidney sunitinib | BRAUN2020 kidney nivolumab,
                     everolimus | GSE91061, GSE78220 melanoma anti-PD-1
  predictive tests   pathway x treatment interaction (Cox OS/PFS; + = treatment benefit grows with pathway)
                     BRAUN2020 CM-025 (randomised nivolumab vs everolimus, PFS)
                     CGGA693 TMZ vs no TMZ (OS; adjusted for grade, IDH, 1p/19q, MGMT, radiotherapy, age)
                     GSE7696 GBM RT+TMZ vs RT (OS; adjusted for MGMT, age)
                     GBM meta: Stouffer of CGGA WHO IV and GSE7696 interaction z

Outputs: robust/results/cohort_biomarkers/{response_pathways.tsv.gz, top10_cohorts.tsv,
         interaction.tsv.gz, top_interactions.tsv, COHORT_BIOMARKERS.md}
Usage: python robust/run_cohort_biomarkers.py [--boot 200]
"""
import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import run_atlas_external as AE  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import curated as CU  # noqa: E402
from run_atlas_cancer_specific import bh  # noqa: E402
from run_atlas_relaxed import logit_score_z  # noqa: E402

OUT = RESULTS / "cohort_biomarkers"
RESPONSE = [("GSE104580", "liver (HCC)", "TACE", None),
            ("EMTAB3267", "kidney (ccRCC)", "SUNITINIB", None),
            ("BRAUN2020_CHECKMATE", "kidney (ccRCC)", "NIVOLUMAB", "NIVOLUMAB"),
            ("BRAUN2020_CHECKMATE", "kidney (ccRCC)", "EVEROLIMUS", "EVEROLIMUS"),
            ("GSE91061", "melanoma", "NIVOLUMAB", None),
            ("GSE78220", "melanoma", "PEMBROLIZUMAB", None)]


def scores(cid):
    expr, clin = CU.load_curated_trial(cid)
    sc = AE.pathway_scores(f"geo_{cid}", expr)
    clin = clin.loc[[s for s in clin.index if s in sc.columns]]
    Z = sc[clin.index].T
    Z = (Z - Z.mean()) / Z.std().replace(0, 1)
    return Z.fillna(0), clin


def response_cohort(cid, cancer, drug, arm, n_boot, rng):
    Z, clin = scores(cid)
    if arm is not None:
        clin = clin[clin["arm"] == arm]
    clin = clin[clin["responder"].notna()]
    y = clin["responder"].astype(float).values
    X = Z.loc[clin.index].values
    z = logit_score_z(X, y)
    hits = np.zeros(X.shape[1])
    used = 0
    for _ in range(n_boot):
        b = rng.integers(0, len(y), len(y))
        if 0 < y[b].sum() < len(b):
            hits[np.argsort(-np.abs(np.nan_to_num(logit_score_z(X[b], y[b]))))[:10]] += 1
            used += 1
    return pd.DataFrame({"cohort": cid, "cancer": cancer, "drug": drug, "n": len(y), "responders": int(y.sum()),
                         "pathway": Z.columns, "z": z, "q": bh(2 * norm.sf(np.abs(z))), "stability": hits / max(used, 1)})


def _cox_int(t, e, x, a, C):
    """z of the pathway x treatment interaction (and main effect) in a Cox model."""
    D = np.column_stack([x, a, x * a] + ([C] if C is not None and C.size else []))
    try:
        r = PHReg(t, D, status=e, ties="efron").fit(disp=0)
        return r.params[2] / r.bse[2], r.params[0] / r.bse[0]
    except Exception:
        return np.nan, np.nan


def interaction(Z, t, e, a, C, jobs=4):
    res = Parallel(n_jobs=jobs)(delayed(_cox_int)(t, e, Z[:, j], a, C) for j in range(Z.shape[1]))
    zi = -np.array([r[0] for r in res])  # sign: + = lower hazard with treatment as pathway increases
    zm = -np.array([r[1] for r in res])
    return zi, zm


def dummies(df, cols):
    out = []
    for c in cols:
        if df[c].dtype.kind in "fi":
            v = df[c].astype(float)
            out.append(v.fillna(v.median()).to_frame(c))
        else:
            out.append(pd.get_dummies(df[c].fillna("NA").astype(str), prefix=c, drop_first=True).astype(float))
    return pd.concat(out, axis=1).values if out else None


def predictive():
    rows = []
    # kidney: CM-025 randomised nivolumab vs everolimus, PFS
    Z, c = scores("BRAUN2020_CHECKMATE")
    c = c[(c["trial"] == "CM-025") & c["pfs_months"].notna() & c["pfs_event"].notna()]
    a = (c["arm"] == "NIVOLUMAB").astype(float).values
    zi, zm = interaction(Z.loc[c.index].values, c["pfs_months"].values, c["pfs_event"].values, a,
                         dummies(c, ["mskcc", "number_of_prior_therapies"]))
    rows.append(pd.DataFrame({"test": "kidney CM-025 nivolumab vs everolimus (PFS)", "n": len(c), "pathway": Z.columns,
                              "interaction_z": zi, "main_z": zm}))
    # glioma: CGGA TMZ vs none
    Z, c = scores("CGGA693")
    c = c[c["arm"].isin(["TMZ", "no TMZ"]) & c["os_months"].notna() & c["os_event"].notna()]
    for label, sub in (("all glioma", c), ("GBM (WHO IV)", c[c["grade"] == "WHO IV"])):
        a = (sub["arm"] == "TMZ").astype(float).values
        cov = ["grade", "idh_status", "codel_1p19q", "mgmt_status", "radio_status", "age"] if label == "all glioma" else \
            ["idh_status", "mgmt_status", "radio_status", "age"]
        zi, zm = interaction(Z.loc[sub.index].values, sub["os_months"].values, sub["os_event"].values, a, dummies(sub, cov))
        rows.append(pd.DataFrame({"test": f"CGGA {label} TMZ vs no TMZ (OS)", "n": len(sub), "pathway": Z.columns,
                                  "interaction_z": zi, "main_z": zm}))
    # GBM: GSE7696 RT+TMZ vs RT
    Z, c = scores("GSE7696")
    c = c[c["os_months"].notna()]
    a = (c["arm"] == "RT+TMZ").astype(float).values
    zi, zm = interaction(Z.loc[c.index].values, c["os_months"].values, c["os_event"].values, a, dummies(c, ["mgmt_status", "age"]))
    rows.append(pd.DataFrame({"test": "GSE7696 GBM RT+TMZ vs RT (OS)", "n": len(c), "pathway": Z.columns,
                              "interaction_z": zi, "main_z": zm}))
    df = pd.concat(rows, ignore_index=True)
    # GBM meta-analysis (pathways present in both)
    g1 = df[df.test == "CGGA GBM (WHO IV) TMZ vs no TMZ (OS)"].set_index("pathway")["interaction_z"]
    g2 = df[df.test == "GSE7696 GBM RT+TMZ vs RT (OS)"].set_index("pathway")["interaction_z"]
    common = g1.index.intersection(g2.index)
    meta = pd.DataFrame({"test": "GBM meta (CGGA WHO IV + GSE7696), Stouffer", "pathway": common,
                         "interaction_z": (g1[common].values + g2[common].values) / np.sqrt(2),
                         "agree_sign": np.sign(g1[common].values) == np.sign(g2[common].values)})
    meta["n"] = int(df[df.test.str.startswith("CGGA GBM")]["n"].iloc[0] + df[df.test.str.startswith("GSE7696")]["n"].iloc[0])
    df = pd.concat([df, meta], ignore_index=True)
    df["q"] = np.nan
    for t, g in df.groupby("test"):
        ok = g["interaction_z"].notna()
        df.loc[g.index[ok], "q"] = bh(2 * norm.sf(np.abs(g.loc[ok, "interaction_z"].values)))
    return df


def name(p):
    return p.replace("REACTOME_", "").replace("_", " ").lower()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=200)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(13)
    resp = pd.concat([response_cohort(*r, a.boot, rng) for r in RESPONSE], ignore_index=True)
    resp.to_csv(OUT / "response_pathways.tsv.gz", sep="\t", index=False)
    tops = []
    for (cid, d), g in resp.groupby(["cohort", "drug"], sort=False):
        tops.append(g.sort_values("z", ascending=False).head(5).assign(direction="response", rank=range(1, 6)))
        tops.append(g.sort_values("z").head(5).assign(direction="resistance", rank=range(1, 6)))
    top = pd.concat(tops, ignore_index=True)
    top.to_csv(OUT / "top10_cohorts.tsv", sep="\t", index=False)
    inter = predictive()
    inter.to_csv(OUT / "interaction.tsv.gz", sep="\t", index=False)
    ti = pd.concat([g.reindex(g["interaction_z"].abs().sort_values(ascending=False).index).head(10)
                    for _, g in inter.groupby("test", sort=False)])
    ti.to_csv(OUT / "top_interactions.tsv", sep="\t", index=False)

    L = ["# Pathway biomarkers in liver, kidney, glioma and melanoma (public treated cohorts)", "",
         "Response cohorts: logistic score-test z (+ = higher in responders); `*` q < 0.1; stab = 200-bootstrap "
         "top-10 frequency. Predictive tests: pathway x treatment interaction z in Cox models (+ = treatment "
         "benefit increases with the pathway). Exploratory.", ""]
    for (cid, d), g in top.groupby(["cohort", "drug"], sort=False):
        h = g.iloc[0]
        L += [f"## {h.cancer} - {d.title()} ({cid}; n = {h.n}, responders = {h.responders})", "",
              "| # | Response pathway | z | stab | Resistance pathway | z | stab |", "|---|---|---|---|---|---|---|"]
        rs, rz = g[g.direction == "response"].reset_index(), g[g.direction == "resistance"].reset_index()
        for i in range(5):
            r1, r2 = rs.iloc[i], rz.iloc[i]
            L.append(f"| {i + 1} | {name(r1.pathway)[:55]}{'*' if r1.q < 0.1 else ''} | {r1.z:+.1f} | {r1.stability:.2f} | "
                     f"{name(r2.pathway)[:55]}{'*' if r2.q < 0.1 else ''} | {r2.z:+.1f} | {r2.stability:.2f} |")
        L.append("")
    L += ["## Predictive (treatment x pathway interaction) tests", ""]
    for t, g in ti.groupby("test", sort=False):
        L += [f"**{t}** (n = {int(g['n'].iloc[0])}; min q = {inter[inter.test == t]['q'].min():.2f})", "",
              "| Pathway | interaction z | q |", "|---|---|---|"]
        L += [f"| {name(r.pathway)[:60]} | {r.interaction_z:+.2f} | {r.q:.2f} |" for r in g.itertuples()]
        L.append("")
    (OUT / "COHORT_BIOMARKERS.md").write_text("\n".join(L) + "\n")
    print(resp.groupby(["cohort", "drug"]).q.min().round(3).to_string())
    print(inter.groupby("test").q.min().round(3).to_string())


if __name__ == "__main__":
    main()
