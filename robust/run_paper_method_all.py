"""The original paper's method (Kong et al. 2020), unchanged, on every organoid-drug-cancer set we can test.

Steps exactly as python/multiple_pathway_prediction.py: ssGSEA pathway scores -> proximal pathways for the drug
-> Ridge / linear SVR / OLS trained on organoid drug response -> top-k pathways (k = 2..10) -> predicted response
in TCGA patients given the drug -> median split -> log-rank test on overall survival (27 tests per set).
Added only: (a) direction of the difference, (b) the same procedure with organoid responses shuffled (500x).
Usage: python robust/run_paper_method_all.py [--perm 500]
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from statsmodels.duration.survfunc import survdiff

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_survlearn as RS  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd.curated import PRECLIN, load_curated_preclinical  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "paper_method_all"
SETS = [  # organoid set, TCGA project, patient drug, organoid drug column, response column
    ("colorectal_vandewetering2015", "COAD", "FLUOROURACIL", "5-FLUOROURACIL", "ln_ic50"),
    ("colorectal_vandewetering2015", "COAD", "OXALIPLATIN", "OXALIPLATIN", "ln_ic50"),
    ("colorectal_vandewetering2015", "COAD", "IRINOTECAN", "CAMPTOTHECIN", "ln_ic50"),
    ("colorectal_vandewetering2015", "COAD", "CETUXIMAB", "CETUXIMAB", "ln_ic50"),
    ("colorectal_vandewetering2015", "READ", "FLUOROURACIL", "5-FLUOROURACIL", "ln_ic50"),
    ("colorectal_vandewetering2015", "READ", "OXALIPLATIN", "OXALIPLATIN", "ln_ic50"),
    ("bladder_lee2018", "BLCA", "CISPLATIN", "CISPLATIN", "response"),
    ("bladder_lee2018", "BLCA", "GEMCITABINE", "GEMCITABINE", "response"),
    ("bladder_lee2018", "BLCA", "DOXORUBICIN", "DOXORUBICIN", "response"),
]


def fit(ml, X, y):
    if ml == "Ridge":
        return RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(X, y).coef_
    if ml == "SVR":
        return SVR(kernel="linear").fit(X, y).coef_[0]
    return LinearRegression().fit(X, y).coef_


def pipeline(X, y, PT, T, E):
    out = {}
    for ml in ("Ridge", "SVR", "LinearRegression"):
        coef = fit(ml, X, y)
        order = np.argsort(-np.abs(coef))
        for k in range(2, min(10, X.shape[1]) + 1):
            idx = order[:k]
            pred = PT[:, idx] @ coef[idx]
            g = (pred > np.median(pred)).astype(int)            # 1 = predicted non-responder
            if g.sum() in (0, len(g)):
                out[(ml, k)] = (1.0, 0.0)
                continue
            _, p = survdiff(T, E, g)
            rate = lambda m: E[m].sum() / max(T[m].sum(), 1e-9)  # noqa: E731  deaths per patient-month
            out[(ml, k)] = (p, np.log(rate(g == 1) / max(rate(g == 0), 1e-9)))
    return out


def scores(tag, expr, gs):
    sc = scoring.cached(tag, scoring.ssgsea, expr, gs)
    return sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=500)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    paper = N.precomputed_proximity()
    long, cdr = R.load_response_labels(), R.load_cdr()
    all_gs = gene_sets(("REACTOME",))
    rows, detail = [], []
    cache = {}
    for oset, proj, drug, ocol, rcol in SETS:
        key = (oset, proj)
        if key not in cache:
            oexpr, _ = load_curated_preclinical(oset)
            texpr = C.tcga_star_fpkm_uq(proj)
            genes = oexpr.index.intersection(texpr.index)        # shared gene universe, as in run_ssGSEA.py
            gs = {k: [g for g in v if g in genes] for k, v in all_gs.items()}
            gs = {k: v for k, v in gs.items() if len(v) >= 5}
            cache[key] = (scores(f"paper_{oset}_{proj}_organoid_ssgsea", oexpr.loc[genes], gs),
                          scores(f"paper_{oset}_{proj}_tcga_ssgsea", texpr.loc[genes], gs))
        osc, tsc = cache[key]
        r = pd.read_csv(PRECLIN / oset / "response.tsv", sep="\t")
        r = r[r["drug"].str.upper() == ocol].groupby("sample")[rcol].median().dropna()
        samples = [s for s in r.index if s in osc.columns]
        if drug in paper.columns:
            z, src = paper[drug].dropna(), "paper table"
        else:
            z, src = RS.proximity(drug), "STRING v12 (drug not in paper table)"
        feats = [p for p in z.index[z <= -1.2816] if p in osc.index and p in tsc.index]
        treated = set(long.loc[(long["project"] == proj) & (long["drug"] == drug), "patient"])
        pats = [p for p in tsc.columns if p in cdr.index and pd.notna(cdr.loc[p, "OS.time"]) and pd.notna(cdr.loc[p, "OS"])]
        P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[feats, pats].T.values), index=pats, columns=feats)
        tp = [p for p in pats if p in treated]
        T, E = cdr.loc[tp, "OS.time"].values / (365 / 12), cdr.loc[tp, "OS"].values.astype(int)
        X = StandardScaler().fit_transform(osc.loc[feats, samples].T.values)
        y = r.loc[samples].values
        obs = pipeline(X, y, P.loc[tp].values, T, E)
        ps = np.array([v[0] for v in obs.values()])
        best = min(obs, key=lambda k: obs[k][0])
        rng = np.random.default_rng(0)
        mins = np.array([min(v[0] for v in pipeline(X, rng.permutation(y), P.loc[tp].values, T, E).values()) for _ in range(a.perm)])
        right = sum(1 for v in obs.values() if v[1] > 0)
        rows.append({"cancer": proj, "drug": drug, "organoids": len(samples), "proximity": src, "pathways": len(feats),
                     "patients": len(tp), "deaths": int(E.sum()), "tests": len(ps), "best_p": ps.min(), "best_setting": f"{best[0]} k={best[1]}",
                     "best_direction_correct": bool(obs[best][1] > 0), "n_p_below_0.05": int((ps < 0.05).sum()),
                     "n_p_below_0.05_correct_direction": int(sum(1 for v in obs.values() if v[0] < 0.05 and v[1] > 0)),
                     "tests_correct_direction": right, "shuffled_as_good": (np.sum(mins <= ps.min()) + 1) / (a.perm + 1)})
        for (ml, k), (p, d) in obs.items():
            detail.append({"cancer": proj, "drug": drug, "model": ml, "k": k, "p": p, "log_rate_ratio": d})
        print(rows[-1], flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(detail).to_csv(OUT / "all_tests.tsv", sep="\t", index=False)
    print(S.drop(columns=["proximity", "best_setting"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
