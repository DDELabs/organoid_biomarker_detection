"""Integrated evidence pipeline: one table of pathway candidates per drug from every source.

Sources (z > 0 always means: higher pathway activity -> better response / more benefit from the drug)
  network    drug-target proximity on STRING (no direction; flag only)            [Kong et al. step 1]
  cells      cancer cell lines (CTRP, PRISM): pathway vs drug sensitivity, general sensitivity removed
  organoids  patient-derived organoids / PDX sets that tested the drug
  notes      TCGA patients given the drug: pathway vs recorded response
  survival   TCGA patients: pathway x treated interaction on progression-free interval
  external   independent trial cohorts: pathway vs response in patients given the drug
  ext_other  the same in patients of those cohorts NOT given the drug (specificity check)

Relaxed tiers (a source "supports" a pathway at |z| >= 1.64, same direction):
  tier 1  >= 2 patient sources (notes, survival, external) + >= 1 lab source (cells, organoids)
  tier 2  >= 2 patient sources
  tier 3  1 patient source at |z| >= 2.58 + >= 1 lab source
Chance expectation: the same rules after shuffling pathway labels independently in every source.
Usage: python robust/run_integrated.py [--null 200]
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import run_survlearn as RS  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.curated import PRECLIN, load_curated_preclinical, load_curated_trial, preclinical_catalog, trial_catalog  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from run_survlearn_external import cohort_scores  # noqa: E402

OUT = RESULTS / "integrated"
ALIAS = {"5-FLUOROURACIL": "FLUOROURACIL", "5-FU": "FLUOROURACIL", "SN-38": "IRINOTECAN"}
CELLS = ("ctrpv2_ccle", "prism_repurposing")
T = 1.64


def resid(A, C):
    """Residualise columns of A on columns of C (with intercept)."""
    C = np.column_stack([np.ones(len(A)), C]) if C is not None and C.size else np.ones((len(A), 1))
    return A - C @ np.linalg.lstsq(C, A, rcond=None)[0]


def corr_z(X, y, k=0):
    """z for the correlation of every column of X with y (rank-based); k = covariates already removed."""
    Xr = np.apply_along_axis(rankdata, 0, X)
    yr = rankdata(y)
    Xr, yr = Xr - Xr.mean(0), yr - yr.mean()
    r = (Xr * yr[:, None]).sum(0) / np.sqrt((Xr ** 2).sum(0) * (yr ** 2).sum() + 1e-12)
    return r * np.sqrt(max(len(y) - 3 - k, 1))


def stouffer(zs, ns):
    zs, w = np.array(zs), np.sqrt(np.array(ns, float))
    return (zs * w[:, None]).sum(0) / np.sqrt((w ** 2).sum())


def lab_scores(name):
    def compute():
        expr, _ = load_curated_preclinical(name)
        gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
        return scoring.rank_score(expr, {k: v for k, v in gs.items() if len(v) >= 5})
    return scoring.cached(f"preclin_{name}_reactome_rank", compute).T


def lab_evidence(drugs, paths):
    cells, orgs = {d: [] for d in drugs}, {d: [] for d in drugs}
    cat = preclinical_catalog()
    for name in cat["set"]:
        r = pd.read_csv(PRECLIN / name / "response.tsv", sep="\t")
        if not r["direction"].astype(str).str.contains("lower = more sensitive").all():
            continue                                             # keep one consistent direction
        r["drug"] = r["drug"].str.upper().map(lambda x: ALIAS.get(x, x))
        R = r.groupby(["sample", "drug"])["response"].median().unstack()
        try:
            sc = lab_scores(name)
        except Exception as exc:  # noqa: BLE001
            print(name, "skipped:", str(exc)[:60])
            continue
        sc = sc.reindex(columns=paths).dropna(axis=1, how="all").fillna(0)
        common = [s for s in R.index if s in sc.index]
        R, sc = R.loc[common], sc.loc[common]
        gen = ((R - R.mean()) / R.std()).mean(axis=1)            # general sensitivity across all drugs
        lin = r.drop_duplicates("sample").set_index("sample")["lineage"].reindex(common) if "lineage" in r else None
        for d in drugs:
            if d not in R or R[d].notna().sum() < (40 if name in CELLS else 8):
                continue
            ok = R[d].notna() & gen.notna()
            C = [gen[ok].values[:, None]]
            if lin is not None and name in CELLS:
                L = pd.get_dummies(lin[ok])
                C.append(L.loc[:, L.sum() >= 5].values.astype(float))
            C = np.hstack(C)
            z = -corr_z(resid(sc[ok].values, C), resid(R.loc[ok, [d]].values, C)[:, 0], C.shape[1])
            (cells if name in CELLS else orgs)[d].append((pd.Series(z, sc.columns).reindex(paths), int(ok.sum()), name))
    return cells, orgs


def score_test(D, Zm):
    """One-step Cox estimates (beta, information) for every pathway given the base covariates."""
    m = SL.StratCox(D["time"], D["event"], D["project"])
    B = D[RS.BASE].values
    b = m.fit(B, 1e-3)
    eta = (B @ b)[m.order]
    X = Zm[m.order]
    mx = np.maximum.reduceat(eta, m.start)[m.sid]
    w = np.exp(eta - mx)
    S0 = m._cum(w)[m.last]
    a = np.zeros(len(w))
    np.add.at(a, m.last, m.e / S0)
    U = X.T @ (m.e - w * m._rcum(a))
    S1, S2 = np.empty_like(X), np.empty_like(X)
    ends = np.r_[m.start[1:], len(w)]
    for s, e in zip(m.start, ends):
        S1[s:e] = np.cumsum(w[s:e, None] * X[s:e], axis=0)
        S2[s:e] = np.cumsum(w[s:e, None] * X[s:e] ** 2, axis=0)
    S1, S2 = S1[m.last], S2[m.last]
    info = (m.e[:, None] * (S2 / S0[:, None] - (S1 / S0[:, None]) ** 2)).sum(0)
    return U / info, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=200)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    Z, cdr, long = RS.load()
    paths = list(Z.columns)
    fold = SL.folds_by_patient(cdr.index, 5, seed=0)
    cdr["g"] = SL.general_risk(Z, cdr, "time", "event", fold)
    cdr["g"] /= cdr["g"].std()
    n = long.drop_duplicates(["patient", "drug"]).groupby("drug").size()
    drugs = list(n[n >= 30].sort_values(ascending=False).index)
    cells, orgs = lab_evidence(drugs, paths)
    trials = {}
    cat = trial_catalog()
    for cid in cat.loc[cat["usable"].astype(str).str.lower().eq("true"), "cohort_id"]:
        try:
            _, clin = load_curated_trial(cid, canonical=False)
            if "responder" in clin and "drugs" in clin:
                trials[cid] = (cohort_scores(cid).reindex(columns=paths).fillna(0), clin)
        except Exception:  # noqa: BLE001
            pass
    rows = []
    for d in drugs:
        ev = pd.DataFrame(index=paths)
        z = RS.proximity(d)
        ev["network_z"] = z.reindex(paths) if z is not None else np.nan
        meta = {}
        for key, src in (("cells", cells), ("organoids", orgs)):
            if src[d]:
                ev[key] = stouffer([s[0].fillna(0).values for s in src[d]], [s[1] for s in src[d]]) if key == "organoids" else \
                    np.mean([s[0].fillna(0).values for s in src[d]], axis=0)
                meta[key] = ";".join(f"{s[2]}({s[1]})" for s in src[d])
        rec = long[long["drug"] == d]
        lab = rec[rec["responder"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")
        if len(lab) >= 30 and lab["responder"].value_counts().min() >= 8:
            P = pd.get_dummies(cdr.loc[lab.index, "project"])
            P = P.loc[:, P.sum() >= 3].values.astype(float)
            ev["notes"] = corr_z(resid(Z.loc[lab.index].values, P), resid(lab[["responder"]].astype(float).values, P)[:, 0], P.shape[1])
            meta["notes"] = f"n={len(lab)}"
        first = set(rec.loc[rec["before_event"], "patient"])
        per = cdr.loc[sorted(first), "project"].value_counts()
        D = cdr[cdr["project"].isin(per[per >= 10].index)].copy()
        D["treated"] = D.index.isin(first).astype(float)
        D["other_systemic"] = D.index.isin(set(long.loc[long["before_event"] & (long["drug"] != d), "patient"])).astype(float)
        et, eu = D.loc[D.treated == 1, "event"].sum(), D.loc[D.treated == 0, "event"].sum()
        if D["treated"].sum() >= 30 and et >= 15 and eu >= 15:
            bt, it = score_test(D[D.treated == 1], Z.loc[D.index[D.treated == 1]].values)
            bu, iu = score_test(D[D.treated == 0], Z.loc[D.index[D.treated == 0]].values)
            ev["survival"] = -(bt - bu) / np.sqrt(1 / it + 1 / iu)
            ev["survival_treated_only"] = -bt * np.sqrt(it)
            meta["survival"] = f"events {int(et)}/{int(eu)}"
        for key, want in (("external", True), ("ext_other", False)):
            zs, ns, names = [], [], []
            for cid, (sc, clin) in trials.items():
                got = clin["drugs"].fillna("").astype(str).str.upper().str.split(";").map(lambda x: d in x)
                if not got.any():
                    continue
                sel = [s for s in clin.index[(got == want) & clin["responder"].notna()] if s in sc.index]
                y = clin.loc[sel, "responder"].astype(float)
                if len(sel) < 15 or y.value_counts().min() < 4:
                    continue
                zs.append(corr_z(sc.loc[sel].values, y.values))
                ns.append(len(sel))
                names.append(f"{cid}({len(sel)})")
            if zs:
                ev[key] = stouffer(zs, ns)
                meta[key] = ";".join(names)
        ev = ev.assign(drug=d).rename_axis("pathway").reset_index()
        for k, v in meta.items():
            ev[f"{k}_data"] = v
        rows.append(ev)
    E = pd.concat(rows, ignore_index=True)
    PAT, LAB = ["notes", "survival", "external"], ["cells", "organoids"]
    for c in PAT + LAB + ["ext_other"]:
        if c not in E:
            E[c] = np.nan

    def tiers(F):
        pz, lz = F[PAT].values, F[LAB].values
        out = np.zeros(len(F), int)
        direction = np.zeros(len(F), int)
        for sgn in (1, -1):
            np_ = np.nansum(sgn * pz >= T, axis=1)
            nl = np.nansum(sgn * lz >= T, axis=1)
            strong = np.nansum(sgn * pz >= 2.58, axis=1)
            t = np.where((np_ >= 2) & (nl >= 1), 1, np.where(np_ >= 2, 2, np.where((strong >= 1) & (nl >= 1), 3, 0)))
            better = (t > 0) & ((out == 0) | (t < out))
            out[better], direction[better] = t[better], sgn
        return out, direction

    E["tier"], E["direction"] = tiers(E)
    E["proximal"] = E["network_z"] <= -1.2816
    E["also_in_untreated"] = (E["direction"] * E["ext_other"] >= T)
    rng = np.random.default_rng(11)
    null = []
    for _ in range(a.null):
        F = E[["drug"] + PAT + LAB].copy()
        for c in PAT + LAB:
            F[c] = F.groupby("drug")[c].transform(lambda v: v.values[rng.permutation(len(v))])
        t, _ = tiers(F)
        null.append([(t == k).sum() for k in (1, 2, 3)])
    null = np.array(null)
    E.to_csv(OUT / "evidence.tsv.gz", sep="\t", index=False)
    C = E[E["tier"] > 0].sort_values(["tier", "drug"])
    C.to_csv(OUT / "candidates.tsv", sep="\t", index=False)
    print("sources available per drug:")
    print(E.groupby("drug")[PAT + LAB + ["ext_other"]].agg(lambda v: int(v.notna().any())).loc[drugs].T.to_string())
    for k in (1, 2, 3):
        obs = int((E.tier == k).sum())
        print(f"tier {k}: observed {obs}; expected by chance {null[:, k-1].mean():.1f} (95% range {np.percentile(null[:, k-1], 2.5):.0f}-{np.percentile(null[:, k-1], 97.5):.0f}); "
              f"proximal among observed: {int(((E.tier == k) & E.proximal).sum())}")
    pd.DataFrame(null, columns=["tier1", "tier2", "tier3"]).to_csv(OUT / "chance_counts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
