"""B1: frozen organoid-only signatures applied to external patient cohorts (no retraining on patients).

Signature = the organoid-only agents' saved Ridge coefficients: rank-1 pathway alone ("rank1") and the top-10
pathways by |coefficient| ("top10", the multi-pathway model). Patient pathway values = within-sample rank scores
(robust/obd/scoring.rank_score, platform-independent; training table + STORM scored by nar/03 loader), each pathway
standardised within cohort; score = sum(coef x z); score standardised within cohort. Higher = predicted resistant
(pre-specified direction: HR or OR for non-response > 1).
Random null: rank1 -> every one of the 673 pathways used alone (exact null); top10 -> 1,000 random 10-pathway sets
with N(0,1) weights. Null p = share of random signatures with |z| >= observed |z| (z of the primary statistic: the
interaction when a comparator arm exists, otherwise the treated-arm effect).
"""
import importlib.util
import sys

import numpy as np
import pandas as pd

from common import DATA, REPO, SEED, TAB, cox, logit, logrank

spec = importlib.util.spec_from_file_location("j_common", REPO / "nar/03_organoid_plus_patient/code/common.py")
J = importlib.util.module_from_spec(spec)
sys.modules["j_common"] = J
spec.loader.exec_module(J)

N_NULL = 1000
NONE = r"^NONE$"
FU = r"FLUOROURACIL|CHEMOTHERAPY"

# signature groups: analysis ids whose signatures are applied to the cohort tests of that group
GROUPS = {
    "CRC_5FU": ["COAD_vandewetering__COAD_authors__FLUOROURACIL", "COAD_vandewetering__COAD__FLUOROURACIL",
                "PDX_gao2015_CRC__COAD__FLUOROURACIL"],
    "CRC_OXALIPLATIN": ["COAD_vandewetering__COAD_authors__OXALIPLATIN", "COAD_vandewetering__COAD__OXALIPLATIN"],
    "CRC_IRINOTECAN": ["COAD_vandewetering__COAD_authors__IRINOTECAN", "COAD_vandewetering__COAD__IRINOTECAN"],
    "BLCA_CISPLATIN": ["BLCA_lee2018_paper__BLCA__CISPLATIN", "BLCA_lee2018_current__BLCA__CISPLATIN"],
    "OV_PLATINUM": ["OV_vias2023__OV__CARBOPLATIN", "OV_vias2023__OV__CISPLATIN"],
    "OV_TAXANE": ["OV_vias2023__OV__PACLITAXEL", "OV_vias2023__OV__DOCETAXEL"],
    "LIHC_SORAFENIB": ["LIHC_ji__LIHC__SORAFENIB", "LIHC_broutier__LIHC__SORAFENIB"],
    "LUNG_PLATINUM": ["LUAD_cellline_substitute__LUAD__CISPLATIN", "LUAD_cellline_substitute__LUAD__CARBOPLATIN"],
}

# tests: (group, cohort, label, treated(arm/drugs regex on 'arm|drugs'), comparator regex or None, endpoint, filter)
T = [
    ("CRC_5FU", "GSE39582", "5-FU-based vs no chemo (stage II-III)", r"\|.*(" + FU + ")", r"\|NONE$", "rfs", "nost4"),
    ("CRC_5FU", "GSE39582", "5-FU-based vs no chemo (stage II-III)", r"\|.*(" + FU + ")", r"\|NONE$", "os", "nost4"),
    ("CRC_5FU", "GSE39582", "5FU alone vs no chemo", r"^5FU\|", r"\|NONE$", "rfs", "nost4"),
    ("CRC_5FU", "GSE39582", "FUFOL vs no chemo", r"^FUFOL\|", r"\|NONE$", "rfs", "nost4"),
    ("CRC_5FU", "GSE39582", "FOLFOX vs no chemo", r"^FOLFOX\|", r"\|NONE$", "rfs", "nost4"),
    ("CRC_5FU", "GSE103479", "5-FU adjuvant vs surgery only", r"FLUOROURACIL", r"\|NONE$", "os", "nost4"),
    ("CRC_5FU", "GSE103479", "5-FU adjuvant vs surgery only", r"FLUOROURACIL", r"\|NONE$", "pfs", "nost4"),
    ("CRC_5FU", "GSE14333", "adjuvant chemo vs none (stage I-III)", r"CHEMOTHERAPY", r"\|NONE$", "dfs", "nost4"),
    ("CRC_5FU", "GSE72970", "all 5-FU-based (metastatic)", r"FLUOROURACIL|CAPECITABINE", None, "pfs", None),
    ("CRC_5FU", "GSE72970", "all 5-FU-based (metastatic)", r"FLUOROURACIL|CAPECITABINE", None, "resp", None),
    ("CRC_5FU", "GSE87211", "5-FU chemoradiotherapy (rectal)", r"FLUOROURACIL", None, "dfs", None),
    ("CRC_5FU", "GSE87211", "5-FU chemoradiotherapy (rectal)", r"FLUOROURACIL", None, "resp", None),
    ("CRC_5FU", "GSE28702", "mFOLFOX6", r"FLUOROURACIL", None, "resp", None),
    ("CRC_5FU", "GSE19860", "mFOLFOX6", r"FLUOROURACIL", None, "resp", None),
    ("CRC_5FU", "GSE62080", "FOLFIRI", r"FLUOROURACIL", None, "resp", None),
    ("CRC_5FU", "GSE104645", "5-FU/tegafur-based (metastatic)", r"FLUOROURACIL|TEGAFUR", None, "resp", None),
    ("CRC_OXALIPLATIN", "GSE72970", "FOLFOX vs FOLFIRI", r"OXALIPLATIN", r"^(?!.*OXALIPLATIN).*IRINOTECAN", "pfs", None),
    ("CRC_OXALIPLATIN", "GSE72970", "FOLFOX vs FOLFIRI", r"OXALIPLATIN", r"^(?!.*OXALIPLATIN).*IRINOTECAN", "resp", None),
    ("CRC_OXALIPLATIN", "GSE39582", "FOLFOX vs no chemo", r"OXALIPLATIN", r"\|NONE$", "rfs", "nost4"),
    ("CRC_OXALIPLATIN", "GSE87211", "5-FU+oxaliplatin+RT vs 5-FU+RT", r"OXALIPLATIN", r"\|FLUOROURACIL;RADIOTHERAPY$", "dfs", None),
    ("CRC_OXALIPLATIN", "GSE87211", "5-FU+oxaliplatin+RT vs 5-FU+RT", r"OXALIPLATIN", r"\|FLUOROURACIL;RADIOTHERAPY$", "resp", None),
    ("CRC_OXALIPLATIN", "GSE28702", "mFOLFOX6", r"OXALIPLATIN", None, "resp", None),
    ("CRC_OXALIPLATIN", "GSE19860", "mFOLFOX6", r"OXALIPLATIN", None, "resp", None),
    ("CRC_OXALIPLATIN", "GSE104645", "oxaliplatin vs irinotecan regimens", r"OXALIPLATIN", r"^(?!.*OXALIPLATIN).*IRINOTECAN", "resp", None),
    ("CRC_IRINOTECAN", "GSE72970", "FOLFIRI vs FOLFOX", r"^(?!.*OXALIPLATIN).*IRINOTECAN", r"^(?!.*IRINOTECAN).*OXALIPLATIN", "pfs", None),
    ("CRC_IRINOTECAN", "GSE72970", "FOLFIRI vs FOLFOX", r"^(?!.*OXALIPLATIN).*IRINOTECAN", r"^(?!.*IRINOTECAN).*OXALIPLATIN", "resp", None),
    ("CRC_IRINOTECAN", "GSE62080", "FOLFIRI", r"IRINOTECAN", None, "resp", None),
    ("CRC_IRINOTECAN", "GSE104645", "irinotecan vs oxaliplatin regimens", r"^(?!.*OXALIPLATIN).*IRINOTECAN", r"^(?!.*IRINOTECAN).*OXALIPLATIN", "resp", None),
    ("CRC_IRINOTECAN", "GSE39582", "FOLFIRI vs no chemo", r"IRINOTECAN", r"\|NONE$", "rfs", "nost4"),
    ("BLCA_CISPLATIN", "GSE169455", "neoadjuvant/induction cisplatin", r"CISPLATIN", None, "resp", None),
    ("BLCA_CISPLATIN", "GSE52219", "neoadjuvant MVAC", r"CISPLATIN", None, "resp", None),
    ("BLCA_CISPLATIN", "GSE48276", "neoadjuvant MVAC vs none", r"CISPLATIN", r"\|NONE$", "os", None),
    ("BLCA_CISPLATIN", "GSE70691", "adjuvant / neoadjuvant MVAC", r"CISPLATIN", None, "os", None),
    ("BLCA_CISPLATIN", "GSE31684", "chemo after cystectomy vs none", r"CHEMOTHERAPY", r"\|NONE$", "os", None),
    ("OV_PLATINUM", "GSE9891", "platinum-based vs no platinum", r"PLATINUM", r"^no platinum\|", "os", None),
    ("OV_PLATINUM", "GSE9891", "platinum-based (all)", r"PLATINUM", None, "rfs", None),
    ("OV_PLATINUM", "GSE63885", "platinum-based", r"PLATINUM", None, "os", None),
    ("OV_PLATINUM", "GSE63885", "platinum-based", r"PLATINUM", None, "resp", None),
    ("OV_PLATINUM", "GSE51373", "carboplatin-paclitaxel", r"CARBOPLATIN", None, "resp", None),
    ("OV_PLATINUM", "GSE15622", "carboplatin vs paclitaxel arm", r"^carboplatin\|", r"^paclitaxel\|", "resp", None),
    ("OV_TAXANE", "GSE9891", "platinum+taxane vs platinum alone", r"PACLITAXEL", r"\|PLATINUM$", "os", None),
    ("OV_TAXANE", "GSE9891", "platinum+taxane vs platinum alone", r"PACLITAXEL", r"\|PLATINUM$", "rfs", None),
    ("OV_TAXANE", "GSE63885", "taxane/platinum vs cyclophosphamide/platinum", r"PACLITAXEL", r"CYCLOPHOSPHAMIDE", "os", None),
    ("OV_TAXANE", "GSE63885", "taxane/platinum vs cyclophosphamide/platinum", r"PACLITAXEL", r"CYCLOPHOSPHAMIDE", "resp", None),
    ("OV_TAXANE", "GSE15622", "paclitaxel vs carboplatin arm", r"^paclitaxel\|", r"^carboplatin\|", "resp", None),
    ("OV_TAXANE", "GSE51373", "carboplatin-paclitaxel", r"PACLITAXEL", None, "resp", None),
    ("LIHC_SORAFENIB", "GSE109211", "sorafenib vs placebo (randomised, STORM)", r"SORAFENIB", r"PLACEBO", "resp", None),
    ("LUNG_PLATINUM", "GSE42127", "adjuvant platinum vs observation", r"PLATINUM", r"\|NONE$", "os", None),
]


def signatures():
    s = pd.concat([pd.read_csv(DATA / f"signatures_{k}.tsv", sep="\t") for k in ("gi", "other")])
    out, seen = [], {}
    for g, ids in GROUPS.items():
        for aid in ids:
            d = s[s.analysis_id == aid].sort_values("rank")
            for kind, dd in (("rank1", d.head(1)), ("top10", d.head(10))):
                key = tuple(zip(dd.pathway, dd.coef.round(10)))
                dup = seen.get((g, key))
                if dup:
                    print("identical signature, skipped:", aid, kind, "=", dup)
                    continue
                seen[(g, key)] = f"{aid}|{kind}"
                out.append({"group": g, "analysis_id": aid, "kind": kind, "sig_id": f"{aid}|{kind}",
                            "pathways": list(dd.pathway), "coefs": dd.coef.values.astype(float)})
    return out


def zcoh(S):
    return (S - S.mean()) / S.std().replace(0, 1)


def test_one(d, score, endpoint, has_comp):
    """d: cohort rows with 'arm1' (1 treated, 0 comparator); score: standardised score (Series aligned)."""
    r = {}
    t = d[d.arm1 == 1]
    st = score.loc[t.index].values
    if endpoint == "resp":
        y = 1 - t["responder"].astype(float).values
        m = logit(y, st, ["score"])
        r.update(n=len(t), events=int(y.sum()), effect="OR non-response per SD")
        if m is not None:
            r.update(est=m.loc["score", "OR"], lo=m.loc["score", "lo"], hi=m.loc["score", "hi"], p=m.loc["score", "p"],
                     z=m.loc["score", "coef"] / m.loc["score", "se"])
        if has_comp:
            c = d[d.arm1 == 0]
            yc = 1 - c["responder"].astype(float).values
            mc = logit(yc, score.loc[c.index].values, ["score"])
            if mc is not None:
                r.update(comp_est=mc.loc["score", "OR"], comp_lo=mc.loc["score", "lo"], comp_hi=mc.loc["score", "hi"],
                         comp_p=mc.loc["score", "p"])
            yy = 1 - d["responder"].astype(float).values
            sc = score.loc[d.index].values
            mi = logit(yy, np.column_stack([sc, d.arm1.values, sc * d.arm1.values]), ["score", "arm", "inter"])
            r.update(comp_n=len(c), comp_events=int(yc.sum()))
            if mi is not None:
                r.update(int_est=mi.loc["inter", "OR"], int_lo=mi.loc["inter", "lo"], int_hi=mi.loc["inter", "hi"],
                         int_p=mi.loc["inter", "p"], int_z=mi.loc["inter", "coef"] / mi.loc["inter", "se"])
        return r
    Tm, Ev = t["T"].values, t["E"].values
    m = cox(Tm, Ev, st, ["score"])
    r.update(n=len(t), events=int(Ev.sum()), effect=f"HR ({endpoint.upper()}) per SD",
             logrank_p_median=logrank(Tm, Ev, (st > np.median(st)).astype(int)))
    if m is not None:
        r.update(est=m.loc["score", "HR"], lo=m.loc["score", "lo"], hi=m.loc["score", "hi"], p=m.loc["score", "p"],
                 z=m.loc["score", "coef"] / m.loc["score", "se"])
    if has_comp:
        c = d[d.arm1 == 0]
        mc = cox(c["T"].values, c["E"].values, score.loc[c.index].values, ["score"])
        r.update(comp_n=len(c), comp_events=int(c["E"].sum()))
        if mc is not None:
            r.update(comp_est=mc.loc["score", "HR"], comp_lo=mc.loc["score", "lo"], comp_hi=mc.loc["score", "hi"],
                     comp_p=mc.loc["score", "p"])
        sc = score.loc[d.index].values
        mi = cox(d["T"].values, d["E"].values, np.column_stack([sc, d.arm1.values, sc * d.arm1.values]),
                 ["score", "arm", "inter"], check=3)
        if mi is not None:
            r.update(int_est=mi.loc["inter", "HR"], int_lo=mi.loc["inter", "lo"], int_hi=mi.loc["inter", "hi"],
                     int_p=mi.loc["inter", "p"], int_z=mi.loc["inter", "coef"] / mi.loc["inter", "se"])
    return r


def cohort_rows(P, cohort, treat, comp, endpoint, filt):
    d = P[P.cohort == cohort].copy()
    key = d["arm"].fillna("NA").astype(str) + "|" + d["drugs"].fillna("NA").astype(str)
    tr = key.str.contains(treat, regex=True)
    cp = key.str.contains(comp, regex=True) & ~tr if comp else pd.Series(False, index=d.index)
    d = d[tr | cp].copy()
    d["arm1"] = tr[d.index].astype(int)
    if filt == "nost4":
        d = d[d["stage"].map(J.stage_num) != 4]
    if endpoint == "resp":
        d = d[d["responder"].notna()]
    else:
        d["T"] = pd.to_numeric(d[f"{endpoint}_months"], errors="coerce")
        d["E"] = pd.to_numeric(d[f"{endpoint}_event"], errors="coerce")
        d = d.dropna(subset=["T", "E"])
        d = d[d["T"] >= 0]
    return d


def main():
    P, S = J.patients()
    sigs = signatures()
    pd.DataFrame([{**{k: v for k, v in s.items() if k not in ("pathways", "coefs")},
                   "pathways": ";".join(s["pathways"]), "coefs": ";".join(f"{c:.5g}" for c in s["coefs"])}
                  for s in sigs]).to_csv(TAB / "B1_frozen_signatures.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED)
    allp = np.array(S.columns)
    rows, scores_out = [], []
    Zc = {}
    for (g, cohort, lab, treat, comp, ep, filt) in T:
        if cohort not in Zc:
            idx = P.index[P.cohort == cohort]
            Zc[cohort] = zcoh(S.loc[idx].astype(float))
        Z = Zc[cohort]
        d = cohort_rows(P, cohort, treat, comp, ep, filt)
        has_comp = comp is not None and (d.arm1 == 0).sum() >= 5
        if (d.arm1 == 1).sum() < 8:
            print("skip (too few treated)", cohort, lab, ep)
            continue
        # null distributions for this test (shared by all signatures of the same size)
        nulls = {}
        for k in (1, 10):
            zs = []
            if k == 1:
                draws = [(np.array([j]), np.array([1.0])) for j in range(len(allp))]
            else:
                draws = [(rng.choice(len(allp), 10, replace=False), rng.standard_normal(10)) for _ in range(N_NULL)]
            for ix, w in draws:
                sc = Z.iloc[:, ix].values @ w
                sc = pd.Series((sc - sc.mean()) / sc.std(), index=Z.index)
                rr = test_one(d, sc, ep, has_comp)
                zs.append(rr.get("int_z" if has_comp else "z", np.nan))
            nulls[k] = np.array(zs, float)
        for s in sigs:
            if s["group"] != g:
                continue
            sc = Z[s["pathways"]].values @ s["coefs"]
            sc = pd.Series((sc - sc.mean()) / sc.std(), index=Z.index)
            r = test_one(d, sc, ep, has_comp)
            zobs = r.get("int_z" if has_comp else "z", np.nan)
            nz = nulls[1 if s["kind"] == "rank1" else 10]
            nz = nz[np.isfinite(nz)]
            r.update(group=g, cohort=cohort, comparison=lab, endpoint=ep.upper(), sig_id=s["sig_id"],
                     analysis_id=s["analysis_id"], kind=s["kind"], has_comparator=has_comp,
                     primary_stat="interaction" if has_comp else "treated-arm effect",
                     null_n=len(nz),
                     null_p_two_sided=float((np.sum(np.abs(nz) >= abs(zobs)) + 1) / (len(nz) + 1)) if np.isfinite(zobs) else np.nan,
                     null_p_directional=float((np.sum(nz >= zobs) + 1) / (len(nz) + 1)) if np.isfinite(zobs) else np.nan)
            rows.append(r)
            if s["kind"] == "rank1":
                for i in d.index:
                    scores_out.append({"cohort": cohort, "comparison": lab, "endpoint": ep, "sig_id": s["sig_id"],
                                       "sample": i, "arm1": d.at[i, "arm1"], "score": sc[i],
                                       "T": d.at[i, "T"] if ep != "resp" else np.nan,
                                       "E": d.at[i, "E"] if ep != "resp" else np.nan,
                                       "responder": d.at[i, "responder"]})
        print(cohort, lab, ep, "n treated", int((d.arm1 == 1).sum()), "comp", int((d.arm1 == 0).sum()), flush=True)
    R = pd.DataFrame(rows)
    R["direction"] = np.where(R["est"] > 1, "expected", "opposite")
    R.to_csv(TAB / "B1_external_results.tsv", sep="\t", index=False)
    pd.DataFrame(scores_out).to_csv(DATA / "b1_patient_scores.tsv.gz", sep="\t", index=False)
    print(R[["sig_id", "cohort", "endpoint", "n", "events", "est", "p", "int_est", "int_p", "null_p_two_sided"]].to_string())


if __name__ == "__main__":
    main()
