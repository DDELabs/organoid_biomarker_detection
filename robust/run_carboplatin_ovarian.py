"""Frozen test of the breast-trained carboplatin score outside breast cancer.

The score is trained only on the carboplatin arms of BrighTNess and I-SPY2 (response at surgery), using the network
shortlist, then frozen and applied to platinum-treated patients elsewhere: TCGA ovarian (survival), other ovarian
cohorts (response, survival) and other TCGA cancers treated with carboplatin. A score that works should give
longer survival / more responses when it is high. Null: 200 random shortlists of the same size, trained the same way.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
import run_carboplatin as RC  # noqa: E402
import run_survlearn as RS  # noqa: E402
import run_train_combined as TC  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import network as N  # noqa: E402
from obd import survlearn as SL  # noqa: E402

OUT = RESULTS / "carboplatin"


def targets(P):
    plat = P["drugs"].fillna("").str.contains("CARBOPLATIN|PLATINUM")
    t = {"TCGA ovarian, carboplatin": P.index[(P["cohort"] == "TCGA-OV") & P["drugs"].fillna("").str.contains("CARBOPLATIN")]}
    for c in ("GSE9891", "GSE63885", "GSE51373", "GSE15622"):
        t[f"{c} ovarian, platinum"] = P.index[(P["cohort"] == c) & plat]
    for c in ("TCGA-UCEC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-HNSC", "TCGA-UCS"):
        t[f"{c}, carboplatin"] = P.index[(P["cohort"] == c) & P["drugs"].fillna("").str.contains("CARBOPLATIN")]
    t["TCGA ovarian, NO carboplatin (comparison)"] = P.index[(P["cohort"] == "TCGA-OV") & ~P["drugs"].fillna("").str.contains("CARBOPLATIN")]
    return t


def evaluate(P, score, idx):
    """Response AUROC and survival hazard ratio per SD of the score (age-adjusted)."""
    d = P.loc[idx].copy()
    d["score"] = (score.loc[idx] - score.loc[idx].mean()) / score.loc[idx].std()
    out = dict(n=len(d))
    y = d["responder"].dropna()
    if (y == 1).sum() >= 8 and (y == 0).sum() >= 8:
        out["n_response"], out["auc"] = len(y), TC.auc(y.astype(int).values, d.loc[y.index, "score"].values)[0]
    for a in ("pfs", "rfs", "os"):
        s = d[d[f"{a}_months"].notna() & d[f"{a}_event"].notna() & (d[f"{a}_months"] > 0)]
        if s[f"{a}_event"].sum() >= 15:
            s = s.assign(time=s[f"{a}_months"], event=s[f"{a}_event"], project="x", age_n=((pd.to_numeric(s["age"], errors="coerce") - 60) / 10).fillna(0))
            r = SL.cox_test(s, ["score", "age_n"])
            out[f"{a}_events"], out[f"{a}_hr"], out[f"{a}_p"] = int(s["event"].sum()), float(np.exp(r.loc["score", "coef"])), float(r.loc["score", "p"])
    return out


def main():
    P, Z = TC.load()
    G = RC.groups(P)
    y = P["responder"].fillna(-1).astype(int)
    train = {"both trials": list(G["bright"][0]) + list(G["ispy"][0]), "BrighTNess only": list(G["bright"][0]), "I-SPY2 only": list(G["ispy"][0]),
             "paclitaxel-only arms (control model)": list(G["bright"][1]) + list(G["ispy"][1])}
    zs, zp = RS.proximity("CARBOPLATIN"), N.precomputed_proximity()["CISPLATIN"].dropna()
    feats = {"network shortlist (STRING v12)": [c for c in Z.columns if c in set(zs.index[zs <= RC.CUT])],
             "paper shortlist": [c for c in Z.columns if c in set(zp.index[zp <= RC.CUT])], "all pathways": list(Z.columns)}
    T = targets(P)
    rng = np.random.default_rng(0)
    rows = []
    for tn, tr in train.items():
        for fn, cols in feats.items():
            m = TC.fit(Z.loc[tr, cols].values, y.loc[tr])
            score = pd.Series(m.decision_function(Z[cols].values), Z.index)
            null = None
            if tn == "both trials" and fn != "all pathways":
                null = []
                for _ in range(200):
                    rc = list(rng.choice(Z.columns, len(cols), replace=False))
                    null.append(pd.Series(TC.fit(Z.loc[tr, rc].values, y.loc[tr]).decision_function(Z[rc].values), Z.index))
            for name, idx in T.items():
                if len(idx) < 15:
                    continue
                r = dict(trained_on=tn, pathways=fn, n_pathways=len(cols), tested_in=name, **evaluate(P, score, idx))
                if null is not None:
                    ne = pd.DataFrame([evaluate(P, s, idx) for s in null])
                    for k, better in (("auc", "ge"), ("pfs_hr", "le"), ("os_hr", "le"), ("rfs_hr", "le")):
                        if k in r and k in ne:
                            r[f"{k}_random_mean"] = ne[k].mean()
                            r[f"{k}_share_random_as_good"] = float((ne[k] >= r[k]).mean() if better == "ge" else (ne[k] <= r[k]).mean())
                rows.append(r)
        print(tn, "done", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "frozen_test_outside_breast.tsv", sep="\t", index=False)
    pd.set_option("display.width", 300, "display.max_columns", 40, "display.max_colwidth", 45, "display.max_rows", 300)
    print(R.round(3).to_string())


if __name__ == "__main__":
    main()
