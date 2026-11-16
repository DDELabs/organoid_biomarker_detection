"""Third test of the carboplatin BENEFIT score, frozen from the two randomised trials.

New data (not used anywhere before):
  GSE106977  triple-negative breast, neoadjuvant chemotherapy with (31) or without (88) carboplatin, complete response.
             Not randomised.
  GSE18864   triple-negative breast, cisplatin alone before surgery (24 patients), Miller-Payne response grade 0-5.
             Single arm: with a platinum given alone, response itself is the drug effect.
The benefit score = pathway x arm coefficients learned in BrighTNess, in I-SPY2, or in both (network shortlist and
all pathways). Nothing is refitted on the new data.
"""
import gzip
import io
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data" / "curated" / "trials"))
warnings.filterwarnings("ignore")
import build_geo_trials as B  # noqa: E402
import build_training_cohorts as BC  # noqa: E402
import run_combination as CB  # noqa: E402
import run_train_combined as TC  # noqa: E402
from build_training_table import score  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd.geo_matrix import read_expr, read_meta  # noqa: E402

OUT = RESULTS / "combination"


def hta2_map():
    p = B.gpl_file("GPL17586")
    txt = open(p, errors="replace").read().split("\n")
    a = next(i for i, l in enumerate(txt) if l.startswith("!platform_table_begin")) + 1
    z = next((i for i, l in enumerate(txt) if l.startswith("!platform_table_end")), len(txt))
    t = pd.read_csv(io.StringIO("\n".join(txt[a:z])), sep="\t", dtype=str, quoting=3, usecols=["ID", "gene_assignment"])
    sym = t["gene_assignment"].str.split(" // ").str[1]
    s = pd.Series(sym.values, index=t["ID"]).dropna().str.strip()
    return s[(s != "") & (s != "---")]


def zscores(expr, name):
    s = score(f"third_{name}_reactome_rank", expr)
    return (s - s.mean()) / s.std().replace(0, 1)


def load_new():
    m1 = read_meta(BC.m("GSE106977"))
    e1 = read_expr(BC.m("GSE106977")); e1.index = e1.index.astype(str)
    e1 = B.collapse(e1.apply(pd.to_numeric, errors="coerce"), hta2_map())
    c1 = pd.DataFrame({"arm": m1["treatment"].str.contains("carboplatin").astype(float),
                       "y": m1["pathological complete response"].map({"yes": 1, "no": 0})}, index=m1.index)
    m2 = read_meta(BC.m("GSE18864"))
    m2 = m2[m2["source_name_ch1"].str.contains("Cisplatin trial")]
    e2 = BC.array_expr("GSE18864", "GPL570")
    mp = pd.to_numeric(m2["miller-payne response"], errors="coerce")
    c2 = pd.DataFrame({"mp": mp, "patient": m2["description"].str.extract(r"patient (\d+)")[0]}, index=m2.index)
    return (zscores(e1, "GSE106977"), c1), (zscores(e2, "GSE18864"), c2, list(m2.index)), len(e1)


def main():
    P, Z = TC.load()
    I = P[(P["cohort"] == "ISPY2") & (P["her2"] == 0)]
    Bt = P[P["cohort"] == "BRIGHTNESS"]
    groups = {"BrighTNess": (Bt.index[Bt["arm"] != "PACLITAXEL"], Bt.index[Bt["arm"] == "PACLITAXEL"], []),
              "I-SPY2": (I.index[I["arm"].str.contains("ABT 888")], I.index[I["arm"] == "Paclitaxel"], ["er"])}
    (Z1, c1), (Z2, c2, cis), ng = load_new()
    print("GSE106977 genes", ng, c1.groupby("arm")["y"].agg(["size", "mean"]).round(3).to_dict(), "| cisplatin samples", len(cis), "patients", c2["patient"].nunique())
    rows = []
    for fname, cols in (("network shortlist", CB.shortlist(["CARBOPLATIN"], Z.columns)), ("all pathways", list(Z.columns))):
        cols = [c for c in cols if c in Z1.columns and c in Z2.columns]
        gam = {}
        for g, (d, c, cov) in groups.items():
            idx = np.array(list(d) + list(c)); arm = np.r_[np.ones(len(d)), np.zeros(len(c))]
            gam[g] = CB.fit_benefit(Z.loc[idx, cols].values, arm, P.loc[idx, cov].fillna(0).values.reshape(len(idx), -1), P.loc[idx, "responder"].astype(int).values)
        gam["both trials (average)"] = (gam["BrighTNess"] / np.linalg.norm(gam["BrighTNess"]) + gam["I-SPY2"] / np.linalg.norm(gam["I-SPY2"])) / 2
        for g, w in gam.items():
            s1 = Z1[cols].fillna(0).values @ w
            t = CB.test(c1["arm"].values, c1["y"].values.astype(int), np.zeros((len(c1), 0)), s1)
            s2 = pd.Series(Z2.loc[cis, cols].fillna(0).values @ w, cis).groupby(c2["patient"]).mean()   # replicate biopsies averaged
            mp = c2.groupby("patient")["mp"].first().loc[s2.index]
            rho, p = spearmanr(s2, mp)
            rows.append(dict(pathways=fname, n_pathways=len(cols), score_learned_in=g, gse106977_gain_high_half=t["gain_high_half"], gse106977_gain_low_half=t["gain_low_half"],
                             gse106977_p_interaction=t["p_interaction"], cisplatin_patients=len(s2), cisplatin_spearman=rho, cisplatin_p=p))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "carboplatin_third_test.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(R.round(3).to_string())


if __name__ == "__main__":
    main()
