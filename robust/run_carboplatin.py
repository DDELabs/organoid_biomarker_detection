"""Follow-up of the carboplatin lead in the two randomised breast trials (BrighTNess, I-SPY2).

A. Single pathways (the paper's idea): for every pathway, test pathway x carboplatin-arm interaction on complete
   response in each trial, then combine the two trials. Compare network-shortlisted pathways with the rest.
   The paper's own bladder cisplatin biomarker is tested as a pre-specified hypothesis.
B. Multi-pathway score restricted to the network shortlist, in the same four tests as run_train_combined part 2,
   against all pathways and against random shortlists of the same size.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
import run_survlearn as RS  # noqa: E402
import run_train_combined as TC  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import network as N  # noqa: E402

OUT = RESULTS / "carboplatin"
CUT = -1.2816
PAPER_BIOMARKER = "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION"


def groups(P):
    B = P[P["cohort"] == "BRIGHTNESS"]
    I = P[(P["cohort"] == "ISPY2") & (P["her2"] == 0)]
    return dict(
        bright=(B.index[B["arm"] != "PACLITAXEL"], B.index[B["arm"] == "PACLITAXEL"]),
        ispy=(I.index[I["arm"].str.contains("ABT 888")], I.index[I["arm"] == "Paclitaxel"]),
        ispy_tn=(I.index[I["arm"].str.contains("ABT 888") & (I["er"] == 0)], I.index[(I["arm"] == "Paclitaxel") & (I["er"] == 0)]))


def one_pathway(P, x, drug, ctrl, cov):
    idx = list(drug) + list(ctrl)
    arm = np.r_[np.ones(len(drug)), np.zeros(len(ctrl))]
    s = x.loc[idx].values; s = (s - s.mean()) / s.std()
    X = pd.DataFrame({"arm": arm, "s": s, "axs": arm * s})
    for c in cov:
        X[c] = P.loc[idx, c].fillna(0).values
    r = sm.Logit(P.loc[idx, "responder"].values.astype(float), sm.add_constant(X)).fit(disp=0)
    return r.params["axs"], r.bse["axs"]


def part_a(P, Z, G, short):
    rows = []
    for p in Z.columns:
        if min(Z.loc[list(G[k][0]) + list(G[k][1]), p].std() for k in ("bright", "ispy")) < 1e-9:
            continue                                   # pathway not measured in one of the trials
        b1, s1 = one_pathway(P, Z[p], *G["bright"], [])
        b2, s2 = one_pathway(P, Z[p], *G["ispy"], ["er"])
        w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
        bm, sm_ = (b1 * w1 + b2 * w2) / (w1 + w2), np.sqrt(1 / (w1 + w2))
        rows.append(dict(pathway=p, coef_brightness=b1, p_brightness=2 * norm.sf(abs(b1 / s1)), coef_ispy2=b2, p_ispy2=2 * norm.sf(abs(b2 / s2)),
                         same_direction=np.sign(b1) == np.sign(b2), coef_combined=bm, p_combined=2 * norm.sf(abs(bm / sm_))))
    A = pd.DataFrame(rows).set_index("pathway")
    from statsmodels.stats.multitest import multipletests
    A["q_combined_all_pathways"] = multipletests(A["p_combined"], method="fdr_bh")[1]
    for name, s in short.items():
        A[f"shortlist_{name}"] = A.index.isin(s)
    A.sort_values("p_combined").to_csv(OUT / "single_pathway_interactions.tsv", sep="\t")
    return A


def part_b(P, Z, G, short, n_random=200, seed=0):
    rng = np.random.default_rng(seed)
    tests = {
        "I-SPY2 carboplatin arm vs paclitaxel alone (cross-fitted)": lambda Zs: TC.drug_specific(P, Zs, *G["ispy"], "", ["er"], seeds=3),
        "BrighTNess carboplatin arms vs paclitaxel alone (cross-fitted)": lambda Zs: TC.drug_specific(P, Zs, *G["bright"], "", [], seeds=3),
        "trained in BrighTNess, tested in I-SPY2 triple-negative": lambda Zs: TC.transfer(P, Zs, G["bright"][0], *G["ispy_tn"], "", []),
        "trained in I-SPY2, tested in BrighTNess": lambda Zs: TC.transfer(P, Zs, G["ispy"][0], *G["bright"], "", []),
    }
    rows = []
    for tname, fn in tests.items():
        base = fn(Z)
        for sname, s in short.items():
            cols = [c for c in Z.columns if c in s]
            r = fn(Z[cols])
            gap = r["auc_in_drug_arm"] - r["auc_in_control_arm"]
            rnd = []
            for _ in range(n_random):
                rr = fn(Z[list(rng.choice(Z.columns, len(cols), replace=False))])
                rnd.append(rr["auc_in_drug_arm"] - rr["auc_in_control_arm"])
            rows.append(dict(test=tname, shortlist=sname, pathways=len(cols), auc_drug_arm=r["auc_in_drug_arm"], auc_control_arm=r["auc_in_control_arm"],
                             gap=gap, p_interaction=r["p_interaction"], gap_all_pathways=base["auc_in_drug_arm"] - base["auc_in_control_arm"],
                             random_shortlist_gap_mean=np.mean(rnd), share_of_random_shortlists_with_gap_at_least_as_big=np.mean(np.array(rnd) >= gap)))
            print(rows[-1], flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "shortlist_score_tests.tsv", sep="\t", index=False)
    return R


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P, Z = TC.load()
    G = groups(P)
    print({k: (len(a), len(b)) for k, (a, b) in G.items()})
    zp = N.precomputed_proximity()["CISPLATIN"].dropna()
    zs = RS.proximity("CARBOPLATIN")
    short = {"paper_table_cisplatin": set(zp.index[zp <= CUT]), "string_v12_cisplatin": set(zs.index[zs <= CUT])}
    A = part_a(P, Z, G, short)
    pd.set_option("display.width", 250, "display.max_colwidth", 70)
    print(A.sort_values("p_combined").head(15).round(4).to_string())
    print("paper biomarker:", {k: (round(float(v), 4) if not isinstance(v, (bool, np.bool_)) else bool(v)) for k, v in A.loc[PAPER_BIOMARKER].items()})
    for name in short:
        s = A[f"shortlist_{name}"]
        print(name, "n", int(s.sum()), "p<0.05 combined in shortlist:", int((A.loc[s, "p_combined"] < 0.05).sum()), f"({(A.loc[s, 'p_combined'] < 0.05).mean():.1%})",
              "| outside:", f"{(A.loc[~s, 'p_combined'] < 0.05).mean():.1%}", "| same direction in shortlist:", f"{A.loc[s, 'same_direction'].mean():.1%}",
              "outside:", f"{A.loc[~s, 'same_direction'].mean():.1%}")
    print("all pathways: p<0.05 combined", int((A["p_combined"] < 0.05).sum()), "of", len(A), "| q<0.1:", int((A["q_combined_all_pathways"] < 0.1).sum()),
          "| same direction:", f"{A['same_direction'].mean():.1%}", "| p<0.05 in both trials same direction:",
          int(((A["p_brightness"] < 0.05) & (A["p_ispy2"] < 0.05) & A["same_direction"]).sum()))
    R = part_b(P, Z, G, short)
    print(R.round(3).to_string())


if __name__ == "__main__":
    main()
