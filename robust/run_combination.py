"""Combination strategy: who gains from ADDING a drug to the backbone, and which combination suits a patient.

Data: the two randomised breast trials in the training table.
  I-SPY2 (HER2-negative patients): paclitaxel alone vs paclitaxel + one added drug (7 arms).
  BrighTNess: paclitaxel alone vs paclitaxel + carboplatin (with or without veliparib).
For each added drug:
  1. network shortlist = pathways close to the ADDED drug's targets (STRING v12 proximity, z <= -1.2816);
  2. benefit model = ridge logistic with pathway main effects and pathway x arm terms, fitted on training folds;
     a patient's benefit score is the pathway x arm part: how much MORE likely response is with the added drug;
  3. test on held-out patients (cross-fitted): does the added drug help more in the high-benefit-score half than
     in the low half? (score x arm interaction), against 200 random shortlists and against all pathways.
Then: carboplatin benefit score learned in one trial and tested in the other; and a recommendation test across
all I-SPY2 arms (do patients who happened to receive the arm ranked best for them respond more often?).
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
import run_survlearn as RS  # noqa: E402
import run_train_combined as TC  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import network as N  # noqa: E402
from obd import survlearn as SL  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "combination"
CUT = -1.2816
LAM = 100.0
TARGETS = {"NERATINIB": ["EGFR", "ERBB2", "ERBB4"], "VELIPARIB": ["PARP1", "PARP2"], "TREBANANIB": ["ANGPT1", "ANGPT2"],
           "GANITUMAB": ["IGF1R"], "GANETESPIB": ["HSP90AA1", "HSP90AB1"], "PEMBROLIZUMAB": ["PDCD1"], "MK-2206": ["AKT1", "AKT2", "AKT3"]}
ARMS = {"Paclitaxel + Neratinib": ["NERATINIB"], "Paclitaxel + ABT 888 + Carboplatin": ["VELIPARIB", "CARBOPLATIN"],
        "Paclitaxel + AMG 386": ["TREBANANIB"], "Paclitaxel + Ganitumab": ["GANITUMAB"], "Paclitaxel + Ganetespib": ["GANETESPIB"],
        "Paclitaxel + Pembrolizumab": ["PEMBROLIZUMAB"], "Paclitaxel + MK-2206": ["MK-2206"]}
_G = {}


def shortlist(drugs, cols):
    keep = set()
    for d in drugs:
        if d == "CARBOPLATIN":
            z = RS.proximity("CARBOPLATIN")
        else:
            def compute(d=d):
                if "g" not in _G:
                    _G["g"] = N.load_string()
                return N.proximity_z(_G["g"], [t for t in TARGETS[d] if t in _G["g"]], gene_sets(("REACTOME",))).to_frame("z")
            z = scoring.cached(f"proximity_STRINGv12_{d}", compute)["z"]
        keep |= set(z.index[z <= CUT])
    return [c for c in cols if c in keep]


def fit_benefit(X, arm, cov, y):
    """Returns the pathway x arm coefficients. Arm and covariates are scaled up so the ridge barely touches them."""
    A = arm[:, None]
    M = np.hstack([10 * A, 10 * cov, X, X * A])
    m = LogisticRegression(C=1.0 / LAM, max_iter=3000).fit(M, y)
    return m.coef_[0][-X.shape[1]:]


def crossfit(P, Z, drug, ctrl, cols, cov_cols, seeds=3):
    idx = np.array(list(drug) + list(ctrl))
    arm = np.r_[np.ones(len(drug)), np.zeros(len(ctrl))]
    y = P.loc[idx, "responder"].astype(int).values
    X = Z.loc[idx, cols].values
    cov = P.loc[idx, cov_cols].fillna(0).values.reshape(len(idx), -1)
    s = np.zeros(len(idx))
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        f = np.empty(len(idx), int)
        for a in (0, 1):                                    # folds balanced within each arm
            i = np.flatnonzero(arm == a)
            f[i] = rng.permutation(len(i)) % 5
        for k in range(5):
            tr, te = f != k, f == k
            s[te] += X[te] @ fit_benefit(X[tr], arm[tr], cov[tr], y[tr]) / seeds
    return idx, arm, y, cov, s


def test(arm, y, cov, s):
    """Added benefit (response rate with the drug minus without) in the high and low halves, and the interaction test."""
    z = (s - s.mean()) / s.std()
    hi = z >= np.median(z)
    def gain(m):
        return y[m & (arm == 1)].mean() - y[m & (arm == 0)].mean()
    D = pd.DataFrame({"arm": arm, "s": z, "axs": arm * z})
    for j in range(cov.shape[1]):
        if np.std(cov[:, j]) > 0:
            D[f"c{j}"] = cov[:, j]
    r = sm.Logit(y.astype(float), sm.add_constant(D)).fit(disp=0)
    return dict(gain_high_half=gain(hi), gain_low_half=gain(~hi), interaction_coef=r.params["axs"], p_interaction=r.pvalues["axs"])


def one(P, Z, label, drug, ctrl, added, cov_cols, rng, n_random=200):
    cols = shortlist(added, Z.columns)
    idx, arm, y, cov, s = crossfit(P, Z, drug, ctrl, cols, cov_cols)
    r = dict(comparison=label, n_drug=len(drug), n_control=len(ctrl), response_drug=y[arm == 1].mean(), response_control=y[arm == 0].mean(),
             shortlist=len(cols), **test(arm, y, cov, s))
    rnd = []
    for _ in range(n_random):
        rc = list(rng.choice(Z.columns, len(cols), replace=False))
        rnd.append(test(arm, y, cov, crossfit(P, Z, drug, ctrl, rc, cov_cols, seeds=1)[4])["interaction_coef"])
    r["random_shortlists_as_good"] = float((np.array(rnd) >= r["interaction_coef"]).mean())
    a = test(arm, y, cov, crossfit(P, Z, drug, ctrl, list(Z.columns), cov_cols)[4])
    r["all_pathways_gain_high"], r["all_pathways_gain_low"], r["all_pathways_p"] = a["gain_high_half"], a["gain_low_half"], a["p_interaction"]
    return r, pd.Series(s, idx)


def transfer(P, Z, label, tr_drug, tr_ctrl, te_drug, te_ctrl, cols, cov_tr, cov_te):
    idx = np.array(list(tr_drug) + list(tr_ctrl)); arm = np.r_[np.ones(len(tr_drug)), np.zeros(len(tr_ctrl))]
    g = fit_benefit(Z.loc[idx, cols].values, arm, P.loc[idx, cov_tr].fillna(0).values.reshape(len(idx), -1), P.loc[idx, "responder"].astype(int).values)
    ti = np.array(list(te_drug) + list(te_ctrl)); ta = np.r_[np.ones(len(te_drug)), np.zeros(len(te_ctrl))]
    return dict(comparison=label, n_drug=len(te_drug), n_control=len(te_ctrl), shortlist=len(cols),
                **test(ta, P.loc[ti, "responder"].astype(int).values, P.loc[ti, cov_te].fillna(0).values.reshape(len(ti), -1), Z.loc[ti, cols].values @ g))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P, Z = TC.load()
    rng = np.random.default_rng(0)
    I = P[(P["cohort"] == "ISPY2") & (P["her2"] == 0)]
    ctrl = I.index[I["arm"] == "Paclitaxel"]
    rows, scores = [], {}
    for arm, added in ARMS.items():
        d = I.index[I["arm"] == arm]
        r, s = one(P, Z, f"I-SPY2: add {' + '.join(added).lower()}", d, ctrl, added, ["er"], rng)
        rows.append(r); scores[arm] = s; print(r, flush=True)
    B = P[P["cohort"] == "BRIGHTNESS"]
    bd, bc = B.index[B["arm"] != "PACLITAXEL"], B.index[B["arm"] == "PACLITAXEL"]
    r, _ = one(P, Z, "BrighTNess: add carboplatin", bd, bc, ["CARBOPLATIN"], [], rng); rows.append(r); print(r, flush=True)
    R = pd.DataFrame(rows); R.to_csv(OUT / "added_drug_benefit.tsv", sep="\t", index=False)

    vc = I.index[I["arm"].str.contains("ABT 888")]
    tn = I[I["er"] == 0]
    tr = []
    for name, cols in (("network shortlist", shortlist(["CARBOPLATIN"], Z.columns)), ("all pathways", list(Z.columns))):
        tr.append(transfer(P, Z, f"carboplatin benefit: learned in BrighTNess, tested in I-SPY2 triple-negative ({name})", bd, bc,
                           tn.index[tn["arm"].str.contains("ABT 888")], tn.index[tn["arm"] == "Paclitaxel"], cols, [], []))
        tr.append(transfer(P, Z, f"carboplatin benefit: learned in I-SPY2, tested in BrighTNess ({name})", vc, ctrl, bd, bc, cols, ["er"], []))
    Tr = pd.DataFrame(tr); Tr.to_csv(OUT / "carboplatin_benefit_transfer.tsv", sep="\t", index=False)

    # recommendation test: each patient's cross-fitted benefit score for the arm they actually received, turned into a
    # within-arm percentile; controls excluded. Does a higher percentile go with response, pooled over arms?
    rec = []
    for arm, s in scores.items():
        d = I.index[I["arm"] == arm]
        rec.append(pd.DataFrame({"arm": arm, "pct": s.loc[d].rank(pct=True), "y": P.loc[d, "responder"].astype(int), "er": P.loc[d, "er"]}))
    rec = pd.concat(rec)
    c_s = pd.concat([s.loc[[i for i in s.index if i in set(ctrl)]].rank(pct=True) for s in scores.values()], axis=1).mean(axis=1)
    top, bot = rec[rec["pct"] >= 0.5], rec[rec["pct"] < 0.5]
    summ = dict(patients_on_added_drug=len(rec), response_if_ranked_high_for_own_arm=top["y"].mean(), response_if_ranked_low_for_own_arm=bot["y"].mean(),
                control_response_high_half=P.loc[c_s.index[c_s >= c_s.median()], "responder"].mean(),
                control_response_low_half=P.loc[c_s.index[c_s < c_s.median()], "responder"].mean())
    pd.Series(summ).to_csv(OUT / "recommendation_summary.tsv", sep="\t", header=False)
    pd.set_option("display.width", 280, "display.max_colwidth", 90, "display.max_columns", 30)
    print(R.round(3).to_string()); print(Tr.round(3).to_string()); print({k: round(float(v), 3) for k, v in summ.items()})


if __name__ == "__main__":
    main()
