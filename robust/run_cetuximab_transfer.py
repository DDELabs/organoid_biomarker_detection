"""Cetuximab: does a pathway signature learned in PDX transfer to other PDX and to patients?

Training: 192 colorectal liver-metastasis PDX treated with cetuximab (Isella 2017 / Bertotti), BRIDGE
domain vector (within-set rank-normalised sensitivity, model-only axes removed against TCGA COAD+HNSC).
Tests (frozen signature; null = signature permuted across pathways, 1000x):
  1. internal 5-fold CV in the training PDX (Spearman, predicted vs observed sensitivity)
  2. cross-lab PDX: Novartis PDXE cetuximab arms (CRC, NSCLC): Spearman with -BestAvgResponse, AUROC CR/PR/SD
  3. patients: ENLIGHT cetuximab cohort (head & neck, cetuximab + platinum/5-FU, n = 40) and TCGA HNSC
     cetuximab RECIST labels
Benchmark: AREG + EREG expression (EGFR-ligand biomarker of cetuximab benefit, Khambata-Ford 2007).

Outputs: robust/results/bridge/cetuximab_transfer.tsv, cetuximab_signature.tsv
Usage: python robust/run_cetuximab_transfer.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402

import run_atlas_external as AE  # noqa: E402
import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import bridge as B  # noqa: E402
from obd import curated as CU  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402

OUT = RESULTS / "bridge"
RNG = np.random.default_rng(21)
N_NULL = 1000


def zcols(Z):
    return ((Z - Z.mean()) / Z.std().replace(0, 1)).fillna(0)


def score(Z, beta):
    f = [p for p in beta.index if p in Z.columns]
    return pd.Series(zcols(Z[f]).values @ beta[f].values, index=Z.index)


def null_p(stat, Z, beta, obs):
    vals = np.array([stat(score(Z, pd.Series(RNG.permutation(beta.values), index=beta.index))) for _ in range(N_NULL)])
    vals = vals[np.isfinite(vals)]
    return float((1 + (vals >= obs).sum()) / (1 + len(vals)))


def ligand_score(expr, samples):
    g = [x for x in ("AREG", "EREG") if x in expr.index]
    e = expr.loc[g, samples].T
    return ((e - e.mean()) / e.std()).mean(1)


def main():
    rows = []
    M, proj = RT.patient_matrix(["COAD", "HNSC"])
    Xref = T.center_within(M, proj)
    pathways = list(Xref.columns)
    v = B.domain_vector("colorectal_pdx_isella2017", "CETUXIMAB", pathways, Xref.values)
    beta = v["beta"]
    beta.sort_values().to_csv(OUT / "cetuximab_signature.tsv", sep="\t", header=["sensitivity_coef"])
    print("training PDX:", v["n"], "dropped axes:", v["dropped_axes"])

    # 1. internal CV
    sc, resp = B.domain_scores("colorectal_pdx_isella2017", pathways)
    r = resp["CETUXIMAB"].dropna()
    r = r[r.index.isin(sc.index)]
    X = B._within(sc.loc[r.index].values.astype(float), np.array(["all"] * len(r)))
    y = -r.values
    pred = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        pred[te] = RidgeCV(alphas=np.logspace(0, 4, 9)).fit(X[tr], y[tr]).predict(X[te])
    rows.append({"test": "Isella PDX 5-fold CV", "n": len(y), "metric": "Spearman", "value": spearmanr(pred, y).correlation,
                 "null_p": np.nan, "AREG_EREG": spearmanr(ligand_score(CU.load_curated_preclinical("colorectal_pdx_isella2017")[0], r.index), y).correlation})

    # 2. cross-lab PDX (PDXE)
    expr, _ = CU.load_curated_preclinical("pdx_novartis_gao2015")
    scx, _ = B.domain_scores("pdx_novartis_gao2015", pathways)
    mr = pd.read_csv(CU.PRECLIN / "pdx_novartis_gao2015" / "mrecist.tsv", sep="\t")
    mr = mr[(mr["drug"] == "CETUXIMAB") & mr["sample"].isin(scx.index)].drop_duplicates("sample").set_index("sample")
    for tt, g in list(mr.groupby("tumour_type")) + [("CRC+NSCLC", mr)]:
        Z = scx.loc[g.index]
        sens = -g["best_avg_response"]
        s = score(Z, beta)
        rho = spearmanr(s, sens).correlation
        benefit = g["mrecist"].isin(["CR", "PR", "SD"]).astype(int)
        auc = roc_auc_score(benefit, s) if benefit.nunique() == 2 else np.nan
        rows.append({"test": f"PDXE cetuximab {tt}", "n": len(g), "metric": "Spearman (sensitivity)", "value": rho,
                     "null_p": null_p(lambda q: spearmanr(q, sens).correlation, Z, beta, rho),
                     "AUROC_benefit": auc, "n_benefit": int(benefit.sum()),
                     "AREG_EREG": spearmanr(ligand_score(expr, g.index), sens).correlation})

    # 3. patients
    for tag, loader in AE.cohorts():
        if "Cetuximab" not in tag:
            continue
        e, c = loader()
        Zc = AE.pathway_scores(tag.replace(":", "_"), e).T
        c = c[c["responder"].notna() & c.index.isin(Zc.index)]
        yv = c["responder"].astype(int)
        Z = Zc.loc[c.index]
        s = score(Z, beta)
        auc = roc_auc_score(yv, s)
        rows.append({"test": f"patients {tag} (HNSC, cetuximab + chemo)", "n": len(c), "metric": "AUROC (response)", "value": auc,
                     "null_p": null_p(lambda q: roc_auc_score(yv, q), Z, beta, auc), "n_benefit": int(yv.sum()),
                     "AREG_EREG": roc_auc_score(yv, ligand_score(e, c.index)) if {"AREG", "EREG"} & set(e.index) else np.nan})
    lab = R.patient_labels("CETUXIMAB")
    lab = lab[lab["responder"].notna() & lab["project"].eq("HNSC") & lab.index.isin(Xref.index)]
    if lab["responder"].nunique() == 2:
        yv = lab["responder"].astype(int)
        Z = Xref.loc[lab.index]
        s = score(Z, beta)
        auc = roc_auc_score(yv, s)
        rows.append({"test": "patients TCGA HNSC cetuximab (RECIST)", "n": len(lab), "metric": "AUROC (response)", "value": auc,
                     "null_p": null_p(lambda q: roc_auc_score(yv, q), Z, beta, auc), "n_benefit": int(yv.sum())})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "cetuximab_transfer.tsv", sep="\t", index=False)
    print(out.round(3).to_string())
    print("top sensitivity pathways:", list(beta.sort_values().index[-8:][::-1]))
    print("top resistance pathways:", list(beta.sort_values().index[:8]))


if __name__ == "__main__":
    main()
