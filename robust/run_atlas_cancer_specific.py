"""Top-10 pathway biomarkers for every drug x cancer pair (TCGA RECIST labels).

For each pair with >= 20 labelled patients and >= 5 in each outcome class:
  local z    per-pathway association with response inside that cancer
             (logistic regression on the within-cancer z-scored pathway, adjusted for
             treatment setting when it varies); BH-FDR q within the pair
  atlas z    the drug's pan-cancer ATLAS total effect z (shared + class + drug layers)
  combined   Stouffer combination (local z + atlas z) / sqrt(2): local evidence,
             stabilised by the pan-cancer prior for the same drug
  stability  fraction of 200 patient bootstraps in which the pathway is among the
             pair's top 10 by |local z|
Top 10 per pair = 5 strongest response (combined z > 0) + 5 strongest resistance.

Outputs: robust/results/atlas/cancer_specific/{top10_by_cancer_drug.tsv, all_pairs_local.tsv.gz,
         TOP10_BY_CANCER.md, top10_by_cancer_drug.xlsx}

Usage: python robust/run_atlas_cancer_specific.py [--boot 200] [--min-n 20]
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

import run_atlas as RA  # noqa: E402
from obd import RESULTS  # noqa: E402

OUT = RESULTS / "atlas" / "cancer_specific"
CANCER = {"ACC": "adrenocortical", "BLCA": "bladder", "BRCA": "breast", "CESC": "cervical", "CHOL": "cholangiocarcinoma",
          "COAD": "colon", "ESCA": "oesophageal", "GBM": "glioblastoma", "HNSC": "head & neck", "KIRC": "kidney clear cell",
          "KIRP": "kidney papillary", "LGG": "lower-grade glioma", "LIHC": "liver", "LUAD": "lung adeno",
          "LUSC": "lung squamous", "MESO": "mesothelioma", "OV": "ovarian", "PAAD": "pancreatic", "READ": "rectal",
          "SARC": "sarcoma", "SKCM": "melanoma", "STAD": "stomach", "TGCT": "testicular germ cell", "THYM": "thymoma",
          "UCEC": "endometrial", "UCS": "uterine carcinosarcoma", "PCPG": "pheochromocytoma", "PRAD": "prostate"}


def local_z(X, y, S):
    """Per-pathway Wald z of a logistic regression y ~ x (+ setting), vectorised by IRLS."""
    n, p = X.shape
    Z = np.full(p, np.nan)
    base = np.ones((n, 1)) if S is None or S.shape[1] == 0 else np.column_stack([np.ones(n), S])
    for j in range(p):
        A = np.column_stack([X[:, j], base])
        beta = np.zeros(A.shape[1])
        for _ in range(25):
            eta = np.clip(A @ beta, -30, 30)
            mu = 1 / (1 + np.exp(-eta))
            W = mu * (1 - mu) + 1e-9
            H = A.T @ (A * W[:, None]) + 1e-6 * np.eye(A.shape[1])
            step = np.linalg.solve(H, A.T @ (y - mu))
            beta += step
            if np.abs(step).max() < 1e-6:
                break
        se = np.sqrt(np.linalg.inv(H)[0, 0])
        Z[j] = beta[0] / se
    return Z


def one_pair(proj, drug, Xp, y, S, az, pathways, n_boot, seed):
    z = local_z(Xp, y, S)
    q = bh(2 * norm.sf(np.abs(z)))
    comb = (z + az) / np.sqrt(2)
    rng = np.random.default_rng(seed)
    hits = np.zeros(len(pathways))  # bootstrap stability of local top-10 membership
    for _ in range(n_boot):
        b = rng.integers(0, len(y), len(y))
        if y[b].min() == y[b].max():
            continue
        zb = local_z(Xp[b], y[b], S[b] if S is not None else None)
        hits[np.argsort(-np.abs(np.nan_to_num(zb)))[:10]] += 1
    df = pd.DataFrame({"cancer": proj, "cancer_name": CANCER.get(proj, proj), "drug": drug, "n": len(y),
                       "responders": int(y.sum()), "pathway": pathways, "local_z": z, "local_q": q,
                       "atlas_z": az, "combined_z": comb, "top10_stability": hits / n_boot})
    resp = df.sort_values("combined_z", ascending=False).head(5).assign(direction="response")
    res = df.sort_values("combined_z").head(5).assign(direction="resistance")
    top = pd.concat([resp, res])
    top["rank"] = list(range(1, 6)) + list(range(1, 6))
    print(proj, drug, len(y), "pts | top response:", resp.pathway.iloc[0], flush=True)
    return df, top


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    r = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--min-n", type=int, default=20)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rec, X, cov, w, pathways, _ = RA.build_records()
    tot = pd.read_csv(RESULTS / "atlas" / "inference_drug_total.tsv.gz", sep="\t")
    atlas_z = tot.pivot(index="pathway", columns="drug", values="z").reindex(pathways)
    setting_cols = cov[:, :-1]  # last covariate column is proliferation
    jobs = []
    for (proj, drug), g in rec.groupby(["project", "drug"]):
        y = g["responder"].astype(float).values
        if len(g) < a.min_n or y.sum() < 5 or (1 - y).sum() < 5:
            continue
        idx = g.index.values
        S = setting_cols[idx]
        S = S[:, S.std(0) > 0] if S.size else None
        az = atlas_z[drug].values if drug in atlas_z else np.zeros(len(pathways))
        jobs.append((proj, drug, X[idx], y, S, az))
    res = Parallel(n_jobs=4)(delayed(one_pair)(*j, pathways, a.boot, k) for k, j in enumerate(jobs))
    all_rows, top_rows = [r[0] for r in res], [r[1] for r in res]

    top = pd.concat(top_rows, ignore_index=True)
    top["pathway_name"] = top["pathway"].str.replace("REACTOME_", "").str.replace("_", " ").str.lower()
    cols = ["cancer", "cancer_name", "drug", "n", "responders", "direction", "rank", "pathway_name",
            "local_z", "local_q", "atlas_z", "combined_z", "top10_stability", "pathway"]
    top = top[cols].sort_values(["cancer", "drug", "direction", "rank"])
    top.to_csv(OUT / "top10_by_cancer_drug.tsv", sep="\t", index=False)
    pd.concat(all_rows, ignore_index=True).to_csv(OUT / "all_pairs_local.tsv.gz", sep="\t", index=False)
    try:
        with pd.ExcelWriter(OUT / "top10_by_cancer_drug.xlsx") as xw:
            for c, gc in top.groupby("cancer"):
                gc.drop(columns=["pathway"]).round(3).to_excel(xw, sheet_name=c[:31], index=False)
    except Exception as exc:
        print("xlsx skipped:", exc)

    # markdown report
    L = ["# Top-10 pathway biomarkers per cancer type and drug (TCGA RECIST response)", "",
         "Each list = 5 response + 5 resistance pathways ranked by the combined z (local evidence in that "
         "cancer, stabilised by the drug's pan-cancer ATLAS effect). `local z` = evidence inside the cancer only; "
         "`stab` = fraction of 200 bootstraps in which the pathway is in that cancer's local top 10. "
         "Local q < 0.1 is marked *. These are exploratory (20-130 patients per pair).", ""]
    for c, gc in top.groupby("cancer"):
        L.append(f"## {c} - {CANCER.get(c, c)}")
        for d, gd in gc.groupby("drug"):
            n, r = gd["n"].iloc[0], gd["responders"].iloc[0]
            L.append(f"\n**{d.title()}** (n = {n}, responders = {r})\n")
            L.append("| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |")
            L.append("|---|---|---|---|---|---|---|---|---|")
            rs = gd[gd.direction == "response"].reset_index(drop=True)
            rz = gd[gd.direction == "resistance"].reset_index(drop=True)
            for i in range(5):
                a1, b1 = rs.iloc[i], rz.iloc[i]
                L.append(f"| {i + 1} | {a1.pathway_name[:55]}{'*' if a1.local_q < 0.1 else ''} | {a1.local_z:+.1f} | "
                         f"{a1.combined_z:+.1f} | {a1.top10_stability:.2f} | {b1.pathway_name[:55]}{'*' if b1.local_q < 0.1 else ''} | "
                         f"{b1.local_z:+.1f} | {b1.combined_z:+.1f} | {b1.top10_stability:.2f} |")
        L.append("")
    (OUT / "TOP10_BY_CANCER.md").write_text("\n".join(L) + "\n")
    print("pairs:", top.groupby(["cancer", "drug"]).ngroups, "cancers:", top.cancer.nunique(), "drugs:", top.drug.nunique())


if __name__ == "__main__":
    main()
