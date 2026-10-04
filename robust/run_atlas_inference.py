"""ATLAS pathway inference from patient data alone (no prior means).

Prior means make drug-layer effects non-zero by construction, so pathway significance is
computed from the prior-free multi-task model (shared + class + drug layers).
Bootstrap over patients (parallel); per layer and per drug total effect:
z = mean / sd, two-sided p, BH-FDR within each output table, sign consistency.

Outputs (robust/results/atlas/):
  inference_layers.tsv.gz        pathway x layer (shared, class:*, drug:*)
  inference_drug_total.tsv.gz    drug x pathway total effect (shared + class + drug)

Usage: python robust/run_atlas_inference.py [--boot 200]
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
from obd import atlas as A  # noqa: E402

OUT = RESULTS / "atlas"


def bh(p):
    """Benjamini-Hochberg q-values (positional, safe with duplicate labels)."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    out = np.empty_like(q)
    out[order] = q
    return pd.Series(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=200)
    a = ap.parse_args()
    rec, X, cov, w, pathways, _ = RA.build_records()
    y = rec["responder"].astype(float).values
    des = RA.variant("atlas", pathways)
    pats = rec["patient"].unique()

    def one(seed):
        rng = np.random.default_rng(seed)
        cnt = pd.Series(rng.choice(pats, len(pats), replace=True)).value_counts()
        ww = w * rec["patient"].map(cnt).fillna(0).values
        keep = ww > 0
        _, b = A.fit(des, X[keep], rec["drug"].values[keep], cov[keep], y[keep], RA.LAM, ww[keep])
        return des.effects(b).values

    cache = OUT / "inference_bootstrap_stack.npy"
    if cache.exists():
        stack = np.load(cache)
    else:
        stack = np.stack(Parallel(n_jobs=4)(delayed(one)(s) for s in range(a.boot)))
        np.save(cache, stack)
    _, b = A.fit(des, X, rec["drug"].values, cov, y, RA.LAM, w)
    eff = des.effects(b)
    cols = list(eff.columns)
    mean, sd = stack.mean(0), stack.std(0) + 1e-12

    long = []
    for j, layer in enumerate(cols):
        z = mean[:, j] / sd[:, j]
        p = pd.Series(2 * norm.sf(np.abs(z)), index=pathways)
        long.append(pd.DataFrame({"pathway": pathways, "layer": layer, "effect": eff[layer].values,
                                  "boot_mean": mean[:, j], "z": z, "p": p.values, "q": bh(p).values,
                                  "sign_consistency": (np.sign(stack[:, :, j]) == np.sign(mean[:, j])).mean(0)}))
    layers = pd.concat(long)
    layers.to_csv(OUT / "inference_layers.tsv.gz", sep="\t", index=False)

    tot = []
    for d in RA.DRUGS:
        idx = [cols.index(c) for c in ("shared", f"class:{A.DRUG_CLASS[d]}", f"drug:{d}")]
        t = stack[:, :, idx].sum(2)
        z = t.mean(0) / (t.std(0) + 1e-12)
        p = pd.Series(2 * norm.sf(np.abs(z)), index=pathways)
        tot.append(pd.DataFrame({"drug": d, "class": A.DRUG_CLASS[d], "pathway": pathways,
                                 "effect": eff[[cols[i] for i in idx]].sum(1).values, "z": z, "p": p.values,
                                 "sign_consistency": (np.sign(t) == np.sign(t.mean(0))).mean(0)}))
    tot = pd.concat(tot, ignore_index=True)
    tot["q"] = bh(tot["p"]).values  # FDR across all drug x pathway tests
    tot.to_csv(OUT / "inference_drug_total.tsv.gz", sep="\t", index=False)
    print("layers q<0.1:", layers[layers.q < 0.1].groupby("layer").size().to_dict())
    print("drug totals q<0.1:", tot[tot.q < 0.1].groupby("drug").size().to_dict())


if __name__ == "__main__":
    main()
