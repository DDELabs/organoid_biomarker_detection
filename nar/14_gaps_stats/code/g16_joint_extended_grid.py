"""G16 part 2: re-run the joint organoid+patient models (nar/03_organoid_plus_patient/code/run_joint.py, network
shortlist features, unchanged data/folds/seed) with the lab weight grid extended to w = 0, 0.05, 0.1, 0.25, 0.5, 0.75.

Output columns in the score files: 'joint' = w and penalty chosen by inner CV over the extended grid;
'joint (w=0.1)' = w fixed at 0.1 (penalty tuned); 'patient-only' = w 0; 'organoid-only' unchanged.
Usage: python g16_joint_extended_grid.py [--jobs 2]
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "_data" / "g16_extended"
OUT.mkdir(parents=True, exist_ok=True)
J = HERE.parents[1] / "03_organoid_plus_patient" / "code"
sys.path.insert(0, str(J))
import pandas as pd  # noqa: E402

import common as C  # noqa: E402  (03_organoid_plus_patient/code/common.py)
import run_joint as RJ  # noqa: E402
from units import UNITS  # noqa: E402

WS_EXT = (0.0, 0.05, 0.1, 0.25, 0.5, 0.75)
C.WS = WS_EXT


def fit_all_ext(Xo, yo, blocks, seed):
    out = {}
    (w0, l0), _ = C.tune(None, None, blocks, ws=(0.0,), seed=seed)
    out["patient-only"] = (C.fit_joint(None, None, blocks, 0.0, l0)[0], dict(w=0.0, lam=l0))
    (wj, lj), resj = C.tune(Xo, yo, blocks, ws=WS_EXT, seed=seed)
    out["joint"] = (C.fit_joint(Xo, yo, blocks, wj, lj)[0], dict(w=wj, lam=lj))
    lf = max((k for k in resj if k[0] == 0.1), key=lambda k: resj[k])[1]
    out["joint (w=0.5)"] = (C.fit_joint(Xo, yo, blocks, 0.1, lf)[0], dict(w=0.1, lam=lf))  # column renamed below
    return out


RJ.fit_all = fit_all_ext


def one(u):
    tag = f"{u['id']}__shortlist_extgrid"
    if (OUT / f"{tag}__scores.tsv.gz").exists():
        return tag
    D, sc, ev, it, fi, meta, extra = RJ.run_unit(u, "shortlist")
    ren = {"joint (w=0.5)": "joint (w=0.1)"}
    sc = sc.rename(columns=ren)
    ev["model"] = ev["model"].replace(ren)
    fi["model"] = fi["model"].replace(ren)
    ev.assign(unit=u["id"]).to_csv(OUT / f"{tag}__eval.tsv", sep="\t", index=False)
    fi.assign(unit=u["id"]).to_csv(OUT / f"{tag}__folds.tsv", sep="\t", index=False)
    out = D[["unit_cohort", "cohort", "group", "time", "event", "endpoint", "responder", "has_surv", "has_resp"]].join(sc)
    out.to_csv(OUT / f"{tag}__scores.tsv.gz", sep="\t")
    print(tag, meta["seconds"], "s", flush=True)
    return tag


if __name__ == "__main__":
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 2
    import multiprocessing as mp
    order = sorted(UNITS, key=lambda u: u["id"] != "CRC_5FU")  # slowest first
    with mp.get_context("fork").Pool(jobs) as pool:          # fork keeps the patched fit_all in the workers
        for t in pool.imap_unordered(one, order):
            print("done", t, flush=True)
