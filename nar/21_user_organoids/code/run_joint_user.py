"""Organoid + patient (joint) method for a user's organoid set, reusing nar/03_organoid_plus_patient/code unchanged.

Usage (called by import_organoids.py, one process):
  python3 run_joint_user.py <set_dir> <set_name> <cancer> <out_dir> [--drugs A,B] [--features shortlist|all|both]
It takes every analysis unit of nar/03 (cancer x drug, with its patient cohorts) that matches the cancer, swaps the lab models of that
unit for the user's organoids, and runs the same code: organoid-only model, patient-only model, and the joint model, evaluated
on held-out patient cohorts.  Everything is written under <out_dir>; nar/03 is only read.
"""
import argparse
import copy
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analysis_common as AC  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("set_dir")
ap.add_argument("set_name")
ap.add_argument("cancer")
ap.add_argument("out_dir")
ap.add_argument("--drugs", default="")
ap.add_argument("--features", default="shortlist", choices=["shortlist", "all", "both"])
a = ap.parse_args()

sys.path.insert(0, str(AC.NAR / "03_organoid_plus_patient" / "code"))
import pandas as pd  # noqa: E402
import common as C  # noqa: E402  (nar/03 code)
import run_joint as RJ  # noqa: E402
from units import UNITS  # noqa: E402
from obd import bridge as BR  # noqa: E402
from obd import curated as CU  # noqa: E402
from obd.reference import common_drug_name  # noqa: E402

out = Path(a.out_dir)
out.mkdir(parents=True, exist_ok=True)
projects = AC.tcga_projects(a.cancer)
group = AC.unit_group(projects)
reg = AC.make_registry(a.set_dir, out, a.set_name)

# ---- redirect writes; reuse nar/03's already computed pathway-score caches for the other sets (read-only symlinks)
cache, tables = out / "cache", out / "tables"
cache.mkdir(exist_ok=True)
tables.mkdir(exist_ok=True)
for p in (AC.NAR / "03_organoid_plus_patient" / "cache").glob("scores_*.tsv.gz"):
    AC.link(p, cache / p.name)
user_scores = pd.read_csv(Path(a.set_dir) / "pathway_scores.tsv.gz", sep="\t", index_col=0)
user_scores.index = user_scores.index.astype(str)
user_scores.to_csv(cache / f"scores_{a.set_name}.tsv.gz", sep="\t")   # identical to what cached_scores would compute
C.CACHE, C.TABLES = cache, tables
C.PRECLIN = CU.PRECLIN = BR.PRECLIN = reg

# ---- units for this cancer, with the user's organoids as the lab models
resp = pd.read_csv(reg / a.set_name / "response.tsv", sep="\t")
resp["drug"] = resp["drug"].map(common_drug_name)
have = {d: g["sample"].nunique() for d, g in resp.groupby("drug")}
want = [d.strip().upper() for d in a.drugs.split(",") if d.strip()]
units, notes = [], []
for u in UNITS:
    if u["cancer"] != group:
        continue
    drug = u["short"].split(":")[-1]
    if want and drug not in want:
        continue
    if have.get(drug, 0) < 4:
        notes.append((u["id"], f"your organoids have {have.get(drug, 0)} responses for {drug} (need >= 4)"))
        continue
    v = copy.deepcopy(u)
    v["lab"] = [(a.set_name, drug, None)]
    units.append(v)
if not group:
    notes.append(("all", f"the repo has no patient cohorts defined for cancer '{a.cancer}' (units exist for colorectal, bladder, pancreas, liver, ovarian, lung, sarcoma)"))
print(f"{len(units)} units for '{a.cancer}'; not run: {notes}", flush=True)
pd.DataFrame(notes, columns=["unit", "reason"]).to_csv(out / "units_not_run.tsv", sep="\t", index=False)

feats = ["shortlist", "all"] if a.features == "both" else [a.features]
for u in units:
    for f in feats:
        t0 = time.time()
        try:
            D, sc, ev, it, fi, meta, extra = RJ.run_unit(u, f)
        except Exception as e:
            print(u["id"], f, "failed:", type(e).__name__, e, flush=True)
            with open(out / "units_not_run.tsv", "a") as fh:
                fh.write(f"{u['id']}\t{f}: error {type(e).__name__}: {e}\n")
            continue
        tag = f"{u['id']}__{f}"
        ev.assign(unit=u["id"], features=f).to_csv(tables / f"{tag}__eval.tsv", sep="\t", index=False)
        it.assign(unit=u["id"], features=f).to_csv(tables / f"{tag}__interaction.tsv", sep="\t", index=False)
        pd.Series(meta).to_csv(tables / f"{tag}__meta.tsv", sep="\t", header=False)
        print(tag, round(time.time() - t0), "s", flush=True)

# ---- one summary table
rows = []
for f in sorted(tables.glob("*__eval.tsv")):
    e = pd.read_csv(f, sep="\t")
    rows.append(e)
if rows:
    pd.concat(rows).to_csv(out / "joint_summary.tsv", sep="\t", index=False)
print("done")
