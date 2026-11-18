"""Paper method (organoid-only) for a user's organoid set, reusing nar/02_organoid_only/code unchanged.

Usage (called by import_organoids.py, one process):
  python3 run_organoid_only_user.py <set_dir> <set_name> <cancer> <out_dir> [--drugs A,B] [--fast]
For every user drug that has a network-proximity entry in the paper (or the platinum proxy used in nar/02) and every
TCGA cohort of the cancer: score pathways with ssGSEA, pick pathways close to the drug targets (z <= -1.28), train
Ridge/SVR/LinearRegression on the organoids, split TCGA patients treated with the drug at the median prediction,
log-rank + Cox.  Everything is written under <out_dir>; nar/02 is only read.
"""
import argparse
import os
import pickle
import sys
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
ap.add_argument("--fast", action="store_true", help="1000 bootstraps / 100 permutations instead of 10000 / 1000")
a = ap.parse_args()
if a.fast:
    os.environ["N_BOOT"], os.environ["N_PERM"] = "1000", "100"

sys.path.insert(0, str(AC.NAR / "02_organoid_only" / "code"))
import pandas as pd  # noqa: E402
import common as C  # noqa: E402  (nar/02 code)
import run_all as RA  # noqa: E402
from obd import scoring  # noqa: E402

out = Path(a.out_dir)
out.mkdir(parents=True, exist_ok=True)
projects = AC.tcga_projects(a.cancer)
if not projects:
    sys.exit(f"No TCGA project is known for cancer '{a.cancer}'.")

# ---- redirect everything this code writes into our own folder (nar/02 is only read)
raw = AC.make_raw_cache()
reg = AC.make_registry(a.set_dir, out, a.set_name)
C.RAW = raw
C.C.EXTERNAL = raw
C.SCORES = out / "scores"
C.PRECLIN = reg
RA.OUT = out / "results"
_ssgsea = scoring.ssgsea
scoring.ssgsea = lambda expr, gs, min_size=2, threads=1: _ssgsea(expr, gs, min_size=min_size, threads=1)   # 1 CPU only

key = a.set_name
C.SETS[key] = (key, a.set_name, projects, "response", f"{a.set_name} ({a.cancer})")
RA.SETS = C.SETS
direction = pd.read_csv(reg / a.set_name / "response.tsv", sep="\t", usecols=["direction"], nrows=1)["direction"].iloc[0]
_set_response = C.set_response


def set_response(k):                                  # the nar/02 code expects lower = more sensitive
    r = _set_response(k)
    return -r if (k == key and str(direction).startswith("higher")) else r


C.set_response = set_response
RA.set_response = set_response

# ---- which (drug, cohort) combinations
resp = C.set_response(key)
want = [d.strip().upper() for d in a.drugs.split(",") if d.strip()]
cohorts = list(projects) + (["COADREAD"] if {"COAD", "READ"} <= set(projects) else [])
specs, skipped = [], []
for d in resp.columns:
    if want and d not in want:
        continue
    pdrug = C.DRUG_MAP.get(d, d)
    z, _ = C.proximity(pdrug)
    if z is None:
        skipped.append((d, "no network-proximity entry in the paper's table"))
        continue
    if resp[d].notna().sum() < 5:
        skipped.append((d, f"only {int(resp[d].notna().sum())} organoids with a response (< 5)"))
        continue
    for c in cohorts:
        specs.append((key, c, d, pdrug))
print(f"{len(specs)} analyses planned for {len({s[2] for s in specs})} drugs x {len(cohorts)} cohorts; {len(skipped)} drugs not analysed", flush=True)

cd = C.cdr()
rows = []
for s in specs:
    try:
        tag, meta = RA.run(s, cd)
    except Exception as e:                           # e.g. a TCGA table cannot be downloaded
        meta = dict(set=key, cohort=s[1], model_drug=s[2], patient_drug=s[3], skip=f"error: {type(e).__name__}: {e}")
        tag = f"{key}__{s[1]}__{s[2]}"
    rows.append(meta)
    print(tag, meta.get("skip", "done"), meta.get("n_treated"), flush=True)
pd.DataFrame(rows).to_csv(out / "run_log.tsv", sep="\t", index=False)
pd.DataFrame(skipped, columns=["drug", "reason"]).to_csv(out / "drugs_not_analysed.tsv", sep="\t", index=False)

# ---- summary table (every number comes from the result pickles)
summ = []
for f in sorted((out / "results").glob("*.pkl")):
    res = pickle.load(open(f, "rb"))
    m, ev = res["meta"], res["ev"]
    r1 = ev[(ev["model"] == "Ridge") & (ev["k"] == 1)].iloc[0]
    loo = res["loo"]["Ridge"]
    summ.append(dict(set=m["set"], cohort=m["cohort"], drug=m["model_drug"], patient_drug=m["patient_drug"], n_organoids=m["n_models"],
                     n_close_pathways=m["n_proximal"], n_treated_patients=m["n_treated"], treated_events=m["treated_events"],
                     n_untreated_patients=m["n_untreated"], ridge_top_pathway=r1["pathways"], HR=r1.get("treated_HR"),
                     HR_low=r1.get("treated_HR_lo"), HR_high=r1.get("treated_HR_hi"), logrank_p=r1["treated_p"],
                     untreated_logrank_p=r1.get("untreated_p"), bootstrap_empirical_p=res["boot_emp_p"], permutation_p=res["perm_rank1_emp_p"],
                     organoid_leave_one_out_spearman=loo[1], organoid_leave_one_out_p=loo[2], borderline_small_n=m.get("borderline")))
pd.DataFrame(summ).to_csv(out / "organoid_only_summary.tsv", sep="\t", index=False)
print("summary rows:", len(summ))
