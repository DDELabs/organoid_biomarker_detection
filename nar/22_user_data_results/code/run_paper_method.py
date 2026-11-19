"""Paper method (organoid-only, nar/02 engine reused unchanged) for the user's sets, with extra network proximity.
Usage: python3 run_paper_method.py <plan_name>      (1 process, 1 CPU)
plans: lung | ovary
Every analysis is saved as a pickle under _data/results/<plan>/ and summarised in tables/paper_method_<plan>.tsv
"""
import os, sys, pickle, time, json
os.environ.setdefault("N_BOOT", "500")
os.environ.setdefault("N_PERM", "300")
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u22
from u22 import C, DATA, OUT, N, pd, np
import run_all as RA
import engine as EN

plan = sys.argv[1]
only = sys.argv[2:]     # optional: restrict to tags containing any of these strings
IMP = DATA / "loader_out" / "imports"
REG = DATA / "registry"
REG.mkdir(parents=True, exist_ok=True)
for s in ("luad_gse165611", "luad_gse165611_listed", "ovary_gse229119", "ovary_gse229119_noMCF7"):
    u22.link(IMP / s, REG / s)
C.PRECLIN = REG
RES = DATA / "results" / plan
RES.mkdir(parents=True, exist_ok=True)
RA.proximity = u22.proximity
C.proximity = u22.proximity
PROX_FOR = {}           # patient drug -> proximity drug override (PARP drug for ovary)
_prox = u22.proximity
RA.proximity = lambda d: _prox(PROX_FOR.get(d, d))

LUNG = ["LUAD", "LUSC"]
PLANS = {
    "lung": dict(set="luad_gse165611", projects=LUNG, label="Lung patient-derived cells (GSE165611), cluster-level response",
        specs=[("GEMCITABINE", "GEMCITABINE"), ("PACLITAXEL", "PACLITAXEL"), ("PACLITAXEL", "DOCETAXEL"), ("ETOPOSIDE", "ETOPOSIDE"),
               ("GEMCITABINE", "CISPLATIN"), ("GEMCITABINE", "CARBOPLATIN"), ("ETOPOSIDE", "CISPLATIN"), ("ETOPOSIDE", "CARBOPLATIN")]),
    "lung_listed": dict(set="luad_gse165611_listed", projects=LUNG, label="Lung (GSE165611), only clusters with a listed value",
        specs=[("GEMCITABINE", "GEMCITABINE"), ("GEMCITABINE", "CISPLATIN"), ("GEMCITABINE", "CARBOPLATIN")]),
    "ovary": dict(set="ovary_gse229119_noMCF7", projects=["OV"], label="Ovarian lines + organoids (GSE229119), 6 lines",
        specs=[("OLAPARIB", "CARBOPLATIN"), ("NIRAPARIB", "CARBOPLATIN"), ("OLAPARIB", "PACLITAXEL"), ("NIRAPARIB", "PACLITAXEL")]),
    "ovary_mcf7": dict(set="ovary_gse229119", projects=["OV"], label="Ovarian lines + organoids + MCF7 (breast), 7 lines",
        specs=[("OLAPARIB", "CARBOPLATIN"), ("NIRAPARIB", "CARBOPLATIN")]),
}
P = PLANS[plan]
key = P["set"]
C.SETS[key] = (key, key, P["projects"], "response", P["label"])
RA.SETS = C.SETS
cd = C.cdr()
rows = []
for cohort_all in P["projects"]:
    pass
cohorts = P["projects"]
for mdrug, pdrug in P["specs"]:
    for cohort in cohorts:
        tag0 = f"{cohort}__{mdrug}__{pdrug}"
        if only and not any(o in tag0 for o in only):
            continue
        out_f = RES / f"{tag0}.pkl"
        if out_f.exists():
            continue
        t0 = time.time()
        # ovary: the PARP drug's own targets give the pathway shortlist, patients get platinum / taxane
        if plan.startswith("ovary"):
            PROX_FOR[pdrug] = mdrug
        else:
            PROX_FOR.pop(pdrug, None)
        RA.OUT = RES / "_tmp"
        try:
            tag, meta = RA.run((key, cohort, mdrug, pdrug), cd)
        except Exception as e:
            meta = dict(set=key, cohort=cohort, model_drug=mdrug, patient_drug=pdrug, skip=f"error {type(e).__name__}: {e}")
            tag = None
        meta["proximity_for"] = PROX_FOR.get(pdrug, C.PROX_PROXY.get(pdrug, pdrug))
        meta["seconds"] = round(time.time() - t0)
        if tag and (RA.OUT / f"{tag}.pkl").exists():
            os.replace(RA.OUT / f"{tag}.pkl", out_f)
        else:
            pickle.dump(dict(meta=meta), open(out_f, "wb"))
        print(tag0, meta.get("skip", "done"), meta.get("n_treated"), meta["seconds"], "s", flush=True)
print("finished", plan)
