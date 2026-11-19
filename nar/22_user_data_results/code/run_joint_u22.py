"""Organoid + patient method (nar/03_organoid_plus_patient code reused unchanged) for the user's lung and ovarian sets, with units defined here.
Patient cohorts come from robust/results/training_table_v2 (20 new lung cohorts included).
Usage: python3 run_joint_u22.py lung|ovary [shortlist|all|both]       (1 process)
Conventions (nar/03): score = predicted benefit (high = longer survival expected); HR per 1 SD < 1 = the expected direction."""
import os, sys, time, copy
os.environ["OMP_NUM_THREADS"] = os.environ["OPENBLAS_NUM_THREADS"] = os.environ["MKL_NUM_THREADS"] = "1"
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
DATA = OUT / "_data"
sys.path.insert(0, str(REPO / "nar" / "03_organoid_plus_patient" / "code"))
import pandas as pd
import common as C
import run_joint as RJ
from obd import bridge as BR, curated as CU

plan = sys.argv[1]
feats_arg = sys.argv[2] if len(sys.argv) > 2 else "shortlist"
SETS = {"lung": "luad_gse165611", "ovary": "ovary_gse229119_noMCF7"}
name = SETS[plan]
run = DATA / "joint" / plan
(run / "tables").mkdir(parents=True, exist_ok=True)
(run / "cache").mkdir(parents=True, exist_ok=True)
reg = DATA / "registry"
for p in (REPO / "data" / "curated" / "preclinical").iterdir():
    if p.is_dir() and not (reg / p.name).exists():
        os.symlink(p, reg / p.name)
for p in (REPO / "nar" / "03_organoid_plus_patient" / "cache").glob("scores_*.tsv.gz"):
    if not (run / "cache" / p.name).exists():
        os.symlink(p, run / "cache" / p.name)
us = pd.read_csv(DATA / "loader_out" / "imports" / name / "pathway_scores.tsv.gz", sep="\t", index_col=0)
us.index = us.index.astype(str)
us.to_csv(run / "cache" / f"scores_{name}.tsv.gz", sep="\t")
C.CACHE, C.TABLES = run / "cache", run / "tables"
C.PRECLIN = CU.PRECLIN = BR.PRECLIN = reg
C.TT = REPO / "robust" / "results" / "training_table_v2"
C._P.clear()

PLAT = r"CISPLATIN|CARBOPLATIN|PLATINUM"
LUNG_PLAT_COHORTS = {"TCGA-LUAD": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy"),
                     "TCGA-LUSC": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy"),
                     "GSE42127": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no adjuvant chemotherapy"),
                     "GSE14814": dict(treated="CISPLATIN", comparison=r"^NONE$", comp_label="observation (randomised, JBR.10)"),
                     "GSE29013": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy"),
                     "GSE68465": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy"),
                     "GSE37745": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy")}
UNITS = []
for d, short in (("GEMCITABINE", "paper:GEMCITABINE"), ("ETOPOSIDE", "ETOPOSIDE"), ("PACLITAXEL", "paper:DOCETAXEL_PACLITAXEL"), ("CAMPTOTHECIN", "paper:CAMPTOTHECIN")):
    UNITS.append(dict(id=f"LUNG_{d}_vs_platinum_proxy", cancer="lung", short="paper:CISPLATIN", lab=[(name, d, None)], align_extra=[], merge={},
                      cohorts=copy.deepcopy(LUNG_PLAT_COHORTS), min_train=(10, 5)))
    if d != "CAMPTOTHECIN":
        UNITS.append(dict(id=f"LUNG_{d}_matched", cancer="lung", short=short, lab=[(name, d, None)], align_extra=[], merge={}, min_train=(10, 5),
                          cohorts={"TCGA-LUAD": dict(treated=d, comparison=r"^NONE$", comp_label="no chemotherapy"),
                                   "TCGA-LUSC": dict(treated=d, comparison=r"^NONE$", comp_label="no chemotherapy")}))
for d in ("OLAPARIB", "NIRAPARIB"):
    UNITS.append(dict(id=f"OV_{d}_vs_platinum", cancer="ovarian", short="paper:CISPLATIN", lab=[(name, d, None)], align_extra=[], merge={}, min_train=(10, 5),
                      cohorts={"TCGA-OV": dict(treated=PLAT, comparison=r"^NONE$", comp_label="no chemotherapy"),
                               "GSE9891": dict(treated="PLATINUM", comparison=None), "GSE63885": dict(treated="PLATINUM", comparison=None),
                               "GSE15622": dict(treated="CARBOPLATIN", comparison="PACLITAXEL", comp_label="paclitaxel alone (randomised)"),
                               "GSE51373": dict(treated="CARBOPLATIN", comparison=None)}))
units = [u for u in UNITS if u["cancer"] == ("lung" if plan == "lung" else "ovarian")]
feats = ["shortlist", "all"] if feats_arg == "both" else [feats_arg]
for u in units:
    for f in feats:
        tag = f"{u['id']}__{f}"
        if (run / "tables" / f"{tag}__eval.tsv").exists():
            continue
        t0 = time.time()
        try:
            D, sc, ev, it, fi, meta, extra = RJ.run_unit(u, f)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(tag, "failed", e, flush=True)
            with open(run / "failed.tsv", "a") as fh:
                fh.write(f"{tag}\t{type(e).__name__}: {e}\n")
            continue
        ev.assign(unit=u["id"], features=f).to_csv(run / "tables" / f"{tag}__eval.tsv", sep="\t", index=False)
        it.assign(unit=u["id"], features=f).to_csv(run / "tables" / f"{tag}__interaction.tsv", sep="\t", index=False)
        fi.assign(unit=u["id"], features=f).to_csv(run / "tables" / f"{tag}__folds.tsv", sep="\t", index=False)
        pd.Series(meta).to_csv(run / "tables" / f"{tag}__meta.tsv", sep="\t", header=False)
        sc.join(D[["unit_cohort", "group", "time", "event", "endpoint"]]).to_csv(run / "cache" / f"{tag}__scores.tsv.gz", sep="\t")
        print(tag, round(time.time() - t0), "s", flush=True)
print("done", plan)
