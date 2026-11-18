"""Audit fix I07: bladder organoids, one line per patient.

The paper's 9 bladder organoid lines with cisplatin response come from 7 patients: SCBO-3_2 and SCBO-11_2 are
recurrence lines of the patients of SCBO-3 and SCBO-11. Sensitivity: rerun the paper-source-data cisplatin analysis
(nar/02 BLCA_OV_LUAD_OTHER run_blca.py, version 'paper') without the two recurrence lines (7 lines, 7 patients).
Fewer permutations / bootstraps than the main run (500 / 2,000) to save CPU.
Writes nar/02_organoid_only/BLCA_lee2018/sensitivity_one_line_per_patient/ (engine outputs) and
nar/02_organoid_only/BLCA_lee2018/tables/BLCA_cisplatin_one_line_per_patient.tsv
"""
import sys
from pathlib import Path
import pandas as pd

CODE = Path(__file__).resolve().parents[2] / "02_organoid_only" / "BLCA_OV_LUAD_OTHER" / "code"
sys.path.insert(0, str(CODE))
import common as K  # noqa: E402
import engine as E  # noqa: E402
import run_blca as B  # noqa: E402

osc, tsc, surv = B.paper_version()
anyd = K.any_drug("BLCA")
z, src = K.proximity("CISPLATIN")
y_all = B.response("CISPLATIN")
rows = []
for lab, drop in (("all 9 lines (7 patients), as main run", []), ("7 lines, one per patient (recurrence lines SCBO-3_2, SCBO-11_2 dropped)", ["SCBO-3_2", "SCBO-11_2"])):
    y = y_all.drop([d for d in drop if d in y_all.index])
    a = E.Analysis(set_id="BLCA_lee2018_paper_onepp" if drop else "BLCA_lee2018_paper_all9", set_label=lab, project="BLCA",
                   drug="CISPLATIN", model_drug="CISPLATIN", osc=osc, tsc=tsc, y=y, surv=surv,
                   treated=K.treated("BLCA", "CISPLATIN"), untreated=set(surv.index) - anyd, z=z, z_source=src,
                   outdir=K.PARENT / "BLCA_lee2018" / "sensitivity_one_line_per_patient" / ("onepp" if drop else "all9"),
                   response_label="AUC of fitted viability curve (CoderData)", n_perm=500, n_boot=2000, isogenic=None,
                   note="audit I07 sensitivity", seed=0)
    r = a.run()
    rows.append({"analysis": lab, "lines": ";".join(a.samples), **{k: r.get(k) for k in (
        "n_models", "n_proximal", "biomarker", "treated_n", "treated_events", "treated_p", "treated_hr", "treated_hr_lo",
        "treated_hr_hi", "untreated_n", "untreated_p", "perm_empirical_p", "boot_empirical_p_rank", "boot_empirical_p_rank1p")}})
    print(rows[-1], flush=True)
T = pd.DataFrame(rows)
out = K.PARENT / "BLCA_lee2018" / "tables" / "BLCA_cisplatin_one_line_per_patient.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
T.to_csv(out, sep="\t", index=False)
print(T.to_string())
