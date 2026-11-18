"""Write tables/lung_data_files.tsv (file, size, sha256) and tables/lung_cohorts_built.tsv (built / skipped per accession)."""
import hashlib, sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "robust"))
import build_training_cohorts as BC
RAW = BC.RAW
files = []
used = {}
for fn in BC.LUNG_BUILDERS:
    g = fn.__name__
    for p in sorted((RAW / g).glob("*")):
        if p.is_file() and (p.name.endswith(("series_matrix.txt.gz", "log2cpm.csv.gz", "exp.tsv.gz", "counts.txt.gz", "pc_rmdupli.txt.gz",
                                             "log2TPM.txt.gz", "metadata.xlsx", "cufflinks.tsv.gz"))):
            files.append(dict(accession=g, file=p.name, bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(files).to_csv(ROOT / "nar/17_lung/tables/lung_data_files.tsv", sep="\t", index=False)
chk = pd.read_csv(ROOT / "nar/17_lung/tables/lung_cohorts_check.tsv", sep="\t").set_index("cohort")
plan = pd.read_csv(ROOT / "nar/16_plan/COHORTS_TO_ADD.tsv", sep="\t")
plan = plan[plan["cancer"].str.contains("lung") | plan["accession"].str.contains("CPTAC|POPLAR")].set_index("accession")
notes = {}
for fn in BC.LUNG_BUILDERS:
    clin = pd.read_csv(BC.OUT / fn.__name__ / "clinical.tsv", sep="\t", nrows=1)
    notes[fn.__name__] = (BC.OUT / fn.__name__ / "NOTE.txt").read_text().strip()
rows = []
for acc, r in plan.iterrows():
    if acc in chk.index:
        c = chk.loc[acc]
        rows.append(dict(accession=acc, priority=r["priority"], status="built", samples_built=c["built_n"], arms_built=c["arms_built"], arms_expected=c["arms_expected"],
                         events_built=c["events_built"], events_expected=c["events_expected"], responders=c["responders"], note=notes[acc]))
    else:
        why = BC.SKIPPED_LUNG.get(acc.split(" ")[0]) or ("data not reachable from this machine (EGA, controlled access)" if acc.startswith("POPLAR") else None)
        rows.append(dict(accession=acc, priority=r["priority"], status="skipped" if why else "not attempted", reason=why or
                         ("priority 4 or excluded in the plan (W1 covers priority 1-3, plus GSE3141, GSE81089 and GSE126044 named in the steps)" if str(r["priority"]) in ("4", "-") else "")))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "nar/17_lung/tables/lung_cohorts_built.tsv", sep="\t", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 70)
print(T[["accession", "priority", "status", "samples_built", "reason"]].to_string())
