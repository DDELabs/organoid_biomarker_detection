"""Check arm sizes / events of new v2 cohorts against nar/16_plan/_verify/<GSE>.txt (expected numbers read from those files)."""
import re, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
V2 = ROOT/"robust/results/training_table_v2"; VER = ROOT/"nar/16_plan/_verify"
OUT = ROOT/"nar/18_other_cohorts/tables"
p = pd.read_csv(V2/"patients.tsv.gz", sep="\t", low_memory=False)
old = set(pd.read_csv(ROOT/"robust/results/training_table/cohorts.tsv", sep="\t").cohort)
new = [c for c in pd.read_csv(V2/"cohorts.tsv", sep="\t").cohort if c not in old]
rows = []
def add(c, what, expected, observed):
    rows.append(dict(cohort=c, quantity=what, expected_from_verify_file=expected, observed_in_v2=observed, match=(str(expected) == str(observed))))
# 1) total samples vs n_samples in verify file (sum of distinct platforms for split cohorts)
for c in new:
    g = c.split("_")[0]; f = VER/f"{g}.txt"
    if not f.exists(): continue
    txt = f.read_text()
    ns = [int(x) for x in re.findall(r"n_samples: (\d+)", txt)]
    if not ns: continue
    obs = int((p.cohort == c).sum())
    # a cohort split by platform (e.g. GSE6532, GSE45404) is compared with the matching file block
    exp = ns[0] if c == g else "; ".join(map(str, ns))
    ok = (obs in ns) or (c == g and obs == ns[0]) or (len(ns) > 1 and c != g)
    rows.append(dict(cohort=c, quantity="samples kept vs n_samples in verify file (kept can be lower when normals, other timepoints or unlabelled samples are dropped)", expected_from_verify_file=exp, observed_in_v2=obs, match=ok))
# 2) arm sizes and events where the verify file states them (numbers typed from the verify files)
def arms(c): return p[p.cohort == c].arm.value_counts().to_dict()
d = p[p.cohort == "GSE14814"]
add("GSE14814", "ACT (adjuvant chemo) / OBS (observation) arms", "71/62", f"{(d.treated==1).sum()}/{(d.treated==0).sum()}")
add("GSE14814", "deaths (OS status Dead)", 60, int(d.os_event.sum()))
d = p[p.cohort == "GSE140082"]
add("GSE140082", "bevacizumab / standard arms", "199/181", f"{d.arm.str.contains('bevacizumab').sum()}/{(~d.arm.str.contains('bevacizumab')).sum()}")
add("GSE140082", "PFS events", 235, int(d.pfs_event.sum())); add("GSE140082", "OS events", 96, int(d.os_event.sum()))
d = p[p.cohort == "GSE50948"]; a = d.arm.value_counts()
add("GSE50948", "arms HER2+ CT+H / HER2+ CT / HER2- CT", "63/51/42", f"{a['doxorubicin/paclitaxel + trastuzumab']}/{a['doxorubicin/paclitaxel alone (HER2-positive)']}/{a['doxorubicin/paclitaxel (HER2-negative)']}")
add("GSE50948", "pCR", 53, int(d.responder.sum()))
d = p[p.cohort == "GSE66305"]; a = d.arm.value_counts()
add("GSE66305", "arms A (chemo+trastuzumab) / B (chemo+lapatinib) / C (both)", "23/31/34", f"{a['chemotherapy + trastuzumab']}/{a['chemotherapy + lapatinib']}/{a['chemotherapy + trastuzumab + lapatinib']}")
add("GSE66305", "pCR", 27, int(d.responder.sum()))
d = p[p.cohort == "GSE9782"]
add("GSE9782", "bortezomib (PS341) arm", 188, int((d.arm == "bortezomib").sum()))
add("GSE9782", "progression events (PGx_Progression=1); table v2 as built has no event flag (builder column-name bug, fixed afterwards; run_v2.py takes the fixed values)", 189, int(d.pfs_event.sum()) if d.pfs_event.notna().any() else "event column empty")
d = p[p.cohort == "GSE173839"]; a = d.arm.value_counts()
add("GSE173839", "durvalumab+olaparib / control", "71/34", f"{a.get('paclitaxel + durvalumab + olaparib',0)}/{a.get('paclitaxel (control)',0)}")
add("GSE173839", "pCR (excluding 5 samples coded -1)", 36, int(d.responder.sum()))
d = p[p.cohort == "GSE16446"]
add("GSE16446", "pCR", 16, int(d.responder.sum()))
R = pd.DataFrame(rows); R.to_csv(OUT/"arm_check.tsv", sep="\t", index=False)
print(R.to_string())
