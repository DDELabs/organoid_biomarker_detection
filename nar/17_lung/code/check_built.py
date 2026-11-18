"""Compare the built lung cohorts with the numbers expected in nar/16_plan/_verify (arm sizes, events)."""
import sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "robust"))
import build_training_cohorts as BC
EXPECT = {  # accession: (n expected in plan, arms expected, events expected in plan)
 "GSE14814": (133, "ACT 71 / OBS 62", "OS 60"), "GSE68465": (443, "chemo yes 89 / no 341 / unknown 13", "OS 236; relapse 205"),
 "GSE37745": (196, "adjuvant yes 29 / no 71 / unknown 96", "OS 145; recurrence 49"), "GSE29013": (55, "adjuvant yes 34 / no 21", "OS 18; progression 28"),
 "GSE31210": (226, "(22 flagged are removed -> 204)", "relapse 64 / OS 35 in 226"), "GSE50081": (181, "not recorded", "OS 75; DFS 51"),
 "GSE8894": (138, "not recorded", "RFS 69"), "GSE3141": (111, "not recorded", "OS (text)"), "GSE30219": (296, "not recorded", "OS 200; DFS 114"),
 "GSE72094": (442, "not recorded", "OS 122"), "GSE41271": (275, "not recorded", "OS 124; recurrence 122"), "GSE157011": (484, "not recorded", "OS 264"),
 "GSE33072": (131, "erlotinib 25 / sorafenib 39 / erl+bex 1", "8-week control, PFS (sparse)"), "GSE161537": (82, "anti-PD-1", "CR1 PR19 SD22 PD34 NA6; OS 61"),
 "GSE162520": (92, "not recorded", "OS 21; PFS 35"), "GSE135222": (27, "anti-PD-1", "PFS 21"), "GSE126044": (16, "anti-PD-1", "responder 5 / non 11"),
 "GSE218989": (355, "PD-1/PD-L1", "responder 168 / non 187"), "GSE207422": (39, "anti-PD-1 + chemo", "bulk pre-treatment 24"), "GSE81089": (199, "not recorded", "dead 93")}
rows = []
for fn in BC.LUNG_BUILDERS:
    g = fn.__name__
    c = pd.read_csv(BC.OUT / g / "clinical.tsv", sep="\t", index_col=0)
    dr = c["drugs"].fillna("(not recorded)").replace({"NONE": "(none)"})
    ev = {k: int(c[k].sum()) for k in ("os_event", "rfs_event", "dfs_event", "pfs_event") if k in c and c[k].notna().any()}
    nev = {k: int(c[k].notna().sum()) for k in ev}
    rows.append(dict(cohort=g, built_n=len(c), plan_n=EXPECT[g][0], arms_built="; ".join(f"{k}: {v}" for k, v in dr.value_counts().items()),
                     arms_expected=EXPECT[g][1], events_built="; ".join(f"{k}={v} (of {nev[k]})" for k, v in ev.items()),
                     events_expected=EXPECT[g][2], responders=f"{int(c['responder'].sum())} of {int(c['responder'].notna().sum())}" if "responder" in c else ""))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "nar/17_lung/tables/lung_cohorts_check.tsv", sep="\t", index=False)
pd.set_option("display.width", 300, "display.max_colwidth", 120)
print(T.to_string())
