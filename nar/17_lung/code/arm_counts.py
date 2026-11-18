"""Arms and events per lung cohort straight from training_table_v2 (feeds the REPORT)."""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[3]
P = pd.read_csv(ROOT / "robust/results/training_table_v2/patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
L = P[P["cancer_group"] == "lung"].copy()
L["arm_label"] = L["drugs"].fillna("(treatment not recorded)").replace({"NONE": "(no drug)"})
g = L.groupby(["cohort", "arm_label"], dropna=False)
T = g.agg(patients=("treated", "size"), with_response=("responder", "count"), responders=("responder", "sum"),
          os_known=("os_event", "count"), os_events=("os_event", "sum"), rfs_known=("rfs_event", "count"), rfs_events=("rfs_event", "sum"),
          dfs_known=("dfs_event", "count"), dfs_events=("dfs_event", "sum"), pfs_known=("pfs_event", "count"), pfs_events=("pfs_event", "sum"),
          flagged_duplicates=("possible_duplicate_of", "count")).reset_index()
T.to_csv(ROOT / "nar/17_lung/tables/lung_arm_counts.tsv", sep="\t", index=False)
print(T[T.cohort.isin(["GSE14814", "GSE33072", "GSE68465", "GSE161537"])].to_string())
