import sys; sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
d, Z, P = breast_data()
y = d["responder"].astype(int)
old = pd.read_csv(ROOT/"robust/results/train_combined/part1_response_leave_one_cohort_out.tsv", sep="\t")
old = old[old.cancer == "breast"].set_index("test_cohort")
rows = []
for c in sorted(d["cohort"].unique()):
    tr, te = d.index[d["cohort"] != c], d.index[d["cohort"] == c]
    a, p = TC.auc(y.loc[te], TC.fit(Z.loc[tr].values, y.loc[tr]).decision_function(Z.loc[te].values))
    rows.append(dict(test_cohort=c, n=len(te), auc_here=a, auc_old=old.loc[c, "auc_pathways"], diff=a - old.loc[c, "auc_pathways"],
                     match_3dp=round(a, 3) == round(old.loc[c, "auc_pathways"], 3)))
R = pd.DataFrame(rows); print(R.round(4).to_string(), len(R), R.match_3dp.all(), Z.shape)
R.to_csv(TAB / "l2_reproduction_check.tsv", sep="\t", index=False)
