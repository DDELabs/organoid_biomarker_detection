"""Baselines (same folds): ER+HER2 only, proliferation score, 200 random pathway sets of shortlist size (L2, lambda=100)."""
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import *
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
PRED = TAB / "preds"
d, Z, P = breast_data()
y = d["responder"].astype(int); coh = d["cohort"]
flags = pd.read_csv(TAB / "pathway_flags.tsv", sep="\t", index_col=0)
cols = list(Z.columns)
SL = [c for c in cols if flags.loc[c, "shortlist_paper"]]; PR = [c for c in cols if flags.loc[c, "proliferation"]]
pd.Series(PR, name="proliferation_pathways").to_csv(TAB / "proliferation_pathways_used.tsv", index=False)
er = d[["er", "her2"]]
o_er, o_pr, hp = [], [], []
for c in sorted(coh.unique()):
    te = coh == c; tr = ~te
    # ER + HER2 only: trained on training rows with both known; scored on test rows with both known
    a, b = er[tr].dropna(), er[te].dropna()
    if len(b) > 20 and b.nunique().max() > 1 and len(a) > 50:
        m = LogisticRegression(C=1e3).fit(a.values, y[a.index])
        o_er.append(pd.DataFrame(dict(cohort=c, score=m.decision_function(b.values), y=y[b.index].values), index=b.index))
    # proliferation score: mean of the z-scored cell-cycle pathways; sign chosen from training cohorts only
    s = Z[PR].mean(axis=1)
    au = [roc_auc_score(y[(coh == k)], s[(coh == k)]) for k in coh[tr].unique()]
    sign = 1.0 if np.mean(au) >= 0.5 else -1.0
    o_pr.append(pd.DataFrame(dict(cohort=c, score=sign * s[te].values, y=y[te].values), index=d.index[te]))
    hp.append(dict(test_cohort=c, proliferation_sign=sign, train_mean_auc=float(np.mean(au))))
pd.concat(o_er).to_csv(PRED / "ER_HER2_only.tsv.gz", sep="\t"); pd.concat(o_pr).to_csv(PRED / "Proliferation_score.tsv.gz", sep="\t")
pd.DataFrame(hp).to_csv(PRED / "Proliferation_score_hp.tsv", sep="\t", index=False)
# random sets
rng = np.random.default_rng(SEED); rows = []
sets = [list(rng.choice(cols, len(SL), replace=False)) for _ in range(200)]
for c in sorted(coh.unique()):
    te = (coh == c).values; tr = ~te
    for i, st in enumerate(sets):
        m = LogisticRegression(C=0.01, max_iter=3000).fit(Z.loc[tr, st].values, y[tr])
        rows.append(dict(test_cohort=c, draw=i, auc=roc_auc_score(y[te], m.decision_function(Z.loc[te, st].values))))
    print(c, flush=True)
pd.DataFrame(rows).to_csv(TAB / "random_pathway_sets_auc.tsv", sep="\t", index=False)
