"""Q2: stability selection. L1 logistic (C=0.1, the value chosen by inner CV in W3) on bootstrap resamples.
(A) 100 bootstraps of all 3,511 patients (resampled inside each cohort). (B) per LOO fold: 20 bootstraps of the 12 training cohorts.
Also L2 weights refit in each leave-one-cohort-out fold."""
import sys
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from common import *
d, Z = breast_data(); y = d["responder"].astype(int).values; coh = d["cohort"].values; cols = list(Z.columns); X = Z.values
cohorts = sorted(set(coh)); C = 0.1

def boot_idx(mask, seed):
    r = np.random.default_rng(seed); out = []
    for c in cohorts:
        i = np.where(mask & (coh == c))[0]
        if len(i): out.append(r.choice(i, len(i), replace=True))
    return np.concatenate(out)

def l1(idx):
    m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=500).fit(X[idx], y[idx])
    return m.coef_[0]

jobs = [("all", 0, boot_idx(np.ones(len(y), bool), SEED + b)) for b in range(100)]
for c in cohorts:
    for b in range(20):
        jobs.append((c, 1, boot_idx(coh != c, SEED + 1000 + 37 * b + len(c))))
res = Parallel(n_jobs=2)(delayed(l1)(j[2]) for j in jobs)
B = np.array(res)
lab = np.array([j[0] for j in jobs])
np.save(TAB / "_l1_boot_coefs.npy", B); pd.Series(lab).to_csv(TAB / "_l1_boot_labels.tsv", sep="\t", index=False)
# L2 weights, full and per fold
w = TC.fit(X, y).coef_[0]
W = np.array([TC.fit(X[coh != c], y[coh != c]).coef_[0] for c in cohorts])
pd.DataFrame(W.T, index=cols, columns=cohorts).assign(full=w).to_csv(TAB / "l2_weights_per_fold.tsv", sep="\t")
print("done", B.shape, (B[lab == "all"] != 0).sum(1).mean())
