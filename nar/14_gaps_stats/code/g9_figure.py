"""G9 figure from the saved permutation results."""
import numpy as np
import pandas as pd
from scipy.stats import binom

from g_common import DATA, OK, TAB, C, save_fig

rt = np.load(DATA / "g9_perm_treated_p.npy")
ru = np.load(DATA / "g9_perm_untreated_p.npy")
S = pd.read_csv(TAB / "G9_permutation_null_summary.tsv", sep="\t")
ct = np.nansum(rt < .05, 1)
cu = np.nansum(ru < .05, 1)
mn = np.nanmin(rt, 1)
plt = C.plot_style()
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
k = np.arange(0, 16)
for a, cnt, obs, ttl in ((ax[0], ct, 7, "treated patients"), (ax[1], cu, 12, "untreated patients (no drug record)")):
    a.bar(k, [np.mean(cnt == i) for i in k], color=OK["sky"], label="permutation null (5,000)")
    a.plot(k, binom.pmf(k, 57, .05), "o-", color=OK["black"], ms=3, lw=.8, label="binomial (independent tests)")
    a.axvline(obs, color=OK["verm"], lw=2, label=f"observed = {obs}")
    a.set_xlabel(f"number of the 57 tests with p < 0.05\n({ttl})")
    a.set_ylabel("probability")
    a.legend(frameon=False, fontsize=7)
C.letter(ax[0], "a"); C.letter(ax[1], "b")
ax[2].hist(-np.log10(mn), bins=40, color=OK["sky"])
ax[2].axvline(-np.log10(0.001098), color=OK["verm"], lw=2, label="observed min p = 0.0011")
ax[2].set_xlabel("-log10 of the smallest of 57 treated p-values")
ax[2].set_ylabel("permutations")
ax[2].legend(frameon=False, fontsize=7)
C.letter(ax[2], "c")
fig.tight_layout()
s = S.set_index("statistic")
pt = s.loc["treated: count p<0.05 (of 57)"]
pu = s.loc["untreated: count p<0.05 (of 57)"]
pm = s.loc["treated: minimum p"]
cap = f"""**Figure G9. How many small p-values would the 57-analysis survey give by chance, taking shared patients into account?**
We shuffled survival (time and death status together) between patients inside each TCGA cohort, 5,000 times. Treated
and untreated patients were shuffled separately. Colon, rectum and colon+rectum analyses share patients, so they were
shuffled together, with one shuffle per cohort used for every analysis. The organoid scores were not changed. This
keeps the overlap between analyses but breaks any real link between score and survival.
(a) Treated patients: blue bars = how often each count of p < 0.05 appears under shuffling; black dots = the simple
binomial count that assumes independent tests; red line = observed 7. Shuffle p = {pt.perm_p_at_least_observed:.3f}.
(b) The same for untreated patients: observed 12, shuffle p = {pu.perm_p_at_least_observed:.4f}.
(c) Smallest of the 57 treated p-values under shuffling; observed 0.0011, shuffle p = {pm.perm_p_at_least_observed:.3f}
(this is the family-wise adjusted p of the best single analysis).
"""
save_fig(fig, "G9_permutation_null", cap)
