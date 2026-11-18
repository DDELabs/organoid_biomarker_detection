"""Audit fix I23: redraw nar/10 SuppFig4 (cross-test matrix) at print size with readable fonts.
Reads nar/10_completeness/tables/SuppFig4_cross_test_matrix.tsv (numbers unchanged); same caption as c04."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
D10 = Path(__file__).resolve().parents[2] / "10_completeness"
sys.path.insert(0, str(D10 / "code"))
import common as K  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

X = pd.read_csv(D10 / "tables" / "SuppFig4_cross_test_matrix.tsv", sep="\t")
bio = list(dict.fromkeys(X.biomarker)); keys = list(dict.fromkeys(X.cohort))
P = X.pivot(index="biomarker", columns="cohort", values="p").reindex(index=bio, columns=keys)
pw = X.drop_duplicates("biomarker").set_index("biomarker")["pathway"]
lab = {k: K.REG[k][0] if k in K.REG else "Bladder, paper's amino acid pathway | cisplatin" for k in bio}
fig, ax = plt.subplots(figsize=(7.2, 7.4))
im = ax.imshow(-np.log10(P.values.astype(float)), cmap="Blues", vmin=0, vmax=3)
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        p = P.values[i, j]
        if np.isfinite(p):
            own = (P.index[i] == P.columns[j]) or (P.index[i] == "BLCA_paper_AMINO" and P.columns[j] == "BLCA_paper_CIS")
            ax.text(j, i, f"{p:.2f}" if p >= 0.01 else f"{p:.0e}".replace("e-0", "e-"), ha="center", va="center", fontsize=6,
                    color="white" if p < 0.01 else "black", fontweight="bold" if own else "normal")
            if own:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=K.OI[3], lw=1.5))
ax.set_xticks(range(len(keys))); ax.set_xticklabels([K.REG[k][0] for k in keys], rotation=70, ha="right", fontsize=7)
ax.set_yticks(range(len(bio))); ax.set_yticklabels([f"{lab[k]}: {K.short(pw[k], 26)}" for k in bio], fontsize=7)
ax.set_xlabel("treated TCGA patients (cohort | drug)", fontsize=8); ax.set_ylabel("biomarker (set | drug: rank-1 pathway)", fontsize=8)
cb = fig.colorbar(im, ax=ax, shrink=0.45); cb.set_label("-log10 log-rank p", fontsize=8); cb.ax.tick_params(labelsize=7)
cap = (D10 / "figures" / "SuppFig4_cross_test_matrix_all_sets.caption.md").read_text()
cap = cap.rstrip() + "\nRedrawn at print size (7.2 in wide) by nar/15_fixes/code/f10_suppfig4.py (audit I23); numbers unchanged.\n" \
    if "Redrawn at print size" not in cap else cap
K.save(fig, "SuppFig4_cross_test_matrix_all_sets", cap)
print("ok")
