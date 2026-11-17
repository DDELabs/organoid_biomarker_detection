"""Network proximity (paper's STRING>700 network from MOESM3) for drugs that are not in the paper's z-score table,
plus cisplatin as a check that our re-computation agrees with the paper's own table."""
import time
import numpy as np, pandas as pd
import common as K

TARGETS = {"CISPLATIN_CHECK": ["A2M", "ATOX1", "MPG", "TF"],
           "PEMETREXED": None, "ETOPOSIDE": None, "VINORELBINE": None,
           "TAMOXIFEN": ["ESR1", "ESR2"], "TRASTUZUMAB": ["ERBB2"], "DACARBAZINE": None,
           "OLAPARIB": ["PARP1", "PARP2", "PARP3"]}
for d, t in TARGETS.items():
    t0 = time.time()
    z, src = K.proximity(d, t if t else K.drugbank_targets()[d])
    print(d, src, len(z), round(time.time() - t0), "s", flush=True)
P = K.paper_proximity()["CISPLATIN"]
z = pd.read_csv(K.CACHE / "proximity_papernet_CISPLATIN_CHECK.tsv", sep="\t", index_col=0)["z"]
c = pd.concat([P, z], axis=1, join="inner").dropna()
from scipy.stats import spearmanr, pearsonr
pd.DataFrame([{"pathways": len(c), "pearson_r": pearsonr(c.iloc[:, 0], c.iloc[:, 1])[0],
               "spearman_r": spearmanr(c.iloc[:, 0], c.iloc[:, 1])[0],
               "proximal_paper": int((c.iloc[:, 0] <= -1.2816).sum()), "proximal_ours": int((c.iloc[:, 1] <= -1.2816).sum()),
               "proximal_both": int(((c.iloc[:, 0] <= -1.2816) & (c.iloc[:, 1] <= -1.2816)).sum())}]).to_csv(
    K.BASE / "SUMMARY" / "tables" / "proximity_recomputation_check_cisplatin.tsv", sep="\t", index=False)
print(open(K.BASE / "SUMMARY" / "tables" / "proximity_recomputation_check_cisplatin.tsv").read())
