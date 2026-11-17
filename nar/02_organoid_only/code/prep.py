"""Compute (once) ssGSEA scores for every lab-model set and matching TCGA project(s)."""
import sys
from common import SETS, scores
keys = sys.argv[1:] or list(SETS)
for k in keys:
    s = scores(k, SETS[k][2])
    print(k, {p: v.shape for p, v in s.items()}, flush=True)
