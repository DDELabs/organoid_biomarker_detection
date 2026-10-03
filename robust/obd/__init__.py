"""Robust organoid-to-patient biomarker discovery (Python 3).

Extends the network-based pipeline of Kong et al. (Nat Commun 2020) with
cohort-independent pathway scoring, organoid/tumour alignment, small-n
resampling, whole-pipeline permutation nulls and adjusted survival tests.
"""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
EXTERNAL = DATA / "external"
RESULTS = REPO / "robust" / "results"
