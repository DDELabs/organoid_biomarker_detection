"""Regression tests: exact reproduction of the paper and scoring invariances.

Run: python -m pytest robust/tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from obd import cohorts as C  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd import survival as S  # noqa: E402
from obd.reference import gene_sets  # noqa: E402
from obd.study import Study  # noqa: E402


def coad_paper_study():
    _, resp = C.coad_organoids()
    feats = N.proximal([N.precomputed_proximity()["FLUOROURACIL"]])
    return Study("COAD", "FLUOROURACIL", C.paper_scores("COAD", "organoid"), resp["FLUOROURACIL"],
                 C.paper_scores("COAD", "TCGA"), C.tcga_biotab_clinical("COAD"),
                 C.tcga_biotab_drugs("COAD")["FLUOROURACIL"], feats)


def test_reproduces_paper_multi_pathway_table():
    st = coad_paper_study()
    assert (len(st.samples), len(st.features), len(st.treated_pats)) == (19, 37, 114)
    ref = pd.read_csv(C.DATA.parent / "python/results/multi_pathway_predictions/multi_pathway_predictions.txt", sep="\t")
    m = st.baseline().merge(ref, left_on=["ML", "k"], right_on=["ML", "num_pathways"])
    assert len(m) == 27
    assert np.allclose(m.logrank_p, m.pvalue, atol=1e-12)
    assert (m.pathways_x == m.pathways_y).all()


def test_rank_score_invariant_to_per_sample_units():
    expr, _ = C.coad_organoids()
    lin = 2 ** expr - 1
    gs = dict(list(gene_sets(("REACTOME",)).items())[:50])
    a = scoring.rank_score(np.log2(lin + 1), gs)
    b = scoring.rank_score(np.log2(lin / lin.quantile(0.75) * 1e3 + 1), gs)  # FPKM-UQ-like rescaling
    c = scoring.rank_score(lin / lin.sum() * 1e6, gs)                         # TPM-like, unlogged
    assert np.allclose(a, b) and np.allclose(a, c)


def test_rank_score_independent_of_cohort():
    expr, _ = C.coad_organoids()
    gs = dict(list(gene_sets(("REACTOME",)).items())[:50])
    full = scoring.rank_score(expr, gs)
    half = scoring.rank_score(expr.iloc[:, :10], gs)
    assert np.allclose(full[half.columns], half)


def test_logrank_and_cindex_sanity():
    t = np.array([1, 2, 3, 4, 5, 6, 7, 8.0])
    e = np.ones(8)
    assert S.logrank(t[:4], e[:4], t[4:], e[4:]) < 0.05
    assert S.c_index(t, e, -t) == 1.0
