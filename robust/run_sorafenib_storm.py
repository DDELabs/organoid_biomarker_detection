"""Liver cancer / sorafenib: frozen organoid signatures validated in the STORM trial.

GSE109211 sorafenib arm (Pinyol et al. 2019 Gut): 67 patients, 21 responders vs 46
non-responders (placebo arm and RFS times are not available from reachable mirrors).
Each ablation stage of the LICOB-trained model is frozen and scored on STORM; a good
signature gives responders a LOWER predicted-resistance score (AUC > 0.5).

Usage: python robust/run_sorafenib_storm.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import run_nextgen as NG  # noqa: E402
from obd import EXTERNAL, RESULTS  # noqa: E402
from obd import network as N  # noqa: E402
from obd import scoring  # noqa: E402
from obd import validation as V  # noqa: E402
from obd.cohorts import to_canonical  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "external_validation"


def storm():
    d = EXTERNAL / "geo" / "GSE109211"
    expr = pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0)
    expr = to_canonical(np.log2(expr.clip(lower=0) + 1))
    clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col=0)
    return expr, clin


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expr, clin = storm()
    gs = {k: [g for g in v if g in expr.index] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    sc = scoring.cached("external_GSE109211_reactome_rank", scoring.rank_score, expr, gs)
    responder = clin["response_binary"].astype(float)

    graph = {}
    ctx = NG.load("LIHC", "SORAFENIB", "licob", None, lambda: graph.setdefault("g", N.load_string()))
    rows = []
    for name, kw in NG.ABLATION:
        if kw.get("transfer") and ctx["prior"] is None:
            continue
        w, _ = NG.make_study(ctx, **kw).signature()
        score = V.signature_score(sc, w)
        r = V.response_test(score, responder)
        null = V.random_null(sc, len(w), lambda s: V.response_test(s, responder).get("AUC"), n=1000)
        r["random_null_p"] = float((1 + (null >= r["AUC"]).sum()) / (1 + len(null)))
        rows.append({"model": name, "pathways": len(w), **r, "signature": ";".join(w.index)})
        print(rows[-1]["model"], {k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()})
    pd.DataFrame(rows).to_csv(OUT / "sorafenib_STORM_GSE109211.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
