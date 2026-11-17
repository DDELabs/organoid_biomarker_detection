"""Extra numbers quoted in REPORT.md: learning-curve endpoints, chosen joint weights, STORM placebo-arm AUROC, file versions."""
import glob
import numpy as np
import pandas as pd
import common as C

lc = pd.read_csv(C.TABLES / "learning_curves.tsv", sep="\t")
g = lc.groupby(["unit", "axis", "model", "fraction"])["concordance"].mean().unstack("fraction")
g.round(3).to_csv(C.TABLES / "learning_curves_mean.tsv", sep="\t")
print(g.round(3).to_string())
F = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(str(C.TABLES / "*__shortlist__folds.tsv"))])
w = F[F["model"] == "joint"].groupby("unit")["w"].agg(lambda v: ";".join(f"{k}:{n}" for k, n in v.value_counts().sort_index().items()))
w.to_csv(C.TABLES / "joint_weight_chosen.tsv", sep="\t")
print(w.to_string())
sc = pd.read_csv(C.CACHE / "LIHC_SORAFENIB__shortlist__scores.tsv.gz", sep="\t", index_col=0)
rows = []
for m in ["organoid-only", "patient-only", "joint"]:
    for grp in ["treated", "comparison"]:
        d = sc[sc["group"] == grp]
        a = C.auc_test(d[m], d["responder"])
        rows.append(dict(model=m, arm="sorafenib" if grp == "treated" else "placebo", n=len(d), responders=int(d["responder"].sum()), AUROC=a["AUROC"], p=a["p"]))
s = pd.DataFrame(rows); s.to_csv(C.TABLES / "storm_auroc_by_arm.tsv", sep="\t", index=False); print(s.round(4).to_string())
files = [C.TT / f for f in ("patients.tsv.gz", "pathway_scores_tcga.tsv.gz", "pathway_scores_geo.tsv.gz")]
files += [C.REPO / "data/curated/trials/GSE109211/expression.tsv.gz", C.REPO / "data/coad_blca_organoid_drugs_zscore_result_reactome.txt",
          C.REPO / "robust/results/survlearn/coefficients.tsv.gz"]
files += sorted(glob.glob(str(C.PRECLIN / "*" / "expression.tsv.gz"))) + sorted(glob.glob(str(C.PRECLIN / "*" / "response.tsv")))
files += sorted(glob.glob(str(C.FOLDER / "data/vias2023_raw/*")))
from pathlib import Path
v = pd.DataFrame([dict(file=str(Path(f).relative_to(C.REPO)), bytes=Path(f).stat().st_size, sha256=C.sha256(f)) for f in files])
v.to_csv(C.TABLES / "data_versions.tsv", sep="\t", index=False)
