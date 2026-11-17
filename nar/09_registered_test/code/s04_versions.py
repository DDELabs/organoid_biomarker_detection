"""Record the input data files used (path, size, sha256)."""
import pandas as pd
import common as C
files = [C.FROZEN, C.TT / "patients.tsv.gz", C.TT / "pathway_scores_tcga.tsv.gz", C.TT / "pathway_scores_geo.tsv.gz", C.PURITY,
         C.DATA / "tcga_response" / "drug_response_long.tsv.gz", C.DATA / "tcga_response" / "survival_cdr.tsv.gz",
         C.REPO / "data" / "curated" / "trials" / "GSE109211" / "expression.tsv.gz",
         C.REPO / "data" / "curated" / "trials" / "GSE109211" / "clinical.tsv", C.HERE / "ANALYSIS_PLAN.md"]
rows = [{"file": str(f.relative_to(C.REPO)), "bytes": f.stat().st_size, "sha256": C.sha256(f)} for f in files]
pd.DataFrame(rows).to_csv(C.TAB / "data_versions.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
