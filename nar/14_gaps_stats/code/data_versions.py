"""Record sha256 + size of every input file used by the G9-G18 scripts."""
import pandas as pd

from g_common import A3, NAR, REPO, TAB, C

files = [A3 / "_data/a3_patients_gi.tsv.gz", A3 / "_data/a3_patients_other.tsv.gz", A3 / "tables/A3_master_results.tsv",
         A3 / "_data/signatures_gi.tsv", A3 / "_data/signatures_other.tsv", A3 / "_data/raw/COADREAD_clinicalMatrix",
         NAR / "02_organoid_only/_data/raw/Survival_SupplementalTable_S1_20171025_xena_sp.tsv",
         REPO / "robust/revamp_audit/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
         REPO / "robust/results/training_table/pathway_scores_tcga.tsv.gz", REPO / "robust/results/training_table/patients.tsv.gz"]
files += sorted((NAR / "03_organoid_plus_patient/cache").glob("*__shortlist__scores.tsv.gz"))
files += sorted((NAR / "02_organoid_only/BLCA_OV_LUAD_OTHER/_data/cache").glob("ssgsea_*_tcga.tsv.gz"))
pd.DataFrame([{"file": str(f.relative_to(REPO)), "bytes": f.stat().st_size, "sha256": C.sha256(f)} for f in files]).to_csv(
    TAB / "data_versions.tsv", sep="\t", index=False)
