# Data builders

Scripts that rebuild the curated external datasets in `data/external/` (git-ignored)
from public mirrors that are reachable when GEO, EGA, journal and Bioconductor hosts
are not (GitHub raw / Git LFS, Xena S3, PyPI). Raw downloads go to `/tmp/obd_raw`.
`provenance/` holds the SOURCE.md of every dataset: exact URLs, sha256 checksums,
column derivations and caveats.

| Script | Builds | Loader |
|---|---|---|
| `build_geo.py` (+ `rda_eset.py`, a pure-Python reader for Bioconductor ExpressionSet .rda) | `data/external/geo/{GSE39582,GSE14333,GSE28702,TCGA-READ}` | `obd.geo_cohorts.load_geo_cohort` |
| `build_sets.py` | `data/external/organoids/{pancreas_tiriac2018,liver_licob}` | `obd.organoid_sets.load_organoid_set` |
| `build_ov.py` | `data/external/organoids/ovarian_vias2023` | `obd.preclinical_sets.load_organoid_set_extra` |
| `build_gdsc.py` | `data/external/gdsc` (GDSC1/2 + CTRPv2 response, RMA expression, COSMIC/DepMap map) | `obd.preclinical_sets.load_gdsc` |

GSE109211 (STORM, sorafenib arm only) was taken from the ENLIGHT mirror
(`PangeaResearch/enlight-data`); see `provenance/geo_GSE109211.md`. **Do not use it for
validation as is:** responders and non-responders separate on the first principal
component (56% of variance, AUC 0.05), a processing batch aligned with the label.
Rebuild it from the original GEO series matrix, which also includes the placebo arm.
