# GSE28702 (Tsuji et al., Br J Cancer 2012; 83 unresectable CRC before mFOLFOX6, GPL570)

Retrieved 2026-10-04 from curatedCRCData (GitHub mirror, release branch). Build: `scratchpad/work/build_geo.py::gse28702`.

| file | URL | sha256 |
|---|---|---|
| GSE28702_eset.rda | https://raw.githubusercontent.com/bioc/curatedCRCData/RELEASE_3_18/data/GSE28702_eset.rda | aaf3e0e0f1143b9811aa7574c0caeb3ed8086697b33881aab667ed69702ebc19 |

## expression.tsv.gz
curatedCRCData gene-level matrix (19,320 x 83), log2 (GEO values were APT PLIER-MM-sketch, log2 per curatedCRCData).
One representative probeset per gene (curatedCRCData collapse), not a median over probes.

## clinical.tsv (83 rows, 83 distinct patients)
- `chemo` = 1, `regimen` = mFOLFOX6, `fu_based` = 1 for all (first-line, unresectable/metastatic).
- `responder` = 1/0 from GEO `mfolfox6: responder/non-responder` (42/41). The paper defines responders as CR+PR by
  RECIST, but per-patient CR/PR/SD/PD is not published, so `response` = NA.
- No survival data published: `os_*`, `rfs_*`, `stage` are NA.
- `lesion`/`sample_type`: 56 primary tumours, 27 metastases (one sample per patient). `set` = authors' training (54) / test (29) split.
