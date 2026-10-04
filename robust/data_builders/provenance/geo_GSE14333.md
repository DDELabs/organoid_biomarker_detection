# GSE14333 (Jorissen et al., Clin Cancer Res 2009; 290 CRC, GPL570)

Retrieved 2026-10-04 from the Bioconductor experiment package curatedCRCData (GitHub mirror of release branch;
GEO/Zenodo are blocked here). Build script: `scratchpad/work/build_geo.py::gse14333`; the .rda ExpressionSet is read
with a small pure-Python reader (`scratchpad/work/rda_eset.py`, uses the `rdata` parser).

| file | URL | sha256 |
|---|---|---|
| GSE14333_eset.rda (ExpressionSet) | https://raw.githubusercontent.com/bioc/curatedCRCData/RELEASE_3_18/data/GSE14333_eset.rda | 1c4109d04a7c9a080ddf5bf308f44b9fb1eb1edb6a9a678f39718d2c1996b012 |

## expression.tsv.gz
curatedCRCData gene-level matrix (19,320 symbols x 290), log2 scale as deposited in the GEO series matrix.
NOTE: curatedCRCData already collapsed probesets to genes by keeping ONE representative probeset per gene
(featureData `probeset`), so this is not a median over probes; no probe-level mirror was found.

## clinical.tsv (290 rows)
Parsed from the authors' GEO characteristics kept in `uncurated_author_metadata`:
`Location, DukesStage, Age_Diag, Gender, DFS_Time, DFS_Cens, AdjXRT, AdjCTX`.
- `rfs_months` = DFS_Time (months; max 142.6). `rfs_event` = 1 - DFS_Cens (authors coded 1 = censored, 0 = event;
  confirmed by curatedCRCData curation script and by event times being shorter). Dukes D patients have no DFS.
- `os_*` NA (not published). `stage` = Dukes A/B/C/D -> I/II/III/IV.
- `chemo` = AdjCTX Y/N (117/172, 1 NA). The regimen is NOT reported per patient: `regimen` = `adjuvant_unspecified`,
  `fu_based` = NA for treated patients. `fu_based_presumed` = 1 for treated patients because adjuvant chemotherapy in
  this Melbourne cohort (operated 1993-2006) was 5-FU-based standard of care; treat this as an assumption.
- Extra: `adjuvant_radiotherapy` (AdjXRT), `location`, `dukes_stage`.

Counts: treated 117 (presumed 5-FU-based); RFS events 50/226 overall; 27 RFS events among the 117 treated.
