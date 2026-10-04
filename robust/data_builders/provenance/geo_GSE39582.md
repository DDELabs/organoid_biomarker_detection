# GSE39582 (Marisa et al., PLoS Med 2013; 566 CRC tumours, Affymetrix HG-U133 Plus 2.0 / GPL570)

GEO/NCBI is not reachable from this environment; all files come from public GitHub mirrors.
Retrieved 2026-10-04. Build script: `scratchpad/work/build_geo.py::gse39582` (not part of the repo).

## Downloads
| file | URL | sha256 |
|---|---|---|
| fRMA probe-level expression (54,675 probesets x 566 tumours, log2) | https://media.githubusercontent.com/media/PoShine/CRC-5-FU/main/1.primary%20data/train/GSE39582.frma.tsv (Git LFS) | a5c56314849eac6747460b04327dd1e08c92915e3421a128021e5a22eda5a8ab |
| PoShine sample table (not used for derived columns, cross-check only) | https://media.githubusercontent.com/media/PoShine/CRC-5-FU/main/1.primary%20data/train/GSE39582.frma.sample.tsv | 0a127c2b25c8e5592f4e5bae6370efe422271a842999f20f773f28c8c8ca8c88 |
| GEO characteristics excerpt incl. `chemotherapy.adjuvant` and `chemotherapy.adjuvant.type` (566 rows) | https://raw.githubusercontent.com/johnnyliu1992/GEV_Linux/master/media/patient_feature_sample.csv | 888138b9cfa377dcab1aa1acafa4de4ebf61ea1741cf3c1e709bd0499f5ffe3e |
| GEO characteristics incl. MMR/CIMP/CIN/CIT subtype/KRAS/BRAF/TP53 (585 rows incl. 19 non-tumour) | https://raw.githubusercontent.com/xmuyulab/scRank-XMBD/master/data/bulkRawClinical/GSE39582_Clinical_data.txt | 35b386d0aa99f0b63348f0e02c9479909c2c282ab51973c7f4ae7e515dd8f3ed |
| GPL570 annotation (GEO platform table GPL570-55999) | https://media.githubusercontent.com/media/Mmynemious/BBB-Gene-Project/main/raw_data/GPL570-55999.txt (Git LFS; sha matches LFS oid) | ecbfa3a4ac145d06b2b42508a05d0c7123fbb4ab681d2b94d6e6330a9b69aa7f |

## expression.tsv.gz
fRMA-normalised log2 intensities (PoShine/CRC-5-FU preprocessing of the GEO CEL files). Probesets mapped to
`Gene Symbol` of GPL570-55999; probesets annotated to several genes (`///`) or to none were dropped; median over
probesets per symbol -> 21,655 genes x 566 samples. Symbols are as in the GEO annotation (not repo-canonicalised;
use `obd.cohorts.to_canonical` if needed).

## clinical.tsv (566 rows)
- `os_months`/`os_event` = GEO `os.delay (months)`/`os.event`; `rfs_months`/`rfs_event` = `rfs.delay`/`rfs.event`.
  Values agree exactly between the two independent mirrors (asserted in the build).
- `stage` = `tnm.stage` 1-4 -> I-IV (4 stage-0 tumours -> NA; original in `tnm_stage`).
- `age` = `age.at.diagnosis (year)`; `sex` 1 = male; `chemo` = `chemotherapy.adjuvant` Y/N -> 1/0 (17 NA).
- `regimen` = `chemotherapy.adjuvant.type` (5FU, FUFOL, FOLFOX, FOLFIRI, other); `unspecified` = chemo Y with type NA;
  `none` = chemo N.
- `fu_based` = 1 for 5FU/FUFOL/FOLFOX/FOLFIRI (164), 0 for no chemo, NA for `other`/unspecified/unknown.
- `response` NA (adjuvant setting). `setting` = adjuvant. Extra: TNM T/N/M, location, MMR, CIMP, CIN, CIT subtype, KRAS/BRAF/TP53.

Counts: chemo 233 Y / 316 N; 5-FU-based 164 (5FU 79, FUFOL 51, FOLFOX 23, FOLFIRI 11); OS events 191/562; RFS events 177/560.
Within 5-FU-treated: OS events 53, RFS events 56.
