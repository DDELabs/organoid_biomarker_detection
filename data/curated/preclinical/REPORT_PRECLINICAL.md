# ATLAS pre-clinical pharmacogenomic data: curation report

Built by `build_preclinical.py` (this directory). Layout per set: `expression.tsv.gz` (HGNC symbol x model,
protein-coding genes, log2 scale), `response.tsv` (long: `sample, drug, response, metric, direction`, plus
set-specific columns), `SOURCE.md` (paper, accessions, URLs, sha256, derivation, ID matching). Summary table:
`CATALOG_PRECLINICAL.tsv`. In every set **lower response = more sensitive**. No set needed a sign flip, because
all metrics are viability AUCs, ln IC50 or % viability.

Gene symbols come from NCBI `Homo_sapiens.gene_info` (protein-coding only; HGNC-authority symbol). Drug names are
generic UPPERCASE names. Source names are used where available. CoderData `improve_drug_id`s are named via the
InChIKey in the repo's DrugBank vocabulary, then by the PRISM name of the same BRD-K ID, then by CoderData's
source name or a synonym.

## Sets obtained

| set | models | n with expression **and** response | drugs | metric | expression units | source |
|---|---|---:|---:|---|---|---|
| `bladder_lee2018` | bladder cancer PDOs | **11** (21 with RNA-seq) | 50 | AUC (fitted, 0-1) | log2(TPM+1) | GEO GSE103990 (NCBI TPM) + Synapse syn64765430 via CoderData 2.1.0 |
| `prism_repurposing` | cancer cell lines (21 lineages) | **476** | 1448 | AUC (PRISM secondary) | log2(TPM+1) | figshare PRISM 19Q4 secondary + DepMap 21Q4 CCLE RNA-seq |
| `ctrpv2_ccle` | cancer cell lines | **819** | 460 | AUC (CoderData refit, 0-1) | log2(TPM+1) | CoderData 2.1.0 CTRPv2 + DepMap 21Q4 CCLE RNA-seq |
| `liver_ji2023` | liver cancer PDOs (HCC, ICC, CHC, HB) | **61** | 73 | AUC (fitted, 0-1) | log2(TPM+1) | Synapse via CoderData 2.2 (figshare 29923646) |
| `pancreas_shi2022` | pancreatic cancer PDOs | **39** (87 with RNA-seq) | 63 | AUC (normalised) | log2(FPKM+1) | GEO GSE194249 + Nat Commun Suppl. Data 7/8 |
| `pancreas_tiriac2018` | pancreatic cancer PDOs | **36** (58 screened) | 5 | AUC (fitted, 0-1) | log2(TPM+1) | CoderData 2.1.0 (Tiriac AACR figshare + HCMI/GDC RNA-seq) |
| `sarcoma_alshihabi2024` | sarcoma PDOs (mixed histologies) | **15** | 33 | viability_pct (single-dose % viability) | log2(TPM+1) | Synapse syn61892224 via CoderData 2.1.0 |
| `liver_broutier2017` | liver tumouroids (HCC, CHC, CC) | **5** | 29 | lnIC50 and AUC | log2(RPKM+1) | Nat Med supplement (PMC5722201) |

Caveats worth knowing before use:
* **Lee 2018 bladder**: the author DESeq2 table on GEO returns 403, so NCBI's uniform TPM re-quantification is used.
  Drug AUCs come from CoderData's refit of the raw Synapse dose-response files. Response is the median over the
  organoid passages screened, and the screened passage can differ from the sequenced one.
* **PRISM**: one curve per line x compound (MTS010 > MTS006 > MTS005 > HTS002), STR-passed lines only. Cell-line
  lineage is in `samples.tsv` and in the `lineage` column of `response.tsv`. It uses DepMap 21Q4
  `CCLE_expression.csv`, because the DepMap portal (Cloudflare challenge) could not be scripted for the newer
  `OmicsExpressionProteinCodingGenesTPMLogp1`.
* **CTRPv2**: `robust/obd/preclinical_sets.py` already pairs CTRPv2 AUC with GDSC RMA arrays (COSMIC IDs). This set
  is a different pairing: CoderData refit on a 0-1 scale plus CCLE RNA-seq on DepMap IDs. Do not pool both. 46/460
  compounds are unnamed CTRP tool compounds and keep BRD/IUPAC names (`drug_map.tsv`).
* **Sarcoma**: the response is a published viability score, not an AUC, so it is labelled `viability_pct`.
* **Broutier**: n = 5. Use for sanity checks or pooled analyses only.
* Sanity checks: ERBB2 vs lapatinib (PRISM) rho = -0.34, EGFR vs erlotinib (PRISM) rho = -0.28, BCL2 vs navitoclax
  (CTRPv2) rho = -0.26. HER2-amplified lines are lapatinib/afatinib/neratinib-sensitive in CTRPv2.

## Skipped / not obtained

| target | reason |
|---|---|
| Lee 2018 Cell supplementary drug tables | PMC5890941 is not open access (Europe PMC refuses the supplement). pmc.ncbi.nlm.nih.gov and www.cell.com are blocked by the egress proxy. The same screen was obtained through Synapse/CoderData instead (above). |
| CoderData colorectal (van de Wetering 2015) | Present in CoderData 2.2, but the released curve fits are degenerate (median fit R^2 ~ 0, all AUCs ~ 0.008, which looks like a dose-unit error). Only 6 organoids have both expression and response. The original Cell 2015 supplement is unreachable (cell.com blocked). |
| Driehuis 2019 PNAS (pancreas) | PMC6936689 is not open access, www.pnas.org is blocked, and the RNA-seq is under EGA controlled access. |
| Kopper 2019 Nat Med (ovarian) | Sequencing is under EGA controlled access. No open supplement with matched expression was found. |
| Yan 2018 Cell Stem Cell (gastric) | Not open access in Europe PMC. Expression is under EGA controlled access. |
| HCMI drug data | No public drug-response data beyond the Tiriac pancreatic PDOs (included). |
| CoderData `novartis` (PDXE) | Patient-derived xenografts, out of scope (organoids and cell lines only). |
| CTRPv2 original CTD2 portal files | ctd2-data.nci.nih.gov is blocked by the egress policy. The CoderData refit was used. |
| Raw Synapse files (Lee bladder, sarcoma, liver) | Synapse file downloads need a login (`Anonymous users have only READ access`), and the Synapse MCP connector returns metadata only. Harmonised CoderData releases on figshare were used instead. |
| api.figshare.com | Blocked (403) by the egress policy. Files were fetched from ndownloader.figshare.com by file ID or as a whole-article zip. |

## Rebuilding

```
python3 data/curated/preclinical/build_preclinical.py            # all sets (~4 min, ~4 GB raw cache)
python3 data/curated/preclinical/build_preclinical.py liver_ji2023 --raw /path/to/cache
```
Raw inputs are cached in `$ATLAS_RAW` (default `~/.cache/atlas_preclinical_raw`), outside the repo. The CoderData
2.2 release zip is 2.8 GB. Expression files are gzip with mtime 0, so rebuilds are byte-reproducible.
