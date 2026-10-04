# Primary liver cancer tumouroids (Broutier et al. 2017)

**Paper**: Broutier L, Mastrogiovanni G, Verstegen MM, et al. *Human primary liver cancer-derived organoid cultures
for disease modeling and drug screening.* Nat Med 2017;23(12):1424-1435. doi:10.1038/nm.4438 (PMC5722201).

**Accessions**: RNA-seq GEO GSE84073 (here taken from the authors' processed RPKM table, Supplementary Dataset 1,
sheet `S1_RPKM_values`); drug screen = Supplementary Dataset 5, sheet `3_Raw data` (29 compounds, 5 tumouroid lines,
4 replicate fits each; GDSC pipeline). Both from the Europe PMC supplementary bundle.

| file | URL | sha256 |
|---|---|---|
| `PMC5722201_supplementary.zip` | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5722201/supplementaryFiles | `6adfceaa1cd35929e34b5e9deaa8528fb529e1538d759e08176f5ccc9a7b3828` |
| `Homo_sapiens.gene_info.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | `0249a78831f1c83f5c76f9c5c5378109559f883154648c6ed9f1d6ca635e47bd` |

## Derivation
* Response: two metrics per line x drug, median over the 4 replicate fits (`n`):
  `lnIC50` = column `IC50` (natural log of IC50 in uM; the sheet defines IC50 (uM) = EXP(IC50)) and
  `AUC` = fitted-curve AUC (fraction viability, 0-1). **Lower = more sensitive** for both. IC50s extrapolated beyond
  `max_conc_uM` are kept as published. Drug names normalised (Gemcitibine -> GEMCITABINE, Dasatanib -> DASATINIB,
  CH5424802 -> ALECTINIB, PD-0332991 -> PALBOCICLIB, BIRB 0796 -> DORAMAPIMOD, AZD8931 -> SAPITINIB,
  EMD-1214063 -> TEPOTINIB, LGK974 -> WNT-974).
* Expression: RPKM -> log2(RPKM + 1); tumouroid columns only (`HCC1_O`, `HCC3_O`, `CHC1_O[a,b]`, `CHC2_O`,
  `CC1_O[a-c]`, `CC2_O`, `CC3_O`); technical/clone replicates averaged on the log scale. Ensembl IDs -> HGNC symbols,
  protein-coding. Units: **log2(RPKM + 1)**.

## Sample-ID matching
Drug-sheet IDs `HCC-1` <-> expression `HCC1_O` (hyphen removed -> `HCC1`). Expression lines: 7;
screened lines: 5; **overlap n = 5**. Genes: 13088; drugs: 29.
Very small n: useful for sanity checks / pooled analyses only.
