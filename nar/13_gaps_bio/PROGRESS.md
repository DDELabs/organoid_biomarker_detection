# Progress (13_gaps_bio)
Scripts in code/ (run from code/, `python3 <script>`; max 2 at a time):
- gcommon.py helpers (fast ssGSEA = closed form, validated vs gseapy: tables/check_ssgsea_fast_vs_gseapy.tsv)
- G3 g3_organoid_bh3.py DONE (figure + tables)
- G2 g2_decompose.py (caches _data/g2_scores.pkl) -> then G5 g5_folfox.py (needs G2 table + G3 table + `python3 g7_bh3_beyond_markers.py g5`)
- G8 g8_ercc2.py then g8_braf.py (draws G8 figure)
- G7 g7_bh3_beyond_markers.py ; G13 g13_odd_pathways.py ; G14 g14_cellline_models.py
- Then REPORT.md (all numbers from tables/)

## Status (updated)
All scripts run: G2, G3, G5, G7, G8 (ERCC2+BRAF), G13, G14 DONE with figure+caption+tables.
Fix made: g7 skipped drugs absent from a screen (WEHI-539 not in GDSC1). Remaining: REPORT.md.
REPORT.md written by code/make_report.py (all numbers read from tables/). Bladder frozen-pathway part of G2 not done.
