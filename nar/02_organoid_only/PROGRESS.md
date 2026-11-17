# Progress (organoid-only, GI sets)

Done: data download + ssGSEA caches (_data/), code (code/), 14 result pickles (_data/results/).
Left, in order:
1. `cd code; python run_all.py 0 2 & python run_all.py 1 2` (resumes; skips finished pickles) until
   `_data/results/run_log_0.tsv` and `run_log_1.tsv` exist.
2. `python plots.py` (all per-combo figures, main composites, proximity heatmaps, cross-set table + forest).
3. `python extras.py` (Fig4b, SuppTable1, regimen split, Fig5a BRAF, SuppFig1, SuppFig10, SuppFig15 GDSC).
4. `python report_tables.py` -> tables/report_numbers.md; then write REPORT_GI.md (with checklist IDs).
5. Delete nothing; `_test/` holds early test figures only.
