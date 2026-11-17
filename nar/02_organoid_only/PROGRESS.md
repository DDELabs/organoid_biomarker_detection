# Progress (organoid-only, GI sets)

Done: data download + ssGSEA caches (_data/), code (code/), 14 result pickles (_data/results/).
Left, in order:
1. `cd code; python run_all.py 0 2 & python run_all.py 1 2` (resumes; skips finished pickles) until
   `_data/results/run_log_0.tsv` and `run_log_1.tsv` exist.
2. `python plots.py` (all per-combo figures, main composites, proximity heatmaps, cross-set table + forest).
3. `python extras.py` (Fig4b, SuppTable1, regimen split, Fig5a BRAF, SuppFig1, SuppFig10, SuppFig15 GDSC).
4. `python report_tables.py` -> tables/report_numbers.md; then write REPORT_GI.md (with checklist IDs).
5. Delete nothing; `_test/` holds early test figures only.

## Finisher agent status (2026-11-17 08:25)
- Workers restarted 08:05 (run2_0.log / run2_1.log); plots.py tested OK on 14 pickles; extras.py all steps ran OK
  (fixed engine.logrank: singular variance -> NaN). Fig4b, Fig5a, SuppFig1, SuppFig15, regimen split done.
- After workers end: `python plots.py; python extras.py clin s10 fig5a; python report_tables.py; python report_sections.py`
  then write REPORT_GI.md.
