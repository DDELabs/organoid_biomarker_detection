# PROGRESS (W3 breast)

- [x] Read brief + PLAN W3. Folder created.
- [x] Reproduced the existing L2 breast leave-one-cohort-out AUCs to 3 decimals, 13/13 cohorts
      (code/repro_l2.py -> tables/l2_reproduction_check.tsv). Max difference 0.0000.
- [x] Pre-computed non-outcome inputs (code/prep.py): proximity shortlist (paper table), pathway overlap graph, proliferation list.
- [x] Models running in 2 background processes (code/bench_breast.py; logs in logs/runA.log, logs/runB.log).
- [x] Baselines (code/baselines.py): ER+HER2 only, proliferation score, 200 random pathway sets.
- [x] Statistics (code/analyze.py) and figures (code/figures.py).
- [x] REPORT_breast.md

Notes
- STRING v12 proximity cache is not in the repo (robust/results/cache has no proximity files). The shortlist therefore uses the
  paper's own proximity table (4 breast drugs in it: paclitaxel, docetaxel, doxorubicin, 5-FU), union = 183 pathways.

All done. Final: REPORT_breast.md (models finished after a restart; RF bug fixed and re-run).

# PROGRESS (W3 positives)  - started 2026-11-19
- [x] Read brief, PLAN W3, breast report, original code of each positive (lung 17_lung/run_lung.py, bladder 02+07, ERCC2/BRAF 10+13, rectal 05, ovarian 03+22).
- [x] code/pzoo.py (model zoo, regression form, inner-CV tuning only inside training data), code/pos_stats.py (metrics, paired bootstrap)
- [x] code/pos_organoid.py: organoid-trained tasks. Bladder paper-data original reproduced (rank-1 log-rank p 0.044; AA pathway p 0.031).
- [ ] running: bladder current / ERCC2 / cetuximab BRAF / rectal / olaparib (logs/pos_orgA.log, pos_orgB.log)
- [ ] lung untreated prognosis, ovarian joint LOCO, analysis, figures, REPORT_positives.md
