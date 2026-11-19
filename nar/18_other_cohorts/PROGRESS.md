# W2 progress
- start: v1 sha256 and lung-only v2 manifest saved in tables/
- builders written in robust/build_training_cohorts.py (OTHER_BUILDERS); build_all.sh builds each cohort separately (log build_other2.log); v1 sha256 OK after restart
- next: build_training_table --others, analysis run_v2.py, report
- all 56 non-lung cohorts built; running build_training_table --others
- v2 table built; v1 sha256 OK; arm_check.tsv written (key trials match); GSE9782 event-flag bug fixed in builder (run_v2 reads fixed values)
- run_v2.py ran clean (31/33 old AUCs reproduced; 2 lung differ as expected); 4 figures + REPORT.md written (code/make_figures.py, make_report.py). DONE; nothing committed.
