# Progress (nar/22_user_data_results)
- Done: input checks (tables/input_checks.tsv), loader import of 4 sets (_data/loader_out), apply_lung (tables/apply_lung.tsv).
- Key data finding: LUAD "IC50" is one log2FC per cluster x drug (not per organoid); LIHC series matrix has no expression; OVARY = 7 lines.
- Running: run_all_chain.sh (one process at a time): apply lung_listed, lung clusters, liver, ovary analysis, apply ovary, joint lung/ovary, paper method lung/ovary.
- To do after: figures, REPORT.md.
- Chain ran; apply ovary / joint lung+ovary failed because training_table_v2 was being rewritten by another agent (truncated gz); rerun after paper_lung/ovary finish (run_rerun.sh).
- Fixed carboplatin->cisplatin proximity proxy; deleting skipped pickles, rerun run_paper_method (lung/ovary) + summarize after run_rerun finishes.
- All steps done; REPORT.md written.
