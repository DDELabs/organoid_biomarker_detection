# POSITIVE_RESULTS

Everything in the project that meets the positive rule, one folder per cancer. Files were **copied** from `nar/` and `robust/results/`; no original was moved or edited (`COPY_MANIFEST.tsv` lists source, destination and checksum of every file).

## The rule
A result is positive if it is **held-out/independent test, p < 0.05 in the pre-stated direction, or BH q < 0.10 in its own family; n and events reported**.
- Held-out or independent means the test patients, cell lines or mutants were not used to train or tune the model (leave-one-cohort-out, a different cohort, TCGA patients scored by a lab-trained model, isogenic cell pairs).
- Pre-stated direction: the direction the method predicts (higher predicted resistance = worse survival, HR above 1; AUC above 0.5; mutant more resistant or sensitive as the drug biology says). Where a direction was not fixed beforehand, or the result goes the other way, it is not positive.
- Every positive is given a tier, **strong** or **nominal** (see `INDEX.md`). Nominal means it meets p < 0.05 but it is one test of many, small, or does not survive family correction. I did not hide those, because the rule asks only for p < 0.05 or q < 0.10; read the tier before you use one.
- Borderline and opposite-direction results are in `BORDERLINE.md`.

## Folders
`BREAST/`, `LUNG/`, `BLADDER/`, `COLORECTAL/`, `RECTAL/`, `OVARIAN/`, `PAN_CANCER_DEPMAP/`. Each has `README.md` (what was tested, data, model, numbers, robustness, caveats, "What I think", sources), `figures/` (PNG + PDF + caption; file names start with the number of the source folder), `tables/` (copies of the source tables), `source_reports/`.
No LIVER, PANCREAS, SARCOMA or melanoma folder: nothing there met the rule in the expected direction (see `BORDERLINE.md`).
Top-level: `INDEX.md` (one table of all positives), `POSITIVE_INDEX.tsv` and `BORDERLINE_INDEX.tsv` (machine-readable; numbers read from the source tables), `COPY_MANIFEST.tsv`, `code/` (`collect.py`, `copy_files.py`, `write_readmes.py`; run in that order, one process, seconds).

## Scan coverage
Reports and tables read: nar/02, 03, 05, 08, 09, 10, 13, 14, 17, 18, 19, 20, 22 and robust/results/train_combined/. nar/09 (registered test: 0 of 7 drugs pass its gates) and nar/14 (shuffle null, meta-analysis, joint-model grid; G9 shuffle p for the 57-test survey 0.030, pooled HR per SD 0.96) contribute only context and borderline rows. The G14 cell-line models, STORM liver and the other drugs and cancers (pancreas, sarcoma, melanoma, liver) produced no result in the expected direction that met the rule outside BORDERLINE.
