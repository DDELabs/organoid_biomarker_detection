**Figure 1. How the test was pre-registered.** Each box is one step, top to bottom. The 35 pathway hypotheses were
fixed (frozen) from cell-line and network data only, before any patient outcome was used. The file's fingerprint
(sha256) was checked by the code at every run. Then the analysis plan was written, still without outcomes. Stage 1
tests each drug's score in TCGA patients, pooled across cancers. Stage 2 tests the same frozen score in outside
patient cohorts that have a comparison group. Orange boxes are the decision gates (G2, G3) and how many drugs passed.
