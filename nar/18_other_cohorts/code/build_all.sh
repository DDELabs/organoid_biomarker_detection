#!/bin/bash
# build each non-lung cohort in its own process; skip those already built
cd /home/user/organoid_biomarker_detection/robust
for c in $(python -c "import build_training_cohorts as B;print(' '.join(f.__name__ for f in B.OTHER_BUILDERS))"); do
  [ -f ../data/external/train_cohorts/$c/clinical.tsv ] && continue
  python build_training_cohorts.py $c 2>&1 | tail -3 | cut -c1-400
done
