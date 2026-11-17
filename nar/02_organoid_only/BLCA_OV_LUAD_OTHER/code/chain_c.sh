#!/bin/bash
# re-render BLCA and OV with the final figure labels (same seeds -> same numbers)
cd "$(dirname "$0")"; L=../_data
python run_blca.py > $L/run_blca2.log 2>&1
for f in ../../OV_vias2023/tables/*_summary.json; do mv "$f" $L/old_ov/; done
python run_sets.py OV > $L/run_OV2.log 2>&1
echo done > $L/chainC.done
