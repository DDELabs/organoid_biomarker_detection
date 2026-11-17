#!/bin/bash
# queue B: waits for the BLCA run, then BLCA extras, olaparib time course, lung PDC expression
cd "$(dirname "$0")"; L=../_data
while pgrep -f "run_blca.py" >/dev/null; do sleep 10; done
python blca_extras.py > $L/blca_extras.log 2>&1
python ov_olaparib.py > $L/ov_olaparib.log 2>&1
python luad_pdc_expression.py > $L/luad_pdc.log 2>&1
echo chainB-done > $L/chainB.done
