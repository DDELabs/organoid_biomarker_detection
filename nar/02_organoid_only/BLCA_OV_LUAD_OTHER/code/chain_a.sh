#!/bin/bash
# queue A: waits for the OV run, then sarcoma, lung substitute, PDX
cd "$(dirname "$0")"; L=../_data
while pgrep -f "run_sets.py OV" >/dev/null; do sleep 10; done
python run_sets.py SARC > $L/run_SARC.log 2>&1
python run_sets.py LUAD > $L/run_LUAD.log 2>&1
python run_sets.py PDX > $L/run_PDX.log 2>&1
echo chainA-done > $L/chainA.done
