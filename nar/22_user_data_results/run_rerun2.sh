#!/bin/bash
cd "$(dirname "$0")/code"; L=../logs
while [ ! -f $L/rerun.done ]; do sleep 30; done
python3 -I run_paper_method.py lung > $L/paper_lung2.log 2>&1
python3 -I run_paper_method.py ovary > $L/paper_ovary2.log 2>&1
python3 -I summarize_paper_method.py lung > $L/summ_lung.log 2>&1
python3 -I summarize_paper_method.py ovary > $L/summ_ovary.log 2>&1
python3 -I make_figures.py > $L/figures2.log 2>&1
echo ok > $L/rerun2.done
