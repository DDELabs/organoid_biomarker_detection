#!/bin/bash
cd "$(dirname "$0")/code"; L=../logs
while [ ! -f $L/chain.done ]; do sleep 30; done
python3 -I apply_model_cohorts.py ovary > $L/apply_ovary.log 2>&1
python3 -I apply_model_cohorts.py ovary_7lines > $L/apply_ovary7.log 2>&1
python3 -I run_joint_u22.py lung shortlist > $L/joint_lung.log 2>&1
python3 -I run_joint_u22.py ovary shortlist > $L/joint_ovary.log 2>&1
python3 -I make_figures.py > $L/figures.log 2>&1
echo ok > $L/rerun.done
