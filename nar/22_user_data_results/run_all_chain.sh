#!/bin/bash
# one process at a time (1 CPU).  Order: cheap steps first, then the heavy ones.
cd "$(dirname "$0")/code"
L=../logs
python3 -I apply_model_cohorts.py lung_listed > $L/apply_lung_listed.log 2>&1
python3 -I lung_clusters.py > $L/lung_clusters.log 2>&1
python3 -I liver_describe.py > $L/liver.log 2>&1
python3 -I ovary_analysis.py > $L/ovary_analysis.log 2>&1
python3 -I apply_model_cohorts.py ovary > $L/apply_ovary.log 2>&1
python3 -I apply_model_cohorts.py ovary_7lines > $L/apply_ovary7.log 2>&1
python3 -I run_joint_u22.py lung shortlist > $L/joint_lung.log 2>&1
python3 -I run_joint_u22.py ovary shortlist > $L/joint_ovary.log 2>&1
python3 -I run_paper_method.py lung > $L/paper_lung.log 2>&1
python3 -I run_paper_method.py ovary > $L/paper_ovary.log 2>&1
python3 -I summarize_paper_method.py lung > $L/summ_lung.log 2>&1
python3 -I summarize_paper_method.py ovary > $L/summ_ovary.log 2>&1
echo chain finished > $L/chain.done
