#!/bin/sh
# Runs the re-analyses that need CPU, one after another (max 2 processes), after the 03 CRC re-run finishes.
N=/home/user/organoid_biomarker_detection/nar
L=$N/15_fixes/logs; mkdir -p $L
until [ -f $N/03_organoid_plus_patient/cache/run_joint_fix15.done ]; do sleep 20; done
cd $N/03_organoid_plus_patient/code
python3 learning_curves.py CRC_5FU CRC_OXALIPLATIN --jobs 2 > $L/lc_crc.log 2>&1
python3 -c "
import pandas as pd,glob
pd.concat([pd.read_csv(f,sep='\t') for f in sorted(glob.glob('../tables/*__learning_curve.tsv'))]).to_csv('../tables/learning_curves.tsv',sep='\t',index=False)" >> $L/lc_crc.log 2>&1
python3 run_joint.py CRC_5FU BLCA_CISPLATIN PAAD_GEMCITABINE LIHC_SORAFENIB OV_PACLITAXEL LUNG_PLATINUM SARC_GEMCITABINE --features shortlist --fold-align --jobs 2 > $L/foldalign.log 2>&1
touch $L/foldalign.done
cd $N/05_A3_B1/code
(python3 b1_external.py > $L/b1.log 2>&1; python3 b1_summary.py > $L/b1_summary.log 2>&1; python3 b1_figures.py > $L/b1_fig.log 2>&1; python3 a3_figures.py > $L/a3_fig.log 2>&1; touch $L/b1.done) &
cd $N/15_fixes/code
python3 f02_bladder_one_line_per_patient.py > $L/bladder.log 2>&1
cd $N/08_depmap/code
python3 -c "
import importlib.util,sys
sys.argv=['x']
spec=importlib.util.spec_from_file_location('pa','02_partA.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
te,ts=m.part_e(); tf,fp,fn,w=m.part_f(); m.fig3(te,ts,tf,fp,fn); print('fig3 ok')" > $L/fig3_08.log 2>&1
wait
touch $L/chain.done
