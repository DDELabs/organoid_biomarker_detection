"""Check that gcommon.ssgsea_fast reproduces gseapy raw ssGSEA ES (used for random gene-set nulls)."""
import sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../../../robust')
import pandas as pd, numpy as np
from gcommon import *
from obd.reference import gene_sets
from obd import cohorts as C, scoring
t=time.time()
e=pd.read_csv(VDW/'expression.tsv.gz',sep='\t',index_col=0); e=C.to_canonical(e); print('load',time.time()-t, e.shape, flush=True)
gs=gene_sets(("REACTOME",)); gs={k:[g for g in v if g in set(e.index)] for k,v in gs.items()}; gs={k:v for k,v in gs.items() if len(v)>=5}
sub={k:gs[k] for k in list(gs)[:5]+[BH3]}
t=time.time(); f=ssgsea_fast(e,sub); print('fast',time.time()-t, flush=True)
t=time.time(); g=scoring.ssgsea(e,sub,threads=1); print('gseapy',time.time()-t, flush=True)
rows=[dict(set=k, r=np.corrcoef(f.loc[k],g.loc[k,f.columns])[0,1], fast0=f.loc[k].iloc[0], gseapy0=g.loc[k,f.columns[0]]) for k in sub]
d=pd.DataFrame(rows); print(d); d.to_csv(TAB/'check_ssgsea_fast_vs_gseapy.tsv',sep='\t',index=False)
