import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
import sys, numpy as np, pandas as pd
sys.path.insert(0,_REPO)
from robust.obd.reference import common_drug_name
R='/tmp/obd_raw/pre_raw/repos/lm687/'
O=_REPO + '/data/external/organoids/ovarian_vias2023/'
jb2org={'JBLAB19902':'119148','JBLAB19904':'54327','JBLAB19905':'118976','JBLAB19906':'119178','JBLAB19907':'119127',
 'JBLAB19917':'23868','JBLAB19921':'119058','JBLAB19925':'54288','JBLAB19937':'54059','JBLAB19938':'151723',
 'JBLAB19939':'151761','JBLAB19940':'54276','JBLAB19941':'32077','JBLAB19942':'151773'}
t=pd.read_csv(R+'RNASeq_DE_resistant_sensitive/files/20191218_ViasM_BJ_orgaBrs_tpm.csv')
t=t[t.gene_name.notna()]
m=t[list(jb2org)].copy(); m.columns=[jb2org[c] for c in m.columns]
m.index=t.gene_name
m=np.log2(m.groupby(level=0).sum()+1).round(4)
m=m[(m>0).any(axis=1)]
m.index.name='gene'
m.to_csv(O+'expression.tsv.gz',sep='\t',compression={'method':'gzip','compresslevel':9,'mtime':0})
a=pd.read_csv(R+'survival_analysis/data/20200419-AUC-organoids.csv')
a['sample']=a['id'].astype(str)
alias={'AZD8186':'AZD8186','AZD2014':'VISTUSERTIB','AZD5363':'CAPIVASERTIB','AZD2281':'OLAPARIB','AZD0156':'AZD0156',
 'AZD6738':'CERALASERTIB','AZD1775':'ADAVOSERTIB','AZD8835':'AZD8835','Elescamol':'ELESCLOMOL','APR-246':'EPRENETAPOPT'}
rows=[]
for c in [c for c in a.columns if c.startswith('auc_ll5.')]:
    d=c.split('.',1)[1]; name=common_drug_name(alias.get(d,d))
    for s,v in zip(a['sample'],a[c]):
        if pd.notna(v): rows.append((s,name,round(float(v),5),'AUC_LL5_viability'))
r=pd.DataFrame(rows,columns=['sample','drug','response','metric'])
r.to_csv(O+'response.tsv',sep='\t',index=False)
print(m.shape, r.shape, r.drug.unique(), sorted(set(r['sample'])&set(m.columns)).__len__())
# sanity: PFI file consistent with AUC file
p=pd.read_csv(R+'survival_analysis/data/20200419-AUC-organoids_PFI.csv'); p['id']=p.id.astype(str)
x=r[r.drug=='PACLITAXEL'].set_index('sample').response
print('paclitaxel agree', np.allclose(x[p.id].values, p['auc_ll5.Paclitaxel'].values, atol=1e-6))
import scipy.stats as st
pf=p.set_index('id'); 
for d in ['PACLITAXEL','OXALIPLATIN']:
    xx=r[r.drug==d].set_index('sample').response
    print(d, 'median AUC resistant PFI', xx[pf.index[pf.PFI=='Resistant']].median(), 'sensitive', xx[pf.index[pf.PFI=='Sensitive']].median())
