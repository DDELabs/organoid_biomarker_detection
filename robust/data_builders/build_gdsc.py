import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
import sys, pandas as pd, numpy as np
sys.path.insert(0, _REPO)
from robust.obd.reference import common_drug_name
R='/tmp/obd_raw/pre_raw/'
O=_REPO + '/data/external/gdsc/'
det=pd.read_excel(R+'Cell_Lines_Details.xlsx',sheet_name='Cell line details')
det.columns=['name','cosmic','wes','cna','expr','meth','drug','tissue1','tissue2','tcga','msi','medium','growth']
det=det[det.cosmic.notna()].copy(); det['cosmic']=det.cosmic.astype(int)
si=pd.read_csv(R+'depmap_sample_info.csv')
md=pd.read_csv(R+'depmap_Model.csv')
sid2cos=dict(zip(si.Sanger_Model_ID.dropna(), si.loc[si.Sanger_Model_ID.notna(),'COSMICID']))
sid2cos.update({k:v for k,v in zip(md.SangerModelID, md.COSMICID) if pd.notna(k) and pd.notna(v)})
norm=lambda s: str(s).upper().replace('-','').replace(' ','').replace('.','').replace('_','').replace('/','').replace(':','')
name2cos={norm(n):c for n,c in zip(det.name,det.cosmic)}
cos2t1=dict(zip(det.cosmic,det.tissue1)); cos2t2=dict(zip(det.cosmic,det.tissue2)); cos2tcga=dict(zip(det.cosmic,det.tcga))
cos2ach=dict(zip(si.COSMICID.dropna().astype(int), si.loc[si.COSMICID.notna(),'DepMap_ID']))
cos2ach.update(dict(zip(md.COSMICID.dropna().astype(int), md.loc[md.COSMICID.notna(),'ModelID'])))
cos2lin=dict(zip(si.COSMICID.dropna().astype(int), si.loc[si.COSMICID.notna(),'lineage']))
cos2lin.update(dict(zip(md.COSMICID.dropna().astype(int), md.loc[md.COSMICID.notna(),'OncotreeLineage'])))
out=[]
for ds in ['GDSC1','GDSC2']:
    d=pd.read_csv(R+ds+'_fitted.csv.gz')
    c=d.SANGER_MODEL_ID.map(sid2cos)
    c=c.fillna(d.CELL_LINE_NAME.map(norm).map(name2cos))
    print(ds, 'rows',len(d),'unmapped', c.isna().sum(), d.loc[c.isna(),'CELL_LINE_NAME'].nunique())
    d['cosmic_id']=c
    d=d[c.notna()].copy(); d['cosmic_id']=d.cosmic_id.astype(int)
    names={n:common_drug_name(n) for n in d.DRUG_NAME.unique()}
    d['drug']=d.DRUG_NAME.map(names)
    out.append(pd.DataFrame({'cosmic_id':d.cosmic_id,'cell_line':d.CELL_LINE_NAME,'sanger_model_id':d.SANGER_MODEL_ID,
        'depmap_id':d.cosmic_id.map(cos2ach),'tissue':d.cosmic_id.map(cos2t1),'tissue_sub':d.cosmic_id.map(cos2t2),
        'tcga_label':d.cosmic_id.map(cos2tcga),'gdsc_cancer_type':d.CANCER_TYPE,'depmap_lineage':d.cosmic_id.map(cos2lin),
        'drug':d.drug,'drug_name_gdsc':d.DRUG_NAME,'drug_id':d.DRUG_ID,'ln_ic50':d.LN_IC50.round(4),'auc':d.AUC.round(4),
        'max_conc_uM':d.MAX_CONC,'dataset':ds}))
# CTRPv2 (AUC over 16-pt curve; different scale from GDSC AUC)
ca=pd.read_csv(R+'ctrp/v22.data.auc_sensitivities.txt',sep='\t')
cl=pd.read_csv(R+'ctrp/v22.meta.per_cell_line.txt',sep='\t')
cp=pd.read_csv(R+'ctrp/v22.meta.per_compound.txt',sep='\t')
strip2cos={norm(a):c for a,c in zip(md.StrippedCellLineName, md.COSMICID) if pd.notna(c)}
strip2cos.update({norm(a):c for a,c in zip(si.stripped_cell_line_name, si.COSMICID) if pd.notna(c) and norm(a) not in strip2cos})
for n,c in name2cos.items(): strip2cos.setdefault(n,c)
cl['cosmic_id']=cl.ccl_name.map(lambda x: strip2cos.get(norm(x)))
ca=ca.merge(cl[['index_ccl','ccl_name','cosmic_id']],on='index_ccl').merge(cp[['index_cpd','cpd_name','top_test_conc_umol']],on='index_cpd')
print('CTRP lines', cl.shape[0], 'mapped to COSMIC', cl.cosmic_id.notna().sum())
ca=ca[ca.cosmic_id.notna()].copy(); ca['cosmic_id']=ca.cosmic_id.astype(int)
names={n:common_drug_name(n) for n in ca.cpd_name.unique()}
out.append(pd.DataFrame({'cosmic_id':ca.cosmic_id,'cell_line':ca.ccl_name,'sanger_model_id':None,
    'depmap_id':ca.cosmic_id.map(cos2ach),'tissue':ca.cosmic_id.map(cos2t1),'tissue_sub':ca.cosmic_id.map(cos2t2),
    'tcga_label':ca.cosmic_id.map(cos2tcga),'gdsc_cancer_type':None,'depmap_lineage':ca.cosmic_id.map(cos2lin),
    'drug':ca.cpd_name.map(names),'drug_name_gdsc':ca.cpd_name,'drug_id':ca.index_cpd,'ln_ic50':float('nan'),'auc':ca.area_under_curve.round(4),
    'max_conc_uM':ca.top_test_conc_umol,'dataset':'CTRPv2'}))
resp=pd.concat(out)
resp.to_csv(O+'response.tsv.gz',sep='\t',index=False)
# expression
if 0:
  e=pd.read_csv(R+'Cell_line_RMA_proc_basalExp.txt',sep='\t')
  e=e[e.GENE_SYMBOLS.notna()].drop(columns='GENE_title')
  e.columns=['gene']+[c.replace('DATA.','') for c in e.columns[1:]]
  e=e.groupby('gene').mean()
  # some columns duplicated as e.g. '1503362.1'
  e.columns=[c.split('.')[0] for c in e.columns]
  e=e.T.groupby(level=0).mean().T.round(3)
  e.to_csv(O+'expression.tsv.gz',sep='\t',compression={'method':'gzip','compresslevel':9})
  print('expr',e.shape)
e=pd.read_csv(O+'expression.tsv.gz',sep='\t',index_col=0,nrows=2)
ex=set(int(c) for c in e.columns)
for drug in ['FLUOROURACIL','CISPLATIN','GEMCITABINE','SORAFENIB','OXALIPLATIN','IRINOTECAN','SN-38','TEMOZOLOMIDE','DOXORUBICIN','PACLITAXEL']:
    for ds in ['GDSC1','GDSC2','CTRPv2']:
        s=resp[(resp.drug==drug)&(resp.dataset==ds)]
        print(drug,ds,s.cosmic_id.nunique(), s.cosmic_id.isin(ex).groupby(s.cosmic_id).first().sum(), s.drug_id.unique())
print(resp.drug.nunique(), resp.cosmic_id.nunique(), resp.cosmic_id.isin(ex).sum())

g=resp[(resp.drug=='TEMOZOLOMIDE')]
gl=g[g.tissue_sub.fillna('').str.contains('glioma',case=False)|g.tcga_label.isin(['GBM','LGG'])|g.depmap_lineage.fillna('').str.contains('CNS|central',case=False)]
print('TMZ glioma/CNS', gl.groupby('dataset').cosmic_id.nunique().to_dict(), 'with expr', gl[gl.cosmic_id.isin(ex)].groupby('dataset').cosmic_id.nunique().to_dict())
