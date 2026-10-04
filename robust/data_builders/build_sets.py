import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..', '..'))
import sys, hashlib, numpy as np, pandas as pd
sys.path.insert(0, _REPO)
from robust.obd.reference import ensembl_to_symbol, common_drug_name
from robust.obd.cohorts import _r_dataframe
RAW = '/tmp/obd_raw/org_raw'
OUT = _REPO + '/data/external/organoids'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()

# ---------------- pancreas (Tiriac 2018) ----------------
d = f'{RAW}/mlmed/data/organoid_pancreas'
for f in ['organoid_pancreas_fpkm.txt', 'organoid_value.tsv', 'file_name_to_organoid_mapping.tsv']:
    print('sha256', f, sha(f'{d}/{f}'))
fp = pd.read_csv(f'{d}/organoid_pancreas_fpkm.txt', sep='\t', index_col=0)
fmap = pd.read_csv(f'{d}/file_name_to_organoid_mapping.tsv', sep='\t').set_index('filename')['organoid']
fp.columns = [fmap[c + '.gz'] for c in fp.columns]
e2s = ensembl_to_symbol()
fp.index = fp.index.str.split('.').str[0]
fp = fp[~fp.index.str.endswith('_PAR_Y')] if False else fp
sym = fp.index.map(e2s)
fp = fp[sym.notna()]; fp.index = sym[sym.notna()]
lg = np.log2(fp + 1)
lg = lg.T.groupby(level=0).mean().T          # organoids with several GDC files -> mean
lg = lg.groupby(level=0).mean()              # duplicate symbols -> mean
lg.index.name = 'gene'
lg = lg[sorted(lg.columns)]
lg.round(4).to_csv(f'{OUT}/pancreas_tiriac2018/expression.tsv.gz', sep='\t', compression='gzip')
r = pd.read_csv(f'{d}/organoid_value.tsv', sep='\t')
r['drug'] = r.inhibitor.map(lambda x: 'SN-38' if x == 'SN-38' else common_drug_name(x))
r = r.rename(columns={'organoid': 'sample', 'value': 'response'})
r['sample'] = r['sample'].str.strip()
r.loc[r.drug == 'DISULFURAM', 'drug'] = 'DISULFIRAM'
print('dup pairs', r.duplicated(['sample','drug']).sum())
r = r.groupby(['sample', 'drug'], as_index=False)['response'].mean()
r['metric'] = 'AUC'
r['direction'] = 'lower=more sensitive'
r = r[['sample', 'drug', 'response', 'metric', 'direction']].sort_values(['sample', 'drug'])
r.to_csv(f'{OUT}/pancreas_tiriac2018/response.tsv', sep='\t', index=False)
ov = sorted(set(lg.columns) & set(r['sample']))
print('PANC expr', lg.shape, 'resp organoids', r['sample'].nunique(), 'drugs', r.drug.nunique(), 'overlap', len(ov))
print('expr-only', sorted(set(lg.columns) - set(r['sample'])), 'resp-only', sorted(set(r['sample']) - set(lg.columns)))
print(sorted(r.drug.unique()))
# sanity: gemcitabine sensitivity vs nothing; check direction later

# ---------------- liver (LICOB) ----------------
import rdata
p = f'{RAW}/wu-yc_iLICOB/data/data_ilicob_org'
print('sha256 data_ilicob_org', sha(p))
full = rdata.parser.parse_file(p).object.value[0].value[3]
ex = _r_dataframe(full.value[0]); auc = _r_dataframe(full.value[5])
ex.index.name = 'gene'
ex.round(4).to_csv(f'{OUT}/liver_licob/expression.tsv.gz', sep='\t', compression='gzip')
long = auc.rename_axis('sample').reset_index().melt(id_vars='sample', var_name='drug_raw', value_name='response').dropna()
long['drug'] = long.drug_raw.map(lambda c: common_drug_name(c.replace('_', '-')))
long['metric'] = 'AUC'; long['direction'] = 'lower=more sensitive'
long[['sample', 'drug', 'response', 'metric', 'direction']].sort_values(['sample', 'drug']).to_csv(f'{OUT}/liver_licob/response.tsv', sep='\t', index=False)
print('LICOB expr', ex.shape, 'resp', auc.shape, 'overlap', len(set(ex.columns) & set(auc.index)), 'AUC range', np.nanmin(auc.values), np.nanmax(auc.values))
print(ex.iloc[:3, :3]); print(ex.describe().T.head(3))
