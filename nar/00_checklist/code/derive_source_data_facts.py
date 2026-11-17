#!/usr/bin/env python3
"""Derive the [derived] facts quoted in the checklist from the paper's source data (MOESM3).
Run: python3 -I nar/00_checklist/code/derive_source_data_facts.py
Writes nar/00_checklist/tables/source_data_facts.tsv
"""
import os, hashlib
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
X = os.path.join(ROOT, '..', '_cache', '41467_2020_19313_MOESM3_ESM.xlsx')
Z = -1.2816
out = []
def add(k, v): out.append((k, v)); print(k, '=', v)
add('MOESM3_sha256', hashlib.sha256(open(X, 'rb').read()).hexdigest())
raw = pd.read_excel(X, sheet_name=None, header=None)
e = raw['Figure 2a-f'].iloc[3:].dropna(how='all')
add('STRING_edge_rows', len(e)); add('STRING_distinct_uniprot', len(set(e[1]) | set(e[3])))
p = pd.read_excel(X, sheet_name='Figure 2a,d', header=2)
add('proximity_pathways', len(p)); add('proximity_drug_columns', p.shape[1] - 1)
n = (p.iloc[:, 1:] <= Z).sum(axis=1)
add('drugs_per_pathway_mean', round(n.mean(), 3)); add('drugs_per_pathway_median', n.median())
add('pathways_proximal_to_>=22_drugs', int((n >= 22).sum())); add('pathways_proximal_to_0_drugs', int((n == 0).sum()))
add('max_drugs_per_pathway', int(n.max()))
def terms(s, h): return set(pd.read_excel(X, sheet_name=s, header=h).iloc[:, 0])
for drug, org, tc, bm in [('FLUOROURACIL', 'Figure 2b', 'Figure 2g', 'REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS'),
                          ('CISPLATIN', 'Figure 2e', 'Figure 2h', 'REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION')]:
    prox = set(p.loc[p[drug] <= Z, 'Pathway']); o = terms(org, 3); t = terms(tc, 2)
    add(f'{drug}_proximal_all', len(prox)); add(f'{drug}_proximal_with_organoid_NES', len(prox & o))
    add(f'{drug}_proximal_with_organoid_and_TCGA_NES', len(prox & o & t))
    add(f'{drug}_biomarker_z', round(float(p.loc[p.Pathway == bm, drug].iloc[0]), 4))
for s in ['Figure 2b', 'Figure 2e', 'Figure 2g', 'Figure 2h', 'Figure 4b', 'Figure 4c']:
    d = raw[s]; h = [i for i in range(6) if str(d.iloc[i, 0]).startswith('Term')][0]
    add(f'{s}_n_pathways', d.shape[0] - h - 1); add(f'{s}_n_sample_columns', int(d.iloc[h, 1:].notna().sum()))
b = pd.read_excel(X, sheet_name='Figure 3e', header=3)
for c in b.columns: add(f'Bolis_{c}', int(b[c].count()))
for s, k in [('Figure 4e', 37), ('Figure 4f', 30)]:
    d = pd.read_excel(X, sheet_name=s, header=2)
    add(f'{s}_iterations', len(d)); add(f'{s}_rank1', int((d.top_rank == 1).sum())); add(f'{s}_rank<=2', int((d.top_rank <= 2).sum()))
    add(f'{s}_NA(top_rank>{k})', int((d.top_rank > k).sum()))
os.makedirs(os.path.join(ROOT, 'tables'), exist_ok=True)
pd.DataFrame(out, columns=['fact', 'value']).to_csv(os.path.join(ROOT, 'tables', 'source_data_facts.tsv'), sep='\t', index=False)
