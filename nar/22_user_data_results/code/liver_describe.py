"""Liver (GSE188541): describe what the files hold and do only what is valid.  Output: tables/liver_*.tsv"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u22
from u22 import C, OUT, pd, np
from scipy.stats import spearmanr, kendalltau
TAB = OUT / "tables"
d = pd.read_csv(TAB / "lihc_ic50_per_organoid.tsv", sep="\t", index_col=0)
rows = []
for a in d.columns:
    for b in d.columns:
        if a < b:
            r = spearmanr(d[a], d[b])
            rows.append(dict(organoid_a=a, organoid_b=b, n_drugs=len(d), spearman_rho=r[0], p=r[1]))
pd.DataFrame(rows).to_csv(TAB / "liver_organoid_pair_rank_agreement.tsv", sep="\t", index=False)
st = pd.DataFrame({"min_IC50": d.min(1), "max_IC50": d.max(1), "fold_range_max_over_min": d.max(1) / d.min(1), "most_sensitive_organoid": d.idxmin(1), "least_sensitive_organoid": d.idxmax(1)})
st.to_csv(TAB / "liver_drug_spread.tsv", sep="\t")
lr = (np.log10(d)).sub(np.log10(d).mean(1), axis=0)
lr.to_csv(TAB / "liver_log10_ic50_centred_per_drug.tsv", sep="\t")
# patients: how many TCGA-LIHC patients got any of these drugs
t = C.treatment_table("LIHC")
vc = t["drug"].value_counts()
names = [x.upper() for x in d.index]
alias = {"FIOXURIDINE": ["FLOXURIDINE", "FUDR"], "NEVARAPINE": ["NEVIRAPINE"]}
cnt = []
for n in names:
    cands = [n] + alias.get(n, [])
    cnt.append(dict(drug=n, tcga_lihc_patients=int(sum(vc.get(c, 0) for c in cands))))
pd.DataFrame(cnt).to_csv(TAB / "liver_drugs_in_tcga_lihc.tsv", sep="\t", index=False)
pd.DataFrame([dict(item="series matrix rows with expression", value=0),
              dict(item="GSE188541 samples", value=6), dict(item="GSE188541 organoid lines", value="1 (iPSC-derived liver organoids on chips; control, APAP, FIAU, tenofovir, tenofovir+inarigivir, plus un-chipped organoids)"),
              dict(item="organoids with IC50 in the uploaded table", value="3 (72.3, 2E, CC3)"),
              dict(item="organoids with both expression and IC50", value=0)]).to_csv(TAB / "liver_data_inventory.tsv", sep="\t", index=False)
print(pd.DataFrame(cnt).to_string()); print(pd.DataFrame(rows))
