"""Ovary: what can be done with 6-7 lines.  (1) IC50 table; (2) correlation of ln IC50 with every Reactome pathway score across the lines,
with an exact permutation test over all label orders (6! = 720); (3) a signature = mean z of the top-10 pathways positively correlated
minus the top-10 negatively correlated; tested in TCGA-OV platinum-treated patients (OS and PFS), GSE9891, GSE63885, and response cohorts;
(4) BRCA1/2 somatic mutation status in TCGA-OV vs the signature.   One CPU.   Output tables: tables/ovary_*.tsv"""
import sys, itertools, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cohort_tools as T
from cohort_tools import pd, np, u22, DATA
from u22 import OUT
from scipy.stats import spearmanr, rankdata, mannwhitneyu
warnings.filterwarnings("ignore")
TAB = OUT / "tables"
IMP = DATA / "loader_out" / "imports"
rng = np.random.default_rng(20260103)

def load(setname):
    S = pd.read_csv(IMP / setname / "pathway_scores.tsv.gz", sep="\t", index_col=0)
    S.index = S.index.astype(str)
    R = pd.read_csv(IMP / setname / "response.tsv", sep="\t")
    y = R.groupby(["sample", "drug"])["response"].median().unstack()
    raw = R.groupby(["sample", "drug"])["value_raw"].median().unstack()
    return S, y, raw

S6, y6, raw6 = load("ovary_gse229119_noMCF7")
S7, y7, raw7 = load("ovary_gse229119")
raw7.to_csv(TAB / "ovary_ic50_raw_7lines.tsv", sep="\t")
rho_lines = pd.DataFrame([dict(lines=n, rho=spearmanr(yy["OLAPARIB"], yy["NIRAPARIB"])[0], p=spearmanr(yy["OLAPARIB"], yy["NIRAPARIB"])[1], n=len(yy))
                          for n, yy in (("6 ovarian lines", y6), ("7 lines incl. MCF7 (breast)", y7))])
rho_lines.to_csv(TAB / "ovary_olaparib_vs_niraparib_rank_agreement.tsv", sep="\t", index=False)

def corr_all(S, y):
    X = S.loc[y.index]
    R = rankdata(X.values, axis=0)
    ry = rankdata(y.values)
    Rc = (R - R.mean(0)) / (R.std(0) + 1e-12)
    yc = (ry - ry.mean()) / ry.std()
    return pd.Series((Rc * yc[:, None]).mean(0), index=X.columns)

def top_sets(rho, k=10):
    r = rho.dropna().sort_values()
    return list(r.index[-k:]), list(r.index[:k])        # positive (higher score = higher IC50 = resistant), negative

rows, sigs = [], {}
exact = {}
for setname, S, y in (("6 ovarian lines", S6, y6), ("7 lines incl. MCF7", S7, y7)):
    for drug in ("OLAPARIB", "NIRAPARIB"):
        yy = y[drug].dropna()
        rho = corr_all(S, yy)
        n = len(yy)
        perms = list(itertools.permutations(range(n)))
        mx = []
        for pm in perms:
            yp = pd.Series(yy.values[list(pm)], index=yy.index)
            mx.append(np.abs(corr_all(S, yp)).max())
        mx = np.array(mx)
        pv = []
        for p_ in rho.index:
            pv.append(np.nan)
        rr = rho.copy()
        # nominal p of Spearman with n points (exact via permutation of the same set): fraction of orders with |rho| >= observed
        obs = rr.abs()
        fw = pd.Series([(np.sum(mx >= o - 1e-12)) / len(mx) for o in obs], index=obs.index)       # family-wise p (max over 671 pathways)
        pos, neg = top_sets(rho, 10)
        sigs[(setname, drug)] = (pos, neg)
        tab = pd.DataFrame({"rho": rho, "familywise_p": fw}).sort_values("rho")
        tab.assign(set=setname, drug=drug).to_csv(TAB / f"ovary_pathway_correlation_{setname.split()[0]}lines_{drug}.tsv", sep="\t")
        rows.append(dict(set=setname, drug=drug, n_lines=n, max_abs_rho=float(rho.abs().max()),
                         smallest_familywise_p=float(fw.min()), n_pathways=len(rho), n_orders=len(perms),
                         median_max_abs_rho_in_shuffles=float(np.median(mx)), share_of_shuffles_with_max_rho_1=float(np.mean(mx >= 1 - 1e-9))))
pd.DataFrame(rows).to_csv(TAB / "ovary_pathway_correlation_summary.tsv", sep="\t", index=False)

# HR / PARP-relevant pathways
kw = ("HOMOLOGY_DIRECTED", "HDR_THROUGH", "DNA_DOUBLE_STRAND", "DOUBLE_STRAND_BREAK_REPAIR", "BRCA", "PARP", "FANCONI", "BASE_EXCISION", "NONHOMOLOGOUS", "G2_M_DNA_DAMAGE", "DNA_REPAIR", "CELL_CYCLE_CHECKPOINTS")
sel = [p for p in S6.columns if any(k in p for k in kw)]
r6 = {d: corr_all(S6, y6[d].dropna()) for d in ("OLAPARIB", "NIRAPARIB")}
pd.DataFrame({"rho_olaparib": r6["OLAPARIB"].reindex(sel), "rho_niraparib": r6["NIRAPARIB"].reindex(sel)}).sort_values("rho_olaparib").to_csv(
    TAB / "ovary_dna_repair_pathways_vs_ic50.tsv", sep="\t")

# (3) signature in patient cohorts
P, SC = T.v2()
cols = [c for c in S6.columns if c in SC.columns]
cohorts = [("TCGA-OV", "OS", PLATREG := r"CISPLATIN|CARBOPLATIN|PLATINUM", r"^NONE$"), ("TCGA-OV", "PFS", PLATREG, r"^NONE$"),
           ("GSE9891", "OS", "PLATINUM", None), ("GSE63885", "OS", "PLATINUM", None)]
NPERM_ORD = 720

def sig_score(Zp, pos, neg):
    return Zp[pos].mean(1) - Zp[neg].mean(1)

out = []
for (setname, drug), (pos, neg) in sigs.items():
    for coh, ep, trx, cmx in cohorts:
        x = P[(P.cohort == coh) & ~P.duplicate_flag]
        Zp = T.zwithin(SC.loc[x.index, cols])
        dr = x["drugs"].fillna("")
        tr = dr.str.contains(trx, regex=True)
        cm = (dr.str.contains(cmx, regex=True) & ~tr) if cmx else pd.Series(False, index=x.index)
        s = sig_score(Zp, pos, neg)
        res = {}
        for lab, m in (("treated", tr), ("comparison", cm)):
            idx = x.index[m.values]
            if len(idx) < 10:
                continue
            Tt, Ee, ok, e_ = T.endpoint(x.loc[idx], (ep.lower(),))
            if Tt is None:
                continue
            ids = idx[ok.values]
            g = pd.Series(T.split(s.loc[idx].values), index=idx).loc[ids].values
            res[lab] = (len(Tt), int(Ee.sum()), T.logrank(Tt, Ee, g), T.cox_binary(Tt, Ee, g), ids, Tt, Ee)
        if "treated" not in res:
            continue
        n, ev, p, hr, ids, Tt, Ee = res["treated"]
        # exact permutation: every order of the line responses, signature rebuilt each time
        yy = (y6 if setname.startswith("6") else y7)[drug].dropna()
        Sx = (S6 if setname.startswith("6") else S7)
        idx_all = x.index[tr.values]
        pp = []
        for pm in itertools.permutations(range(len(yy))):
            yp = pd.Series(yy.values[list(pm)], index=yy.index)
            ps, ng = top_sets(corr_all(Sx, yp), 10)
            s2 = sig_score(Zp, ps, ng)
            g2 = pd.Series(T.split(s2.loc[idx_all].values), index=idx_all).loc[ids].values
            pp.append(T.logrank(Tt, Ee, g2))
        pp = np.array(pp)
        row = dict(set=setname, drug=drug, cohort=coh, endpoint=ep, n_treated=n, events_treated=ev, logrank_p=p, HR=hr[0], HR_lo=hr[1], HR_hi=hr[2],
                   permutation_p=float((np.sum(pp <= p + 1e-12)) / len(pp)), n_orders=len(pp))
        if "comparison" in res:
            n2, e2, p2, h2, ids2, T2, E2 = res["comparison"]
            row.update(n_comparison=n2, events_comparison=e2, logrank_p_comparison=p2, HR_comparison=h2[0], HR_comparison_lo=h2[1], HR_comparison_hi=h2[2])
        out.append(row)
        print(setname, drug, coh, ep, round(p, 3), flush=True)
pd.DataFrame(out).to_csv(TAB / "ovary_signature_survival.tsv", sep="\t", index=False)

# response cohorts
resp_rows = []
for (setname, drug), (pos, neg) in sigs.items():
    for coh in ("GSE15622", "GSE51373"):
        x = P[(P.cohort == coh) & ~P.duplicate_flag & P["responder"].notna() & (P["drugs"].fillna("").str.contains("CARBOPLATIN"))]
        if len(x) < 10:
            continue
        Zp = T.zwithin(SC.loc[x.index, cols])
        s = sig_score(Zp, pos, neg)
        r = x["responder"]
        a = s[r == 0].values; b = s[r == 1].values
        auc = float(np.mean([(ai > bj) + .5 * (ai == bj) for ai in a for bj in b]))
        resp_rows.append(dict(set=setname, drug=drug, cohort=coh, n=len(x), responders=int(r.sum()), AUROC_nonresponse=auc, mannwhitney_p=mannwhitneyu(a, b).pvalue))
pd.DataFrame(resp_rows).to_csv(TAB / "ovary_signature_response.tsv", sep="\t", index=False)

# (4) BRCA1/2 somatic status in TCGA-OV
m = pd.read_csv(REPO_MUT := u22.NAR / "02_organoid_only/BLCA_OV_LUAD_OTHER/_data/raw/TCGA-OV.somaticmutation_wxs.tsv.gz", sep="\t")
hit = m[m["gene"].isin(["BRCA1", "BRCA2"]) & ~m["effect"].str.contains("synonymous|intron|UTR|upstream|downstream|intergenic", na=False)]
brca_pat = set(hit["sample"].str[:12])
mut_pat = set(m["sample"].str[:12])
x = P[(P.cohort == "TCGA-OV") & ~P.duplicate_flag]
x = x[x["patient"].isin(mut_pat)]
x = x.assign(BRCA=x["patient"].isin(brca_pat).astype(int))
Zp = T.zwithin(SC.loc[x.index, cols])
b_rows = []
for (setname, drug), (pos, neg) in sigs.items():
    s = sig_score(Zp, pos, neg)
    a = s[x["BRCA"] == 1]; b = s[x["BRCA"] == 0]
    b_rows.append(dict(set=setname, drug=drug, n_patients_with_mutation_data=len(x), n_BRCA1_2_mutated=int(x["BRCA"].sum()),
                       mean_score_mutated=a.mean(), mean_score_wildtype=b.mean(), mannwhitney_p_two_sided=mannwhitneyu(a, b).pvalue,
                       AUROC_mutated_has_lower_score=float(np.mean([(bj > ai) + .5 * (bj == ai) for ai in a for bj in b]))))
pd.DataFrame(b_rows).to_csv(TAB / "ovary_signature_vs_BRCA_somatic.tsv", sep="\t", index=False)
# survival by BRCA status among platinum-treated
xt = P[(P.cohort == "TCGA-OV") & ~P.duplicate_flag]
xt = xt[xt["patient"].isin(mut_pat) & xt["drugs"].fillna("").str.contains(PLATREG)]
xt = xt.assign(BRCA=xt["patient"].isin(brca_pat).astype(int))
srows = []
for ep in ("os", "pfs"):
    Tt, Ee, ok, e_ = T.endpoint(xt, (ep,))
    g = xt.loc[ok.values, "BRCA"].values
    hr = T.cox_binary(Tt, Ee, g)
    srows.append(dict(endpoint=ep.upper(), n=len(Tt), events=int(Ee.sum()), n_BRCA_mut=int(g.sum()), logrank_p=T.logrank(Tt, Ee, g), HR_mut_vs_wt=hr[0], HR_lo=hr[1], HR_hi=hr[2]))
pd.DataFrame(srows).to_csv(TAB / "ovary_BRCA_somatic_survival_platinum.tsv", sep="\t", index=False)
print("ovary done")
