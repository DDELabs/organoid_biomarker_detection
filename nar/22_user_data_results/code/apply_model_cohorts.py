"""Organoid-only (paper-style) model applied to the patient cohorts of training_table_v2 (rank-scored pathways, z-scored inside each cohort).
 plans:  lung | ovary
Score = Ridge (RidgeCV, as the paper) on the pathways close to the drug's targets; 'rank1' = the paper's single top pathway,
'all' = every close pathway.  High score = predicted LESS sensitive.  Patients are split at the median of their own arm.
Permutation: the organoid response is shuffled 500 times, the model refitted (including the choice of the top pathway) and the treated-arm
log-rank p recomputed; p_perm = (number of shuffles with p <= observed + 1) / (500 + 1).
Output: tables/apply_<plan>.tsv"""
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cohort_tools as T
from cohort_tools import pd, np, u22, DATA
from u22 import OUT
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")
plan = sys.argv[1]
NPERM = 500
rng = np.random.default_rng(20260102)
TAB = OUT / "tables"
IMP = DATA / "loader_out" / "imports"

PLATLUNG = r"CISPLATIN|CARBOPLATIN|PLATINUM"
CFG = {
    "lung": dict(set="luad_gse165611", resp_col="response",
        # (organoid drug, patient regimen label, proximity drug)
        drugs=[("GEMCITABINE", "GEMCITABINE", "GEMCITABINE"), ("PACLITAXEL", "PACLITAXEL", "PACLITAXEL"), ("ETOPOSIDE", "ETOPOSIDE", "ETOPOSIDE"),
               ("GEMCITABINE", "PLATINUM(proxy)", "CISPLATIN"), ("ETOPOSIDE", "PLATINUM(proxy)", "CISPLATIN"),
               ("PACLITAXEL", "PLATINUM(proxy)", "CISPLATIN"), ("CAMPTOTHECIN", "PLATINUM(proxy)", "CISPLATIN")],
        cohorts=[  # cohort, regimen label, treated regex, comparison regex, kind
            ("TCGA-LUAD", "GEMCITABINE", "GEMCITABINE", "^NONE$", "surv"), ("TCGA-LUSC", "GEMCITABINE", "GEMCITABINE", "^NONE$", "surv"),
            ("TCGA-LUAD", "PACLITAXEL", "PACLITAXEL", "^NONE$", "surv"), ("TCGA-LUSC", "PACLITAXEL", "PACLITAXEL", "^NONE$", "surv"),
            ("TCGA-LUAD", "ETOPOSIDE", "ETOPOSIDE", "^NONE$", "surv"), ("TCGA-LUSC", "ETOPOSIDE", "ETOPOSIDE", "^NONE$", "surv"),
            ("TCGA-LUAD", "PLATINUM(proxy)", PLATLUNG, "^NONE$", "surv"), ("TCGA-LUSC", "PLATINUM(proxy)", PLATLUNG, "^NONE$", "surv"),
            ("GSE14814", "PLATINUM(proxy)", "CISPLATIN", "^NONE$", "surv"), ("GSE29013", "PLATINUM(proxy)", "PLATINUM", "^NONE$", "surv"),
            ("GSE68465", "PLATINUM(proxy)", "PLATINUM", "^NONE$", "surv"), ("GSE37745", "PLATINUM(proxy)", "PLATINUM", "^NONE$", "surv"),
            ("GSE42127", "PLATINUM(proxy)", "PLATINUM", "^NONE$", "surv"),
            ("GSE207422", "PLATINUM(proxy)", "CARBOPLATIN", None, "resp")]),
    "ovary": dict(set="ovary_gse229119_noMCF7", resp_col="response",
        drugs=[("OLAPARIB", "PLATINUM", "OLAPARIB"), ("NIRAPARIB", "PLATINUM", "NIRAPARIB"),
               ("OLAPARIB", "PLATINUM(proxy prox)", "CISPLATIN"), ("NIRAPARIB", "PLATINUM(proxy prox)", "CISPLATIN"),
               ("OLAPARIB", "PACLITAXEL", "OLAPARIB"), ("NIRAPARIB", "PACLITAXEL", "NIRAPARIB")],
        cohorts=[("TCGA-OV", "PLATINUM", PLATLUNG, "^NONE$", "surv"), ("GSE9891", "PLATINUM", "PLATINUM", None, "surv"),
                 ("GSE63885", "PLATINUM", "PLATINUM", None, "surv"),
                 ("TCGA-OV", "PLATINUM(proxy prox)", PLATLUNG, "^NONE$", "surv"), ("GSE9891", "PLATINUM(proxy prox)", "PLATINUM", None, "surv"),
                 ("GSE63885", "PLATINUM(proxy prox)", "PLATINUM", None, "surv"),
                 ("TCGA-OV", "PACLITAXEL", "PACLITAXEL", "^NONE$", "surv"), ("GSE9891", "PACLITAXEL", "PACLITAXEL", "^PLATINUM$", "surv"),
                 ("GSE63885", "PACLITAXEL", "PACLITAXEL", "CYCLOPHOSPHAMIDE", "surv"),
                 ("GSE15622", "PLATINUM", "CARBOPLATIN", None, "resp"), ("GSE51373", "PLATINUM", "CARBOPLATIN", None, "resp")]),
}
CFG["lung_listed"] = dict(CFG["lung"], set="luad_gse165611_listed",
    drugs=[("GEMCITABINE", "GEMCITABINE", "GEMCITABINE"), ("GEMCITABINE", "PLATINUM(proxy)", "CISPLATIN"), ("PACLITAXEL", "PACLITAXEL", "PACLITAXEL")])
CFG["ovary_7lines"] = dict(CFG["ovary"], set="ovary_gse229119", drugs=[("OLAPARIB", "PLATINUM", "OLAPARIB"), ("NIRAPARIB", "PLATINUM", "NIRAPARIB")])
cfg = CFG[plan]
S = pd.read_csv(IMP / cfg["set"] / "pathway_scores.tsv.gz", sep="\t", index_col=0)
S.index = S.index.astype(str)
R = pd.read_csv(IMP / cfg["set"] / "response.tsv", sep="\t")
R["sample"] = R["sample"].astype(str)
Pv, SC = T.v2()
cols = [c for c in S.columns if c in SC.columns]
rows = []
for odrug, regimen, proxdrug in cfg["drugs"]:
    y = R[R["drug"] == odrug].groupby("sample")["response"].median()
    y = y[y.index.isin(S.index)]
    if y.nunique() < 2:
        continue
    feats, pcol = T.shortlist(proxdrug, cols)
    if len(feats) < 2:
        continue
    coef = T.organoid_fit(S[cols], y, feats)
    top = coef.abs().idxmax()
    # null: shuffled response
    nulls = [T.organoid_fit(S[cols], pd.Series(rng.permutation(y.values), index=y.index), feats) for _ in range(NPERM)]
    for coh, reg, trx, cmx, kind in cfg["cohorts"]:
        if reg != regimen:
            continue
        x = Pv[(Pv["cohort"] == coh) & ~Pv["duplicate_flag"]]
        if not len(x):
            continue
        Zp = T.zwithin(SC.loc[x.index, cols])
        dr = x["drugs"].fillna("")
        tr = dr.str.contains(trx, regex=True)
        cm = (dr.str.contains(cmx, regex=True) & ~tr) if cmx else pd.Series(False, index=x.index)
        for how in ("rank1", "all"):
            sc_fn = (lambda c: T.rank1_score(Zp, c)[0]) if how == "rank1" else (lambda c: T.multi_score(Zp, c))
            sc = sc_fn(coef)
            row = dict(set=cfg["set"], organoid_drug=odrug, regimen=regimen, proximity=pcol, n_organoids=len(y), n_close_pathways=len(feats),
                       top_pathway=top, score=how, cohort=coh)
            if kind == "surv":
                res = {}
                for lab, mask in (("treated", tr), ("comparison", cm)):
                    idx = x.index[mask.values]
                    if len(idx) < 10:
                        continue
                    pref = ("os",) if coh.startswith("TCGA") or coh in ("GSE9891",) else ("os",)
                    Tt, Ee, ok, ep = T.endpoint(x.loc[idx], pref)
                    if Tt is None:
                        continue
                    ids = idx[ok.values]
                    g = pd.Series(T.split(sc.loc[idx].values), index=idx).loc[ids].values
                    hr = T.cox_binary(Tt, Ee, g)
                    res[lab] = (len(Tt), int(Ee.sum()), T.logrank(Tt, Ee, g), hr, ids, Tt, Ee)
                if "treated" not in res:
                    continue
                n, ev, p, hr, ids, Tt, Ee = res["treated"]
                idx_all = x.index[tr.values]
                pp = []
                for cn in nulls:
                    s2 = sc_fn(cn)
                    g2 = pd.Series(T.split(s2.loc[idx_all].values), index=idx_all).loc[ids].values
                    pp.append(T.logrank(Tt, Ee, g2))
                row.update(outcome="OS", n_treated=n, events_treated=ev, logrank_p=p, HR=hr[0], HR_lo=hr[1], HR_hi=hr[2],
                           permutation_p=float((np.sum(np.array(pp) <= p + 1e-12) + 1) / (NPERM + 1)), n_shuffles=NPERM)
                if "comparison" in res:
                    n2, e2, p2, h2, ids2, T2, E2 = res["comparison"]
                    both = x.loc[list(ids) + list(ids2)]
                    Tb, Eb, okb, _ = T.endpoint(both, ("os",))
                    med = sc.loc[list(idx_all) + list(x.index[cm.values])].median()
                    gb = (sc.loc[both.index[okb.values]] > med).astype(int).values
                    ar = both.index[okb.values].isin(ids).astype(int)
                    it = T.interaction(Tb, Eb, gb, ar)
                    row.update(n_comparison=n2, events_comparison=e2, logrank_p_comparison=p2, HR_comparison=h2[0], HR_comparison_lo=h2[1],
                               HR_comparison_hi=h2[2], interaction_HR=it[0], interaction_lo=it[1], interaction_hi=it[2], interaction_p=it[3])
            else:      # response cohorts: AUROC of the 'resistance' score for non-response
                xt = x.loc[tr.values & x["responder"].notna().values]
                if len(xt) < 10 or xt["responder"].sum() < 3 or (1 - xt["responder"]).sum() < 3:
                    continue
                s = sc.loc[xt.index]
                def auc(s_):
                    a = s_[xt["responder"] == 0].values
                    b = s_[xt["responder"] == 1].values
                    return float(np.mean([(ai > bj) + 0.5 * (ai == bj) for ai in a for bj in b]))
                a0 = auc(s)
                mw = mannwhitneyu(s[xt["responder"] == 0], s[xt["responder"] == 1]).pvalue
                nul = [auc(sc_fn(cn).loc[xt.index]) for cn in nulls]
                pperm = float((np.sum(np.abs(np.array(nul) - 0.5) >= abs(a0 - 0.5) - 1e-12) + 1) / (NPERM + 1))
                row.update(outcome="response", n_treated=len(xt), events_treated=int(xt["responder"].sum()), AUROC_nonresponse=a0,
                           logrank_p=mw, permutation_p=pperm, n_shuffles=NPERM)
            rows.append(row)
    print(odrug, regimen, flush=True)
pd.DataFrame(rows).to_csv(TAB / f"apply_{plan}.tsv", sep="\t", index=False)
print("rows", len(rows))
