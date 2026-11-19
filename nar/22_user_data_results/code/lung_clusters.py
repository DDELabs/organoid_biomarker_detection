"""Lung: the drug table holds ONE value per cluster and drug (C1-C4).  This script
 (1) shows what the clusters look like (size, histology, pathway scores),
 (2) tests whether pathway scores can tell the clusters apart (cross-validated, with label shuffles),
 (3) predicts each patient's cluster in TCGA-LUAD / TCGA-LUSC and the lung GEO cohorts, turns it into an expected drug value
     (sum over clusters of probability x cluster value) and tests it against survival in treated patients and in the comparison group.
Output tables: tables/lung_cluster_*.tsv.  One CPU."""
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cohort_tools as T
from cohort_tools import pd, np, u22, DATA
from u22 import OUT
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from scipy.stats import kruskal
warnings.filterwarnings("ignore")
TAB = OUT / "tables"
rng = np.random.default_rng(20260101)
imp = DATA / "loader_out" / "imports" / "luad_gse165611"
S = pd.read_csv(imp / "pathway_scores.tsv.gz", sep="\t", index_col=0)
S.index = S.index.astype(str).str.upper()
cl = pd.read_csv(TAB / "luad_cluster_labels.tsv", sep="\t", index_col=0)["cluster"]
cl.index = cl.index.str.upper()
hist = pd.read_csv(TAB / "luad_histology.tsv", sep="\t", index_col=0)["histology"]
hist.index = hist.index.str.upper()
cl = cl.loc[S.index]
vals = pd.read_csv(TAB / "luad_cluster_by_drug_value_listed.tsv", sep="\t", index_col=0)      # listed values only
vals0 = vals.reindex(columns=["C1", "C2", "C3", "C4"]).fillna(0.0)                             # user's Excel: 0 where nothing is listed

# (1) clusters
comp = pd.crosstab(cl, hist.loc[cl.index]).assign(n=lambda d: d.sum(1))
comp.to_csv(TAB / "lung_cluster_histology.tsv", sep="\t")

# (2) pathway scores vs cluster
Z = (S - S.mean()) / S.std().replace(0, np.nan)
Z = Z.fillna(0.0)
kw = []
for p in Z.columns:
    g = [Z.loc[cl == c, p].values for c in ("C1", "C2", "C3", "C4")]
    h, pv = kruskal(*g)
    kw.append((p, h, pv, *[x.mean() for x in g]))
kw = pd.DataFrame(kw, columns=["pathway", "H", "p", "mean_z_C1", "mean_z_C2", "mean_z_C3", "mean_z_C4"]).sort_values("p")
kw["q_BH"] = np.minimum.accumulate((kw["p"] * len(kw) / np.arange(1, len(kw) + 1))[::-1])[::-1].clip(upper=1)
kw.to_csv(TAB / "lung_cluster_pathway_kruskal.tsv", sep="\t", index=False)


def cv_acc(Zm, y, seed, C_=0.05):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    pred = np.empty(len(y), dtype=object)
    for tr, te in skf.split(Zm, y):
        m = LogisticRegression(C=C_, max_iter=300).fit(Zm[tr], y[tr])
        pred[te] = m.predict(Zm[te])
    cls = np.unique(y)
    recalls = [np.mean(pred[y == c] == c) for c in cls]
    return float(np.mean(pred == y)), float(np.mean(recalls))


y = cl.values.astype(str)
acc, bal = cv_acc(Z.values, y, 1)
accs = [cv_acc(Z.values, y, s) for s in range(2, 6)]
perm = []
for i in range(100):
    yp = rng.permutation(y)
    perm.append(cv_acc(Z.values, yp, 1)[1])
perm = np.array(perm)
pd.DataFrame([dict(model="multinomial logistic (L2, C=0.05), 5-fold CV, 5 repeats", n=len(y), accuracy=np.mean([acc] + [a[0] for a in accs]),
                   balanced_accuracy=np.mean([bal] + [a[1] for a in accs]), chance_balanced_accuracy=0.25,
                   shuffled_balanced_accuracy_mean=perm.mean(), shuffled_balanced_accuracy_p95=np.quantile(perm, .95),
                   permutation_p=float((np.sum(perm >= np.mean([bal] + [a[1] for a in accs])) + 1) / (len(perm) + 1)), n_shuffles=len(perm))]).to_csv(
    TAB / "lung_cluster_classifier_cv.tsv", sep="\t", index=False)
final = LogisticRegression(C=0.05, max_iter=500).fit(Z.values, y)
classes = list(final.classes_)

# (3) patients
P, SC = T.v2()
rows = []
COH = {  # cohort -> (treated regex, comparison regex, label)
    "TCGA-LUAD": [("GEMCITABINE", "GEMCITABINE", "^NONE$"), ("PACLITAXEL", "PACLITAXEL", "^NONE$"), ("ETOPOSIDE", "ETOPOSIDE", "^NONE$"),
                  ("PLATINUM(proxy)", "CISPLATIN|CARBOPLATIN|PLATINUM", "^NONE$")],
    "TCGA-LUSC": [("GEMCITABINE", "GEMCITABINE", "^NONE$"), ("PACLITAXEL", "PACLITAXEL", "^NONE$"), ("ETOPOSIDE", "ETOPOSIDE", "^NONE$"),
                  ("PLATINUM(proxy)", "CISPLATIN|CARBOPLATIN|PLATINUM", "^NONE$")],
    "GSE14814": [("PLATINUM(proxy)", "CISPLATIN", "^NONE$")],
    "GSE29013": [("PLATINUM(proxy)", "PLATINUM", "^NONE$")], "GSE68465": [("PLATINUM(proxy)", "PLATINUM", "^NONE$")],
    "GSE37745": [("PLATINUM(proxy)", "PLATINUM", "^NONE$")], "GSE42127": [("PLATINUM(proxy)", "PLATINUM", "^NONE$")],
}
# which cluster-valued drug stands for which regimen
MODEL_DRUGS = {"GEMCITABINE": ["Gemcitabine"], "PACLITAXEL": ["Paclitaxel"], "ETOPOSIDE": ["Etoposide"],
               "PLATINUM(proxy)": ["Gemcitabine", "Etoposide", "Paclitaxel", "Camptothecin"]}
PERM = 500
for coh, defs in COH.items():
    x = P[(P["cohort"] == coh) & ~P["duplicate_flag"]]
    Sx = SC.loc[x.index]
    Zp = T.zwithin(Sx)
    Zp = Zp.reindex(columns=Z.columns).fillna(0.0)
    prob = pd.DataFrame(final.predict_proba(Zp.values), index=Zp.index, columns=classes)
    prob.groupby(lambda i: coh).mean().assign(cohort=coh, n=len(prob)).to_csv(TAB / "lung_cluster_predicted_share.tsv", sep="\t", mode="a", header=not (TAB / "lung_cluster_predicted_share.tsv").exists() or coh == "TCGA-LUAD")
    drugs_s = x["drugs"].fillna("")
    for reg, trx, cmx in defs:
        tr = drugs_s.str.contains(trx, regex=True)
        if reg == "PLATINUM(proxy)":
            pass
        cmpm = drugs_s.str.contains(cmx, regex=True) & ~tr
        for dn in MODEL_DRUGS[reg]:
            vv = vals0.loc[dn]
            # shuffles of the organoid labels (cluster -> drug value link kept, expression link broken)
            def score_from(model):
                pr = pd.DataFrame(model.predict_proba(Zp.values), index=Zp.index, columns=list(model.classes_))
                return pr.reindex(columns=["C1", "C2", "C3", "C4"]).fillna(0.0).values @ vv.values
            sc = pd.Series(score_from(final), index=Zp.index)
            res = {}
            for lab, mask in (("treated", tr), ("comparison", cmpm)):
                idx = x.index[mask.values]
                if len(idx) < 10:
                    continue
                Tt, Ee, ok, ep = T.endpoint(x.loc[idx], ("os",))
                if Tt is None:
                    continue
                ids = idx[ok.values]
                g = T.split(sc.loc[idx].values)
                g = pd.Series(g, index=idx).loc[ids].values
                hr = T.cox_binary(Tt, Ee, g)
                res[lab] = dict(n=len(Tt), events=int(Ee.sum()), logrank_p=T.logrank(Tt, Ee, g), HR=hr[0], HR_lo=hr[1], HR_hi=hr[2], ids=ids, Tt=Tt, Ee=Ee)
            if "treated" not in res:
                continue
            obs = res["treated"]["logrank_p"]
            ids, Tt, Ee = res["treated"]["ids"], res["treated"]["Tt"], res["treated"]["Ee"]
            idx_all = x.index[tr.values]
            pp = []
            for _ in range(PERM):
                m = LogisticRegression(C=0.05, max_iter=200).fit(Z.values, rng.permutation(y))
                s2 = pd.Series(score_from(m), index=Zp.index)
                g2 = pd.Series(T.split(s2.loc[idx_all].values), index=idx_all).loc[ids].values
                pp.append(T.logrank(Tt, Ee, g2))
            pperm = float((np.sum(np.array(pp) <= obs + 1e-12) + 1) / (PERM + 1))
            inter = (np.nan,) * 4
            if "comparison" in res:
                both = x.loc[list(res["treated"]["ids"]) + list(res["comparison"]["ids"])]
                Tb, Eb, okb, _ = T.endpoint(both, ("os",))
                gb = (sc.loc[both.index[okb.values]] > sc.loc[list(idx_all) + list(x.index[cmpm.values])].median()).astype(int).values
                armb = both.index[okb.values].isin(res["treated"]["ids"]).astype(int)
                inter = T.interaction(Tb, Eb, gb, armb)
            row = dict(cohort=coh, regimen=reg, organoid_drug=dn.upper(), endpoint="OS",
                       n_treated=res["treated"]["n"], events_treated=res["treated"]["events"], logrank_p_treated=obs,
                       HR_treated=res["treated"]["HR"], HR_lo=res["treated"]["HR_lo"], HR_hi=res["treated"]["HR_hi"],
                       permutation_p=pperm, n_shuffles=PERM)
            if "comparison" in res:
                row.update(n_comparison=res["comparison"]["n"], events_comparison=res["comparison"]["events"],
                           logrank_p_comparison=res["comparison"]["logrank_p"], HR_comparison=res["comparison"]["HR"],
                           HR_comparison_lo=res["comparison"]["HR_lo"], HR_comparison_hi=res["comparison"]["HR_hi"],
                           interaction_HR=inter[0], interaction_lo=inter[1], interaction_hi=inter[2], interaction_p=inter[3])
            rows.append(row)
            print(row["cohort"], reg, dn, row["n_treated"], round(obs, 3), pperm, flush=True)
pd.DataFrame(rows).to_csv(TAB / "lung_cluster_expected_value_survival.tsv", sep="\t", index=False)
