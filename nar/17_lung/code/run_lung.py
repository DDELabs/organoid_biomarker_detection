"""W1 step 6: lung analyses a-d on robust/results/training_table_v2.

a. response, leave-one-cohort-out, lung only (immunotherapy group; chemotherapy group TCGA-LUAD/LUSC), plus the frozen
   breast model and pan-cancer models applied to each lung response cohort;
b. adjuvant-chemo benefit: treated x score interaction (Cox, stratified by cohort), learned on the other lung cohorts and
   tested on the held-out cohort (GSE14814 is the randomised one);
c. the same score in patients who did not get chemotherapy (is it prognostic rather than drug-specific?);
d. (b) and (c) again with the paper's network-proximity shortlist for cisplatin instead of all pathways.
Every number in REPORT.md comes from the tables written here.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "robust"))
warnings.filterwarnings("ignore")
import run_train_combined as TC  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import network as N  # noqa: E402
from obd import survlearn as SL  # noqa: E402

TC.T = RESULTS / "training_table_v2"
TAB = ROOT / "nar/17_lung/tables"
TAB.mkdir(parents=True, exist_ok=True)
SEED = 20260101
np.random.seed(SEED)
BASE = ["age_n", "stage_n"]
STAGE = {"I": 1, "II": 2, "III": 3, "IV": 4}
N_RANDOM = 100          # random pathway sets of the same size as the shortlist


def ci_hr(coef, se):
    return np.exp(coef), np.exp(coef - 1.96 * se), np.exp(coef + 1.96 * se)


# ------------------------------------------------------------------------------------------------ a. response
def a_response(P, Z):
    R = P[P["has_response"] & P["treated"].eq(1)]
    lung = R[R["cancer_group"] == "lung"]
    group = {}
    for c, d in lung.groupby("cohort"):
        if d["drugs"].str.contains("ICI").all():
            group[c] = "immunotherapy (anti-PD-1/PD-L1)"
        elif c.startswith("TCGA-LU"):
            group[c] = "chemotherapy (TCGA lung)"
        else:
            group[c] = "targeted (BATTLE: erlotinib/sorafenib)"
    rows = []

    def one(model, grp, c, tr_idx):
        te = lung.index[lung["cohort"] == c]
        y = lung.loc[te, "responder"].astype(int)
        ytr = R.loc[tr_idx, "responder"].astype(int)
        if ytr.nunique() < 2 or ytr.value_counts().min() < 5 or y.value_counts().min() < 5 or len(y.unique()) < 2:
            return
        s = TC.fit(Z.loc[tr_idx].values, ytr).decision_function(Z.loc[te].values)
        a, se = TC.auc_se(y.values, s)
        rows.append(dict(group=grp, test_cohort=c, model=model, n=len(y), responders=int(y.sum()), train_n=len(tr_idx),
                         train_cohorts=";".join(sorted(R.loc[tr_idx, "cohort"].unique())) if R.loc[tr_idx, "cohort"].nunique() < 6
                         else f"{R.loc[tr_idx, 'cohort'].nunique()} cohorts", auc=a, ci_lo=a - 1.96 * se, ci_hi=a + 1.96 * se,
                         p=TC.auc(y.values, s)[1]))
    for c, grp in group.items():
        peers = [k for k, g in group.items() if g == grp and k != c]
        if peers:
            one("same group, leave one cohort out", grp, c, R.index[R["cohort"].isin(peers)])
        one("breast model (trained on all breast response cohorts)", grp, c, R.index[R["cancer_group"] == "breast"])
        one("pan-cancer model (all other cohorts, lung included)", grp, c, R.index[R["cohort"] != c])
        one("pan-cancer model without any lung cohort", grp, c, R.index[R["cancer_group"] != "lung"])
    T = pd.DataFrame(rows)
    T.to_csv(TAB / "lung_loco_response.tsv", sep="\t", index=False)
    return T


# ------------------------------------------------------------------------------------------------ survival frames
def frame(P, cohorts=None, treated_known=True, endpoint="os"):
    d = P[P["cancer_group"] == "lung"].copy()
    if cohorts is not None:
        d = d[d["cohort"].isin(cohorts)]
    if treated_known:
        d = d[d["treated"].notna()]
    ep = np.where(d["cohort"].eq("GSE8894"), "rfs", endpoint)       # GSE8894 only has relapse-free survival
    d["time"] = np.where(ep == "rfs", d["rfs_months"], d[f"{endpoint}_months"])
    d["event"] = np.where(ep == "rfs", d["rfs_event"], d[f"{endpoint}_event"])
    d["endpoint"] = np.where(ep == "rfs", "relapse-free survival", "overall survival")
    d = d[(d["time"] > 0) & d["event"].notna()]
    st = d["stage"].astype(str).str.upper().str.replace("STAGE", "").str.strip().str.extract(r"^(IV|III|II|I)")[0].map(STAGE)
    d = d[~(st == 4)]                                   # stage IV: not an adjuvant setting
    st = st.loc[d.index]
    d["stage_n"] = (st - 2).fillna(0.0)
    d["age_n"] = ((pd.to_numeric(d["age"], errors="coerce") - 65) / 10).fillna(0)
    d["project"] = d["cohort"]
    return d


def cox_row(df, cols, term):
    r = SL.cox_test(df, cols)
    hr, lo, hi = ci_hr(r.loc[term, "coef"], r.loc[term, "se"])
    return dict(n=r.attrs["n"], events=r.attrs["events"], hr=hr, ci_lo=lo, ci_hi=hi, p=r.loc[term, "p"],
                coef=r.loc[term, "coef"], se=r.loc[term, "se"])


def learn(tr, Z, cols):
    alpha, gamma = SL.fit_interaction(tr, Z.loc[tr.index, cols], BASE, 300.0, 300.0)
    return alpha, gamma


def std(s):
    return (s - s.mean()) / s.std()


def pooled(rows):
    """Fixed-effect (inverse-variance) pool of log hazard ratios."""
    r = pd.DataFrame(rows)
    w = 1 / r["se"] ** 2
    b = (w * r["coef"]).sum() / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(n=int(r["n"].sum()), events=int(r["events"].sum()), hr=np.exp(b), ci_lo=np.exp(b - 1.96 * se), ci_hi=np.exp(b + 1.96 * se),
                p=2 * norm.sf(abs(b / se)), coef=b, se=se)


# ------------------------------------------------------------------------------------------------ b, d. adjuvant interaction
def b_adjuvant(P, Z, sets):
    ADJ = [c for c, d in P[(P["cancer_group"] == "lung") & P["treated"].notna()].groupby("cohort")
           if d["treated"].nunique() == 2 and d["os_event"].notna().any()]
    allfr = frame(P, ADJ)
    rows, scores, rng = [], {}, np.random.default_rng(SEED)
    for set_name, cols in sets.items():
        per = []
        for c in ADJ:
            te = allfr[allfr["cohort"] == c]
            tr = allfr[allfr["cohort"] != c]
            n_t, n_u, ev = int(te["treated"].sum()), int((1 - te["treated"]).sum()), int(te["event"].sum())
            if min(n_t, n_u, ev) < 15:
                rows.append(dict(setting="adjuvant chemotherapy", pathways=set_name, test_cohort=c, n=len(te), treated=n_t, untreated=n_u,
                                 events=ev, note="not tested: needs at least 15 treated, 15 untreated and 15 events"))
                continue
            alpha, gamma = learn(tr, Z, cols)
            te = te.assign(score=-(Z.loc[te.index, cols].values @ gamma.values))     # high = more benefit predicted
            te["score"] = std(te["score"]); te["tx_score"] = te["treated"] * te["score"]
            r = cox_row(te, BASE + ["treated", "score", "tx_score"], "tx_score")
            hi = te[te["score"] >= te["score"].median()]; lo = te[te["score"] < te["score"].median()]
            h = cox_row(hi, BASE + ["treated"], "treated"); l = cox_row(lo, BASE + ["treated"], "treated")
            allt = cox_row(te, BASE + ["treated"], "treated")
            row = dict(setting="adjuvant chemotherapy", pathways=set_name, test_cohort=c, randomised=bool(P.loc[te.index, "randomised"].iloc[0]),
                       endpoint="overall survival", n=r["n"], treated=n_t, untreated=n_u, events=r["events"], train_n=len(tr),
                       interaction_hr=r["hr"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p"], coef=r["coef"], se=r["se"],
                       treated_vs_untreated_hr_all=allt["hr"], all_lo=allt["ci_lo"], all_hi=allt["ci_hi"], all_p=allt["p"],
                       hr_high_benefit_half=h["hr"], high_lo=h["ci_lo"], high_hi=h["ci_hi"], high_p=h["p"], high_n=h["n"], high_events=h["events"],
                       hr_low_benefit_half=l["hr"], low_lo=l["ci_lo"], low_hi=l["ci_hi"], low_p=l["p"], low_n=l["n"], low_events=l["events"])
            if set_name.startswith("cisplatin"):                       # compare with random pathway sets of the same size
                nul = []
                for _ in range(N_RANDOM):
                    rc = list(rng.choice(Z.columns, len(cols), replace=False))
                    _, g2 = learn(tr, Z, rc)
                    t2 = te.assign(score=-(Z.loc[te.index, rc].values @ g2.values)); t2["score"] = std(t2["score"]); t2["tx_score"] = t2["treated"] * t2["score"]
                    nul.append(SL.cox_test(t2, BASE + ["treated", "score", "tx_score"]).loc["tx_score", "coef"])
                nul = np.array(nul)
                row["random_sets"] = N_RANDOM
                row["share_random_sets_with_more_benefit_coef"] = float((nul <= r["coef"]).mean())
            rows.append(row)
            per.append(row)
            scores[(set_name, c)] = te["score"]
            print("b/d", set_name, c, round(row["interaction_hr"], 3), round(row["p"], 4), flush=True)
        if per:
            pl = pooled(per)
            rows.append(dict(setting="adjuvant chemotherapy", pathways=set_name, test_cohort="ALL TESTED COHORTS (inverse-variance pool)",
                             randomised=False, endpoint="overall survival", n=pl["n"], events=pl["events"], interaction_hr=pl["hr"], ci_lo=pl["ci_lo"],
                             ci_hi=pl["ci_hi"], p=pl["p"], coef=pl["coef"], se=pl["se"],
                             note="pool of held-out tests; the training sets overlap, so treat the pooled p-value as approximate"))
    T = pd.DataFrame(rows)
    T.to_csv(TAB / "lung_adjuvant_interaction.tsv", sep="\t", index=False)
    return T, scores, ADJ, allfr


# ------------------------------------------------------------------------------------------------ c. untreated prognosis
def c_untreated(P, Z, sets, ADJ, allfr):
    lungc = P[P["cancer_group"] == "lung"]
    cand = []
    for c, d in lungc.groupby("cohort"):
        if d["os_event"].notna().sum() + d["rfs_event"].notna().sum() < 15:
            continue
        if d["drugs"].fillna("").str.contains("ICI|ERLOTINIB|SORAFENIB").any():
            continue                                    # immunotherapy / targeted cohorts are not chemo-free
        cand.append(c)
    un = frame(P, cand, treated_known=False)
    un = un[~un["treated"].eq(1)]                       # untreated (0) or not recorded (NaN)
    rec = un.groupby("cohort")["treated"].transform(lambda s: s.notna().any())
    un = un[~(rec & un["treated"].isna())]              # where treatment is recorded, drop patients with 'unknown'
    rows = []
    for set_name, cols in sets.items():
        per = {"benefit score (drug-specific part)": [], "prognostic score (main-effect part)": []}
        for c in sorted(un["cohort"].unique()):
            te = un[un["cohort"] == c]
            kind = "no chemotherapy (recorded)" if te["treated"].notna().all() else "treatment not recorded (surgery cohort)"
            tr = allfr[allfr["cohort"] != c]
            ev = int(te["event"].sum())
            if ev < 15:
                rows.append(dict(pathways=set_name, test_cohort=c, group=kind, n=len(te), events=ev, note="not tested: fewer than 15 events"))
                continue
            alpha, gamma = learn(tr, Z, cols)
            Zt = Z.loc[te.index, cols].values
            for name, s in (("benefit score (drug-specific part)", -(Zt @ gamma.values)), ("prognostic score (main-effect part)", Zt @ alpha.values)):
                t = te.assign(score=std(pd.Series(s, te.index)))
                cols_c = [x for x in BASE if (x != "stage_n" or (t["stage_n"] != 0).mean() > 0.8)]
                if "age_n" in cols_c and (t["age_n"] != 0).mean() < 0.8:
                    cols_c.remove("age_n")
                r = cox_row(t, cols_c + ["score"], "score")
                row = dict(pathways=set_name, test_cohort=c, group=kind, endpoint=te["endpoint"].iloc[0], score=name, n=r["n"], events=r["events"],
                           hr_per_sd=r["hr"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p"], coef=r["coef"], se=r["se"], adjusted_for=",".join(cols_c) or "nothing")
                rows.append(row); per[name].append(row)
        for name, lst in per.items():
            if lst:
                pl = pooled(lst)
                rows.append(dict(pathways=set_name, test_cohort="ALL TESTED COHORTS (inverse-variance pool)", group="", score=name, n=pl["n"], events=pl["events"],
                                 hr_per_sd=pl["hr"], ci_lo=pl["ci_lo"], ci_hi=pl["ci_hi"], p=pl["p"], coef=pl["coef"], se=pl["se"]))
                print("c", set_name, name, round(pl["hr"], 3), round(pl["p"], 4), flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(TAB / "lung_untreated_prognosis.tsv", sep="\t", index=False)
    return T


# ------------------------------------------------------------------------------------------------ duplicate check
def dup_table():
    T = RESULTS / "training_table_v2"
    P = pd.read_csv(T / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    S = pd.concat([pd.read_csv(T / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")]).loc[P.index]
    old = ["GSE42127", "TCGA-LUAD", "TCGA-LUSC"]
    new = [c for c in P.loc[P["cancer_group"] == "lung", "cohort"].unique() if c not in old]
    rows = []
    for c in new:
        a, b = S[P["cohort"].eq(c).values].dropna(axis=1), S[P["cohort"].isin(old).values]
        b = b[a.columns].dropna(axis=1); a = a[b.columns]
        r = np.corrcoef(((a - a.mean()) / a.std()).values, ((b - b.mean()) / b.std()).values)[:len(a), len(a):]
        rows.append(dict(new_cohort=c, samples=len(a), pathways_compared=a.shape[1], highest_correlation=r.max(),
                         samples_above_0_95=int((r.max(axis=1) > 0.95).sum()), compared_with="GSE42127, TCGA-LUAD, TCGA-LUSC"))
    D = pd.DataFrame(rows)
    D.to_csv(TAB / "lung_duplicate_check.tsv", sep="\t", index=False)
    return D


def main():
    P, Z = TC.load()
    print("table", P.shape, Z.shape, flush=True)
    dup_table()
    A = a_response(P, Z); print(A.round(3).to_string(), flush=True)
    prox = N.precomputed_proximity()["CISPLATIN"].dropna()
    short = [c for c in Z.columns if c in set(prox.index[prox <= -1.2816])]
    sets = {"all pathways (n=%d)" % Z.shape[1]: list(Z.columns), "cisplatin network shortlist (n=%d)" % len(short): short}
    pd.Series(short, name="pathway").to_csv(TAB / "lung_cisplatin_shortlist.tsv", index=False)
    B, scores, ADJ, allfr = b_adjuvant(P, Z, sets)
    print(B.drop(columns=[c for c in B.columns if c.startswith(("all_", "high_", "low_"))]).round(3).to_string(), flush=True)
    C = c_untreated(P, Z, sets, ADJ, allfr)
    print(C.round(3).to_string(), flush=True)
    pd.to_pickle({k: v for k, v in scores.items()}, ROOT / "nar/17_lung/tables/_scores_gse14814.pkl")


if __name__ == "__main__":
    main()
