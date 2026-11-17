"""Stage 1: pooled TCGA interaction test of each drug's frozen hypothesis score (see ANALYSIS_PLAN.md §2)."""
import json
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

import common as C
import s00_fetch_tcga_drug as F

warnings.filterwarnings("ignore")
R = F.R
from run_atlas_relaxed import SUPPORTIVE, THERAPY  # noqa: E402

N_NULL, N_BOOT = 1000, 200


# ------------------------------------------------------------------ data
def build():
    H = C.frozen()
    P, S = C.load_table()
    P = P[P["source"].eq("TCGA")]
    S = S.loc[P.index]
    cdr = R.load_cdr()
    P = P[P["patient"].isin(cdr.index)]
    S = S.loc[P.index]
    proj = P["cohort"].str.replace("TCGA-", "")
    Zm = C.zwithin(S, proj)
    d = pd.DataFrame(index=P.index)
    d["patient"], d["project"] = P["patient"], proj
    c = cdr.loc[d["patient"]]
    for k in ("PFI", "OS"):
        d[f"{k}_time"], d[f"{k}_event"] = c[f"{k}.time"].values, c[k].values
    d["age"] = c["age"].values
    d["stage"] = [C.stage_num(s) for s in c["stage"].values]
    d["male"] = c["gender"].str.upper().eq("MALE").astype(float).values
    pur = pd.read_csv(C.PURITY, sep="\t")
    pur = pur.assign(patient=pur["array"].str[:12]).drop_duplicates("patient").set_index("patient")["purity"]
    d["purity"] = pur.reindex(d["patient"]).values
    for drug, w in H.items():
        d[f"score_{drug}"] = C.signed_score(Zm, w, proj).values
    d["score_PROLIF"] = C.signed_score(Zm, {p: 1.0 for p in C.PROLIF}, proj).values
    g = d.groupby("project")["purity"]
    d["score_PURITY"] = ((d["purity"] - g.transform("mean")) / g.transform("std")).values
    # drug records
    long = R.load_response_labels()
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    long = long[long["patient"].isin(set(d["patient"]))].copy()
    long["start"] = long["start_days"].fillna(0).clip(lower=0)
    first = long.groupby(["patient", "drug"])["start"].min().unstack()
    for drug in DRUG_UNITS:
        names = DRUG_UNITS[drug]
        cols = [n for n in names if n in first.columns]
        d[f"start_{drug}"] = first[cols].min(axis=1).reindex(d["patient"]).values
    anyst = long.groupby("patient")["start"].min()
    d["any_systemic"] = d["patient"].isin(anyst.index).astype(float)
    d = d[(d["PFI_time"] > 0) & d["PFI_event"].notna() & d[[f"score_{x}" for x in H]].notna().all(1)]
    return d, Zm.loc[d.index], long


DRUG_UNITS = {**{x: [x] for x in C.DRUGS}, "FLUOROURACIL_CAPE": ["FLUOROURACIL", "CAPECITABINE"]}


def included_projects(d, unit, min_n=10):
    n = d[d[f"start_{unit}"].notna()].groupby("project").size()
    return sorted(n[n >= min_n].index)


# ------------------------------------------------------------------ design
def covariates(d, purity=True):
    X = pd.DataFrame(index=d.index)
    g = d.groupby("project")
    X["age10"] = ((d["age"] - g["age"].transform("mean")) / 10).fillna(0)
    X["stage"] = (d["stage"] - g["stage"].transform("mean")).fillna(0)
    X["stage_missing"] = d["stage"].isna().astype(float)
    X["male"] = d["male"] - g["male"].transform("mean")
    if purity:
        X["purity"] = (d["purity"] - g["purity"].transform("mean")).fillna(0)
        X["purity_missing"] = d["purity"].isna().astype(float)
    return X


def split_rows(d, unit, ep):
    """Start-stop rows with time-varying exposure."""
    t, e, s = d[f"{ep}_time"].values, d[f"{ep}_event"].values, d[f"start_{unit}"].values
    idx, en, ex, ev, ex_flag = [], [], [], [], []
    for i in range(len(d)):
        si = s[i]
        if np.isnan(si) or si >= t[i]:
            idx.append(i); en.append(0.0); ex.append(t[i]); ev.append(e[i]); ex_flag.append(0.0)
        elif si <= 0:
            idx.append(i); en.append(0.0); ex.append(t[i]); ev.append(e[i]); ex_flag.append(1.0)
        else:
            idx += [i, i]; en += [0.0, si]; ex += [si, t[i]]; ev += [0.0, e[i]]; ex_flag += [0.0, 1.0]
    return np.array(idx), np.array(en), np.array(ex), np.array(ev), np.array(ex_flag)


class Design:
    """Holds the split rows for one drug unit x endpoint so many scores can be tested quickly."""

    def __init__(self, d, unit, ep, purity=True, extra=None, weights=None):
        self.d = d
        self.idx, self.entry, self.exit, self.event, self.E = split_rows(d, unit, ep)
        Xc = covariates(d, purity)
        if extra is not None:
            for k, v in extra.items():
                Xc[k] = v
        self.Xc = Xc.values[self.idx]
        self.cov_names = list(Xc.columns)
        self.strata = d["project"].values[self.idx]
        self.w = None if weights is None else np.asarray(weights)[self.idx]
        self.D1 = float(self.event[self.E == 1].sum())
        self.D0 = float(self.event[self.E == 0].sum())
        self.n = len(d)
        self.n_exposed = int(len(set(self.idx[self.E == 1])))

    def fit(self, score, extra_int=None):
        s = np.asarray(score, float)[self.idx]
        cols = [s, self.E, s * self.E]
        if extra_int is not None:
            xi = np.asarray(extra_int, float)[self.idx]
            cols += [xi * self.E]
        X = np.column_stack(cols + [self.Xc])
        keep = X.std(0) > 0
        b, se, cov = C.cox_fit(self.exit, self.event, X[:, keep], entry=self.entry, strata=self.strata, weights=self.w)
        bb, ss = np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan)
        bb[keep], ss[keep] = b, se
        ci = np.where(keep)[0]
        c_main_int = cov[np.searchsorted(ci, 0), np.searchsorted(ci, 2)]
        return {"b_score": bb[0], "se_score": ss[0], "b_exp": bb[1], "se_exp": ss[1], "b_int": bb[2], "se_int": ss[2],
                "se_score_exposed": float(np.sqrt(ss[0] ** 2 + ss[2] ** 2 + 2 * c_main_int))}


def summarise(r, prefix=""):
    z, p = C.wald(r["b_int"], r["se_int"])
    return {f"{prefix}HR_int": np.exp(r["b_int"]), f"{prefix}lo": np.exp(r["b_int"] - 1.96 * r["se_int"]),
            f"{prefix}hi": np.exp(r["b_int"] + 1.96 * r["se_int"]), f"{prefix}z": z, f"{prefix}p": p}


# ------------------------------------------------------------------ null worker
_G = {}


def _null_init(d, unit, Zv, groups):
    _G["des"] = Design(d, unit, "PFI")
    _G["Z"], _G["groups"] = Zv, groups


def _null_chunk(Wchunk):
    Zv, groups = _G["Z"], _G["groups"]
    S = Zv @ Wchunk / 5.0
    S = pd.DataFrame(S)
    g = S.groupby(groups)
    S = ((S - g.transform("mean")) / g.transform("std")).values
    out = []
    for j in range(S.shape[1]):
        r = _G["des"].fit(S[:, j])
        out.append(r["b_int"] / r["se_int"])
    return out


def random_weights(n_path, n=N_NULL, k=5, seed=C.SEED):
    rng = np.random.default_rng(seed)
    W = np.zeros((n_path, n))
    for j in range(n):
        W[rng.choice(n_path, k, replace=False), j] = rng.choice([-1.0, 1.0], k)
    return W


# ------------------------------------------------------------------ IPW
def ipw_weights(d, unit):
    y = d[f"start_{unit}"].notna() & (d[f"start_{unit}"] < d["PFI_time"])
    X = covariates(d, True)
    X = pd.concat([X, pd.get_dummies(d["project"], drop_first=True).astype(float)], axis=1).values
    m = LogisticRegression(C=1e6, max_iter=2000).fit(X, y)
    ps = np.clip(m.predict_proba(X)[:, 1], 1e-3, 1 - 1e-3)
    pt = y.mean()
    w = np.where(y, pt / ps, (1 - pt) / (1 - ps))
    lo, hi = np.percentile(w, [1, 99])
    return np.clip(w, lo, hi)


def main():
    d_all, Zm, long = build()
    d_all.to_csv(C.DATA / "tcga_stage1_patients.tsv.gz", sep="\t")
    rng = np.random.default_rng(C.SEED)
    prim, sens, nulls, pc, spec, powr, ipwrows = [], [], {}, [], [], [], []
    W = random_weights(Zm.shape[1])
    pd.DataFrame(W, index=Zm.columns).to_csv(C.DATA / "random_signed_sets.tsv.gz", sep="\t")
    incl = {u: included_projects(d_all, u) for u in DRUG_UNITS}
    for drug in C.DRUGS:
        projs = incl[drug]
        d = d_all[d_all["project"].isin(projs)].copy()
        sc = d[f"score_{drug}"].values
        row = {"drug": drug, "projects": ";".join(projs), "n_projects": len(projs)}
        for ep in ("PFI", "OS"):
            dd = d[d[f"{ep}_time"] > 0]
            des = Design(dd, drug, ep)
            r = des.fit(dd[f"score_{drug}"].values)
            s = summarise(r)
            row.update({f"{ep}_{k}": v for k, v in s.items()})
            row.update({f"{ep}_n": des.n, f"{ep}_n_exposed": des.n_exposed, f"{ep}_events_exposed": des.D1,
                        f"{ep}_events_unexposed": des.D0,
                        f"{ep}_HR_score_unexposed": np.exp(r["b_score"]),
                        f"{ep}_HR_score_unexposed_lo": np.exp(r["b_score"] - 1.96 * r["se_score"]),
                        f"{ep}_HR_score_unexposed_hi": np.exp(r["b_score"] + 1.96 * r["se_score"]),
                        f"{ep}_HR_score_exposed": np.exp(r["b_score"] + r["b_int"]),
                        f"{ep}_HR_score_exposed_lo": np.exp(r["b_score"] + r["b_int"] - 1.96 * r["se_score_exposed"]),
                        f"{ep}_HR_score_exposed_hi": np.exp(r["b_score"] + r["b_int"] + 1.96 * r["se_score_exposed"]),
                        f"{ep}_se_int": r["se_int"]})
            if ep == "PFI":
                des_pfi = des
                se_exp = np.sqrt(1 / des.D1 + 1 / des.D0)
                powr.append({"drug": drug, "endpoint": ep, "test": "Stage 1 primary interaction", "n": des.n,
                             "events_exposed": des.D1, "events_unexposed": des.D0, "expected_se": se_exp,
                             "observed_se": r["se_int"],
                             "MDHR_alpha05_expected": C.mdhr(se_exp), "MDHR_holm_first_step_expected": C.mdhr(se_exp, 0.05 / 7),
                             "MDHR_alpha05_observed": C.mdhr(r["se_int"]),
                             "power_HR1.25_expected": C.power_two_sided(np.log(1.25), se_exp),
                             "power_HR1.43_expected": C.power_two_sided(np.log(1 / 0.7), se_exp),
                             "power_HR1.25_observed": C.power_two_sided(np.log(1.25), r["se_int"]),
                             "power_HR1.43_observed": C.power_two_sided(np.log(1 / 0.7), r["se_int"]),
                             "power_HR1.43_holm_first_step_observed": C.power_two_sided(np.log(1 / 0.7), r["se_int"], 0.05 / 7)})
        prim.append(row)
        print(drug, {k: round(v, 4) if isinstance(v, float) else v for k, v in row.items() if k.startswith("PFI_")}, flush=True)

        # sensitivity analyses (PFI)
        def add(name, r, des_):
            s = summarise(r)
            sens.append({"drug": drug, "analysis": name, "n": des_.n, "events_exposed": des_.D1,
                         "events_unexposed": des_.D0, **s})
        des = Design(d, drug, "PFI", purity=False)
        add("S1 no purity adjustment", des.fit(sc), des)
        # S2 IPW with bootstrap
        w = ipw_weights(d, drug)
        des = Design(d, drug, "PFI", weights=w)
        r_ipw = des.fit(sc)
        boots = []
        pats = np.arange(len(d))
        for bt in range(N_BOOT):
            ii = rng.choice(pats, len(pats), replace=True)
            db = d.iloc[ii].copy()
            db.index = [f"{x}#{k}" for k, x in enumerate(db.index)]
            try:
                wb = ipw_weights(db, drug)
                rb = Design(db, drug, "PFI", weights=wb).fit(db[f"score_{drug}"].values)
                boots.append(rb["b_int"])
            except Exception:
                pass
        boots = np.array(boots)
        se_b = boots.std(ddof=1)
        z = r_ipw["b_int"] / se_b
        sens.append({"drug": drug, "analysis": "S2 IPW (bootstrap SE)", "n": des.n, "events_exposed": des.D1,
                     "events_unexposed": des.D0, "HR_int": np.exp(r_ipw["b_int"]),
                     "lo": np.exp(np.percentile(boots, 2.5)), "hi": np.exp(np.percentile(boots, 97.5)), "z": z,
                     "p": 2 * C.stats.norm.sf(abs(z)), "n_boot": len(boots)})
        ipwrows.append({"drug": drug, "w_min": w.min(), "w_max": w.max(), "w_mean": w.mean(),
                        "ess": w.sum() ** 2 / (w ** 2).sum(), "n": len(w)})
        # S3 active comparator
        da = d[(d["any_systemic"] == 1)]
        des = Design(da, drug, "PFI")
        add("S3 active comparator (other systemic therapy)", des.fit(da[f"score_{drug}"].values), des)
        # S4 proliferation-adjusted
        des = Design(d, drug, "PFI", extra={"prolif": d["score_PROLIF"].values})
        add("S4 + proliferation score and its interaction", des.fit(sc, extra_int=d["score_PROLIF"].values), des)
        if drug == "FLUOROURACIL":
            u = "FLUOROURACIL_CAPE"
            dc = d_all[d_all["project"].isin(incl[u])]
            des = Design(dc, u, "PFI")
            add("S5 5-FU or capecitabine", des.fit(dc[f"score_{drug}"].values), des)

        # negative controls
        r = des_pfi.fit(d["score_PROLIF"].values)
        pc.append({"drug": drug, "control": "proliferation", "n": des_pfi.n, **summarise(r)})
        dp = d[d["purity"].notna()]
        des = Design(dp, drug, "PFI", purity=False)
        r = des.fit(dp["score_PURITY"].values)
        pc.append({"drug": drug, "control": "tumour purity", "n": des.n, **summarise(r)})

        # specificity matrix: other drugs' scores on this drug's exposure
        for other in C.DRUGS:
            r = des_pfi.fit(d[f"score_{other}"].values)
            spec.append({"exposure": drug, "score": other, **summarise(r)})

        # per-cancer
        for pj in projs:
            dj = d[d["project"] == pj]
            dj_des = Design(dj, drug, "PFI")
            if dj_des.D1 < 5 or dj_des.D0 < 5:
                prim_pc = {"drug": drug, "project": pj, "n": dj_des.n, "n_exposed": dj_des.n_exposed,
                           "events_exposed": dj_des.D1, "events_unexposed": dj_des.D0, "fitted": False}
            else:
                r = dj_des.fit(dj[f"score_{drug}"].values)
                prim_pc = {"drug": drug, "project": pj, "n": dj_des.n, "n_exposed": dj_des.n_exposed,
                           "events_exposed": dj_des.D1, "events_unexposed": dj_des.D0, "fitted": True,
                           "b_int": r["b_int"], "se_int": r["se_int"], **summarise(r)}
            nulls.setdefault("_percancer", []).append(prim_pc)

        # random-signature null (2 processes)
        Zv = Zm.loc[d.index].values
        chunks = np.array_split(W, 20, axis=1)
        with Pool(2, initializer=_null_init, initargs=(d, drug, Zv, d["project"].values)) as pool:
            zs = np.concatenate(pool.map(_null_chunk, chunks))
        nulls[drug] = zs
        print(drug, "null done", flush=True)

    T = pd.DataFrame(prim)
    T["PFI_holm_p"] = C.holm(T["PFI_p"].values)
    T["OS_holm_p_secondary"] = C.holm(T["OS_p"].values)
    T["PFI_p_one_sided"] = C.stats.norm.sf(T["PFI_z"])
    for drug in C.DRUGS:
        zs = nulls[drug]
        zs = zs[np.isfinite(zs)]
        T.loc[T["drug"] == drug, "random_null_n_valid"] = len(zs)
        zo = T.loc[T["drug"] == drug, "PFI_z"].iloc[0]
        T.loc[T["drug"] == drug, "random_null_p_one_sided"] = (1 + (zs >= zo).sum()) / (1 + len(zs))
        T.loc[T["drug"] == drug, "random_null_pct_rank"] = (zs < zo).mean()
    T.to_csv(C.TAB / "stage1_primary.tsv", sep="\t", index=False)
    Sx = pd.DataFrame(sens)
    Sx["holm_within_analysis"] = Sx.groupby("analysis")["p"].transform(lambda p: C.holm(p.values))
    Sx.to_csv(C.TAB / "stage1_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(ipwrows).to_csv(C.TAB / "stage1_ipw_weights.tsv", sep="\t", index=False)
    NC = pd.DataFrame(pc)
    NC["holm_within_control"] = NC.groupby("control")["p"].transform(lambda p: C.holm(p.values))
    NC.to_csv(C.TAB / "stage1_negative_controls.tsv", sep="\t", index=False)
    pd.DataFrame(spec).to_csv(C.TAB / "stage1_specificity_matrix.tsv", sep="\t", index=False)
    pd.DataFrame(powr).to_csv(C.TAB / "stage1_power.tsv", sep="\t", index=False)
    PC = pd.DataFrame(nulls.pop("_percancer"))
    PC.to_csv(C.TAB / "stage1_per_cancer.tsv", sep="\t", index=False)
    pd.DataFrame(nulls).to_csv(C.TAB / "stage1_random_null_z.tsv.gz", sep="\t", index=False)
    meta = []
    for drug, g in PC[PC["fitted"]].groupby("drug"):
        m = C.dersimonian_laird(g["b_int"], g["se_int"])
        meta.append({"drug": drug, "k": m["k"], "HR_re": np.exp(m["re"]), "lo_re": np.exp(m["re"] - 1.96 * m["se_re"]),
                     "hi_re": np.exp(m["re"] + 1.96 * m["se_re"]), "p_re": m["p_re"], "HR_fe": np.exp(m["fe"]),
                     "p_fe": m["p_fe"], "tau2": m["tau2"], "I2": m["I2"], "Q": m["Q"]})
    pd.DataFrame(meta).to_csv(C.TAB / "stage1_meta_analysis.tsv", sep="\t", index=False)
    # gate G2
    G = []
    for _, r in T.iterrows():
        drug = r["drug"]
        s = Sx[Sx["drug"] == drug].set_index("analysis")
        nc = NC[NC["drug"] == drug].set_index("control")
        c1 = (r["PFI_holm_p"] < 0.05) and (r["PFI_HR_int"] > 1)
        c2 = all((s.loc[a, "HR_int"] > 1) and (s.loc[a, "p"] < 0.05) for a in ("S1 no purity adjustment", "S2 IPW (bootstrap SE)"))
        c3 = bool((nc["p"] >= 0.05).all())
        c4 = abs(np.log(r["PFI_HR_score_unexposed"])) < abs(np.log(r["PFI_HR_score_exposed"]))
        G.append({"drug": drug, "G2_1_holm_and_direction": c1, "G2_2_purity_ipw_robust": c2, "G2_3_negative_controls_null": c3,
                  "G2_4_weaker_in_comparator": c4, "G2_pass": c1 and c2 and c3 and c4})
    pd.DataFrame(G).to_csv(C.TAB / "stage1_gate_G2.tsv", sep="\t", index=False)
    print(T[["drug", "PFI_HR_int", "PFI_lo", "PFI_hi", "PFI_p", "PFI_holm_p", "random_null_p_one_sided"]].to_string())


if __name__ == "__main__":
    main()
