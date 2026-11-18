"""G9: dependence-aware null for the 57-analysis survey.

Within each TCGA survival pool (COAD+READ+COADREAD share patients, so they are one pool; the paper's BLCA biotab
survival and the 2019 COAD files are their own pools), the (time, event) pairs are shuffled across patients, separately
inside the 'treated with any drug in this survey' and 'no drug record' strata. The same shuffle is used for every
analysis in the pool, so analyses sharing patients stay correlated exactly as in the real data. The frozen scores and
median splits are not touched. All 57 treated and 57 untreated log-rank tests are recomputed per permutation.
"""
import numpy as np
import pandas as pd

from g_common import DATA, SEED, TAB, load_patients, logrank_fast, master

NPERM = 5000


def clusters(p):
    """Connected components of treated patient sets with overlap >= 50% of the smaller set."""
    t = p[p.group == "treated"]
    sets = {a: set(d.patient) for a, d in t.groupby("analysis_id")}
    ids = sorted(sets)
    par = {a: a for a in ids}

    def f(a):
        while par[a] != a:
            a = par[a]
        return a
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if len(sets[a] & sets[b]) >= 0.5 * min(len(sets[a]), len(sets[b])):
                par[f(a)] = f(b)
    roots = {a: f(a) for a in ids}
    lab = {r: k for k, r in enumerate(sorted(set(roots.values())))}
    return {a: lab[roots[a]] for a in ids}


def main():
    rng = np.random.default_rng(SEED)
    p = load_patients()
    M = master()
    aids = list(M.analysis_id)
    # canonical survival per pool x patient; stratum = ever treated in the pool
    canon = p.drop_duplicates(["pool", "patient"]).set_index(["pool", "patient"])[["T", "E"]]
    ever_t = p[p.group == "treated"].groupby("pool").patient.apply(set).to_dict()
    pools = {}
    for pool, d in p.groupby("pool"):
        pats = np.array(sorted(d.patient.unique()))
        strat = np.array([x in ever_t.get(pool, set()) for x in pats])
        pools[pool] = (pats, strat)
    # per analysis arrays (index into pool patients)
    an = {}
    for (aid, grp), d in p.groupby(["analysis_id", "group"]):
        pool = d.pool.iloc[0]
        pats = pools[pool][0]
        pos = np.searchsorted(pats, d.patient.values)
        an[(aid, grp)] = (pool, pos, d.split.values, d["T"].values, d.E.values)
    cl = clusters(p)
    obs = {}
    for k, (pool, pos, g, T, E) in an.items():
        obs[k] = logrank_fast(T, E, g)
    # also observed direction (events in predicted non-responders > expected): use observed hazard sign via O-E
    res_t = np.full((NPERM, len(aids)), np.nan)
    res_u = np.full((NPERM, len(aids)), np.nan)
    dir_t = np.zeros((NPERM, len(aids)), bool)
    for b in range(NPERM):
        perm = {}
        for pool, (pats, strat) in pools.items():
            src = np.arange(len(pats))
            for s in (True, False):
                ix = np.where(strat == s)[0]
                src[ix] = rng.permutation(ix)
            perm[pool] = src
        for j, aid in enumerate(aids):
            for grp, R in (("treated", res_t), ("untreated", res_u)):
                if (aid, grp) not in an:
                    continue
                pool, pos, g, T, E = an[(aid, grp)]
                ct = CAN[pool]
                donor = perm[pool][pos]
                Tp, Ep = ct[0][donor], ct[1][donor]
                R[b, j] = logrank_fast(Tp, Ep, g)
                if grp == "treated":
                    # expected direction = more events than expected in split==1 (predicted non-responders)
                    dir_t[b, j] = Ep[g == 1].sum() / max(g.sum(), 1) > Ep[g == 0].sum() / max((1 - g).sum(), 1)
        if b % 500 == 0:
            print(b, flush=True)
    np.save(DATA / "g9_perm_treated_p.npy", res_t)
    np.save(DATA / "g9_perm_untreated_p.npy", res_u)
    np.save(DATA / "g9_perm_treated_dir.npy", dir_t)
    ot = np.array([obs[(a, "treated")] for a in aids])
    ou = np.array([obs.get((a, "untreated"), np.nan) for a in aids])
    cid = np.array([cl[a] for a in aids])
    expdir = M.direction.str.startswith("expected").values

    def summ(name, obsv, null):
        return {"statistic": name, "observed": obsv, "null_mean": float(np.mean(null)), "null_median": float(np.median(null)),
                "null_2.5pct": float(np.percentile(null, 2.5)), "null_97.5pct": float(np.percentile(null, 97.5)),
                "null_95pct": float(np.percentile(null, 95)),
                "perm_p_at_least_observed": float((np.sum(null >= obsv) + 1) / (len(null) + 1))}

    def summ_low(name, obsv, null):
        r = summ(name, obsv, null)
        r["perm_p_at_least_observed"] = float((np.sum(null <= obsv) + 1) / (len(null) + 1))
        r["note"] = "p = share of permutations with a value as small or smaller"
        return r
    cnt_t = np.nansum(res_t < 0.05, 1)
    cnt_u = np.nansum(res_u < 0.05, 1)
    cnt_te = np.sum((res_t < 0.05) & dir_t, 1)
    ncl = len(set(cid))
    clhit = np.array([len(set(cid[res_t[b] < 0.05])) for b in range(NPERM)])
    rows = [summ("treated: count p<0.05 (of 57)", int(np.sum(ot < 0.05)), cnt_t),
            summ("treated: count p<0.05 in expected direction", int(np.sum((ot < 0.05) & expdir)), cnt_te),
            summ("treated: count p<0.01", int(np.sum(ot < 0.01)), np.sum(res_t < 0.01, 1)),
            summ_low("treated: minimum p", float(np.min(ot)), np.nanmin(res_t, 1)),
            summ("untreated: count p<0.05 (of %d)" % np.isfinite(ou).sum(), int(np.nansum(ou < 0.05)), cnt_u),
            summ_low("untreated: minimum p", float(np.nanmin(ou)), np.nanmin(res_u, 1)),
            summ(f"treated: independent patient clusters with >=1 p<0.05 (of {ncl})", len(set(cid[ot < 0.05])), clhit)]
    S = pd.DataFrame(rows)
    S["n_permutations"] = NPERM
    S.to_csv(TAB / "G9_permutation_null_summary.tsv", sep="\t", index=False)
    # per analysis permutation p (single test, for reference) and cluster ids
    per = pd.DataFrame({"analysis_id": aids, "cluster": cid, "treated_p_obs": ot, "untreated_p_obs": ou,
                        "treated_perm_p_single": [(np.sum(res_t[:, j] <= ot[j]) + 1) / (NPERM + 1) for j in range(len(aids))],
                        "familywise_minP_adjusted_p": [(np.sum(np.nanmin(res_t, 1) <= ot[j]) + 1) / (NPERM + 1) for j in range(len(aids))],
                        "null_share_p_below_0.05": np.mean(res_t < 0.05, 0)})
    per.to_csv(TAB / "G9_per_analysis_permutation.tsv", sep="\t", index=False)
    # check: binomial independent reference
    from scipy.stats import binom
    ref = pd.DataFrame({"k": np.arange(0, 21), "perm_null_prob": [np.mean(cnt_t == k) for k in range(21)],
                        "binomial_57_0.05_prob": binom.pmf(np.arange(21), 57, 0.05),
                        "perm_null_untreated_prob": [np.mean(cnt_u == k) for k in range(21)]})
    ref.to_csv(TAB / "G9_null_count_distribution.tsv", sep="\t", index=False)
    print(S.to_string())
    print(per.sort_values("treated_p_obs").head(8).to_string())


if __name__ == "__main__":
    import g_common
    p0 = load_patients()
    canon0 = p0.drop_duplicates(["pool", "patient"])
    CAN = {}
    for pool, d in canon0.groupby("pool"):
        d = d.sort_values("patient")
        CAN[pool] = (d["T"].values.astype(float), d.E.values.astype(int))
    main()
