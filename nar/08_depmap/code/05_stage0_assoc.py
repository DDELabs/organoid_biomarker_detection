"""ANCHOR Stage 0, label-free evidence sources 1 and 2 (no patient data are read here).

Source 1 - expression association: for every drug in every screen (GDSC1, GDSC2, CTRPv2, PRISM), partial Spearman
  correlation of each Reactome pathway rank score with drug response, removing cancer type (lineage dummies) and a
  general drug-sensitivity factor (leave-this-drug-out mean of z-scored responses of the line). Converted to
  z = atanh(r) * sqrt(n - 3 - k). GDSC1 and GDSC2 share lines and expression, so their z are AVERAGED into one
  "Sanger" source; Sanger, CTRP and PRISM are combined with Stouffer's method.
Source 2 - CRISPR co-dependency: the same adjusted drug response against Project Score knockout effects of every
  gene (lineage-adjusted), per-gene z combined across screens as above, then a pathway z = (mean gene z - null mean)
  / null SD over 2,000 random gene sets of the same size.

python -I nar/08_depmap/code/05_stage0_assoc.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from common import DERIVED, SEED, TAB, lineage_dummies  # noqa: E402
from data import SCREENS, crispr, screen  # noqa: E402
from stats import bh  # noqa: E402

from obd.reference import gene_sets  # noqa: E402

# target drug -> {screen: drug name inside that screen}
SOURCES = {
    "FLUOROURACIL": {s: "FLUOROURACIL" for s in SCREENS},
    "OXALIPLATIN": {"GDSC2": "OXALIPLATIN", "CTRP": "OXALIPLATIN", "PRISM": "OXALIPLATIN"},
    "IRINOTECAN": {"GDSC1": "SN-38", "GDSC2": "IRINOTECAN", "CTRP": "SN-38", "PRISM": "IRINOTECAN"},
    "CISPLATIN": {"GDSC1": "CISPLATIN", "GDSC2": "CISPLATIN", "PRISM": "CISPLATIN"},
    "CARBOPLATIN": {"CTRP": "CARBOPLATIN", "PRISM": "CARBOPLATIN"},
    "GEMCITABINE": {s: "GEMCITABINE" for s in SCREENS},
    "PACLITAXEL": {s: "PACLITAXEL" for s in SCREENS},
    "DOXORUBICIN": {"GDSC1": "DOXORUBICIN", "PRISM": "DOXORUBICIN"},
    "PEMETREXED": {"GDSC1": "PEMETREXED"},
    "SORAFENIB": {s: "SORAFENIB" for s in SCREENS},
}
FAMILIES = {
    "antimetabolite": ["FLUOROURACIL", "GEMCITABINE", "PEMETREXED", "METHOTREXATE", "CYTARABINE", "FLOXURIDINE",
                       "DECITABINE", "CLOFARABINE", "FLUDARABINE", "NELARABINE", "RALTITREXED", "AZACITIDINE"],
    "platinum": ["CISPLATIN", "CARBOPLATIN", "OXALIPLATIN"],
    "topoisomerase I": ["IRINOTECAN", "SN-38", "TOPOTECAN", "CAMPTOTHECIN", "BELOTECAN"],
    "topoisomerase II / anthracycline": ["DOXORUBICIN", "EPIRUBICIN", "DAUNORUBICIN", "IDARUBICIN", "MITOXANTRONE",
                                         "ETOPOSIDE", "TENIPOSIDE", "AMRUBICIN", "PIXANTRONE"],
    "microtubule": ["PACLITAXEL", "DOCETAXEL", "CABAZITAXEL", "VINORELBINE", "VINBLASTINE", "VINCRISTINE", "ERIBULIN",
                    "IXABEPILONE"],
    "multi-kinase (VEGFR/RAF)": ["SORAFENIB", "SUNITINIB", "PAZOPANIB", "REGORAFENIB", "AXITINIB", "CABOZANTINIB",
                                 "LENVATINIB", "VANDETANIB"],
}
FAMILY_OF = {"FLUOROURACIL": "antimetabolite", "PEMETREXED": "antimetabolite", "GEMCITABINE": "antimetabolite",
             "OXALIPLATIN": "platinum", "CISPLATIN": "platinum", "CARBOPLATIN": "platinum",
             "IRINOTECAN": "topoisomerase I", "DOXORUBICIN": "topoisomerase II / anthracycline",
             "PACLITAXEL": "microtubule", "SORAFENIB": "multi-kinase (VEGFR/RAF)"}
SAME = {"IRINOTECAN": {"IRINOTECAN", "SN-38"}}
MIN_LINES = 100


def rz(M):
    """Column-wise rank -> normal scores (NaN-safe for complete columns)."""
    r = M.rank(axis=0)
    return stats.norm.ppf((r - 0.5) / r.count())


def resid(Y, C):
    X = np.column_stack([np.ones(len(C)), C])
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X @ b


def general_sensitivity(R):
    Z = (R - R.mean()) / R.std()
    s, c = Z.sum(axis=1, skipna=True), Z.notna().sum(axis=1)
    return Z, s, c


def screen_assoc(name):
    """drugs x pathways matrix of adjusted association z (Source 1) for every drug with >= MIN_LINES lines."""
    S = screen(name)
    lines = S.resp.index.intersection(S.scores.columns)
    R = S.resp.loc[lines]
    R = R.loc[:, R.notna().sum() >= MIN_LINES]
    Z, s, c = general_sensitivity(R)
    P = S.scores[lines].T
    out, ns, gsl = {}, {}, {}
    for d in R.columns:
        L = R[d].dropna().index
        gs_loo = ((s[L] - Z.loc[L, d]) / (c[L] - 1)).values
        lin = lineage_dummies(S.lineage.loc[L].fillna("other").values)
        C = np.column_stack([lin.values, gs_loo])
        y = resid(stats.norm.ppf((R.loc[L, d].rank() - 0.5) / len(L)).reshape(-1, 1), C)[:, 0]
        X = resid(rz(P.loc[L]), C)
        r = (X * y[:, None]).sum(0) / np.sqrt((X ** 2).sum(0) * (y ** 2).sum())
        dof = len(L) - 3 - C.shape[1]
        out[d] = np.arctanh(np.clip(r, -0.999, 0.999)) * np.sqrt(max(dof, 1))
        ns[d] = len(L)
        gsl[d] = (pd.Series(y, index=L), C, L)
    M = pd.DataFrame(out, index=P.columns).T
    return M, pd.Series(ns), gsl


def combine(zs):
    """zs: {screen: Series}. Sanger = mean(GDSC1, GDSC2); Stouffer over Sanger/CTRP/PRISM."""
    parts = []
    sanger = [zs[k] for k in ("GDSC1", "GDSC2") if k in zs]
    if sanger:
        parts.append(sum(sanger) / len(sanger))
    parts += [zs[k] for k in ("CTRP", "PRISM") if k in zs]
    return sum(parts) / np.sqrt(len(parts)), len(parts)


def crispr_codep(gsl_by_screen, drug):
    K = crispr()
    K = K.T.fillna(K.T.mean())          # lines x genes
    zs = {}
    for scr, dname in SOURCES[drug].items():
        S = screen(scr)
        y, C, L = gsl_by_screen[scr][dname]
        sid = S.sidm.reindex(L)
        keep = [i for i, v in zip(L, sid) if isinstance(v, str) and v in K.index]
        if len(keep) < 50:
            continue
        pos = [list(L).index(i) for i in keep]
        Ck, yk = C[pos], y.loc[keep].values
        G = K.loc[[S.sidm[i] for i in keep]]
        G = G.loc[:, G.std() > 0]
        X = resid(rz(G), Ck)
        yk = resid(yk.reshape(-1, 1), Ck)[:, 0]
        r = (X * yk[:, None]).sum(0) / np.sqrt((X ** 2).sum(0) * (yk ** 2).sum())
        dof = len(keep) - 3 - Ck.shape[1]
        zs[scr] = pd.Series(np.arctanh(np.clip(r, -0.999, 0.999)) * np.sqrt(dof), index=G.columns)
    if not zs:
        return None, 0
    common = sorted(set.intersection(*(set(v.index) for v in zs.values())))
    zg, m = combine({k: v[common] for k, v in zs.items()})
    return zg, m


def pathway_z(zg, rng, n_null=2000):
    gs = gene_sets(("REACTOME",))
    genes = np.array(zg.index)
    vals = zg.values
    idx = {g: i for i, g in enumerate(genes)}
    null_cache, rows = {}, {}
    for pw, members in gs.items():
        ii = [idx[g] for g in members if g in idx]
        k = len(ii)
        if k < 5:
            continue
        if k not in null_cache:
            draws = np.array([vals[rng.choice(len(vals), k, replace=False)].mean() for _ in range(n_null)])
            null_cache[k] = (draws.mean(), draws.std())
        mu, sd = null_cache[k]
        rows[pw] = ((vals[ii].mean() - mu) / sd, vals[ii].mean(), k)
    t = pd.DataFrame(rows, index=["crispr_z", "mean_gene_z", "genes_in_screen"]).T
    t["crispr_p"] = 2 * stats.norm.sf(np.abs(t["crispr_z"]))
    t["crispr_q"] = bh(t["crispr_p"].values)
    return t


def main():
    rng = np.random.default_rng(SEED)
    mats, gsl_all = {}, {}
    for scr in SCREENS:
        M, ns, gsl = screen_assoc(scr)
        mats[scr] = M
        gsl_all[scr] = gsl
        M.to_pickle(DERIVED / f"stage0_assoc_{scr}.pkl")
        print(scr, M.shape, flush=True)
    rows, crows = [], []
    for drug, src in SOURCES.items():
        zs = {s: mats[s].loc[d] for s, d in src.items() if d in mats[s].index}
        meta, m = combine(zs)
        p = 2 * stats.norm.sf(np.abs(meta))
        q = bh(p)
        signs = pd.DataFrame({s: np.sign(v) for s, v in zs.items()})
        same = (signs.values == np.sign(meta.values)[:, None]).sum(1)
        # specificity vs all drugs and vs related drugs (directional z, hypothesis direction = sign of meta)
        fam = [x for x in FAMILIES[FAMILY_OF[drug]] if x not in SAME.get(drug, {drug})]
        pct, rel, relused = [], [], set()
        for s, d in src.items():
            if d not in mats[s].index:
                continue
            A = mats[s]
            direc = np.sign(meta.values)
            Ad = A.values * direc[None, :]
            own = A.loc[d].values * direc
            pct.append((Ad < own[None, :]).sum(0) / (len(A) - 1))
            rr = [x for x in fam if x in A.index]
            relused |= set(rr)
            if rr:
                rel.append((A.loc[rr].values * direc[None, :]).mean(0))
        t = pd.DataFrame({"drug": drug, "pathway": meta.index, "expr_meta_z": meta.values, "expr_p": p, "expr_q": q,
                          "expr_sources": m, "screens_used": ",".join(f"{s}:{src[s]}" for s in zs),
                          "screens_same_sign": same, "screens_total": len(zs),
                          "pct_vs_all_drugs": np.mean(pct, axis=0),
                          "related_mean_dir_z": np.mean(rel, axis=0) if rel else np.nan,
                          "related_drugs": ",".join(sorted(relused))})
        for s, v in zs.items():
            t[f"z_{s}"] = v.values
        t["own_dir_z_mean"] = np.mean([zs[s].values * np.sign(meta.values) for s in zs], axis=0)
        rows.append(t)
        zg, mc = crispr_codep(gsl_all, drug)
        if zg is not None:
            ct = pathway_z(zg, rng)
            ct["drug"] = drug
            ct["crispr_sources"] = mc
            crows.append(ct.reset_index().rename(columns={"index": "pathway"}))
            zg.rename("gene_z").to_frame().assign(drug=drug).to_csv(DERIVED / f"crispr_gene_z_{drug}.tsv.gz", sep="\t")
        print(drug, "done", flush=True)
    E = pd.concat(rows, ignore_index=True)
    E.to_csv(TAB / "B_source1_expression_association.tsv.gz", sep="\t", index=False)
    Cc = pd.concat(crows, ignore_index=True)
    Cc.to_csv(TAB / "B_source2_crispr_codependency.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
