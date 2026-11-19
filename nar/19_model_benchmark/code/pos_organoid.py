"""W3 positives, organoid-trained tasks: bladder cisplatin (log-rank) + ERCC2, colon cetuximab vs BRAF V600E,
rectal GSE87211 5-FU signature, ovarian olaparib signature in GSE63885 (+ TCGA-OV, GSE9891 as context).

Usage: python pos_organoid.py TASK [TASK ...]      (one process per call; at most 2 at a time)
Every model is trained on the lab models only (organoids / cell lines), hyper-parameters chosen by leave-one-out
(<= 20 lab models) or 5-fold CV INSIDE the lab models. The patients are never used for training or tuning.
Outputs: tables/preds_pos/<task>.pkl
"""
import os, sys, pickle, time, json
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
import pzoo as Z
import pos_stats as S

HERE = Z.HERE; ROOT = Z.ROOT; NAR = ROOT / "nar"
OUT = HERE / "tables" / "preds_pos"; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(NAR / "07_B3_models/code"))
import b3common as B3                                   # bundles, mutation-free helpers
CUT = -1.2816
AA = "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION"
KW = ["CELL_CYCLE", "MITOTIC", "G1_S", "G2_M", "E2F", "DNA_REPLICATION", "PRE_REPLICATIVE", "CHECKPOINT"]
N_RANDOM = 200
SILENT = {"synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "upstream_gene_variant",
          "downstream_gene_variant", "intergenic_region", "non_coding_transcript_exon_variant", "splice_region_variant"}


def zy(y):
    y = np.asarray(y, float); return (y - y.mean()) / (y.std() or 1)


def ridge_paper(X, y):
    return RidgeCV(cv=min(3, len(y)), alphas=np.arange(0.1, 1, 0.1)).fit(X, y).coef_


# ------------------------------------------------------------------------------------------------ task definitions
def mut_sets(files, gene, rule):
    m = pd.concat([pd.read_csv(f, sep="\t", usecols=["sample", "gene", "Amino_Acid_Change", "effect"]) if "Amino_Acid_Change" in pd.read_csv(f, sep="\t", nrows=1).columns
                   else pd.read_csv(f, sep="\t", usecols=["sample", "gene", "effect"]).assign(Amino_Acid_Change="") for f in files])
    m = m[m["sample"].str[13:15] == "01"].copy(); m["patient"] = m["sample"].str[:12]
    s = m[(m.gene == gene) & ~m.effect.isin(SILENT)]
    return set(s.loc[s.Amino_Acid_Change.fillna("").map(rule), "patient"]), set(m.patient)


TASKS = {
    "bladder_logrank_paperdata":   dict(kind="surv", bundle="BLCA_paper", drug="CISPLATIN", orig="Paper_AA_pathway", label="Bladder cisplatin (paper data): organoids -> TCGA-BLCA treated, survival"),
    "bladder_logrank_currentdata": dict(kind="surv", bundle="BLCA_current", drug="CISPLATIN", orig="Paper_AA_pathway", label="Bladder cisplatin (current data): organoids -> TCGA-BLCA treated, survival"),
    "bladder_ERCC2_paperdata":     dict(kind="bin", bundle="BLCA_paper", drug="CISPLATIN", orig="Paper_AA_pathway", gene="ERCC2", sign=-1, label="Bladder ERCC2 mutants vs wild type (paper data)"),
    "bladder_ERCC2_currentdata":   dict(kind="bin", bundle="BLCA_current", drug="CISPLATIN", orig="Paper_AA_pathway", gene="ERCC2", sign=-1, label="Bladder ERCC2 mutants vs wild type (current data)"),
    "colon_cetuximab_BRAF":        dict(kind="bin", bundle="COADREAD_vdW", drug="CETUXIMAB", orig="Paper_rank1", gene="BRAF", sign=+1, label="Colon cetuximab: BRAF V600E mutants vs wild type"),
    "rectal_GSE87211_5FU":         dict(kind="surv", bundle="COADREAD_vdW", drug="FLUOROURACIL", orig="Paper_top10_frozen", label="Rectal GSE87211 5-FU chemoradiotherapy, DFS (organoid signature)"),
    "ovarian_olaparib_GSE63885":   dict(kind="surv", orig="Corr_signature", cohort="GSE63885", label="Ovarian olaparib signature, GSE63885 platinum-treated, OS"),
    "ovarian_olaparib_TCGAOV":     dict(kind="surv", orig="Corr_signature", cohort="TCGA-OV", label="Ovarian olaparib signature, TCGA-OV platinum-treated, OS (context)"),
    "ovarian_olaparib_GSE9891":    dict(kind="surv", orig="Corr_signature", cohort="GSE9891", label="Ovarian olaparib signature, GSE9891 platinum-treated, OS (context)"),
}


def load_bundle_task(t):
    b = B3.bundle(t["bundle"]); d = b["drugs"][t["drug"]]
    y = d["y"]; Sm = [s for s in y.index if s in b["osc"].columns]
    U = sorted(set(b["osc"].index) & set(b["tsc"].index))
    Xo = b["osc"].loc[U, Sm].T.values.astype(float); yv = y.loc[Sm].values.astype(float)
    z = d["z"].reindex(U)
    out = dict(feats=U, Xo=Xo, y=yv, z=z, org=Sm)
    if t["kind"] == "bin":
        sets = {"ERCC2": (NAR / "02_organoid_only/BLCA_OV_LUAD_OTHER/_data/raw/TCGA-BLCA.somaticmutation_wxs.tsv.gz",),
                "BRAF": (NAR / "10_completeness/_data/raw/xena_mut/TCGA-COAD.somaticmutation_wxs.tsv.gz", NAR / "10_completeness/_data/raw/xena_mut/TCGA-READ.somaticmutation_wxs.tsv.gz")}[t["gene"]]
        rule = (lambda a: True) if t["gene"] == "ERCC2" else (lambda a: a in ("p.V600E", "p.V640E"))
        mut, prof = mut_sets(list(sets), t["gene"], rule)
        pats = [p for p in b["tsc"].columns if p in prof]
        Xp = B3.zcols(b["tsc"].loc[U, :].T.values)[0]                                # z over all tumours (as B3)
        Xp = pd.DataFrame(Xp, index=b["tsc"].columns)
        Xp = Xp.loc[pats].values
        lab = np.array([1 if p in mut else 0 for p in pats])
        out.update(Xp=Xp, ev=pd.DataFrame({"label": lab}, index=pats), n_mut=int(lab.sum()))
    else:
        pats = [p for p in b["tsc"].columns if p in b["gp"].columns and p in b["surv"].index]
        Xp_all = B3.zcols(b["tsc"].loc[U, pats].T.values)[0]
        tr = [i for i, p in enumerate(pats) if p in set(d["treated"])]
        if t["bundle"] == "COADREAD_vdW":                                             # rectal GSE87211: patients from the training table
            P, Sx = rect_patients()
            cols = [c for c in U if c in Sx.columns]
            assert len(cols) == len(U), (len(cols), len(U))
            Zc = ((Sx[U] - Sx[U].mean()) / Sx[U].std().replace(0, 1)).fillna(0)
            ev = P[["T", "E"]]
            Xp = Zc.loc[ev.index].values
            out.update(Xp=Xp, ev=ev, Zc=Zc)
        else:
            sv = b["surv"].loc[pats]
            ev = pd.DataFrame({"T": np.clip(sv["T"].values, 0, None), "E": sv["E"].values.astype(int)}, index=pats).iloc[tr]
            out.update(Xp=Xp_all[tr], ev=ev)
    return out


def rect_patients():
    TT = ROOT / "robust" / "results" / "training_table"
    P = pd.read_csv(TT / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    Sx = pd.concat([pd.read_csv(TT / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")]).loc[P.index]
    idx = P.index[P.cohort == "GSE87211"]
    d = P.loc[idx].copy()
    key = d["arm"].fillna("NA").astype(str) + "|" + d["drugs"].fillna("NA").astype(str)
    d = d[key.str.contains("FLUOROURACIL", regex=True)]
    d["T"] = pd.to_numeric(d["dfs_months"], errors="coerce"); d["E"] = pd.to_numeric(d["dfs_event"], errors="coerce")
    d = d.dropna(subset=["T", "E"]); d = d[d["T"] >= 0]; d["E"] = d["E"].astype(int)
    return d, Sx.loc[idx]            # Sx = all samples of the cohort (z-scoring within cohort as in B1)


def load_olaparib_task(t):
    IMP = NAR / "22_user_data_results/_data/loader_out/imports/ovary_gse229119_noMCF7"
    Sc = pd.read_csv(IMP / "pathway_scores.tsv.gz", sep="\t", index_col=0); Sc.index = Sc.index.astype(str)
    R = pd.read_csv(IMP / "response.tsv", sep="\t"); R["sample"] = R["sample"].astype(str)
    y = R[R["drug"] == "OLAPARIB"].groupby("sample")["response"].median()
    TT2 = ROOT / "robust" / "results" / "training_table_v2"
    P = pd.read_csv(TT2 / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    SC = pd.concat([pd.read_csv(TT2 / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")]).loc[P.index]
    P["dup"] = P["possible_duplicate_of"].notna()
    cols = [c for c in Sc.columns if c in SC.columns]
    x = P[(P.cohort == t["cohort"]) & ~P["dup"]]
    Zp = SC.loc[x.index, cols]; Zp = ((Zp - Zp.mean()) / Zp.std().replace(0, np.nan)).fillna(0.0)
    tr = x["drugs"].fillna("").str.contains(r"CISPLATIN|CARBOPLATIN|PLATINUM", regex=True)
    ok = x["os_months"].notna() & x["os_event"].notna() & (x["os_months"] > 0)
    keep = x.index[tr.values & ok.values]
    ev = pd.DataFrame({"T": x.loc[keep, "os_months"].astype(float).values, "E": x.loc[keep, "os_event"].astype(int).values}, index=keep)
    prox = pd.read_csv(NAR / "22_user_data_results/_data/cache/prox_OLAPARIB.tsv", sep="\t", index_col=0).iloc[:, 0]
    # NOTE: original used the median of the whole treated arm for the split (same as here: arm = treated)
    return dict(feats=cols, Xo=Sc.loc[y.index, cols].values.astype(float), y=y.values.astype(float), z=prox.reindex(cols),
                org=list(y.index), Xp=Zp.loc[keep].values, ev=ev, ZpAll=Zp, Sc=Sc[cols])


# ------------------------------------------------------------------------------------------------ the models
def fit_all_models(D, t):
    feats, Xo, y = D["feats"], D["Xo"], D["y"]
    F = [i for i, p in enumerate(feats) if pd.notna(D["z"].iloc[i]) and D["z"].iloc[i] <= CUT]
    Xz, mu, sd = B3.zcols(Xo); yz = zy(y)
    Xp = D["Xp"]
    prior = (-D["z"].fillna(0).values).astype(float)
    scores, hps = {}, []
    # --- the original (paper) methods
    if len(F) >= 2:
        c = ridge_paper(Xz[:, F], yz)
        o = np.argsort(-np.abs(c), kind="stable")
        scores["Paper_rank1"] = Xp[:, F][:, o[:1]] @ c[o[:1]]
        scores["Paper_top10"] = Xp[:, F][:, o[:10]] @ c[o[:10]]
        D["rank1_pathway"] = feats[F[o[0]]]
        D["top10_pathways"] = [feats[F[i]] for i in o[:10]]
        if AA in feats and feats.index(AA) in F:
            k = F.index(feats.index(AA))
            scores["Paper_AA_pathway"] = Xp[:, feats.index(AA)] * c[k]
            D["AA_rank_in_shortlist"] = int(np.where(o == k)[0][0]) + 1
    # --- proliferation baseline: mean of cell-cycle pathways, sign learned from the lab models
    pr = [i for i, p in enumerate(feats) if any(k in p for k in KW)]
    if pr:
        from scipy.stats import spearmanr
        r = spearmanr(Xz[:, pr].mean(1), yz)[0]
        sg = -1.0 if (np.isfinite(r) and r < 0) else 1.0
        scores["Proliferation_score"] = sg * Xp[:, pr].mean(1)
        D["prolif_sign"] = sg; D["n_prolif"] = len(pr)
    # --- the zoo
    for name in Z.MODEL_NAMES:
        if name.endswith("shortlist") and len(F) < 3:
            continue
        f, hp, best = Z.fit_model(name, Xz, yz, feats, F, groups=None, prior=prior)
        scores[name] = f(Xp)
        hps.append(dict(model=name, inner_cv_spearman=best, **hp))
        print(f"   {name} {hp} inner={best:.2f}", flush=True)
    return pd.DataFrame(scores), pd.DataFrame(hps), F, Xz, yz, prior


def metric_value(kind, ev, s, sign=1):
    if kind == "bin":
        return S.auc_fast(ev["label"].values, sign * np.asarray(s))
    return S.cindex(ev["T"].values, ev["E"].values, s)


def random_sets(D, t, F, Xz, yz, prior, ev):
    rng = np.random.default_rng(Z.SEED)
    feats = D["feats"]; rows = []
    k = max(len(F), 2); sign = t.get("sign", 1)
    for i in range(N_RANDOM):
        ix = np.sort(rng.choice(len(feats), k, replace=False))
        f, hp, _ = Z.fit_model("L2_shortlist", Xz, yz, feats, ix, groups=None, prior=prior)
        s = f(D["Xp"])
        row = dict(draw=i, metric=metric_value(t["kind"], ev, s, sign))
        if t["kind"] == "bin":
            row["p"] = S.auc_stats(ev["label"].values, sign * s)["p"]
        else:
            row["p"] = S.cox(ev["T"].values, ev["E"].values, s)["p"]; row["lr_p"] = S.split_stats(ev["T"].values, ev["E"].values, s)["lr_p"]
        rows.append(row)
    return pd.DataFrame(rows)


def corr_signature(D, Xz_or_raw):
    """original ovarian signature: mean z of the 10 pathways most positively correlated with IC50 minus the 10 most negative (6 lines)."""
    from scipy.stats import rankdata
    X = D["Xo"]; y = D["y"]
    R = rankdata(X, axis=0); ry = rankdata(y)
    Rc = (R - R.mean(0)) / (R.std(0) + 1e-12); yc = (ry - ry.mean()) / ry.std()
    rho = pd.Series((Rc * yc[:, None]).mean(0), index=D["feats"]).dropna().sort_values()
    pos, neg = list(rho.index[-10:]), list(rho.index[:10])
    return pos, neg


def run(task):
    t0 = time.time(); t = TASKS[task]
    D = load_olaparib_task(t) if task.startswith("ovarian_olaparib") else load_bundle_task(t)
    ev = D["ev"]
    print(task, "organoids", len(D["y"]), "features", len(D["feats"]), "eval n", len(ev), flush=True)
    scores, hps, F, Xz, yz, prior = fit_all_models(D, t)
    scores.index = ev.index
    if task.startswith("ovarian_olaparib"):
        pos, neg = corr_signature(D, None)
        Zp = D["ZpAll"]; ix = ev.index
        scores["Corr_signature"] = (Zp.loc[ix, pos].mean(1) - Zp.loc[ix, neg].mean(1)).values
        D["corr_pos"], D["corr_neg"] = pos, neg
    if task == "rectal_GSE87211_5FU":
        fr = pd.read_csv(NAR / "05_A3_B1/tables/B1_frozen_signatures.tsv", sep="\t")
        r = fr[fr.sig_id == "COAD_vandewetering__COAD__FLUOROURACIL|top10"].iloc[0]
        pw = r.pathways.split(";"); co = np.array([float(x) for x in r.coefs.split(";")])
        Zc = D["Zc"]; sc = Zc.loc[ev.index, pw].values @ co
        scores["Paper_top10_frozen"] = sc
        D["same_top10_as_frozen"] = set(pw) == set(D["top10_pathways"])
    rnd = random_sets(D, t, F, Xz, yz, prior, ev)
    res = dict(task=task, kind=t["kind"], label=t["label"], orig=t["orig"], sign=t.get("sign", 1), eval=ev, scores=scores, hps=hps,
               random=rnd, n_org=len(D["y"]), n_short=len(F), n_feat=len(D["feats"]),
               meta={k: D[k] for k in ("rank1_pathway", "top10_pathways", "AA_rank_in_shortlist", "prolif_sign", "n_prolif", "n_mut", "same_top10_as_frozen", "corr_pos", "corr_neg") if k in D})
    pickle.dump(res, open(OUT / f"{task}.pkl", "wb"))
    print("saved", task, f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    for task in sys.argv[1:]:
        run(task)
