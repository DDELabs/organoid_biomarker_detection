"""Combine every patient cohort with tumour gene activity, treatment and outcome into one training table.

Sources: TCGA (33 cancers), the curated GEO cohorts under data/curated/trials, and the new cohorts prepared by
build_training_cohorts.py. Output (robust/results/training_table):
  patients.tsv.gz         one row per tumour sample: cohort, cancer, drugs, outcome, covariates
  pathway_scores_*.tsv.gz Reactome pathway activity (within-sample rank score), same row ids
  cohorts.tsv             one row per cohort with usable counts
  drug_counts.tsv         patients with an outcome per cancer group and drug
  skipped.tsv             downloaded cohorts that could not be used, and why
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import build_training_cohorts as BC  # noqa: E402
import run_survlearn as RS  # noqa: E402
import run_triad as RT  # noqa: E402
from obd import RESULTS, scoring  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd.cohorts import to_canonical  # noqa: E402
from obd.curated import load_curated_trial, trial_catalog  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

OUT = RESULTS / "training_table"
COLS = ["cohort", "source", "cancer", "cancer_group", "setting", "arm", "drugs", "treated", "randomised", "control_arm",
        "responder", "os_months", "os_event", "pfs_months", "pfs_event", "rfs_months", "rfs_event", "dfs_months", "dfs_event",
        "age", "sex", "stage", "er", "her2", "patient"]
GROUP = [("breast", "breast"), ("colo", "colorectal"), ("rectal", "colorectal"), ("colon", "colorectal"), ("bladder", "bladder"),
         ("ovar", "ovarian"), ("hepato", "liver"), ("glio", "brain"), ("kidney", "kidney"), ("renal", "kidney"),
         ("melanoma", "melanoma"), ("nsclc", "lung"), ("oesoph", "oesophagus"), ("lung", "lung")]
TCGA_GROUP = {"BRCA": "breast", "COAD": "colorectal", "READ": "colorectal", "BLCA": "bladder", "OV": "ovarian", "LIHC": "liver",
              "GBM": "brain", "LGG": "brain", "KIRC": "kidney", "KIRP": "kidney", "KICH": "kidney", "SKCM": "melanoma",
              "LUAD": "lung", "LUSC": "lung", "ESCA": "oesophagus", "STAD": "stomach", "PAAD": "pancreas", "HNSC": "head and neck",
              "SARC": "sarcoma", "UCEC": "uterus", "UCS": "uterus", "CESC": "cervix", "PRAD": "prostate", "TGCT": "testis"}


# cohorts from the same hospital: a new sample whose pathway profile correlates > 0.95 with an earlier one is flagged
DUP_CHECK = [(["GSE25066"], ["GSE20194", "GSE20271", "GSE22093", "GSE23988"]), (["GSE48276", "GSE70691"], ["GSE52219"])]


def group(c):
    c = str(c).lower()
    return next((g for k, g in GROUP if k in c), c)


def score(name, expr):
    def compute():
        e = to_canonical(expr)
        gs = {k: [g for g in v if g in e.index] for k, v in gene_sets(("REACTOME",)).items()}
        return scoring.rank_score(e, {k: v for k, v in gs.items() if len(v) >= 5})
    return scoring.cached(name, compute).T


def tidy(clin, cohort, source, cols=COLS):
    c = clin.copy()
    c["cohort"], c["source"] = cohort, source
    if "cancer" not in c:
        c["cancer"] = np.nan
    c["cancer_group"] = c["cancer"].map(group)
    if "drugs" not in c:
        c["drugs"] = np.nan
    for x in cols:
        if x not in c:
            c[x] = np.nan
    none = c["drugs"].isin(["", "NONE"]) | (c["drugs"].isna() & c["control_arm"].eq(True))
    c["drugs"] = c["drugs"].mask(none, "NONE")   # NONE = known to have had no drug; empty = not recorded
    c["treated"] = np.where(c["drugs"].isna(), np.nan, (~none).astype(float))
    c.index = [f"{cohort}:{i}" for i in c.index]
    return c[cols]


def tcga():
    M, proj = RT.patient_matrix(R.PROJECTS)
    _, cdr, long = RS.load()
    cdr = cdr.loc[[p for p in cdr.index if p in M.index]]
    first = long[long["before_event"]]
    drugs = first.groupby("patient")["drug"].apply(lambda s: ";".join(sorted(set(s))))
    lab = first[first["responder"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")["responder"]
    raw = R.load_cdr().loc[cdr.index]
    c = pd.DataFrame(index=cdr.index)
    c["cancer"] = "TCGA-" + cdr["project"]
    c["setting"], c["arm"] = "as recorded by the hospital", np.nan
    c["drugs"] = drugs.reindex(c.index).fillna("")
    c["randomised"], c["control_arm"] = False, c["drugs"].eq("")
    c["responder"] = lab.reindex(c.index).astype(float)
    c["os_months"], c["os_event"] = raw["OS.time"] / 30.44, raw["OS"]
    c["pfs_months"], c["pfs_event"] = raw["PFI.time"] / 30.44, raw["PFI"]
    c["age"], c["sex"], c["stage"] = raw["age"], raw["gender"].str.capitalize(), raw["stage"]
    c["patient"] = c.index
    out = []
    for p, idx in c.groupby(cdr["project"]).groups.items():
        t = tidy(c.loc[idx], f"TCGA-{p}", "TCGA")
        t["cancer_group"] = TCGA_GROUP.get(p, p)
        out.append(t)
    S = M.loc[c.index].copy()
    S.index = [f"TCGA-{cdr.loc[i, 'project']}:{i}" for i in S.index]
    return pd.concat(out), S


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P, S = tcga()
    paths = list(S.columns)
    clin_all, sc_all = [P], [S]
    cat = trial_catalog()
    for _, r in cat[cat["usable"]].iterrows():
        cid = r["cohort_id"]
        expr, clin = load_curated_trial(cid, canonical=False)
        s = score(f"curated_{cid}_reactome_rank", expr)
        if "cancer" not in clin:
            clin["cancer"] = r["cancer"]
        clin_all.append(tidy(clin.loc[[i for i in clin.index if i in s.index]], cid, "GEO (already in use)"))
        s = s.reindex(columns=paths); s.index = [f"{cid}:{i}" for i in s.index]
        sc_all.append(s.loc[clin_all[-1].index])
    for fn in BC.BUILDERS:
        cid = fn.__name__
        d = BC.OUT / cid
        clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col=0, low_memory=False)
        s = score(f"train_{cid}_reactome_rank", pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0))
        clin_all.append(tidy(clin.loc[[i for i in clin.index if i in s.index]], cid, "GEO (new)"))
        s = s.reindex(columns=paths); s.index = [f"{cid}:{i}" for i in s.index]
        sc_all.append(s.loc[clin_all[-1].index])
        print(cid, s.shape, "pathways missing:", int(s.isna().all().sum()), flush=True)
    P, S = pd.concat(clin_all), pd.concat(sc_all)
    P.index.name = S.index.name = "sample_id"
    P["possible_duplicate_of"] = pd.Series(index=P.index, dtype=object)
    for new, old in DUP_CHECK:
        a, b = S[P["cohort"].isin(new).values].dropna(axis=1), S[P["cohort"].isin(old).values]
        b = b[a.columns].dropna(axis=1); a = a[b.columns]
        r = np.corrcoef(((a - a.mean()) / a.std()).values, ((b - b.mean()) / b.std()).values)[:len(a), len(a):]
        for i in np.where(r.max(axis=1) > 0.95)[0]:
            P.loc[a.index[i], "possible_duplicate_of"] = b.index[r[i].argmax()]
    ev = P[["os_event", "pfs_event", "rfs_event", "dfs_event"]]
    P["has_survival"], P["has_response"] = ev.notna().any(axis=1), P["responder"].notna()
    P["has_outcome"] = P["has_survival"] | P["has_response"]
    P.to_csv(OUT / "patients.tsv.gz", sep="\t")
    tc = P["source"].eq("TCGA")
    S[tc.values].round(3).to_csv(OUT / "pathway_scores_tcga.tsv.gz", sep="\t")
    S[~tc.values].round(3).to_csv(OUT / "pathway_scores_geo.tsv.gz", sep="\t")
    g = P.groupby(["source", "cohort"], sort=False)
    C = pd.DataFrame({"cancer_group": g["cancer_group"].agg(lambda s: s.mode().iat[0]), "samples": g.size(),
                      "treated": g["treated"].sum().astype(int), "untreated_comparison": g["treated"].apply(lambda s: int((s == 0).sum())),
                      "with_response": g["has_response"].sum(), "responders": g["responder"].sum().astype(int),
                      "with_survival": g["has_survival"].sum(), "possible_duplicates": g["possible_duplicate_of"].count(),
                      "events": ev.groupby([P["source"], P["cohort"]], sort=False).max().max(axis=1).pipe(lambda _: ev.fillna(0).max(axis=1).groupby([P["source"], P["cohort"]], sort=False).sum().astype(int)),
                      "randomised": g["randomised"].first()})
    C.to_csv(OUT / "cohorts.tsv", sep="\t")
    L = P[P["treated"].eq(1) & P["has_outcome"] & P["possible_duplicate_of"].isna()].assign(drug=lambda d: d["drugs"].str.split(";")).explode("drug")
    D = L.groupby(["cancer_group", "drug"]).agg(patients=("cohort", "size"), cohorts=("cohort", "nunique"),
                                                with_response=("has_response", "sum"), with_survival=("has_survival", "sum"),
                                                tcga=("source", lambda s: int((s == "TCGA").sum())))
    D = D[D["patients"] >= 20].sort_values(["cancer_group", "patients"], ascending=[True, False])
    D.to_csv(OUT / "drug_counts.tsv", sep="\t")
    pd.Series(BC.SKIPPED, name="reason").rename_axis("accession").to_csv(OUT / "skipped.tsv", sep="\t")
    print(P.shape, S.shape)
    print(P.groupby("source")[["has_outcome", "has_response", "has_survival"]].agg(["size", "sum"]).to_string())


# ---------------------------------------------------------------------------------------------------------------
# Version 2 (added for the lung harvest, nar/17_lung). Starts from the finished `training_table`, keeps every old row
# exactly as it was, and appends the cohorts in BC.LUNG_BUILDERS (plus a `histology` column). Existing code above is unchanged.
COLS_V2 = COLS + ["histology"]
# new lung cohorts are compared with the lung cohorts already in the table (same hospitals are possible)
DUP_CHECK_V2 = [([fn.__name__ for fn in BC.LUNG_BUILDERS], ["GSE42127", "TCGA-LUAD", "TCGA-LUSC"])]


def _flag_duplicates(P, S, checks):
    for new, old in checks:
        a, b = S[P["cohort"].isin(new).values].dropna(axis=1), S[P["cohort"].isin(old).values]
        b = b[a.columns].dropna(axis=1); a = a[b.columns]
        r = np.corrcoef(((a - a.mean()) / a.std()).values, ((b - b.mean()) / b.std()).values)[:len(a), len(a):]   # same test as DUP_CHECK
        for i in np.where(r.max(axis=1) > 0.95)[0]:
            P.loc[a.index[i], "possible_duplicate_of"] = b.index[r[i].argmax()]
    return P


def write_tables(P, S, out, skipped):
    """Same summary tables as main(), written to `out`."""
    out.mkdir(parents=True, exist_ok=True)
    ev = P[["os_event", "pfs_event", "rfs_event", "dfs_event"]]
    P["has_survival"], P["has_response"] = ev.notna().any(axis=1), P["responder"].notna()
    P["has_outcome"] = P["has_survival"] | P["has_response"]
    P.to_csv(out / "patients.tsv.gz", sep="\t")
    tc = P["source"].eq("TCGA")
    S[tc.values].round(3).to_csv(out / "pathway_scores_tcga.tsv.gz", sep="\t")
    S[~tc.values].round(3).to_csv(out / "pathway_scores_geo.tsv.gz", sep="\t")
    g = P.groupby(["source", "cohort"], sort=False)
    C = pd.DataFrame({"cancer_group": g["cancer_group"].agg(lambda s: s.mode().iat[0]), "samples": g.size(),
                      "treated": g["treated"].sum().astype(int), "untreated_comparison": g["treated"].apply(lambda s: int((s == 0).sum())),
                      "with_response": g["has_response"].sum(), "responders": g["responder"].sum().astype(int),
                      "with_survival": g["has_survival"].sum(), "possible_duplicates": g["possible_duplicate_of"].count(),
                      "events": ev.fillna(0).max(axis=1).groupby([P["source"], P["cohort"]], sort=False).sum().astype(int),
                      "randomised": g["randomised"].first()})
    C.to_csv(out / "cohorts.tsv", sep="\t")
    L = P[P["treated"].eq(1) & P["has_outcome"] & P["possible_duplicate_of"].isna()].assign(drug=lambda d: d["drugs"].str.split(";")).explode("drug")
    D = L.groupby(["cancer_group", "drug"]).agg(patients=("cohort", "size"), cohorts=("cohort", "nunique"),
                                                with_response=("has_response", "sum"), with_survival=("has_survival", "sum"),
                                                tcga=("source", lambda s: int((s == "TCGA").sum())))
    D = D[D["patients"] >= 20].sort_values(["cancer_group", "patients"], ascending=[True, False])
    D.to_csv(out / "drug_counts.tsv", sep="\t")
    pd.Series(skipped, name="reason").rename_axis("accession").to_csv(out / "skipped.tsv", sep="\t")


def main_v2(out_name="training_table_v2"):
    out = RESULTS / out_name
    base = OUT                                         # the finished version-1 table, read only
    P = pd.read_csv(base / "patients.tsv.gz", sep="\t", index_col=0, low_memory=False)
    P = P.drop(columns=["has_survival", "has_response", "has_outcome"])
    for c in ("histology",):
        P[c] = pd.Series(np.nan, index=P.index, dtype=object)
    P.loc[P["cohort"].eq("TCGA-LUAD"), "histology"] = "adenocarcinoma"
    P.loc[P["cohort"].eq("TCGA-LUSC"), "histology"] = "squamous"
    S = pd.concat([pd.read_csv(base / f"pathway_scores_{k}.tsv.gz", sep="\t", index_col=0) for k in ("tcga", "geo")])
    paths = list(S.columns)
    clin_all, sc_all = [P], [S]
    for fn in BC.LUNG_BUILDERS:
        cid = fn.__name__
        d = BC.OUT / cid
        clin = pd.read_csv(d / "clinical.tsv", sep="\t", index_col=0, low_memory=False)
        s = score(f"train_{cid}_reactome_rank", pd.read_csv(d / "expression.tsv.gz", sep="\t", index_col=0))
        t = tidy(clin.loc[[i for i in clin.index if i in s.index]], cid, "GEO (new)", COLS_V2)
        s = s.reindex(columns=paths); s.index = [f"{cid}:{i}" for i in s.index]
        clin_all.append(t); sc_all.append(s.loc[t.index])
        print(cid, s.shape, "pathways missing:", int(s.isna().all().sum()), flush=True)
    P, S = pd.concat(clin_all), pd.concat(sc_all)
    P.index.name = S.index.name = "sample_id"
    P = _flag_duplicates(P, S, DUP_CHECK_V2)
    skipped = {**BC.SKIPPED, **BC.SKIPPED_LUNG}
    write_tables(P, S, out, skipped)
    print(P.shape, S.shape)


if __name__ == "__main__":
    if "--out" in sys.argv:        # e.g. --out training_table_v2 : extend version 1 with the lung cohorts
        main_v2(sys.argv[sys.argv.index("--out") + 1])
    else:
        main()
