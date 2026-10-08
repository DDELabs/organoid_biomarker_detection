"""Relaxed-criteria top-10 pathway biomarkers per cancer x drug, reaching cancers with few RECIST labels.

The standard analysis (run_atlas_cancer_specific.py) needs >= 20 RECIST-labelled patients and >= 5 per
class, which leaves out liver, glioblastoma, kidney, ovarian and others. Here every TCGA cancer x drug
pair is tried with the strongest endpoint it supports, in this order:

  A  RECIST response (CR/PR vs SD/PD), any anti-cancer drug, n >= 10 and >= 3 per class
  B  first-course treatment outcome (TCGA-CDR treatment_outcome_first_course: CR/PR/no measurable
     tumour vs SD/PD/persistent) among patients who received the drug first-line, n >= 15, >= 4 per class
  C  progression-free interval among patients who received the drug (Cox), n >= 20 and >= 8 events;
     also reported: the same pathway in untreated patients of that cancer and the difference
     (treated - untreated log HR) as a crude predictive-vs-prognostic check

Statistics are per-pathway score tests (logistic / Cox) on within-cancer standardised Reactome rank
scores, fast enough for bootstraps. local z > 0 = higher pathway -> better outcome. For drugs in ATLAS
the combined z = (local z + ATLAS z) / sqrt(2); for other drugs combined z = local z (no prior).
Pairs already covered by the standard analysis are kept from it (tier 'A20').

Outputs (robust/results/atlas/cancer_specific/): relaxed_top10.tsv, relaxed_all_pairs.tsv.gz,
relaxed_pairs_summary.tsv, TOP10_RELAXED.md

Usage: python robust/run_atlas_relaxed.py [--boot 100]
"""
import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from run_atlas_cancer_specific import CANCER, bh  # noqa: E402

OUT = RESULTS / "atlas" / "cancer_specific"
THERAPY = {"Chemotherapy", "Targeted Molecular therapy", "Immunotherapy", "Hormone Therapy", "Chemotherapy|Hormone Therapy"}
SUPPORTIVE = {"DEXAMETHASONE", "LEVOTHYROXINE", "LIOTHYRONINE", "LEVOTHROID", "ARMOUR THYROID", "PREDNISONE",
              "LEUCOVORIN", "MESNA", "ONDANSETRON", "FILGRASTIM", "PEGFILGRASTIM", "ZOLEDRONIC ACID",
              "DENOSUMAB", "HYDROCORTISONE", "METHYLPREDNISOLONE", "LYMPHOCYTE INFUSION"}
CR_LIKE = {"Complete Remission/Response", "Partial Remission/Response", "No Measureable Tumor or Tumor Markers"}
NR_LIKE = {"Progressive Disease", "Stable Disease", "Persistent Disease"}


def logit_score_z(X, y):
    """Per-column logistic score-test z (no covariates)."""
    yc = y - y.mean()
    Xc = X - X.mean(0)
    den = np.sqrt(y.mean() * (1 - y.mean()) * (Xc ** 2).sum(0)) + 1e-12
    return (Xc * yc[:, None]).sum(0) / den


def cox_score(X, t, e):
    """Per-column Cox score statistic U and information I at beta = 0 (Breslow, ties approximate)."""
    o = np.argsort(-t, kind="stable")
    X, e = X[o], e[o]
    n_risk = np.arange(1, len(t) + 1)[:, None]
    m1 = np.cumsum(X, 0) / n_risk
    m2 = np.cumsum(X ** 2, 0) / n_risk
    ev = e.astype(bool)
    U = (X[ev] - m1[ev]).sum(0)
    I = (m2[ev] - m1[ev] ** 2).sum(0) + 1e-12
    return U, I


def cox_score_z(X, t, e):
    U, I = cox_score(X, t, e)
    return -U / np.sqrt(I)  # > 0: higher pathway -> longer PFI


def stability(stat, n, boot, rng, ok=lambda b: True):
    hits, used = None, 0
    for _ in range(boot):
        b = rng.integers(0, n, n)
        if not ok(b):
            continue
        z = np.nan_to_num(stat(b))
        hits = np.zeros_like(z) if hits is None else hits
        hits[np.argsort(-np.abs(z))[:10]] += 1
        used += 1
    return hits / max(used, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=100)
    a = ap.parse_args()
    long = R.load_response_labels()
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    cdr = R.load_cdr()
    M, proj = RT.patient_matrix(sorted(long["project"].unique()))
    Xc = T.center_within(M, proj)
    Xc = Xc / Xc.groupby(proj).transform("std").replace(0, 1)
    pathways = list(Xc.columns)
    tot = pd.read_csv(RESULTS / "atlas" / "inference_drug_total.tsv.gz", sep="\t")
    atlas_z = tot.pivot(index="pathway", columns="drug", values="z").reindex(pathways)
    std = pd.read_csv(OUT / "all_pairs_local.tsv.gz", sep="\t")
    std_pairs = set(zip(std["cancer"], std["drug"]))
    rng = np.random.default_rng(11)
    rows, summ = [], []
    systemic = set(long["patient"])

    for (p, d), g in long.groupby(["project", "drug"]):
        if (p, d) in std_pairs:
            continue
        pats = [x for x in g["patient"].unique() if x in Xc.index]
        if len(pats) < 10:
            continue
        res = None
        # tier A: RECIST
        lab = g[g["recist"].notna()].sort_values("start_days").drop_duplicates("patient").set_index("patient")
        lab = lab[lab.index.isin(Xc.index)]
        y = lab["responder"].astype(float).values
        if len(lab) >= 10 and y.sum() >= 3 and (1 - y).sum() >= 3:
            X = Xc.loc[lab.index].values
            res = ("A", "RECIST response", logit_score_z(X, y), len(y), int(y.sum()),
                   stability(lambda b: logit_score_z(X[b], y[b]), len(y), a.boot, rng, lambda b: 0 < y[b].sum() < len(b)))
        # tier B: first-course outcome among first-line recipients
        if res is None:
            first = g[g["setting_proxy"].eq("early")]["patient"].unique()
            oc = cdr.reindex([x for x in first if x in Xc.index])["treatment_outcome_first_course"]
            oc = oc[oc.isin(CR_LIKE | NR_LIKE)]
            y = oc.isin(CR_LIKE).astype(float).values
            if len(oc) >= 15 and y.sum() >= 4 and (1 - y).sum() >= 4:
                X = Xc.loc[oc.index].values
                res = ("B", "first-course outcome", logit_score_z(X, y), len(y), int(y.sum()),
                       stability(lambda b: logit_score_z(X[b], y[b]), len(y), a.boot, rng, lambda b: 0 < y[b].sum() < len(b)))
        # tier C: PFI among treated
        extra = {}
        if res is None:
            s = cdr.reindex(pats)[["PFI", "PFI.time"]].dropna()
            s = s[s["PFI.time"] > 0]
            if len(s) >= 20 and s["PFI"].sum() >= 8:
                X, t, e = Xc.loc[s.index].values, s["PFI.time"].values, s["PFI"].values
                z = cox_score_z(X, t, e)
                res = ("C", "PFI among treated (Cox)", z, len(s), int(e.sum()),
                       stability(lambda b: cox_score_z(X[b], t[b], e[b]), len(t), a.boot, rng, lambda b: e[b].sum() >= 3))
                un = cdr[(cdr["project"] == p) & ~cdr.index.isin(systemic) & cdr.index.isin(Xc.index)][["PFI", "PFI.time"]].dropna()
                un = un[un["PFI.time"] > 0]
                if len(un) >= 20 and un["PFI"].sum() >= 8:
                    Ut, It = cox_score(X, t, e)
                    Uu, Iu = cox_score(Xc.loc[un.index].values, un["PFI.time"].values, un["PFI"].values)
                    extra = {"untreated_z": -Uu / np.sqrt(Iu),
                             "treated_minus_untreated_z": -(Ut / It - Uu / Iu) / np.sqrt(1 / It + 1 / Iu),
                             "n_untreated": len(un)}
        if res is None:
            continue
        tier, endpoint, z, n, k, stab = res
        az = atlas_z[d].values if d in atlas_z else np.full(len(pathways), np.nan)
        comb = np.where(np.isnan(az), z, (z + np.nan_to_num(az)) / np.sqrt(2))
        df = pd.DataFrame({"cancer": p, "cancer_name": CANCER.get(p, p), "drug": d, "tier": tier, "endpoint": endpoint,
                           "n": n, "responders_or_events": k, "pathway": pathways, "local_z": z,
                           "local_q": bh(2 * norm.sf(np.abs(z))), "atlas_z": az, "combined_z": comb,
                           "top10_stability": stab})
        for kk, v in extra.items():
            df[kk] = v
        rows.append(df)
        summ.append({"cancer": p, "drug": d, "tier": tier, "endpoint": endpoint, "n": n, "responders_or_events": k,
                     "min_local_q": float(df["local_q"].min()), "max_stability": float(stab.max())})
        print(p, d, tier, n, flush=True)

    allp = pd.concat(rows, ignore_index=True)
    allp.to_csv(OUT / "relaxed_all_pairs.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summ).sort_values(["tier", "cancer", "drug"]).to_csv(OUT / "relaxed_pairs_summary.tsv", sep="\t", index=False)
    tops = []
    for _, g in allp.groupby(["cancer", "drug"]):
        r = g.sort_values("combined_z", ascending=False).head(5).assign(direction="response", rank=range(1, 6))
        s = g.sort_values("combined_z").head(5).assign(direction="resistance", rank=range(1, 6))
        tops += [r, s]
    top = pd.concat(tops, ignore_index=True)
    top["pathway_name"] = top["pathway"].str.replace("REACTOME_", "").str.replace("_", " ").str.lower()
    top.to_csv(OUT / "relaxed_top10.tsv", sep="\t", index=False)

    L = ["# Top-10 pathway biomarkers per cancer x drug - relaxed criteria", "",
         "Pairs not covered by the standard analysis (TOP10_BY_CANCER.md). Tier A = RECIST (n >= 10, >= 3 per class); "
         "B = TCGA first-course treatment outcome among first-line recipients; C = progression-free interval among "
         "treated patients (Cox; positive z = longer PFI). For tier C, `pred` = treated minus untreated effect z "
         "(crude predictive check; untreated = no recorded systemic therapy). `*` = local q < 0.1. "
         "Combined z uses the ATLAS pan-cancer prior when the drug is in ATLAS. Exploratory.", ""]
    for c, gc in top.groupby("cancer"):
        L.append(f"## {c} - {CANCER.get(c, c)}")
        for d, gd in gc.groupby("drug"):
            h = gd.iloc[0]
            L.append(f"\n**{d.title()}** - tier {h.tier} ({h.endpoint}), n = {h.n}, "
                     f"{'events' if h.tier == 'C' else 'responders'} = {h.responders_or_events}\n")
            pred = "treated_minus_untreated_z" in gd and gd["treated_minus_untreated_z"].notna().any()
            L.append("| # | Better-outcome pathway | local z | comb z | stab" + (" | pred" if pred else "") +
                     " | Worse-outcome pathway | local z | comb z | stab" + (" | pred" if pred else "") + " |")
            L.append("|---" * (10 if pred else 9) + "|")
            rs = gd[gd.direction == "response"].reset_index(drop=True)
            rz = gd[gd.direction == "resistance"].reset_index(drop=True)

            def cell(r):
                s = (f"{r.pathway_name[:55]}{'*' if r.local_q < 0.1 else ''} | {r.local_z:+.1f} | "
                     f"{r.combined_z:+.1f} | {r.top10_stability:.2f}")
                return s + (f" | {r.treated_minus_untreated_z:+.1f}" if pred else "")
            for i in range(5):
                L.append(f"| {i + 1} | {cell(rs.iloc[i])} | {cell(rz.iloc[i])} |")
        L.append("")
    (OUT / "TOP10_RELAXED.md").write_text("\n".join(L) + "\n")
    print("relaxed pairs:", len(summ), "cancers:", allp.cancer.nunique(), "tiers:",
          pd.DataFrame(summ)["tier"].value_counts().to_dict())


if __name__ == "__main__":
    main()
