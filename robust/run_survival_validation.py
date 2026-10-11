"""Survival validation of the per-cancer x drug top-10 pathway biomarkers in TCGA patients NOT used for discovery.

Discovery used RECIST labels (tier A20 / A) or first-course outcome (B) or PFI (C). Here, for every
cancer x drug pair, the signature score = mean z of its 5 response pathways - mean z of its 5 resistance
pathways (within-cancer standardised Reactome rank scores) is tested on:

  held-out treated  patients who received the drug but were not in the discovery set (for tier C
                    every treated patient was used for PFI, so only OS is reported and it is not independent)
  untreated         patients of the same cancer with no recorded systemic therapy (prognostic control)
  interaction       score x treated (held-out treated vs untreated): is the survival benefit of a high
                    score larger under the drug? (predictive test)

Survival: TCGA-CDR overall survival (OS) and progression-free interval (PFI). Statistics:
  paper-style      median split of treated patients, log-rank p (as in python/*_pathway_prediction.py)
  Cox              HR per SD of the score, adjusted for age and stage
Expected direction if the biomarker holds: HR < 1 in treated patients, and interaction HR < 1.

Outputs: robust/results/survival_validation/{pairs.tsv, pathways.tsv.gz, SURVIVAL_VALIDATION.md}
Usage: python robust/run_survival_validation.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest, norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import survival as SV  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.cohorts import _stage  # noqa: E402
from run_atlas_relaxed import CR_LIKE, NR_LIKE, SUPPORTIVE, THERAPY  # noqa: E402

OUT = RESULTS / "survival_validation"
D = RESULTS / "atlas" / "cancer_specific"
MIN_N, MIN_EV = 20, 8


def cox(df, cols):
    """HR, p for the first column of a Cox model with the listed columns (complete cases)."""
    d = df.dropna(subset=["time", "event"] + cols)
    d = d[d["time"] > 0]
    if len(d) < MIN_N or d["event"].sum() < MIN_EV:
        return np.nan, np.nan, len(d), int(d["event"].sum()) if len(d) else 0
    X = d[cols].astype(float).values
    keep = X.std(0) > 0
    keep[0] = True
    try:
        r = PHReg(d["time"].values, X[:, keep], status=d["event"].values, ties="efron").fit(disp=0)
        return float(np.exp(r.params[0])), float(r.pvalues[0]), len(d), int(d["event"].sum())
    except Exception:
        return np.nan, np.nan, len(d), int(d["event"].sum())


def covars(cdr, pats):
    c = cdr.reindex(pats)
    st = c["stage"].map(_stage)
    out = pd.DataFrame({"age": c["age"].astype(float)}, index=pats)
    out["age"] = out["age"].fillna(out["age"].median())
    for s in ("II", "III", "IV"):
        out[f"stage_{s}"] = (st == s).astype(float)
    return out


def endpoint(cdr, pats, ep):
    c = cdr.reindex(pats)
    return pd.DataFrame({"time": c[f"{ep}.time"].astype(float) / 30.44, "event": c[ep].astype(float)}, index=pats)


def discovery_patients(long, cdr, Xc, cancer, drug, tier):
    g = long[(long["project"] == cancer) & (long["drug"] == drug)]
    if tier in ("A20", "A"):
        return set(g[g["recist"].notna()]["patient"])
    if tier == "B":
        first = g[g["setting_proxy"].eq("early")]["patient"].unique()
        oc = cdr.reindex(first)["treatment_outcome_first_course"]
        return set(oc[oc.isin(CR_LIKE | NR_LIKE)].index)
    return set(g["patient"])  # tier C used every treated patient


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    long = R.load_response_labels()
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    cdr = R.load_cdr()
    top = pd.concat([pd.read_csv(D / "top10_by_cancer_drug.tsv", sep="\t").assign(tier="A20"),
                     pd.read_csv(D / "relaxed_top10.tsv", sep="\t")], ignore_index=True)
    M, proj = RT.patient_matrix(sorted(top["cancer"].unique()))
    Xc = T.center_within(M, proj)
    systemic = set(long["patient"])
    pair_rows, pw_rows = [], []
    for (cancer, drug), g in top.groupby(["cancer", "drug"]):
        tier = g["tier"].iloc[0]
        inC = proj.index[proj == cancer]
        treated = set(long[(long["project"] == cancer) & (long["drug"] == drug)]["patient"]) & set(inC)
        disc = discovery_patients(long, cdr, Xc, cancer, drug, tier)
        held = sorted(treated - disc) if tier != "C" else sorted(treated)
        untreated = sorted(set(inC) - systemic)
        resp = g[g["direction"] == "response"]["pathway"].tolist()
        res = g[g["direction"] == "resistance"]["pathway"].tolist()
        Z = Xc.loc[inC]
        score = Z[resp].mean(1) - Z[res].mean(1)
        score = (score - score.mean()) / score.std()
        row = {"cancer": cancer, "drug": drug, "tier": tier, "n_treated": len(treated), "n_discovery": len(disc & treated),
               "n_heldout": len(held), "n_untreated": len(untreated), "independent": tier != "C"}
        for ep in ("OS", "PFI"):
            if tier == "C" and ep == "PFI":
                continue
            E = endpoint(cdr, list(inC), ep)
            C = covars(cdr, list(inC))
            dfh = pd.concat([E.loc[held], C.loc[held]], axis=1).assign(score=score.loc[held])
            hr, p, n, ev = cox(dfh, ["score", "age", "stage_II", "stage_III", "stage_IV"])
            row.update({f"{ep}_treated_HR": hr, f"{ep}_treated_p": p, f"{ep}_treated_n": n, f"{ep}_treated_events": ev})
            # paper-style median split log-rank among held-out treated
            dd = dfh.dropna(subset=["time", "event"])
            dd = dd[dd["time"] > 0]
            if len(dd) >= MIN_N and dd["event"].sum() >= MIN_EV:
                hi = dd["score"] > dd["score"].median()
                row[f"{ep}_treated_logrank_p"] = SV.logrank(dd.loc[hi, "time"], dd.loc[hi, "event"], dd.loc[~hi, "time"], dd.loc[~hi, "event"])
            dfu = pd.concat([E.loc[untreated], C.loc[untreated]], axis=1).assign(score=score.loc[untreated])
            hr, p, n, ev = cox(dfu, ["score", "age", "stage_II", "stage_III", "stage_IV"])
            row.update({f"{ep}_untreated_HR": hr, f"{ep}_untreated_p": p, f"{ep}_untreated_n": n})
            both = pd.concat([dfh.assign(treated=1.0), dfu.assign(treated=0.0)])
            both["sxt"] = both["score"] * both["treated"]
            hr, p, _, _ = cox(both, ["sxt", "score", "treated", "age", "stage_II", "stage_III", "stage_IV"])
            row.update({f"{ep}_interaction_HR": hr, f"{ep}_interaction_p": p})
            # single pathways in held-out treated
            for pw, sign in [(x, 1) for x in resp] + [(x, -1) for x in res]:
                d1 = dfh.assign(score=Z.loc[held, pw].values)
                h1, p1, n1, _ = cox(d1, ["score", "age", "stage_II", "stage_III", "stage_IV"])
                pw_rows.append({"cancer": cancer, "drug": drug, "tier": tier, "endpoint": ep, "pathway": pw,
                                "direction": "response" if sign > 0 else "resistance", "n": n1, "HR_per_SD": h1, "p": p1,
                                "concordant": (h1 < 1) if sign > 0 else (h1 > 1)})
        pair_rows.append(row)
        print(cancer, drug, tier, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()
                                   if k in ("n_heldout", "OS_treated_HR", "OS_treated_p", "OS_interaction_HR")}, flush=True)
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT / "pairs.tsv", sep="\t", index=False)
    pw = pd.DataFrame(pw_rows)
    pw.to_csv(OUT / "pathways.tsv.gz", sep="\t", index=False)
    report(pairs, pw)


def report(pairs, pw):
    L = ["# Survival validation of the top-10 pathway biomarkers (TCGA, held-out patients)", "",
         "Signature score per cancer x drug = mean of its 5 response pathways - mean of its 5 resistance pathways. "
         "Held-out treated = received the drug but not used for discovery. HR per SD, Cox adjusted for age and stage "
         "(HR < 1 = higher score, longer survival). Interaction HR < 1 = the survival advantage of a high score is larger "
         "in treated than in untreated patients (predictive). Logrank = the original paper's median split.", ""]
    for ep in ("OS", "PFI"):
        ok = pairs[pairs["independent"] & pairs[f"{ep}_treated_HR"].notna()]
        k = int((ok[f"{ep}_treated_HR"] < 1).sum())
        bt = binomtest(k, len(ok), 0.5, alternative="greater").pvalue if len(ok) else np.nan
        sig = ok[(ok[f"{ep}_treated_p"] < 0.05)]
        z = norm.isf(ok[f"{ep}_treated_p"].clip(1e-12, 1) / 2) * np.sign(1 - ok[f"{ep}_treated_HR"])
        L += [f"## {ep}: independent pairs", "",
              f"- evaluable pairs: {len(ok)}; HR < 1 (expected direction): {k}/{len(ok)} (binomial p = {bt:.3f})",
              f"- nominal p < 0.05: {len(sig)} ({int((sig[f'{ep}_treated_HR'] < 1).sum())} in the expected direction)",
              f"- Stouffer z across pairs (expected direction positive): {z.sum() / np.sqrt(len(z)):.2f}" if len(z) else "",
              f"- median-split log-rank p < 0.05 (paper style): {int((ok.get(f'{ep}_treated_logrank_p', pd.Series(dtype=float)) < 0.05).sum())}",
              f"- untreated HR < 1: {int((ok[f'{ep}_untreated_HR'] < 1).sum())}/{int(ok[f'{ep}_untreated_HR'].notna().sum())}; "
              f"interaction HR < 1: {int((ok[f'{ep}_interaction_HR'] < 1).sum())}/{int(ok[f'{ep}_interaction_HR'].notna().sum())}, "
              f"of which p < 0.05: {int(((ok[f'{ep}_interaction_HR'] < 1) & (ok[f'{ep}_interaction_p'] < 0.05)).sum())}", ""]
        cols = ["cancer", "drug", "tier", "n_heldout", f"{ep}_treated_events", f"{ep}_treated_HR", f"{ep}_treated_p",
                f"{ep}_treated_logrank_p", f"{ep}_untreated_HR", f"{ep}_interaction_HR", f"{ep}_interaction_p"]
        cols = [c for c in cols if c in ok]
        t = ok.sort_values(f"{ep}_treated_p")[cols].head(15)
        L += ["Top 15 pairs by treated p:", "", "| " + " | ".join(c.replace(f"{ep}_", "") for c in cols) + " |",
              "|" + "---|" * len(cols)]
        for r in t.itertuples(index=False):
            L.append("| " + " | ".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in r) + " |")
        L.append("")
        q = pw[(pw["endpoint"] == ep) & pw["HR_per_SD"].notna() & pw["tier"].ne("C")]
        if len(q):
            kk = int(q["concordant"].sum())
            L += [f"Single pathways ({ep}, held-out treated): {kk}/{len(q)} in the expected direction "
                  f"(binomial p = {binomtest(kk, len(q), 0.5, alternative='greater').pvalue:.3f}; pathways overlap, so optimistic); "
                  f"p < 0.05 and concordant: {int(((q['p'] < 0.05) & q['concordant']).sum())}, discordant: "
                  f"{int(((q['p'] < 0.05) & ~q['concordant']).sum())}", ""]
    (OUT / "SURVIVAL_VALIDATION.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
