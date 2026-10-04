"""Run the full robust analysis for one drug across data versions and write a report."""
import dataclasses
import json
import warnings

import numpy as np
import pandas as pd

from . import plots
from . import survival as S
from .study import Study

warnings.filterwarnings("ignore")


def jackknife(study, top_k=7, n_boot=100):
    """Leave one organoid out: signature membership and patient HR each time."""
    member, hrs = {}, []
    for i, s in enumerate(study.samples):
        keep = np.arange(len(study.y)) != i
        w, _ = study.signature(y=study.y[keep], X=study.X[keep], top_k=top_k, n_boot=n_boot, seed=i)
        for f in w.index:
            member[f] = member.get(f, 0) + 1
        hrs.append(study.evaluate(w, interaction=False).get("adj_HR", np.nan))
    n = len(study.samples)
    return pd.Series(member) / n, np.array(hrs, float)


def patient_bootstrap(study, weights, n=200, seed=3):
    """Resample treated patients; distribution of the adjusted HR."""
    rng = np.random.default_rng(seed)
    tp = np.array(study.treated_pats)
    score = study.patient_score(weights, list(tp))
    base = S.covariate_frame(study.clinical, list(tp))
    base["score"] = (score - score.mean()) / score.std()
    base["months"] = study.clinical.loc[tp, "months"].values
    base["event"] = study.clinical.loc[tp, "event"].values
    out = []
    for _ in range(n):
        d = base.iloc[rng.integers(0, len(base), len(base))].reset_index(drop=True)
        try:
            out.append(S.cox(d, ["score", "stage_III", "stage_IV", "age", "sex"]).loc["score", "HR"])
        except Exception:
            continue
    out = np.array(out, float)
    return out[np.isfinite(out)]


def pathway_patient_hr(study, pathways):
    """Per-pathway adjusted Cox HR (per SD) in treated patients."""
    rows = {}
    tp = study.treated_pats
    for p in pathways:
        d = S.covariate_frame(study.clinical, tp)
        d["x"] = study.P.loc[tp, p].values
        d["months"] = study.clinical.loc[tp, "months"].values
        d["event"] = study.clinical.loc[tp, "event"].values
        try:
            r = S.cox(d, ["x", "stage_III", "stage_IV", "age", "sex"]).loc["x"]
            rows[p] = (r.HR, r.p)
        except Exception:
            rows[p] = (np.nan, np.nan)
    return pd.DataFrame(rows, index=["patient_adj_HR", "patient_adj_p"]).T


def run(versions, proximity, outdir, n_perm=500, cutoffs=(-1.0, -1.2816, -1.645), top_k=7, n_jobs=4, extra_z=None):
    """versions: {name: Study}; the first one is the primary analysis."""
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(versions)
    primary = versions[names[0]]
    drug, cancer = primary.drug, primary.cancer
    report = {"cancer": cancer, "drug": drug, "model": primary.notes.get("model", "patient-derived organoids"),
              "notes": primary.notes, "versions": {}}

    # ---------------------------------------------------------- every data version
    sigs, tables, evals = {}, {}, {}
    for name, st in versions.items():
        print(f"[{name}] {st.summary()}")
        base = st.baseline()
        base.to_csv(outdir / f"baseline_{name}.tsv", sep="\t", index=False)
        w, tab = st.signature(top_k=top_k)
        wa, _ = st.signature(top_k=top_k, align=True)
        ev, eva = st.evaluate(w), st.evaluate(wa)
        sigs[name], tables[name], evals[name] = w, tab, ev
        tab.to_csv(outdir / f"stability_{name}.tsv", sep="\t")
        tp = st.treated_pats
        score = st.patient_score(w, tp)
        r = score <= np.median(score)
        plots.km_plot(st.clinical.loc[tp, "months"].values, st.clinical.loc[tp, "event"].values, r,
                      f"{cancer} {drug.lower()} | robust signature | {name}", outdir / f"km_robust_{name}.png")
        plots.stability_plot(tab, f"Stability selection | {name}", outdir / f"stability_{name}.png")
        report["versions"][name] = {
            "summary": st.summary(),
            "baseline_best": base.sort_values("logrank_p").iloc[0][["ML", "k", "logrank_p"]].to_dict(),
            "baseline_k7_ridge_p": float(base[(base.ML == "Ridge") & (base.k == 7)].logrank_p.iloc[0]),
            "signature": w.round(4).to_dict(),
            "robust": ev, "robust_precise_aligned": eva,
        }

    # ---------------------------------------------------------- primary: deeper checks
    st = primary
    print("[primary] organoid LOOCV")
    report["organoid_cv"] = st.organoid_cv(top_k=top_k)
    print("[primary] jackknife")
    jk_member, jk_hr = jackknife(st, top_k=top_k)
    print("[primary] patient bootstrap")
    boot = patient_bootstrap(st, sigs[names[0]])
    print("[primary] cut-off sensitivity")
    cut = {}
    for c in cutoffs:
        feats = [p for p in proximity.index if proximity[p] <= c]
        sc = dataclasses.replace(st, features=feats, label=f"z<={c}")  # keeps landmark/covariates/priors
        if len(sc.features) < 3:
            cut[str(c)] = {"n_pathways": len(sc.features), "signature": [], "adj_HR": np.nan, "adj_p": np.nan}
            continue
        wc, _ = sc.signature(top_k=top_k)
        ec = sc.evaluate(wc, interaction=False)
        cut[str(c)] = {"n_pathways": len(sc.features), "signature": list(wc.index), "adj_HR": ec.get("adj_HR"), "adj_p": ec.get("adj_p")}
    print(f"[primary] permutation null ({n_perm})")
    perm = st.permutation_null(n_perm=n_perm, top_k=top_k, n_jobs=n_jobs)
    perm.to_csv(outdir / "permutation_null.tsv", sep="\t", index=False)
    obs = evals[names[0]]
    base0 = pd.read_csv(outdir / f"baseline_{names[0]}.tsv", sep="\t")
    report["permutation"] = {
        "n": n_perm,
        "robust_adjHR_empirical_p": (float((1 + (perm.log_adj_HR >= np.log(obs["adj_HR"])).sum()) / (1 + n_perm))
                                     if np.isfinite(obs.get("adj_HR", np.nan)) else float("nan")),
        "robust_logrank_empirical_p": float((1 + (perm.robust_logrank_p <= obs["logrank_p"]).sum()) / (1 + n_perm)),
        "baseline_min_p_observed": float(base0.logrank_p.min()),
        "baseline_min_p_empirical_p": float((1 + (perm.baseline_min_p <= base0.logrank_p.min()).sum()) / (1 + n_perm)),
    }
    report["fragility"] = {
        "jackknife_HR_gt1_fraction": float(np.nanmean(jk_hr > 1)),
        "jackknife_HR_range": [float(np.nanmin(jk_hr)), float(np.nanmax(jk_hr))],
        "bootstrap_HR_CI95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "bootstrap_HR_gt1_fraction": float(np.mean(boot > 1)),
        "cutoff_sensitivity": cut,
    }

    # ---------------------------------------------------------- pathway scorecard
    card = tables[names[0]].copy()
    card.insert(0, "proximity_z", proximity.reindex(card.index))
    for k, z in (extra_z or {}).items():
        card.insert(1, k, z.reindex(card.index))
    card["jackknife_retention"] = jk_member.reindex(card.index).fillna(0)
    for name in names[1:]:
        card[f"freq_{name}"] = tables[name]["selection_freq"].reindex(card.index)
        card[f"coef_{name}"] = tables[name]["mean_coef"].reindex(card.index)
    card = card.join(pathway_patient_hr(st, list(card.index)))
    # organoid coefficient > 0 means resistance, so the patient HR should be > 1
    card["patient_direction_agrees"] = np.sign(card["mean_coef"]) == np.sign(np.log(card["patient_adj_HR"]))
    version_ok = np.ones(len(card), bool)
    for name in names[1:]:
        version_ok &= (card[f"freq_{name}"] >= 0.4) & (np.sign(card[f"coef_{name}"]) == np.sign(card["mean_coef"]))
    card["consistent_across_versions"] = version_ok
    checks = pd.DataFrame({
        "stable": card.selection_freq >= 0.5,
        "sign": card.sign_agreement >= 0.9,
        "jackknife": card.jackknife_retention >= 0.7,
        "versions": card.consistent_across_versions,
        "patient_direction": card.patient_direction_agrees,
    })
    card["checks_passed"] = checks.sum(axis=1)
    card["tier"] = np.select([checks.all(axis=1), card.checks_passed >= 3], ["robust", "supported"], "fragile")
    card = card.sort_values(["checks_passed", "selection_freq"], ascending=False)
    card.to_csv(outdir / "pathway_scorecard.tsv", sep="\t")
    report["robust_pathways"] = list(card.index[card.tier == "robust"])
    report["supported_pathways"] = list(card.index[card.tier == "supported"])

    # ---------------------------------------------------------- signature verdict
    p_perm = report["permutation"]["robust_adjHR_empirical_p"]
    hr_dirs = [evals[n].get("adj_HR", np.nan) > 1 for n in names]
    inter = obs.get("interaction_p", np.nan)
    report["verdict"] = {
        "permutation_significant": p_perm < 0.05,
        "direction_consistent_all_versions": bool(all(hr_dirs)),
        "predictive_interaction_p_lt_0.1": bool(inter < 0.1),
        "events_in_treated": int(obs["events"]),
        "adequately_powered": int(obs["events"]) >= 30,
    }
    hrs = np.array([evals[n].get("adj_HR", np.nan) for n in names], float)
    report["verdict"]["direction"] = ("as expected (resistant score -> worse survival)" if np.all(hrs > 1) else
                                      "reversed (resistant score -> better survival)" if np.all(hrs < 1) else "inconsistent")
    v = report["verdict"]
    v["signature_status"] = ("validated" if v["permutation_significant"] and v["direction_consistent_all_versions"]
                             and v["predictive_interaction_p_lt_0.1"] else
                             "promising, not validated" if v["direction_consistent_all_versions"] else "not supported")

    plots.forest_plot([(n, evals[n].get("adj_HR", np.nan), evals[n].get("adj_HR_low", np.nan),
                        evals[n].get("adj_HR_high", np.nan)) for n in names],
                      f"{cancer} {drug.lower()}: signature across data versions", outdir / "forest_versions.png")
    with open(outdir / "report.json", "w") as f:
        json.dump(report, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    write_markdown(report, card, outdir)
    return report


def _fmt(x, d=3):
    try:
        return f"{float(x):.{d}g}"
    except (TypeError, ValueError):
        return str(x)


def write_markdown(rep, card, outdir):
    L = [f"# {rep['cancer']} / {rep['drug']}: robust biomarker report", "",
         f"Pre-clinical model: **{rep['model']}**. Drug targets: {', '.join(rep['notes'].get('targets', [])) or 'see drug_drugTarget.txt'}. "
         f"Network: {rep['notes'].get('network', 'STRING > 700 (original precomputed proximity)')}.", ""]
    v = rep["verdict"]
    L += [f"**Signature status: {v['signature_status']}** "
          f"(permutation p = {_fmt(rep['permutation']['robust_adjHR_empirical_p'])}, "
          f"direction: {v['direction']}, "
          f"treatment-interaction p = {_fmt(rep['versions'][next(iter(rep['versions']))]['robust'].get('interaction_p'))}, "
          f"events in treated patients = {v['events_in_treated']}"
          f"{'' if v['adequately_powered'] else ' - underpowered, fewer than 30 events'}).", ""]
    L += ["## Data versions", "", "| version | organoids | pathways | treated pts | events | paper-style best p | robust log-rank p | adj HR (95% CI) | adj p | C-index | interaction p | PRECISE adj HR |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, r in rep["versions"].items():
        s, e, a = r["summary"], r["robust"], r["robust_precise_aligned"]
        L.append(f"| {n} | {s['organoids']} | {s['candidate_pathways']} | {s['treated_patients']} | {e['events']} | "
                 f"{_fmt(r['baseline_best']['logrank_p'])} ({r['baseline_best']['ML']}, k={r['baseline_best']['k']}) | "
                 f"{_fmt(e['logrank_p'])} | {_fmt(e.get('adj_HR'))} ({_fmt(e.get('adj_HR_low'))}-{_fmt(e.get('adj_HR_high'))}) | "
                 f"{_fmt(e.get('adj_p'))} | {_fmt(e['c_index'])} | {_fmt(e.get('interaction_p'))} | {_fmt(a.get('adj_HR'))} |")
    p = rep["permutation"]
    L += ["", "## Honest significance (whole-pipeline permutation of organoid IC50)", "",
          f"- Robust signature, adjusted Cox HR: empirical p = {_fmt(p['robust_adjHR_empirical_p'])} ({p['n']} permutations)",
          f"- Original approach (best of 3 models x k = 2..10): observed min p = {_fmt(p['baseline_min_p_observed'])}, "
          f"empirical p = {_fmt(p['baseline_min_p_empirical_p'])}", ""]
    cv = rep["organoid_cv"]
    f = rep["fragility"]
    L += ["## Organoid cross-validation (leave one organoid out, selection inside folds)", "",
          f"- Robust learner: Spearman rho = {_fmt(cv['robust_loocv_spearman'])} (p = {_fmt(cv['robust_loocv_p'])})",
          f"- Original Ridge top-7: Spearman rho = {_fmt(cv['baseline_loocv_spearman'])} (p = {_fmt(cv['baseline_loocv_p'])})", "",
          "## Fragility", "",
          f"- Leave-one-organoid-out: HR > 1 in {_fmt(100 * f['jackknife_HR_gt1_fraction'])}% of refits (range {_fmt(f['jackknife_HR_range'][0])}-{_fmt(f['jackknife_HR_range'][1])})",
          f"- Patient bootstrap: adjusted HR 95% interval {_fmt(f['bootstrap_HR_CI95'][0])}-{_fmt(f['bootstrap_HR_CI95'][1])}, HR > 1 in {_fmt(100 * f['bootstrap_HR_gt1_fraction'])}%",
          "- Proximity cut-off sensitivity:"]
    for c, r in f["cutoff_sensitivity"].items():
        L.append(f"  - z <= {c}: {r['n_pathways']} pathways, adj HR {_fmt(r['adj_HR'])} (p {_fmt(r['adj_p'])}), signature: {', '.join(x.replace('REACTOME_', '') for x in r['signature'])}")
    L += ["", "## Pathway scorecard (top 15)", "",
          "| pathway | tier | checks | proximity z | selection freq | sign agree | jackknife | organoid coef | patient adj HR |",
          "|---|---|---|---|---|---|---|---|---|"]
    for pw, r in card.head(15).iterrows():
        L.append(f"| {pw.replace('REACTOME_', '')} | {r.tier} | {r.checks_passed}/5 | {_fmt(r.proximity_z)} | {_fmt(r.selection_freq)} | "
                 f"{_fmt(r.sign_agreement)} | {_fmt(r.jackknife_retention)} | {_fmt(r.mean_coef)} | {_fmt(r.patient_adj_HR)} |")
    L += ["", "Checks: stable (selection freq >= 0.5), sign agreement >= 0.9, jackknife retention >= 0.7, "
          "consistent in every data version, patient HR direction matches organoid coefficient. "
          "Tier robust = 5/5, supported = 3-4/5.", "",
          "Figures: `km_robust_*.png`, `stability_*.png`, `forest_versions.png`."]
    (outdir / "REPORT.md").write_text("\n".join(L) + "\n")
