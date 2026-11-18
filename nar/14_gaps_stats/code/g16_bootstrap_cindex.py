"""G16 part 1+3: paired bootstrap 95% CI of the held-out C-index difference between models, per unit.

Uses the held-out (cross-fitted) patient scores saved by nar/03_organoid_plus_patient (cache/<unit>__shortlist__scores)
and by g16_joint_extended_grid.py (extended lab-weight grid). Same cohorts and inclusion rules as the original
evaluation (treated patients; survival cohorts with >= 10 patients and >= 5 events -> Harrell's C; response cohorts with
>= 10 patients and >= 3 per class -> AUROC). Pooled C = concordant pairs / comparable pairs summed within cohorts
(pairs never cross cohorts), exactly as the original 'POOLED' row. Bootstrap: patients resampled with replacement
inside each cohort, the same resample for all models (paired), 2,000 replicates.
"""
import numpy as np
import pandas as pd

from g_common import DATA, NAR, OK, SEED, TAB, C, save_fig

J = NAR / "03_organoid_plus_patient"
B = 2000


def pair_mats(t, e, s_list):
    t = np.asarray(t, float); e = np.asarray(e, int)
    P = (e[:, None] == 1) & (t[None, :] > t[:, None])
    As = []
    for s in s_list:
        s = np.asarray(s, float)
        As.append(P * ((s[:, None] < s[None, :]) + 0.5 * (s[:, None] == s[None, :])))
    return P.astype(float), As


def auc_mats(y, s_list):
    y = np.asarray(y, int)
    P = (y[:, None] == 1) & (y[None, :] == 0)
    As = []
    for s in s_list:
        s = np.asarray(s, float)
        As.append(P * ((s[:, None] > s[None, :]) + 0.5 * (s[:, None] == s[None, :])))
    return P.astype(float), As


def run(scores, models, tag, rng):
    rows = []
    for kind in ("survival", "response"):
        mats = []
        for c, x in scores[scores.group == "treated"].groupby("unit_cohort"):
            if kind == "survival":
                x = x[x.has_surv.astype(bool)]
                if len(x) < 10 or x.event.sum() < 5:
                    continue
                x = x.dropna(subset=models)
                P, As = pair_mats(x.time, x.event, [x[m] for m in models])
            else:
                x = x[x.has_resp.astype(bool)]
                if len(x) < 10 or x.responder.sum() < 3 or (1 - x.responder).sum() < 3:
                    continue
                x = x.dropna(subset=models)
                P, As = auc_mats(x.responder, [x[m] for m in models])
            mats.append((c, len(x), P, As))
        if not mats:
            continue
        def stat(ws):
            den = sum(w @ P @ w for (c, n, P, As), w in zip(mats, ws))
            return np.array([sum(w @ A[k] @ w for (c, n, P, A), w in zip(mats, ws)) / den for k in range(len(models))])
        obs = stat([np.ones(n) for c, n, P, A in mats])
        boot = np.array([stat([rng.multinomial(n, np.ones(n) / n).astype(float) for c, n, P, A in mats]) for _ in range(B)])
        base = {"unit": tag, "outcome": kind, "n_cohorts": len(mats), "n_patients": sum(n for c, n, P, A in mats),
                "cohorts": ";".join(c for c, n, P, A in mats)}
        for k, m in enumerate(models):
            base[f"C_{m}"] = obs[k]
        for a, b in (("joint", "patient-only"), ("joint", "organoid-only"), ("patient-only", "organoid-only"),
                     ("joint (w=0.1)", "patient-only"), ("joint (w=0.5)", "patient-only"), ("joint", "joint (w=0.5)")):
            if a in models and b in models:
                ia, ib = models.index(a), models.index(b)
                d = boot[:, ia] - boot[:, ib]
                rows.append({**base, "comparison": f"{a} minus {b}", "diff": obs[ia] - obs[ib],
                             "lo95": np.percentile(d, 2.5), "hi95": np.percentile(d, 97.5),
                             "share_boot_diff_le_0": float(np.mean(d <= 0))})
    return rows


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for f in sorted((J / "cache").glob("*__shortlist__scores.tsv.gz")):
        unit = f.name.split("__")[0]
        S = pd.read_csv(f, sep="\t", index_col=0)
        rows += [dict(r, grid="original (w 0.25/0.5/0.75)") for r in run(S, ["organoid-only", "patient-only", "joint", "joint (w=0.5)"], unit, rng)]
        g = DATA / "g16_extended" / f"{unit}__shortlist_extgrid__scores.tsv.gz"
        if g.exists():
            S2 = pd.read_csv(g, sep="\t", index_col=0)
            rows += [dict(r, grid="extended (w 0/0.05/0.1/0.25/0.5/0.75)") for r in run(S2, ["organoid-only", "patient-only", "joint", "joint (w=0.1)"], unit, rng)]
    R = pd.DataFrame(rows)
    R["CI_above_0"] = R.lo95 > 0
    R["CI_below_0"] = R.hi95 < 0
    R.to_csv(TAB / "G16_paired_bootstrap_cindex_diff.tsv", sep="\t", index=False)
    # check against the original POOLED rows
    chk = []
    for f in sorted((J / "tables").glob("*__shortlist__eval.tsv")):
        e = pd.read_csv(f, sep="\t")
        e = e[e.cohort == "POOLED"]
        for _, r in e.iterrows():
            chk.append({"unit": f.name.split("__")[0], "outcome": r.outcome, "model": r.model, "C_original_table": r.C})
    ck = pd.DataFrame(chk)
    o = R[R.grid.str.startswith("original")].drop_duplicates(["unit", "outcome"])
    m = []
    for _, r in ck.iterrows():
        x = o[(o.unit == r.unit) & (o.outcome == r.outcome)]
        col = f"C_{r.model}"
        m.append(x[col].iloc[0] if len(x) and col in x else np.nan)
    ck["C_recomputed"] = m
    ck["abs_diff"] = (ck.C_original_table - ck.C_recomputed).abs()
    ck.to_csv(TAB / "G16_check_pooled_C_vs_original.tsv", sep="\t", index=False)
    print("max abs diff vs original pooled C:", ck.abs_diff.max())
    # chosen weights in the extended grid
    wrows = []
    for f in sorted((DATA / "g16_extended").glob("*__folds.tsv")):
        fo = pd.read_csv(f, sep="\t")
        j = fo[fo.model == "joint"]
        wrows.append({"unit": f.name.split("__")[0], "n_fits": len(j),
                      **{f"chosen_w_{w}": int((j.w.round(2) == w).sum()) for w in (0.0, 0.05, 0.1, 0.25, 0.5, 0.75)}})
    W = pd.DataFrame(wrows)
    if len(W):
        W.loc[len(W)] = ["ALL"] + list(W.iloc[:, 1:].sum())
    W.to_csv(TAB / "G16_extended_grid_chosen_weight.tsv", sep="\t", index=False)
    print(W.to_string())
    print(R[["unit", "grid", "outcome", "comparison", "diff", "lo95", "hi95"]].to_string())
    figure(R)


def figure(R):
    plt = C.plot_style()
    fig, axs = plt.subplots(1, 3, figsize=(14, 5.5), sharey=True)
    o = R[R.grid.str.startswith("original")]
    e = R[R.grid.str.startswith("extended")]
    panels = [(o, "joint minus patient-only", "a  joint - patient-only (original grid)"),
              (o, "joint minus organoid-only", "b  joint - organoid-only (original grid)"),
              (e, "joint minus patient-only", "c  joint - patient-only (extended grid, w down to 0)")]
    labs = sorted(set(zip(R.unit, R.outcome)))
    ypos = {k: i for i, k in enumerate(labs)}
    for ax, (d, comp, title) in zip(axs, panels):
        d = d[d.comparison == comp]
        for _, r in d.iterrows():
            y = ypos[(r.unit, r.outcome)]
            col = OK["green"] if r.lo95 > 0 else (OK["verm"] if r.hi95 < 0 else OK["blue"])
            ax.plot([r.lo95, r.hi95], [y, y], color=col, lw=1.5)
            ax.plot(r["diff"], y, "o", color=col, ms=5)
        ax.axvline(0, color="k", lw=.7, ls="--")
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
        ax.set_xlabel("difference in held-out C-index / AUROC\n(green: CI above 0; red: CI below 0)")
    axs[0].set_yticks(range(len(labs)))
    axs[0].set_yticklabels([f"{u} ({k})" for u, k in labs], fontsize=7)
    axs[0].invert_yaxis()
    fig.tight_layout()
    cap = """**Figure G16. Does adding organoid data to a patient model improve prediction in unseen patients?**
Each row is one cancer x drug unit and outcome type (survival: C-index; response: AUROC; both 0.5 = chance). Dot = the
difference between two models in held-out patients; line = 95% paired bootstrap interval (patients resampled within
each cohort, the same resample for both models, 2,000 times). (a) joint organoid+patient model minus patient-only
model; (b) joint model minus organoid-only model; (c) as (a) after letting the inner cross-validation choose the lab
weight from a wider grid (w = 0, 0.05, 0.1, 0.25, 0.5, 0.75; w = weight of the lab-model loss, 0 = patients only).
A line that crosses the dashed zero line means the data cannot tell the two models apart.
"""
    save_fig(fig, "G16_joint_model_bootstrap", cap)


if __name__ == "__main__":
    main()
