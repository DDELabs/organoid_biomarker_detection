"""Markdown tables for nar/03 REPORT.md in the harmonised convention (audit I05: HR > 1 = expected direction).

Reads nar/03_organoid_plus_patient/tables/summary_main.tsv and interaction_all.tsv (written by make_figures.py,
which adds HR_risk / ratio_risk columns) and the fold-aligned sensitivity runs (I06) if present.
Writes nar/15_fixes/tables/03_pooled_table.md, 03_interaction_table.md, 03_foldalign_comparison.tsv/.md.
"""
from pathlib import Path
import numpy as np
import pandas as pd

F = Path(__file__).resolve().parents[2]
T3 = F / "03_organoid_plus_patient" / "tables"
OUT = F / "15_fixes" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def fp(p):
    return "<0.001" if p < 0.001 else (f"{p:.2g}" if p < 0.01 else f"{p:.2f}" if p >= 0.1 else f"{p:.3f}")


M = pd.read_csv(T3 / "summary_main.tsv", sep="\t")
M = M[M.features == "shortlist"]
lines = ["| Unit | Model | Survival n / events | C-index | HR per SD of risk score (95% CI), p | Response: n / responders, AUROC |",
         "|---|---|---|---|---|---|"]
units = list(dict.fromkeys(M.unit))
for u in units:
    first = True
    for m in ("organoid-only", "patient-only", "joint"):
        s = M[(M.unit == u) & (M.model == m) & (M.outcome == "survival")]
        r = M[(M.unit == u) & (M.model == m) & (M.outcome == "response")]
        lab = m if not (u.startswith(("OV_PLATINUM", "LUNG")) and m == "organoid-only") else "lab-only (cell lines)"
        if len(s):
            s = s.iloc[0]
            sv = (f"{int(s.n)} / {int(s.events)}" if first else "", f"{s.C:.3f}",
                  f"{s.HR_risk:.2f} ({s.HR_risk_low:.2f}-{s.HR_risk_high:.2f}), {fp(s.p)}")
        else:
            sv = ("no public survival times" if first else "", "", "")
        rv = ""
        if len(r):
            r = r.iloc[0]
            rv = f"{int(r.n)} / {int(r.events)}, {r.C:.3f}" if first else f"{r.C:.3f}"
        elif first:
            rv = "—"
        lines.append(f"| {u if first else ''} | {lab} | {sv[0]} | {sv[1]} | {sv[2]} | {rv} |")
        first = False
(OUT / "03_pooled_table.md").write_text("\n".join(lines) + "\n")

I = pd.read_csv(T3 / "interaction_all.tsv", sep="\t")
I = I[(I.features == "shortlist") & (I.model.isin(["joint", "organoid-only"]))]
rows = ["| Cohort | Unit | Outcome | Treated / comparison (events or responders) | Joint interaction (95% CI), p | Organoid/lab-only interaction, p |",
        "|---|---|---|---|---|---|"]
for (u, c, o), g in I.groupby(["unit", "cohort", "outcome"], sort=False):
    j = g[g.model == "joint"].iloc[0]
    k = g[g.model == "organoid-only"]
    meas = "HR" if j.measure == "interaction HR" else "OR"
    kk = f"{k.iloc[0].ratio_risk:.2f}, {fp(k.iloc[0].p)}" if len(k) else ""
    rows.append(f"| {c} | {u} | {o} | {int(j.n_treated)} / {int(j.n_comparison)} ({int(j.events_treated)} / {int(j.events_comparison)}) | "
                f"{meas} {j.ratio_risk:.2f} ({j.low_risk:.2f}-{j.high_risk:.2f}), {fp(j.p)} | {kk} |")
(OUT / "03_interaction_table.md").write_text("\n".join(rows) + "\n")

# I06: fold-aligned vs transductive alignment
fa = sorted(T3.glob("*__shortlist_foldalign__eval.tsv"))
if fa:
    rr = []
    for f in fa:
        u = f.name.split("__")[0]
        a = pd.read_csv(f, sep="\t")
        b = pd.read_csv(T3 / f"{u}__shortlist__eval.tsv", sep="\t")
        for m in ("organoid-only", "patient-only", "joint"):
            for kind in ("survival", "response"):
                x = a[(a.cohort == "POOLED") & (a.model == m) & (a.outcome == kind)]
                y = b[(b.cohort == "POOLED") & (b.model == m) & (b.outcome == kind)]
                if not len(x) or not len(y):
                    continue
                x, y = x.iloc[0], y.iloc[0]
                d = dict(unit=u, model=m, outcome=kind, n=int(y.n), C_transductive=y.C, C_fold_aligned=x.C, C_diff=x.C - y.C)
                if kind == "survival":
                    d.update(HR_risk_transductive=1 / y.HR, HR_risk_fold_aligned=1 / x.HR, p_transductive=y.p, p_fold_aligned=x.p)
                rr.append(d)
    R = pd.DataFrame(rr)
    R.to_csv(OUT / "03_foldalign_comparison.tsv", sep="\t", index=False)
    L = ["| Unit | Model | Outcome | C / AUROC transductive | C / AUROC fold-aligned | HR (risk) transductive, p | HR (risk) fold-aligned, p |",
         "|---|---|---|---|---|---|---|"]
    for r in R.itertuples():
        h1 = f"{r.HR_risk_transductive:.2f}, {fp(r.p_transductive)}" if r.outcome == "survival" else ""
        h2 = f"{r.HR_risk_fold_aligned:.2f}, {fp(r.p_fold_aligned)}" if r.outcome == "survival" else ""
        L.append(f"| {r.unit} | {r.model} | {r.outcome} | {r.C_transductive:.3f} | {r.C_fold_aligned:.3f} | {h1} | {h2} |")
    (OUT / "03_foldalign_comparison.md").write_text("\n".join(L) + "\n")
    print(R.round(3).to_string())
print("\n".join(lines)); print("\n".join(rows))
