"""ANCHOR Stage 0: triangulate the three label-free sources and FREEZE <= 5 signed pathways per drug.

No patient outcome data are read anywhere in this agent. Rules (fixed before looking at the result table):
  Source 1 (expression association) passes: BH q < 0.05 for the meta z AND same sign in >= min(2, screens) screens.
  Source 2 (CRISPR co-dependency) passes:   BH q < 0.05 for the pathway co-dependency z.
  Source 3 (network proximity) passes:      STRING v12 proximity z <= -1.2816 (paper cut-off, one-sided 10%).
  Signed:   direction = sign of the Source-1 meta z; needs nominal Source-1 p < 0.05.
  Specific: percentile of the drug's directional z among all drugs in the same screens >= 0.75 AND larger than
            the mean directional z of related drugs (same mechanism class).
  Eligible: >= 2 sources passed AND signed AND specific.
  Ranking:  sources passed (desc), then |meta z| (desc); a pathway sharing > 50% genes (Jaccard) with one already
            chosen for the drug is skipped; keep at most 5.

python -I nar/08_depmap/code/06_freeze.py
"""
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import OI, TAB, panel, savefig, style  # noqa: E402
from obd.reference import gene_sets  # noqa: E402

CUT = -1.2816
DRUGS = ["FLUOROURACIL", "OXALIPLATIN", "IRINOTECAN", "CISPLATIN", "CARBOPLATIN", "GEMCITABINE", "PACLITAXEL",
         "DOXORUBICIN", "PEMETREXED", "SORAFENIB"]


def jacc(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(1, len(a | b))


def main():
    E = pd.read_csv(TAB / "B_source1_expression_association.tsv.gz", sep="\t")
    C = pd.read_csv(TAB / "B_source2_crispr_codependency.tsv.gz", sep="\t")
    Z = pd.read_csv(TAB / "B_network_proximity_z.tsv", sep="\t", index_col=0)
    gs = gene_sets(("REACTOME",))
    D = E.merge(C, on=["drug", "pathway"], how="left")
    D["proximity_z_v12"] = [Z.at[p, d] if p in Z.index else np.nan for p, d in zip(D.pathway, D.drug)]
    D["proximity_z_paper"] = [Z.at[p, f"{d}__paper_table"] if p in Z.index else np.nan for p, d in zip(D.pathway, D.drug)]
    need = np.minimum(2, D["screens_total"])
    D["S1_expression"] = (D.expr_q < 0.05) & (D.screens_same_sign >= need)
    D["S2_crispr"] = D.crispr_q < 0.05
    D["S3_network"] = D.proximity_z_v12 <= CUT
    D["n_sources"] = D[["S1_expression", "S2_crispr", "S3_network"]].sum(1)
    D["signed"] = D.expr_p < 0.05
    D["specific"] = (D.pct_vs_all_drugs >= 0.75) & ((D.own_dir_z_mean > D.related_mean_dir_z) | D.related_mean_dir_z.isna())
    D["eligible"] = (D.n_sources >= 2) & D.signed & D.specific
    D["direction"] = np.where(D.expr_meta_z > 0, "+ (higher score -> more RESISTANT)", "- (higher score -> more SENSITIVE)")
    D["crispr_direction"] = np.where(D.crispr_z > 0, "dependency on pathway genes -> more sensitive",
                                     np.where(D.crispr_z < 0, "dependency on pathway genes -> more resistant", ""))
    cand = D[D.n_sources >= 2].sort_values(["drug", "n_sources", "expr_meta_z"], ascending=[True, False, True])
    cand.to_csv(TAB / "B_candidates_two_or_more_sources.tsv", sep="\t", index=False)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    frozen, summary = [], []
    for d in DRUGS:
        g = D[(D.drug == d) & D.eligible].copy()
        g["absz"] = g.expr_meta_z.abs()
        g = g.sort_values(["n_sources", "absz"], ascending=[False, False])
        chosen = []
        for r in g.itertuples():
            if any(jacc(gs[r.pathway], gs[c.pathway]) > 0.5 for c in chosen):
                continue
            chosen.append(r)
            if len(chosen) == 5:
                break
        for rank, r in enumerate(chosen, 1):
            frozen.append({"frozen_at_utc": stamp, "drug": d, "rank": rank, "pathway": r.pathway,
                           "sign": "+" if r.expr_meta_z > 0 else "-", "direction": r.direction,
                           "n_sources": int(r.n_sources),
                           "sources": ";".join(k for k, v in (("expression", r.S1_expression), ("crispr", r.S2_crispr),
                                                              ("network", r.S3_network)) if v),
                           "expr_meta_z": round(r.expr_meta_z, 3), "expr_q": float(f"{r.expr_q:.3g}"),
                           "screens_same_sign": f"{int(r.screens_same_sign)}/{int(r.screens_total)}",
                           "crispr_z": round(r.crispr_z, 3) if pd.notna(r.crispr_z) else np.nan,
                           "crispr_q": float(f"{r.crispr_q:.3g}") if pd.notna(r.crispr_q) else np.nan,
                           "proximity_z_STRINGv12": round(r.proximity_z_v12, 3) if pd.notna(r.proximity_z_v12) else np.nan,
                           "pct_vs_all_drugs": round(r.pct_vs_all_drugs, 3),
                           "own_dir_z": round(r.own_dir_z_mean, 3), "related_mean_dir_z": round(r.related_mean_dir_z, 3)
                           if pd.notna(r.related_mean_dir_z) else np.nan,
                           "related_drugs": r.related_drugs, "screens_used": r.screens_used,
                           "patient_outcomes_used": "NO"})
        dd = D[D.drug == d]
        summary.append({"drug": d, "screens": dd.screens_used.iloc[0], "pathways_tested": len(dd),
                        "S1_pass": int(dd.S1_expression.sum()), "S2_pass": int(dd.S2_crispr.sum()),
                        "S3_pass": int(dd.S3_network.sum()), "two_plus_sources": int((dd.n_sources >= 2).sum()),
                        "two_plus_signed_specific": int(dd.eligible.sum()), "frozen": len(chosen)})
    F = pd.DataFrame(frozen)
    out = TAB / "frozen_hypotheses.tsv"
    if out.exists():
        print("frozen_hypotheses.tsv already exists - pre-registration is not overwritten; writing frozen_hypotheses.rerun.tsv")
        out = TAB / "frozen_hypotheses.rerun.tsv"
    F.to_csv(out, sep="\t", index=False)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    (TAB / (out.stem + ".sha256")).write_text(f"{sha}  {out.name}\nfrozen_at_utc {stamp}\n")
    S = pd.DataFrame(summary)
    S.to_csv(TAB / "B_stage0_summary.tsv", sep="\t", index=False)
    print(S.to_string())
    print(F[["drug", "rank", "pathway", "sign", "n_sources", "sources", "expr_meta_z", "crispr_z", "proximity_z_STRINGv12"]].to_string())
    heatmap(F, D, S)


def heatmap(F, D, S):
    plt = style()
    pws = list(dict.fromkeys(F.sort_values(["drug", "rank"]).pathway))
    drugs = DRUGS
    M = np.full((len(drugs), len(pws)), np.nan)
    N = np.zeros_like(M)
    for i, d in enumerate(drugs):
        for j, p in enumerate(pws):
            r = D[(D.drug == d) & (D.pathway == p)]
            if r.empty:
                continue
            r = r.iloc[0]
            N[i, j] = r.n_sources
            frozen = ((F.drug == d) & (F.pathway == p)).any()
            M[i, j] = np.sign(r.expr_meta_z) * (r.n_sources if frozen else 0.35 * r.n_sources)
    fig, axes = plt.subplots(1, 2, figsize=(max(9, 0.33 * len(pws) + 4), 7.2), gridspec_kw={"width_ratios": [len(pws) or 1, 3]})
    ax = axes[0]
    panel(ax, "a")
    if pws:
        from matplotlib.colors import ListedColormap, BoundaryNorm
        cmap = ListedColormap(["#08519c", "#6baed6", "#c6dbef", "#f7f7f7", "#fdd0a2", "#fd8d3c", "#a63603"])
        norm = BoundaryNorm([-3.5, -2.5, -1.5, -0.05, 0.05, 1.5, 2.5, 3.5], cmap.N)
        ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
        for i in range(len(drugs)):
            for j in range(len(pws)):
                if ((F.drug == drugs[i]) & (F.pathway == pws[j])).any():
                    ax.text(j, i, f"{int(N[i, j])}", ha="center", va="center", fontsize=7, fontweight="bold",
                            color="white" if N[i, j] >= 3 else "black")
        ax.set_xticks(range(len(pws)))
        ax.set_xticklabels([p.replace("REACTOME_", "").replace("_", " ").lower()[:45] for p in pws], rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([d.title() for d in drugs])
    ax.set_title("Frozen Stage-0 hypotheses (number = evidence sources; red = higher score -> resistant,\n"
                 "blue = higher score -> sensitive; pale = same pathway, not frozen for that drug)", fontsize=8.5)
    ax = axes[1]
    panel(ax, "b")
    y = np.arange(len(S))
    for k, (col, c) in enumerate((("S1_pass", OI["blue"]), ("S2_pass", OI["orange"]), ("S3_pass", OI["green"]),
                                  ("two_plus_sources", OI["grey"]), ("frozen", OI["red"]))):
        ax.barh(y + k * 0.16 - 0.32, S[col], height=0.16, color=c, label=col.replace("_", " "))
    ax.set_yticks(y)
    ax.set_yticklabels([d.title() for d in S.drug])
    ax.invert_yaxis()
    ax.set_xscale("symlog")
    ax.set_xlabel("Number of pathways")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.set_title("Evidence funnel per drug", fontsize=9)
    fig.tight_layout()
    savefig(fig, "fig5_anchor_stage0_frozen_hypotheses", """
**Figure 5. ANCHOR Stage 0: pre-registered, outcome-free mechanism hypotheses.**
(a) Each row is a drug; each column a Reactome pathway that was frozen for at least one drug. Coloured cells are
frozen hypotheses: red = a HIGHER pathway score is expected to go with RESISTANCE, blue = with SENSITIVITY; the
number is how many of the three independent label-free sources supported it (1: expression association across
GDSC1/GDSC2, CTRPv2 and PRISM after removing cancer type and general drug sensitivity; 2: CRISPR co-dependency in
Project Score; 3: STRING network proximity to the drug's targets). Pale cells show the same pathway in another drug
where it was not frozen. Every frozen hypothesis also passed a specificity rule against all drugs and against drugs
of the same mechanism class. No patient data were used. (b) How many pathways passed each source and each step.
Numbers: `tables/frozen_hypotheses.tsv`, `tables/B_stage0_summary.tsv`, `tables/B_candidates_two_or_more_sources.tsv`.
""")
    plt.close(fig)


if __name__ == "__main__":
    main()
