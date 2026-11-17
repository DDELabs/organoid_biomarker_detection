"""Cross-set summary table and forest plot for every (set, drug) run by run_blca.py and run_sets.py."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import common as K  # noqa: E402
import engine as E  # noqa: E402

DIRS = ["BLCA_lee2018/paper", "BLCA_lee2018/current", "OV_vias2023", "SARC_alshihabi2024", "LUAD_GSE165611_substitute",
        "PDX_novartis_gao2015"]
OUT = K.BASE / "SUMMARY"


def main():
    rows = []
    for d in DIRS:
        for f in sorted((K.PARENT / d / "tables").glob("*_summary.json")):
            r = json.load(open(f))
            r["folder"] = d
            rows.append(r)
    S = pd.DataFrame(rows)
    keep = ["folder", "set", "set_label", "project", "drug", "model_drug", "n_models", "n_proximal", "biomarker",
            "same_top_all_models", "treated_n", "treated_events", "treated_p", "treated_hr", "treated_hr_lo", "treated_hr_hi",
            "treated_surv5_resp", "treated_surv5_nonresp", "untreated_n", "untreated_events", "untreated_p", "untreated_hr",
            "untreated_hr_lo", "untreated_hr_hi", "loo_Ridge_spearman", "loo_Ridge_spearman_p", "multik_min_p",
            "multik_min_setting", "multik_n_p_below_0.05", "ridge_k2_p", "perm_empirical_p", "boot_observed_rank",
            "boot_empirical_p_rank", "boot_empirical_p_rank1p", "cut-2.1925_p", "cut-1.1824_p", "proximity_source", "note", "skipped"]
    S = S[[c for c in keep if c in S.columns]]
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    S.to_csv(OUT / "tables" / "cross_set_summary.tsv", sep="\t", index=False)
    S = S[S.treated_p.notna()].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 0.36 * len(S) + 1.8))
    lab = []
    for i, r in S.iterrows():
        ax.errorbar(r.treated_hr, i - 0.15, xerr=[[r.treated_hr - r.treated_hr_lo], [r.treated_hr_hi - r.treated_hr]], fmt="o",
                    color=E.OK["verm"], capsize=2, ms=4)
        if pd.notna(r.get("untreated_hr")) and np.isfinite(r.untreated_hr_hi):
            ax.errorbar(r.untreated_hr, i + 0.15, xerr=[[r.untreated_hr - r.untreated_hr_lo], [r.untreated_hr_hi - r.untreated_hr]],
                        fmt="s", color=E.OK["blue"], mfc="white", capsize=2, ms=4)
        lab.append(f"{r.set.replace('_', ' ')} | {r.project} {r.drug.lower()} (n={int(r.treated_n)}, p={r.treated_p:.2g}, perm p={r.perm_empirical_p:.2g})")
    ax.set_yticks(range(len(S)))
    ax.set_yticklabels(lab, fontsize=7.5)
    ax.axvline(1, color="k", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.08, 15)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("Hazard ratio, predicted non-responders vs responders (95% CI, log scale)")
    ax.invert_yaxis()
    ax.plot([], [], "o", color=E.OK["verm"], label="patients given the drug")
    ax.plot([], [], "s", color=E.OK["blue"], mfc="white", label="patients with no drug record")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.set_title("Organoid-only strategy (Kong et al. 2020), rank-1 Ridge pathway, every set and drug")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "cross_set_forest.png", dpi=300)
    fig.savefig(OUT / "figures" / "cross_set_forest.pdf")
    n05 = int((S.treated_p < 0.05).sum())
    nperm = int((S.perm_empirical_p < 0.05).sum())
    nright = int(((S.treated_p < 0.05) & (S.treated_hr > 1)).sum())
    open(OUT / "figures" / "cross_set_forest.caption.md", "w").write(
        "# All sets and drugs in one plot\n\nEach row is one lab-model set and one drug. The paper's procedure picked one pathway "
        "from the lab models; TCGA patients were split at the median of that pathway. Red circles: patients given the drug; "
        "blue squares: patients with no drug record. HR > 1 means predicted non-responders died sooner (the expected direction). "
        f"Rows: {len(S)}; log-rank p < 0.05 in treated patients: {n05} (of which {nright} in the expected direction); "
        f"permutation empirical p < 0.05: {nperm}. Numbers in `SUMMARY/tables/cross_set_summary.tsv`.\n")
    print(S[["set", "project", "drug", "n_models", "treated_n", "treated_p", "treated_hr", "untreated_p", "perm_empirical_p"]].to_string())


if __name__ == "__main__":
    main()
