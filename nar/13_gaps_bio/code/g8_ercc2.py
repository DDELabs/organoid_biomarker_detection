"""G8 part 1: ERCC2-mutant bladder tumours and the cisplatin organoid model -- random-pathway nulls and TMB adjustment.

Reuses nar/02_organoid_only/BLCA_OV_LUAD_OTHER/code/blca_extras.setup() (paper source data and current data).
Score = Ridge coefficient x scaled pathway score (higher = predicted more cisplatin-resistant); ERCC2 mutants are
expected to score LOWER (one-sided Mann-Whitney 'less'), as in the paper's Fig. 5b.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from gcommon import O2, TAB, version_rows
from g8_nulls import SILENT, adjust, all_pathway_null, mwu, pipeline_null, rz, tmb

sys.path.insert(0, str(O2 / "BLCA_OV_LUAD_OTHER" / "code"))
import blca_extras as BX  # noqa: E402
import common as K  # noqa: E402

AA = BX.AA


def main():
    mutf = K.RAW / "TCGA-BLCA.somaticmutation_wxs.tsv.gz"
    T, m = tmb(mutf)
    prof = set(T.index)
    mut = set(m.loc[(m.gene == "ERCC2") & ~m.effect.isin(SILENT), "patient"])
    cd = K.cdr()
    st = cd.reindex(T.index)["stage"].astype(str).str.upper().str.replace("STAGE", "").str.strip()
    hi_stage = pd.Series(np.where(st.str.startswith(("III", "IV")), 1.0, np.where(st.str.startswith("I"), 0.0, np.nan)),
                         index=T.index)
    cov_t = pd.DataFrame({"logTMB_rz": rz(np.log1p(T.values))}, index=T.index)
    cov_ts = cov_t.join(pd.DataFrame({"stage_III_IV": hi_stage}))
    rows, nulls = [], {}
    for v in ("paper", "current"):
        d = BX.setup(v)
        Pz = pd.DataFrame(StandardScaler().fit_transform(d["tsc"].loc[d["uni"]].T), index=d["tsc"].columns, columns=d["uni"])
        c = d["c"]
        r1 = c.abs().idxmax()
        for name, pw in (("paper biomarker (amino acid synthesis)", AA), ("rank-1 pathway", r1)):
            if pw not in Pz.columns:
                continue
            s = c[pw] * Pz[pw]
            s = s[s.index.isin(prof)]
            p = mwu(s, mut, "less")
            sh_all, n_all, ps_all = all_pathway_null(Pz, prof, mut, "less", p)
            Xo = d["X"] if set(d["uni"]) <= set(d["X"].columns) else None
            # organoid matrix with every pathway of the universe (scaled across organoids), for the pipeline null
            Xall = pd.DataFrame(StandardScaler().fit_transform(d["osc"].loc[d["uni"], d["S"]].T), index=d["S"], columns=d["uni"])
            sh_pipe, ps_pipe = pipeline_null(Xall, d["y"].values.astype(float), Pz, list(d["uni"]), len(d["F"]), prof, mut,
                                             "less", p, n=2000, seed=1)
            nulls[(v, name)] = (p, ps_all, ps_pipe)
            a_t = adjust(s, mut, cov_t, -1)
            a_ts = adjust(s, mut, cov_ts, -1)
            rho_tmb = pd.Series(s).corr(np.log1p(T.reindex(s.index)), method="spearman")
            rows.append(dict(version=v, model=name, pathway=pw.replace("REACTOME_", ""), coef=c[pw],
                             ridge_rank=int((c.abs() > abs(c[pw])).sum() + 1), n_proximal=len(d["F"]),
                             n_mut=int(s.index.isin(mut).sum()), n_wt=int((~s.index.isin(mut)).sum()), p_one_sided=p,
                             share_all_pathways_both_signs=sh_all, n_null_tests=n_all,
                             share_random_proximal_sized_sets=sh_pipe, n_random_sets=2000,
                             spearman_score_vs_logTMB=rho_tmb,
                             median_logTMB_mut=float(np.median(np.log1p(T.reindex(list(mut & set(s.index)))))),
                             median_logTMB_wt=float(np.median(np.log1p(T.reindex([x for x in s.index if x not in mut])))),
                             **{f"adjTMB_{k}": val for k, val in a_t.items()},
                             **{f"adjTMBstage_{k}": val for k, val in a_ts.items()}))
            print(v, name, p, sh_all, sh_pipe, a_t.get("ols_p_one_sided"), flush=True)
    R = pd.DataFrame(rows)
    R.insert(0, "test", "ERCC2 (bladder, cisplatin)")
    R.to_csv(TAB / "G8_ERCC2_null_and_adjustment.tsv", sep="\t", index=False)
    pd.to_pickle(nulls, HERE.parent / "_data" / "g8_ercc2_nulls.pkl")
    pd.DataFrame(version_rows([mutf])).to_csv(TAB / "G8_ERCC2_data_versions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
