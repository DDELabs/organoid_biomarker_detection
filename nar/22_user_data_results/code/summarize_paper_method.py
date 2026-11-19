"""Read the pickles written by run_paper_method.py and write tables/paper_method_<plan>.tsv  (every number comes from the pickles)."""
import sys, pickle
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u22
from u22 import DATA, OUT, pd, np
import engine as EN
plan = sys.argv[1]
rows = []
for f in sorted((DATA / "results" / plan).glob("*.pkl")):
    r = pickle.load(open(f, "rb"))
    m = r["meta"]
    base = dict(plan=plan, set=m.get("set"), cohort=m.get("cohort"), organoid_drug=m.get("model_drug"), patient_drug=m.get("patient_drug"),
                proximity_for=m.get("proximity_for"), proximity_source=m.get("proximity_column"))
    if "ev" not in r:
        rows.append({**base, "n_organoids": m.get("n_models"), "n_treated": m.get("n_treated"), "note": m.get("skip")})
        continue
    ev = r["ev"]
    r1 = ev[(ev.model == "Ridge") & (ev.k == 1)].iloc[0]
    kmax = ev[ev.model == "Ridge"].k.max()
    rk = ev[(ev.model == "Ridge") & (ev.k == kmax)].iloc[0]
    T_u, E_u, g_u, _ = r["km"][("untreated", 1)]
    su = EN.surv_stats(T_u, E_u, g_u, True)
    coef1 = float(r["coefs"]["Ridge"].loc[r1["pathways"]])
    rows.append({**base, "n_organoids": m["n_models"], "n_close_pathways": m["n_proximal"], "top_pathway": r1["pathways"], "top_pathway_ridge_coef": coef1,
                 "n_treated": r1["treated_n"], "events_treated": r1["treated_events"], "logrank_p": r1["treated_p"], "HR": r1["treated_HR"],
                 "HR_lo": r1["treated_HR_lo"], "HR_hi": r1["treated_HR_hi"],
                 "n_untreated": su["n"], "events_untreated": su["events"], "logrank_p_untreated": su["p"], "HR_untreated": su["HR"],
                 "HR_untreated_lo": su["HR_lo"], "HR_untreated_hi": su["HR_hi"],
                 "permutation_p_rank1": r["perm_rank1_emp_p"], "permutation_p_best_of_27": r["perm_best_emp_p"], "bootstrap_p": r["boot_emp_p"],
                 "organoid_LOO_spearman": r["loo"]["Ridge"][1], "organoid_LOO_p": r["loo"]["Ridge"][2],
                 f"multi_k{kmax}_logrank_p": rk["treated_p"], f"multi_k{kmax}_HR": rk["treated_HR"], "multi_k": kmax,
                 "smallest_p_of_27_tests": ev["treated_p"].min(), "borderline_small_n": m.get("borderline"), "note": ""})
out = pd.DataFrame(rows)
out.to_csv(OUT / "tables" / f"paper_method_{plan}.tsv", sep="\t", index=False)
print(out[["cohort", "organoid_drug", "patient_drug", "n_treated", "events_treated", "logrank_p", "HR", "permutation_p_rank1"]].to_string())
