"""Audit 4: power available in the curated external cohorts (data/curated/trials/CATALOG_GEO.tsv and
data/external/trials/CATALOG.tsv).

For each cohort (and each arm when arms exist): labelled n, responders, OS/PFS/RFS/DRFS events,
power for AUROC 0.60 / 0.65 (binary endpoints), power for a main-effect HR 0.7/SD, and, where two arms
exist, power for a score x arm interaction HR 0.7 (var = 1/D_arm1 + 1/D_arm2 for time-to-event;
1/n11+1/n10+1/n21+1/n20 on the log-OR scale for binary endpoints).

Output: robust/revamp_audit/out/external_cohort_power.tsv, external_power_summary.json
"""
import json

import numpy as np
import pandas as pd

from common import OUT, REPO, auc_to_delta, power_auc, power_from_var, power_main

SRC = [(REPO / "data/curated/trials/CATALOG_GEO.tsv", REPO / "data/curated/trials"),
       (REPO / "data/external/trials/CATALOG.tsv", REPO / "data/external/trials")]
TTE = [("os", "os_event"), ("pfs", "pfs_event"), ("rfs", "rfs_event"), ("drfs", "drfs_event")]


def main():
    rows = []
    for cat_path, root in SRC:
        if not cat_path.exists():
            continue
        cat = pd.read_csv(cat_path, sep="\t")
        for _, c in cat.iterrows():
            cid = c["cohort_id"]
            f = root / cid / "clinical.tsv"
            if not f.exists():
                continue
            cl = pd.read_csv(f, sep="\t")
            base = {"cohort": cid, "source": cat_path.parent.name, "cancer": c.get("cancer"), "drugs": c.get("drugs"),
                    "setting": c.get("setting"), "randomised": c.get("randomised", c.get("has_control_arm")),
                    "usable": c.get("usable", True), "qc_flag": c.get("qc_flag", False), "n": len(cl)}
            y = cl["responder"] if "responder" in cl else pd.Series(dtype=float)
            n1, n0 = int((y == 1).sum()), int((y == 0).sum())
            base.update({"n_labelled": n1 + n0, "responders": n1, "nonresponders": n0,
                         "power_auc060": power_auc(n1, n0, 0.60), "power_auc065": power_auc(n1, n0, 0.65)})
            for ep, ev in TTE:
                if ev in cl:
                    D = int(pd.to_numeric(cl[ev], errors="coerce").fillna(0).sum())
                    base[f"{ep}_events"] = D
            evs = {k: base.get(f"{k}_events", 0) for k, _ in TTE}
            best_ep = max(evs, key=evs.get)
            Dbest = evs[best_ep]
            base["best_tte"] = best_ep if Dbest else ""
            base["best_tte_events"] = Dbest
            base["power_main_HR0.7"] = power_main(Dbest) if Dbest else np.nan
            # two-arm interaction
            arms = cl["arm"].value_counts() if "arm" in cl else pd.Series(dtype=int)
            if len(arms) >= 2:
                a1, a2 = arms.index[:2]
                base["arm1"], base["arm2"] = a1, a2
                if Dbest:
                    evc = dict(TTE)[best_ep]
                    d1 = pd.to_numeric(cl.loc[cl["arm"] == a1, evc], errors="coerce").sum()
                    d2 = pd.to_numeric(cl.loc[cl["arm"] == a2, evc], errors="coerce").sum()
                    base["arm_events"] = f"{int(d1)}/{int(d2)}"
                    base["power_int_tte_HR0.7"] = power_from_var(np.log(0.7), 1 / d1 + 1 / d2) if d1 > 0 and d2 > 0 else np.nan
                if n1 + n0:
                    cc = [((cl["arm"] == a) & (y == v)).sum() for a in (a1, a2) for v in (1, 0)]
                    if min(cc) > 0:
                        # interaction on logit scale: difference of two per-SD log-ORs (delta for AUC 0.65 vs 0.5 in one arm)
                        var = sum(1 / x for x in cc)
                        base["power_int_binary_auc065_vs_050"] = power_from_var(auc_to_delta(0.65), var)
            rows.append(base)
    E = pd.DataFrame(rows)
    E.to_csv(OUT / "external_cohort_power.tsv", sep="\t", index=False)
    s = {"cohorts": int(len(E)),
         "cohorts_power_auc065_ge_0.8": int((E["power_auc065"] >= 0.8).sum()),
         "cohorts_power_auc060_ge_0.8": int((E["power_auc060"] >= 0.8).sum()),
         "cohorts_power_main_hr07_ge_0.8": int((E["power_main_HR0.7"] >= 0.8).sum()),
         "cohorts_power_int_tte_ge_0.8": int((E.get("power_int_tte_HR0.7", pd.Series(dtype=float)) >= 0.8).sum()),
         "cohorts_power_int_binary_ge_0.8": int((E.get("power_int_binary_auc065_vs_050", pd.Series(dtype=float)) >= 0.8).sum()),
         "list_auc065_ge_0.8": E.loc[E["power_auc065"] >= 0.8, "cohort"].tolist(),
         "list_main_ge_0.8": E.loc[E["power_main_HR0.7"] >= 0.8, "cohort"].tolist(),
         "list_int_tte_ge_0.5": E.loc[E.get("power_int_tte_HR0.7", pd.Series(0, index=E.index)) >= 0.5, "cohort"].tolist()}
    json.dump(s, open(OUT / "external_power_summary.json", "w"), indent=1)
    print(json.dumps(s, indent=1))
    cols = ["cohort", "cancer", "n_labelled", "responders", "power_auc065", "best_tte", "best_tte_events", "power_main_HR0.7",
            "arm1", "arm2", "arm_events", "power_int_tte_HR0.7", "power_int_binary_auc065_vs_050", "qc_flag"]
    pd.set_option("display.width", 250)
    print(E[[c for c in cols if c in E]].round(2).to_string())


if __name__ == "__main__":
    main()
