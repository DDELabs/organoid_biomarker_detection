"""Audit 1: TCGA response-label quality and regimen confounding.

Outputs (robust/revamp_audit/out/):
  label_vs_pfi_by_category.tsv    12-month progression-free / alive rate by RECIST category and setting
  label_vs_pfi_agreement.tsv      responder label vs 12-month PFI agreement (kappa, PPV/NPV, C-index)
  label_vs_firstcourse.tsv        drug-record label vs CDR treatment_outcome_first_course
  regimen_confounding.tsv         per drug: multi-drug fraction, top partner, label sharing
  regimen_identical_pairs.tsv     cancer x drug pairs whose labelled patients are (nearly) the same set
  label_summary.json              headline numbers
"""
import json

import numpy as np
import pandas as pd

from common import CR_LIKE, NR_LIKE, OUT, R, systemic_long


def cohen_kappa(a, b):
    a, b = np.asarray(a, int), np.asarray(b, int)
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


def harrell_c(score, time, event):
    """C-index: higher score = longer time."""
    s, t, e = map(lambda x: np.asarray(x, float), (score, time, event))
    num = den = 0.0
    for i in range(len(t)):
        if not e[i]:
            continue
        m = t > t[i]
        den += m.sum()
        num += (s[m] > s[i]).sum() + 0.5 * (s[m] == s[i]).sum()
    return num / den if den else np.nan


def main():
    long = systemic_long()
    cdr = R.load_cdr()
    lab = long[long["recist"].notna()].copy()
    # one label per patient x drug: earliest record (same rule as patient_labels(how='first'))
    order = {"CR": 0, "PR": 1, "SD": 2, "PD": 3}
    lab["_r"] = lab["recist"].map(order)
    lab = lab.sort_values(["start_days", "_r"], na_position="last").drop_duplicates(["patient", "drug"])
    c = cdr.reindex(lab["patient"])
    lab["PFI"] = c["PFI"].values
    lab["PFI.time"] = c["PFI.time"].values
    lab["OS"] = c["OS"].values
    lab["OS.time"] = c["OS.time"].values
    lab["first_course"] = c["treatment_outcome_first_course"].values
    summ = {"n_labelled_patient_drug": int(len(lab)), "n_labelled_patients": int(lab["patient"].nunique())}

    # ---------------- 1. RECIST label vs 12-month outcome measured from drug start
    rows, agree = [], []
    for setting in ["early", "post_progression", "late_no_event", "all"]:
        d = lab if setting == "all" else lab[lab["setting_proxy"] == setting]
        d = d[d["start_days"].notna()]
        if setting == "post_progression":  # PFI event already happened before the drug -> use OS
            ep, t = "OS", d["OS.time"] - d["start_days"]
            ev = d["OS"]
        else:
            ep, t = "PFI", d["PFI.time"] - d["start_days"]
            ev = d["PFI"]
        ok = t.notna() & ev.notna() & (t >= 0)
        d, t, ev = d[ok], t[ok], ev[ok]
        # event-free at 12 m: known only if event after 365 d or censored after 365 d
        known = (ev.eq(1) & (t <= 365)) | (t > 365)
        free12 = (t > 365).astype(int)
        for cat in ["CR", "PR", "SD", "PD"]:
            m = known & d["recist"].eq(cat)
            rows.append({"setting": setting, "endpoint_from_drug_start": f"{ep}-free at 12 m", "recist": cat,
                         "n_known": int(m.sum()), "event_free_12m": float(free12[m].mean()) if m.sum() else np.nan,
                         "n_all": int(d["recist"].eq(cat).sum())})
        k = known
        y, f = d.loc[k, "responder"].astype(int), free12[k]
        if len(y) > 20:
            agree.append({
                "setting": setting, "endpoint": f"{ep}-free at 12 m from drug start", "n": int(k.sum()),
                "label_responder_rate": float(y.mean()), "event_free_rate": float(f.mean()),
                "raw_agreement": float((y == f).mean()), "kappa": float(cohen_kappa(y, f)),
                "PPV_responder_eventfree": float(f[y == 1].mean()), "NPV_nonresponder_event": float(1 - f[y == 0].mean()),
                "C_index_responder_vs_time": float(harrell_c(d["responder"], t, ev)),
            })
    pd.DataFrame(rows).to_csv(OUT / "label_vs_pfi_by_category.tsv", sep="\t", index=False)
    A = pd.DataFrame(agree)
    A.to_csv(OUT / "label_vs_pfi_agreement.tsv", sep="\t", index=False)
    summ["agreement"] = A.round(3).to_dict(orient="records")

    # CR share among responders by setting (adjuvant CR = 'no recurrence', not shrinkage)
    rr = lab[lab["responder"] == 1].groupby("setting_proxy")["recist"].apply(lambda s: float((s == "CR").mean()))
    summ["CR_share_among_responders_by_setting"] = rr.round(3).to_dict()
    summ["setting_share_of_labels"] = lab["setting_proxy"].value_counts(normalize=True).round(3).to_dict()

    # ---------------- 2. drug-record label vs CDR first-course outcome (early setting)
    e = lab[(lab["setting_proxy"] == "early") & lab["first_course"].isin(CR_LIKE | NR_LIKE)]
    fc = e["first_course"].isin(CR_LIKE).astype(int)
    fcrow = {"n": int(len(e)), "raw_agreement": float((fc == e["responder"]).mean()), "kappa": float(cohen_kappa(e["responder"], fc)),
             "label_responder_rate": float(e["responder"].mean()), "firstcourse_responder_rate": float(fc.mean())}
    pd.DataFrame([fcrow]).to_csv(OUT / "label_vs_firstcourse.tsv", sep="\t", index=False)
    summ["label_vs_firstcourse_early"] = fcrow

    # ---------------- 3. regimen confounding
    L = long.copy()
    L["regkey"] = L["patient"] + "|" + L["start_days"].astype(str)
    # own count of distinct systemic (non-supportive) drugs started the same day (start known)
    nd = L[L["start_days"].notna()].groupby("regkey")["drug"].nunique()
    L["n_systemic_same_day"] = L["regkey"].map(nd)
    lab_all = L[L["recist"].notna()]
    top = lab.groupby("drug")["patient"].nunique().sort_values(ascending=False)
    top = top[top >= 30].index
    reg = []
    for drug in top:
        d = lab_all[lab_all["drug"] == drug].drop_duplicates(["patient"])
        multi = (d["n_drugs_in_regimen"] > 1).mean()
        dk = d[d["start_days"].notna()]
        multi_sys = float((dk["n_systemic_same_day"] > 1).mean()) if len(dk) else np.nan
        partners = L[L["regkey"].isin(d["regkey"]) & (L["drug"] != drug) & d["start_days"].notna().any()]
        partners = partners[partners["regkey"].isin(d.loc[d["start_days"].notna(), "regkey"])]
        pc = partners.groupby("drug")["patient"].nunique().sort_values(ascending=False)
        p1 = pc.index[0] if len(pc) else ""
        # label concordance with the partner drug's own label in the same regimen
        conc = np.nan
        if p1:
            a = lab_all[lab_all["drug"] == drug].drop_duplicates("regkey").set_index("regkey")["recist"]
            b = lab_all[lab_all["drug"] == p1].drop_duplicates("regkey").set_index("regkey")["recist"]
            j = a.index.intersection(b.index)
            conc = float((a[j] == b[j]).mean()) if len(j) >= 10 else np.nan
        reg.append({"drug": drug, "n_labelled_patients": int(len(d)), "frac_multidrug_regimen": float(multi), "frac_multi_systemic_same_day": multi_sys,
                    "top_partner": p1, "frac_with_top_partner": float(pc.iloc[0] / max(1, d["start_days"].notna().sum())) if len(pc) else 0.0,
                    "label_identical_to_partner_label": conc, "median_n_drugs": float(d["n_drugs_in_regimen"].median())})
    G = pd.DataFrame(reg)
    G.to_csv(OUT / "regimen_confounding.tsv", sep="\t", index=False)
    w = lab_all.drop_duplicates(["patient", "drug"])
    summ["frac_labelled_records_multidrug"] = float((w["n_drugs_in_regimen"] > 1).mean())
    summ["frac_labelled_records_multidrug_weighted_top_drugs"] = float(np.average(G["frac_multidrug_regimen"], weights=G["n_labelled_patients"]))
    summ["median_label_identity_with_partner"] = float(G["label_identical_to_partner_label"].median())

    # identical patient sets within a cancer (non-identifiable drug effects)
    ident = []
    for proj, g in lab.groupby("project"):
        sets = {d: set(x["patient"]) for d, x in g.groupby("drug") if len(x) >= 10}
        ds = sorted(sets)
        for i, a in enumerate(ds):
            for b in ds[i + 1:]:
                ja = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
                if ja >= 0.6:
                    ident.append({"project": proj, "drug_a": a, "drug_b": b, "n_a": len(sets[a]), "n_b": len(sets[b]), "jaccard": round(ja, 3)})
    I = pd.DataFrame(ident).sort_values("jaccard", ascending=False)
    I.to_csv(OUT / "regimen_identical_pairs.tsv", sep="\t", index=False)
    summ["n_cancer_drug_pairs_jaccard_ge_0.6"] = int(len(I))
    summ["n_cancer_drug_pairs_jaccard_ge_0.9"] = int((I["jaccard"] >= 0.9).sum()) if len(I) else 0

    # class imbalance in the most-used labels
    imb = lab.groupby(["project", "drug"]).agg(n=("responder", "size"), resp=("responder", "mean"))
    imb = imb[imb["n"] >= 20]
    summ["pairs_n_ge_20"] = int(len(imb))
    summ["pairs_responder_rate_gt_0.85_or_lt_0.15"] = int(((imb["resp"] > 0.85) | (imb["resp"] < 0.15)).sum())
    json.dump(summ, open(OUT / "label_summary.json", "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float))
    print(G.to_string())
    print(I.head(30).to_string())


if __name__ == "__main__":
    main()
