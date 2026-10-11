"""Audit 6: how much power does pooling by mechanism class (within cancer) or by drug (across cancers) buy?

Output: out/tcga_events_by_cancer_class.tsv (cancer x drug class; events among treated / untreated, power).
"""
import pandas as pd

from common import OUT, R, events_main, expr_patients, mdhr, power_interaction, power_main, systemic_long

CLASSES = {
    "platinum": {"CISPLATIN", "CARBOPLATIN", "OXALIPLATIN"},
    "taxane": {"PACLITAXEL", "DOCETAXEL"},
    "anthracycline": {"DOXORUBICIN", "EPIRUBICIN"},
    "fluoropyrimidine": {"FLUOROURACIL", "CAPECITABINE"},
    "gemcitabine/pemetrexed": {"GEMCITABINE", "PEMETREXED"},
    "alkylating": {"TEMOZOLOMIDE", "CYCLOPHOSPHAMIDE", "DACARBAZINE", "IFOSFAMIDE", "LOMUSTINE"},
}


def main():
    long = systemic_long()
    cdr = R.load_cdr()
    systemic = set(long["patient"])
    rows = []
    for proj in sorted(long["project"].unique()):
        pts = expr_patients(proj)
        c = cdr.reindex(sorted(pts))
        c = c[c["PFI.time"].notna()]
        if c.empty:
            continue
        un = c.index.difference(list(systemic))
        Du = c.loc[un, "PFI"].sum()
        g = long[long["project"] == proj]
        for cls, drugs in CLASSES.items():
            tr = sorted(set(g[g["drug"].isin(drugs)]["patient"]) & set(c.index))
            if len(tr) < 10:
                continue
            Dt = c.loc[tr, "PFI"].sum()
            rows.append({"project": proj, "class": cls, "n_treated": len(tr), "PFI_events_treated": int(Dt), "PFI_events_untreated": int(Du),
                         "power_main_HR0.7": power_main(Dt), "power_int_HR0.7": power_interaction(Dt, Du), "MDHR": mdhr(Dt),
                         "main_testable": Dt >= events_main(0.7), "int_testable": power_interaction(Dt, Du) >= 0.8})
    D = pd.DataFrame(rows).sort_values("PFI_events_treated", ascending=False)
    D.to_csv(OUT / "tcga_events_by_cancer_class.tsv", sep="\t", index=False)
    print(D.round(2).head(30).to_string())
    print("main testable:", int(D["main_testable"].sum()), "of", len(D), "; interaction testable:", int(D["int_testable"].sum()))


if __name__ == "__main__":
    main()
