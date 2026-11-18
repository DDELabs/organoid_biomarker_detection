"""Re-derive the treated / untreated median-split log-rank p of every headline biomarker from the bundles and the
coefficient tables, and compare with the numbers in the organoid-only summary tables (verification step)."""
from common import *

g = pd.read_csv(O2 / "tables/cross_set_summary.tsv", sep="\t")
o = pd.read_csv(O2 / "BLCA_OV_LUAD_OTHER/SUMMARY/tables/cross_set_summary.tsv", sep="\t")
rows = []
for k, (lab, bn, drug, cf, grp) in REG.items():
    b = bundle(bn)
    pw, c = rank1(k)
    if pw not in b["tsc"].index:
        rows.append(dict(key=k, label=lab, pathway=pw, note="pathway not in patient scores")); continue
    sc = b["tsc"].loc[pw] * np.sign(c)
    t, _ = treated_split(b, drug, sc); u, _ = untreated_split(b, sc)
    rep = np.nan
    if "COAD_vandewetering" in cf or "PAAD" in cf or "LIHC" in cf or "isella" in cf:
        coh = {"COAD_paper": "COAD_authors", "COADREAD_vdW": "COADREAD"}.get(bn, b["meta"]["project"])
        folder = cf.split("/")[0]
        m = g[(g.folder == folder) & (g.cohort == coh) & (g.patient_drug == drug)]
        rep = float(m.treated_p.iloc[0]) if len(m) else np.nan
    else:
        m = o[(o.folder == cf.rsplit("/tables", 1)[0]) & (o.drug == drug)]
        if bn.startswith("PDX"):
            m = o[(o.folder == "PDX_novartis_gao2015") & (o.drug == drug) & (o.project == b["meta"]["project"])]
        rep = float(m.treated_p.iloc[0]) if len(m) else np.nan
    rows.append(dict(key=k, label=lab, bundle=bn, drug=drug, pathway=pw, ridge_coef=c, treated_n=t["n"],
                     treated_events=t["events"], treated_p=t["p"], treated_HR=t["HR"], treated_HR_lo=t["lo"],
                     treated_HR_hi=t["hi"], untreated_n=u["n"], untreated_p=u["p"], reported_treated_p=rep))
df = pd.DataFrame(rows)
df["match"] = (df.treated_p - df.reported_treated_p).abs() < 1e-6
df.to_csv(TAB / "verify_headline_biomarkers.tsv", sep="\t", index=False)
print(df[["key", "pathway", "treated_n", "treated_p", "reported_treated_p", "match"]].to_string())
