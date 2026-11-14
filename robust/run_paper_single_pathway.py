"""The paper's HEADLINE analysis (Kong et al. 2020, Fig. 2): single-pathway biomarker.

Procedure of python/single_pathway_prediction.py: proximal pathways -> ridge on organoid IC50 -> the pathway with the
largest |coefficient| is the biomarker -> patients given the drug are split at the median of coef * pathway score ->
log-rank test on overall survival. Paper: COAD 5-FU P = 0.014 (114 patients; untreated control P = 0.16);
BLCA cisplatin P = 0.01 (77 patients; untreated control P = 0.066).

Part A  authors' own committed data (colon, 5-FU): exact reproduction + checks.
Part B  current public data, all testable organoid-drug-cancer sets, same procedure.
Checks added: direction, untreated control, every proximal pathway tested alone, shuffled organoid response (2000x).
"""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from statsmodels.duration.survfunc import survdiff

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

from obd import RESULTS  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = RESULTS / "paper_single_pathway"


def fit(ml, X, y):
    if ml == "Ridge":
        return RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(X, y).coef_
    if ml == "SVR":
        return SVR(kernel="linear").fit(X, y).coef_[0]
    return LinearRegression().fit(X, y).coef_


def split_test(pred, T, E):
    """Median split as in the paper; returns (p, deaths in predicted non-responders, deaths in responders)."""
    g = (pred > np.median(pred)).astype(int)
    if g.sum() in (0, len(g)):
        return 1.0, 0, 0
    return survdiff(T, E, g)[1], int(E[g == 1].sum()), int(E[g == 0].sum())


def analyse(tag, X, y, feats, Pt, Tt, Et, Pu, Tu, Eu, n_perm, rng):
    rows = []
    coefs = {ml: fit(ml, X, y) for ml in ("Ridge", "SVR", "LinearRegression")}
    top = {ml: int(np.argmax(np.abs(c))) for ml, c in coefs.items()}
    i = top["Ridge"]
    c = coefs["Ridge"][i]
    p, dn, dr = split_test(c * Pt[:, i], Tt, Et)
    pu, un, ur = split_test(c * Pu[:, i], Tu, Eu) if len(Tu) else (np.nan, 0, 0)
    each = [split_test(coefs["Ridge"][j] * Pt[:, j], Tt, Et) for j in range(len(feats))]
    ep = np.array([e[0] for e in each])
    right = np.array([e[1] > e[2] for e in each])
    sh = []
    for _ in range(n_perm):
        cc = fit("Ridge", X, rng.permutation(y))
        j = int(np.argmax(np.abs(cc)))
        sh.append(split_test(cc[j] * Pt[:, j], Tt, Et)[0])
    sh = np.array(sh)
    return {"set": tag, "organoids": len(y), "proximal_pathways": len(feats), "biomarker_ridge": feats[i].replace("REACTOME_", ""),
            "same_top_pathway_all_3_models": len(set(top.values())) == 1,
            "patients": len(Tt), "deaths": int(Et.sum()), "p_treated": p, "deaths_nonresponders": dn, "deaths_responders": dr,
            "direction_correct": dn > dr, "untreated_patients": len(Tu), "untreated_deaths": int(Eu.sum()) if len(Tu) else 0,
            "p_untreated": pu, "untreated_deaths_nonresp": un, "untreated_deaths_resp": ur,
            "pathways_alone_p_le_observed": int((ep <= p + 1e-12).sum()), "pathways_alone_p_below_0.05": int((ep < 0.05).sum()),
            "pathways_alone_p_below_0.05_right_direction": int(((ep < 0.05) & right).sum()),
            "shuffled_as_good": (np.sum(sh <= p + 1e-12) + 1) / (n_perm + 1), "shuffled_p_below_0.05": float((sh < 0.05).mean())}


def authors_data():
    nes = lambda p: pd.read_csv(p, sep="\t", comment="#", index_col=0)  # noqa: E731
    org = nes(REPO / "python/results/COAD/organoid/reactome_ssgsea_result/gseapy.samples.normalized.es.txt")
    pat = nes(REPO / "python/results/COAD/TCGA/reactome_ssgsea_result/gseapy.samples.normalized.es.txt")
    org.columns = [c.upper() for c in org.columns]
    resp = pd.read_csv(REPO / "data/organoid_COAD/drug_response/median_response.txt", sep="\t")
    resp = resp[resp["Drug Name"].str.upper() == "5-FLUOROURACIL"]
    resp = dict(zip(resp["Organoid"].str.upper(), resp["median_IC50"]))
    z = pd.read_csv(REPO / "data/coad_blca_organoid_drugs_zscore_result_reactome.txt", sep="\t")
    prox = set(z.loc[z["FLUOROURACIL"] <= -1.2816, "Pathway"])
    syn = {}
    for line in csv.reader(open(REPO / "data/drugbank_vocabulary.csv", encoding="utf-8", errors="replace")):
        if "DrugBank" in line[0]:
            continue
        for s in line[5].split(" | "):
            s = s.upper()
            if s:
                syn[s] = syn[s.replace("-", "")] = syn[s.replace(" ", "")] = line[2].upper()
    d = pd.read_csv(REPO / "data/TCGA_COAD/clinical_drug/dae92461-d098-40a3-997e-1a2d3e5b8dfe/nationwidechildrens.org_clinical_drug_coad.txt", sep="\t")
    treated, anydrug = set(), set()
    for p, dn in zip(d["bcr_patient_barcode"], d["pharmaceutical_therapy_drug_name"]):
        if "TCGA" not in str(p):
            continue
        anydrug.add(p)
        for t in str(dn).upper().replace("-", "").replace(" ", "").split("AND"):
            if "UNKNOWN" not in t and "NOTAVAILABLE" not in t and syn.get(t, t) == "FLUOROURACIL":
                treated.add(p)
    c = pd.read_csv(REPO / "data/TCGA_COAD/clinical_patient/4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt", sep="\t")
    sur = {}
    for p, st, lc, dd in zip(c["bcr_patient_barcode"], c["vital_status"], c["last_contact_days_to"], c["death_days_to"]):
        st = str(st).upper()
        days = lc if st == "ALIVE" else dd if st == "DEAD" else None
        if "TCGA" in str(p) and days is not None and "[" not in str(days):
            sur[p] = (float(days) / (365. / 12.), 1 if st == "DEAD" else 0)
    feats = sorted(set(org.index) & set(pat.index) & prox)
    samples = sorted(set(org.columns) & set(resp))
    X = StandardScaler().fit_transform(org.loc[feats, samples].T.values)
    y = np.array([resp[s] for s in samples])
    pats = sorted(set(sur) & set(pat.columns))
    P = pd.DataFrame(StandardScaler().fit_transform(pat.loc[feats, pats].T.values), index=pats, columns=feats)
    tp, up = [p for p in pats if p in treated], [p for p in pats if p not in anydrug]
    te = lambda ps: (np.array([sur[p][0] for p in ps]), np.array([sur[p][1] for p in ps]))  # noqa: E731
    return X, y, feats, P.loc[tp].values, *te(tp), P.loc[up].values, *te(up)


def current_data(oset, proj, drug, ocol, rcol):
    import run_paper_method_all as M
    from obd import network as N, tcga_response as R
    from obd.curated import PRECLIN
    import run_survlearn as RS
    osc = pd.read_csv(RESULTS / "cache" / f"paper_{oset}_{proj}_organoid_ssgsea.tsv.gz", sep="\t", index_col=0)
    tsc = pd.read_csv(RESULTS / "cache" / f"paper_{oset}_{proj}_tcga_ssgsea.tsv.gz", sep="\t", index_col=0)
    paper, long, cdr = N.precomputed_proximity(), R.load_response_labels(), R.load_cdr()
    r = pd.read_csv(PRECLIN / oset / "response.tsv", sep="\t")
    r = r[r["drug"].str.upper() == ocol].groupby("sample")[rcol].median().dropna()
    samples = [s for s in r.index if s in osc.columns]
    z = paper[drug].dropna() if drug in paper.columns else RS.proximity(drug)
    feats = [p for p in z.index[z <= -1.2816] if p in osc.index and p in tsc.index]
    pats = [p for p in tsc.columns if p in cdr.index and pd.notna(cdr.loc[p, "OS.time"]) and pd.notna(cdr.loc[p, "OS"])]
    P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[feats, pats].T.values), index=pats, columns=feats)
    lp = long[long["project"] == proj]
    treated, anydrug = set(lp.loc[lp["drug"] == drug, "patient"]), set(lp["patient"])
    tp, up = [p for p in pats if p in treated], [p for p in pats if p not in anydrug]
    te = lambda ps: (cdr.loc[ps, "OS.time"].values / (365 / 12), cdr.loc[ps, "OS"].values.astype(int))  # noqa: E731
    X = StandardScaler().fit_transform(osc.loc[feats, samples].T.values)
    return X, r.loc[samples].values, feats, P.loc[tp].values, *te(tp), P.loc[up].values, *te(up)


def main():
    import run_paper_method_all as M
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    rows = [analyse("A. authors' data: COAD FLUOROURACIL", *authors_data(), 2000, rng)]
    print(rows[-1], flush=True)
    for oset, proj, drug, ocol, rcol in M.SETS:
        rows.append(analyse(f"B. current data: {proj} {drug}", *current_data(oset, proj, drug, ocol, rcol), 2000, rng))
        print(rows[-1]["set"], round(rows[-1]["p_treated"], 4), flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 300)
    print(S.drop(columns=["biomarker_ridge"]).round(3).T.to_string())
    print(S[["set", "biomarker_ridge"]].to_string(index=False))


if __name__ == "__main__":
    main()
