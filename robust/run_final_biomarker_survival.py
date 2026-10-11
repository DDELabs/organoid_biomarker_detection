"""Survival validation of every final cancer x drug pathway biomarker (one consolidated table).

For each biomarker row (the final summary table: response pathways up, resistance pathways down), the score
= mean z(response pathways) - mean z(resistance pathways), standardised within cohort/cancer; HR < 1 means a
higher score goes with longer survival.

TCGA (TCGA-CDR OS and PFI; Cox adjusted for age and stage; HR per SD):
  held-out treated   received the drug, not used for discovery (RECIST- or first-course-labelled patients
                     removed). If < 20 patients / 8 events remain, the all-treated estimate is also given and
                     flagged as not independent.
  untreated          same cancer, no recorded systemic therapy (prognostic control)
  interaction        score x treated (held-out treated vs untreated): HR < 1 = predictive of drug benefit
  log-rank           the original paper's test: median split of treated patients
External cohorts (curated GEO / trial data): treated-only Cox and, where a control arm exists, the
score x arm interaction.

Verdict:
  validated-predictive  independent treated HR < 1 (p < 0.05) and interaction HR < 1 (p < 0.05)
  supported             independent treated HR < 1 (p < 0.05) somewhere, no contrary significant result
  trend                 independent treated or interaction HR < 1 with p < 0.2, nothing contrary
  contradicted          independent treated HR > 1 with p < 0.05
  not supported         testable but none of the above
  not testable          no independent set with >= 20 patients and >= 8 events

Outputs: robust/results/survival_validation/{FINAL_BIOMARKER_SURVIVAL.tsv, .xlsx, .md, final_long.tsv}
Usage: python robust/run_final_biomarker_survival.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import run_atlas_external as AE  # noqa: E402
import run_triad as RT  # noqa: E402
from obd import RESULTS  # noqa: E402
from obd import survival as SV  # noqa: E402
from obd import tcga_response as R  # noqa: E402
from obd import triad as T  # noqa: E402
from obd.cohorts import _stage  # noqa: E402
from run_atlas_relaxed import CR_LIKE, NR_LIKE, SUPPORTIVE, THERAPY  # noqa: E402

OUT = RESULTS / "survival_validation"
MIN_N, MIN_EV = 20, 8
P = "REACTOME_"
BCR, PROST, FA, RFA = P + "SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR", P + "PROSTANOID_LIGAND_RECEPTORS", P + "FANCONI_ANEMIA_PATHWAY", P + "REGULATION_OF_THE_FANCONI_ANEMIA_PATHWAY"
GAPD, FGFR3, RIGI, CC = P + "GAP_JUNCTION_DEGRADATION", P + "SIGNALING_BY_FGFR3_MUTANTS", P + "RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS", P + "CELL_CYCLE"
AQP, RNA, SEMA = P + "REGULATION_OF_WATER_BALANCE_BY_RENAL_AQUAPORINS", P + "METABOLISM_OF_RNA", P + "SEMAPHORIN_INTERACTIONS"
UTER = ([P + "SIGNALING_BY_FGFR1_MUTANTS", P + "IL_3_5_AND_GM_CSF_SIGNALING"], [PROST])
LUNG = ([CC, RNA], [RIGI])
SARC = ([RNA, RFA], [AQP])
MESO = ([CC, P + "BASE_EXCISION_REPAIR"], [PROST, AQP])
OVAR = ([P + "CELL_CELL_JUNCTION_ORGANIZATION"], [PROST, P + "SIGNALING_BY_BMP"])
ANTHR = ([FA, P + "DNA_REPAIR"], [GAPD])
FU = ([FA], [P + "PHASE_II_CONJUGATION"])
OXA = ([P + "IRAK2_MEDIATED_ACTIVATION_OF_TAK1_COMPLEX_UPON_TLR7_8_OR_9_STIMULATION"], [P + "TRANSCRIPTION_COUPLED_NER_TC_NER"])

# external cohorts: (tag, treated selector, control selector or None, label)
arm = lambda pat: (lambda c: c["arm"].astype(str).str.contains(pat, case=False, regex=False))  # noqa: E731
LUNG_EXT = [("trial:GSE14814_jbr10", arm("adjuvant chemo"), arm("observation"), "JBR.10 cisplatin/vinorelbine vs observation"),
            ("geo:GSE42127", arm("adjuvant chemo"), arm("observation"), "GSE42127 adjuvant platinum vs observation"),
            ("trial:GSE68465_lung_adj", arm("adjuvant chemo"), arm("observation"), "GSE68465 adjuvant chemo vs observation"),
            ("trial:GSE37745_lung_adj", arm("adjuvant chemo"), arm("observation"), "GSE37745 adjuvant chemo vs observation")]
FU_EXT = [("geo:GSE103479", arm("adjuvant chemo"), arm("surgery only"), "GSE103479 adjuvant 5-FU vs surgery"),
          ("geo:GSE72970", lambda c: pd.Series(True, index=c.index), None, "GSE72970 metastatic FOLFOX/FOLFIRI (all 5-FU)")]
OXA_EXT = [("geo:GSE72970", arm("FOLFOX"), arm("FOLFIRI"), "GSE72970 FOLFOX vs FOLFIRI")]
OV_PLAT = [("geo:GSE9891", lambda c: c["arm"].isin(["platinum+taxane", "platinum"]), arm("no platinum"), "GSE9891 platinum vs no platinum"),
           ("geo:GSE63885", lambda c: pd.Series(True, index=c.index), None, "GSE63885 platinum-based (all)")]
OV_TAX = [("geo:GSE9891", lambda c: c["arm"] == "platinum+taxane", lambda c: c["arm"] == "platinum", "GSE9891 platinum+taxane vs platinum"),
          ("geo:GSE63885", arm("taxane"), arm("cyclophosphamide"), "GSE63885 taxane/platinum vs platinum/cyclophosphamide")]

ROWS = [
    # cancer, drug, label, response, resistance, evidence, discovery rule, external
    ("BRCA", "PACLITAXEL", "breast", [BCR], [PROST], "strong", "recist", []),
    ("BRCA", "DOXORUBICIN", "breast", *ANTHR, "moderate", "recist", []),
    ("BRCA", "EPIRUBICIN", "breast", *ANTHR, "strong", "recist", []),
    ("BRCA", "CYCLOPHOSPHAMIDE", "breast", *ANTHR, "moderate", "recist", []),
    ("BLCA", "CISPLATIN", "bladder", [FGFR3, P + "CTLA4_INHIBITORY_SIGNALING"], [P + "MICRORNA_MIRNA_BIOGENESIS"], "moderate", "recist", []),
    ("BLCA", "GEMCITABINE", "bladder", [P + "BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY"], [SEMA], "moderate", "recist", []),
    ("HNSC", "CISPLATIN", "head & neck", [FGFR3], [RIGI], "moderate-strong", "recist", []),
    ("HNSC", "CETUXIMAB", "head & neck", ["GENE:AREG", "GENE:EREG"], [], "moderate", "gene", []),
    ("CESC", "CISPLATIN", "cervical", [CC], [RIGI, AQP], "moderate", "recist", []),
    ("LUAD", "CISPLATIN", "lung adeno", *LUNG, "moderate", "recist", LUNG_EXT),
    ("LUAD", "CARBOPLATIN", "lung adeno", *LUNG, "moderate", "recist", LUNG_EXT),
    ("LUSC", "CISPLATIN", "lung squamous", *LUNG, "moderate", "recist", LUNG_EXT),
    ("LUSC", "CARBOPLATIN", "lung squamous", *LUNG, "moderate", "recist", LUNG_EXT),
    ("LUAD", "PACLITAXEL", "lung adeno", [P + "GLUCOSE_METABOLISM"], [PROST], "weak-moderate", "recist", []),
    ("LUSC", "PACLITAXEL", "lung squamous", [P + "GLUCOSE_METABOLISM"], [PROST], "weak-moderate", "recist", []),
    ("LUAD", "PEMETREXED", "lung adeno", [P + "PHOSPHOLIPID_METABOLISM"], [P + "DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR"], "weak", "recist", []),
    ("STAD", "CISPLATIN", "stomach", [CC, FA], [PROST], "moderate", "recist", []),
    ("STAD", "EPIRUBICIN", "stomach", [RFA], [GAPD], "moderate", "recist", []),
    ("COAD", "FLUOROURACIL", "colon", *FU, "moderate", "recist", FU_EXT),
    ("READ", "FLUOROURACIL", "rectal", *FU, "moderate", "recist", []),
    ("COAD", "OXALIPLATIN", "colon", *OXA, "moderate", "recist", OXA_EXT),
    ("READ", "OXALIPLATIN", "rectal", *OXA, "moderate", "recist", []),
    ("COAD", "CETUXIMAB", "colon", ["SIG:cetuximab"], [], "model-only (PDX)", "recist", []),
    ("PAAD", "GEMCITABINE", "pancreatic", [P + "INFLAMMASOMES", P + "E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION", P + "PERK_REGULATED_GENE_EXPRESSION"],
     [SEMA, P + "TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"], "moderate", "recist",
     [("trial:PRINCE_paad", lambda c: pd.Series(True, index=c.index), None, "PRINCE gemcitabine/nab-paclitaxel + IO (all)")]),
    ("UCEC", "CARBOPLATIN", "endometrial", *UTER, "weak-moderate", "recist", []),
    ("UCEC", "PACLITAXEL", "endometrial", *UTER, "weak-moderate", "recist", []),
    ("UCS", "CARBOPLATIN", "uterine carcinosarcoma", *UTER, "weak-moderate", "recist", []),
    ("UCS", "PACLITAXEL", "uterine carcinosarcoma", *UTER, "weak-moderate", "recist", []),
    ("MESO", "CISPLATIN", "mesothelioma", *MESO, "weak", "recist", []),
    ("MESO", "PEMETREXED", "mesothelioma", *MESO, "weak", "recist", []),
    ("SARC", "GEMCITABINE", "sarcoma", *SARC, "weak", "recist", []),
    ("SARC", "DOCETAXEL", "sarcoma", *SARC, "weak", "recist", []),
    ("SARC", "DOXORUBICIN", "sarcoma", *SARC, "weak", "recist", []),
    ("OV", "CARBOPLATIN", "ovarian", *OVAR, "weak", "firstcourse", OV_PLAT),
    ("OV", "PACLITAXEL", "ovarian", *OVAR, "weak", "firstcourse", OV_TAX),
    ("LIHC", "TACE", "liver", [P + "XENOBIOTICS", P + "PHASE1_FUNCTIONALIZATION_OF_COMPOUNDS"],
     [P + "REGULATION_OF_HYPOXIA_INDUCIBLE_FACTOR_HIF_BY_OXYGEN", P + "GLUCOSE_TRANSPORT"], "moderate (cohort caveat)", "none", []),
    ("LIHC", "SORAFENIB", "liver", [P + "PRE_NOTCH_TRANSCRIPTION_AND_TRANSLATION", P + "NOTCH_HLH_TRANSCRIPTION_PATHWAY"],
     [P + "METABOLISM_OF_STEROID_HORMONES_AND_VITAMINS_A_AND_D"], "weak", "recist", []),
    ("KIRC", "SUNITINIB", "kidney", [], [], "none identified", "none", []),
    ("GBM", "TEMOZOLOMIDE", "glioblastoma", [], [P + "DNA_REPAIR", P + "MITOTIC_G2_G2_M_PHASES"], "weak", "none",
     [("geo:CGGA693", lambda c: (c["arm"] == "TMZ") & (c["grade"] == "WHO IV"), lambda c: (c["arm"] == "no TMZ") & (c["grade"] == "WHO IV"),
       "CGGA GBM TMZ vs none (discovery cohort)"),
      ("geo:GSE7696", arm("RT+TMZ"), lambda c: c["arm"] == "RT", "GSE7696 RT+TMZ vs RT (discovery cohort)")]),
    ("LGG", "TEMOZOLOMIDE", "lower-grade glioma", "TOP10", "TOP10", "survival-derived", "recist",
     [("geo:CGGA693", lambda c: (c["arm"] == "TMZ") & c["grade"].isin(["WHO II", "WHO III"]),
       lambda c: (c["arm"] == "no TMZ") & c["grade"].isin(["WHO II", "WHO III"]), "CGGA grade II/III TMZ vs none")]),
    ("SKCM", "DACARBAZINE", "melanoma", [P + "AMINE_DERIVED_HORMONES"], [P + "CIRCADIAN_CLOCK"], "weak-moderate", "recist", []),
    ("SKCM", "ANTI-PD-1", "melanoma", [], [], "none identified", "none", []),
    ("DLBC", "RITUXIMAB", "lymphoma (DLBCL)", [P + "G2_M_DNA_DAMAGE_CHECKPOINT"],
     [P + "TRAF6_MEDIATED_NFKB_ACTIVATION", P + "RIP_MEDIATED_NFKB_ACTIVATION_VIA_DAI"], "moderate", "firstcourse", []),
]
DISCOVERY_COHORTS = {"geo:CGGA693|GBM", "geo:GSE7696|GBM"}


def _fit(d, cols):
    X = d[cols].astype(float).fillna(d[cols].astype(float).median()).values
    r = PHReg(d["time"].values, X, status=d["event"].values, ties="efron").fit(disp=0)
    if not np.isfinite(r.params[0]) or not np.isfinite(r.bse[0]):
        raise ValueError("non-finite fit")
    ci = r.conf_int()[0]
    return float(np.exp(r.params[0])), float(np.exp(ci[0])), float(np.exp(ci[1])), float(r.pvalues[0])


def cox(df, cols):
    """HR per SD of cols[0]; covariates dropped when sparse (binary minority < 5) and, if the fit still fails,
    adjustment is reduced step-wise (all -> age only -> none). 'adj' records what was used."""
    d = df.dropna(subset=["time", "event", cols[0]])
    d = d[d["time"] > 0]
    n, ev = len(d), int(d["event"].sum()) if len(d) else 0
    nan = dict(HR=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=n, events=ev, adj="")
    if n < MIN_N or ev < MIN_EV:
        return nan
    keep = [cols[0]]
    for c in cols[1:]:
        v = d[c].astype(float)
        if v.nunique() <= 2:
            if min((v == v.max()).sum(), (v != v.max()).sum()) >= 5 and d.loc[v == v.max(), "event"].sum() >= 1:
                keep.append(c)
        elif v.std() > 0:
            keep.append(c)
    tries = [keep, [c for c in keep if c in (cols[0], "age", "treated", "score", "arm")], [c for c in keep if not c.startswith("stage")][:3], [cols[0]]]
    seen = []
    for t in tries:
        if t in seen:
            continue
        seen.append(t)
        try:
            hr, lo, hi, p = _fit(d, t)
            return dict(HR=hr, lo=lo, hi=hi, p=p, n=n, events=ev, adj="+".join(t[1:]) or "none")
        except Exception:
            continue
    return nan


def logrank_split(df):
    d = df.dropna(subset=["time", "event", "score"])
    d = d[d["time"] > 0]
    if len(d) < MIN_N or d["event"].sum() < MIN_EV:
        return np.nan
    hi = d["score"] > d["score"].median()
    return SV.logrank(d.loc[hi, "time"], d.loc[hi, "event"], d.loc[~hi, "time"], d.loc[~hi, "event"])


def make_score(Z, resp, res, weights=None):
    if weights is not None:
        f = [p for p in weights.index if p in Z.columns]
        s = Z[f].values @ weights[f].values
        s = pd.Series(s, index=Z.index)
    else:
        r = [p for p in resp if p in Z.columns]
        x = [p for p in res if p in Z.columns]
        s = (Z[r].mean(1) if r else 0) - (Z[x].mean(1) if x else 0)
        if not isinstance(s, pd.Series):
            return None
    return (s - s.mean()) / (s.std() or 1)


def tcga_covars(cdr, pats):
    c = cdr.reindex(pats)
    st = c["stage"].map(_stage)
    out = pd.DataFrame({"age": c["age"].astype(float)}, index=pats)
    for s in ("II", "III", "IV"):
        out[f"stage_{s}"] = (st == s).astype(float)
    return out


def tcga_tests(row, long, cdr, Xc, proj, systemic, top):
    cancer, drug, _, resp, res, _, rule, _ = row
    out = {}
    if cancer not in set(proj):
        return out
    inC = list(proj.index[proj == cancer])
    Z = Xc.loc[inC]
    weights = None
    if resp == "TOP10":
        g = top[(top["cancer"] == cancer) & (top["drug"] == drug)]
        resp, res = list(g[g.direction == "response"].pathway), list(g[g.direction == "resistance"].pathway)
    if resp and str(resp[0]).startswith("SIG:"):
        weights = pd.read_csv(RESULTS / "bridge" / "cetuximab_signature.tsv", sep="\t", index_col=0).iloc[:, 0]
    if resp and str(resp[0]).startswith("GENE:"):
        out["note"] = "gene-level biomarker; TCGA expression not cached (pathway scores only)"
        return out
    if not resp and not res:
        return out
    score = make_score(Z, resp, res, weights)
    g = long[(long["project"] == cancer) & (long["drug"] == drug)]
    treated = sorted(set(g["patient"]) & set(inC))
    if rule == "recist":
        disc = set(g[g["recist"].notna()]["patient"])
    elif rule == "firstcourse":
        first = g[g["setting_proxy"].eq("early")]["patient"].unique()
        oc = cdr.reindex(first)["treatment_outcome_first_course"]
        disc = set(oc[oc.isin(CR_LIKE | NR_LIKE)].index)
    else:
        disc = set()
    held = sorted(set(treated) - disc)
    untreated = sorted(set(inC) - systemic)
    out.update({"tcga_n_treated": len(treated), "tcga_n_discovery": len(disc & set(treated)), "tcga_n_heldout": len(held),
                "tcga_n_untreated": len(untreated)})
    C = tcga_covars(cdr, inC)
    for ep in ("OS", "PFI"):
        c = cdr.reindex(inC)
        E = pd.DataFrame({"time": c[f"{ep}.time"].astype(float) / 30.44, "event": c[ep].astype(float)}, index=inC)
        base = pd.concat([E, C], axis=1).assign(score=score)
        adj = ["score", "age", "stage_II", "stage_III", "stage_IV"]
        h = cox(base.loc[held], adj) if held else cox(base.iloc[:0], adj)
        out.update({f"{ep}_heldout_{k}": v for k, v in h.items()})
        out[f"{ep}_heldout_logrank_p"] = logrank_split(base.loc[held]) if held else np.nan
        if (h["n"] < MIN_N or h["events"] < MIN_EV) and treated:
            a = cox(base.loc[treated], adj)
            out.update({f"{ep}_alltreated_{k}": v for k, v in a.items()})
        u = cox(base.loc[untreated], adj)
        out.update({f"{ep}_untreated_{k}": v for k, v in u.items()})
        both = pd.concat([base.loc[held].assign(treated=1.0), base.loc[untreated].assign(treated=0.0)])
        both["sxt"] = both["score"] * both["treated"]
        i = cox(both, ["sxt", "score", "treated", "age", "stage_II", "stage_III", "stage_IV"])
        out.update({f"{ep}_interaction_HR": i["HR"], f"{ep}_interaction_p": i["p"]})
    return out


def external_tests(row, cache, top):
    cancer, drug, _, resp, res, _, _, ext = row
    rows = []
    for tag, tsel, csel, label in ext:
        if tag not in cache:
            continue
        expr, clin, sc = cache[tag]
        if resp == "TOP10":
            g = top[(top["cancer"] == cancer) & (top["drug"] == drug)]
            r_, x_ = list(g[g.direction == "response"].pathway), list(g[g.direction == "resistance"].pathway)
        else:
            r_, x_ = resp, res
        Z = sc.T.loc[clin.index]
        Z = ((Z - Z.mean()) / Z.std().replace(0, 1)).fillna(0)
        s = make_score(Z, r_, x_)
        t = tsel(clin).fillna(False).astype(bool)
        ctl = csel(clin).fillna(False).astype(bool) if csel is not None else pd.Series(False, index=clin.index)
        for ep in ("os", "pfs", "rfs", "dfs"):
            if f"{ep}_months" not in clin or f"{ep}_event" not in clin or clin[f"{ep}_months"].notna().sum() < MIN_N:
                continue
            d = pd.DataFrame({"time": clin[f"{ep}_months"].astype(float), "event": clin[f"{ep}_event"].astype(float), "score": s,
                              "age": clin["age"].astype(float) if "age" in clin else np.nan})
            cols = ["score"] + (["age"] if d["age"].notna().sum() > MIN_N else [])
            tr = cox(d[t], cols)
            rec = {"cohort": label, "endpoint": ep.upper(), "independent": f"{tag}|{cancer}" not in DISCOVERY_COHORTS,
                   **{f"treated_{k}": v for k, v in tr.items()}, "treated_logrank_p": logrank_split(d[t])}
            if ctl.sum() >= 10:
                b = d[t | ctl].assign(arm=t[t | ctl].astype(float))
                b["sxa"] = b["score"] * b["arm"]
                ii = cox(b, ["sxa", "score", "arm"] + cols[1:])
                rec.update({"interaction_HR": ii["HR"], "interaction_p": ii["p"], "control_n": int(ctl.sum())})
            rows.append(rec)
    return rows


def _classify(tests):
    if not tests:
        return "not testable"
    if any(h > 1 and p < 0.05 for h, p, k in tests if k == "treated"):
        return "contradicted"
    good_t = any(h < 1 and p < 0.05 for h, p, k in tests if k == "treated")
    good_i = any(h < 1 and p < 0.05 for h, p, k in tests if k == "interaction")
    if good_t and good_i:
        return "predictive"
    if good_i:
        return "predictive (interaction only)"
    if good_t:
        return "survival-associated"
    if any(h < 1 and p < 0.2 for h, p, _ in tests):
        return "trend"
    return "not supported"


def verdict(r, ext):
    """(TCGA held-out verdict, external verdict, combined)."""
    tt = []
    for ep in ("OS", "PFI"):
        if pd.notna(r.get(f"{ep}_heldout_HR")):
            tt.append((r[f"{ep}_heldout_HR"], r[f"{ep}_heldout_p"], "treated"))
            if pd.notna(r.get(f"{ep}_interaction_HR")):
                tt.append((r[f"{ep}_interaction_HR"], r[f"{ep}_interaction_p"], "interaction"))
    et = []
    for e in ext:
        if not e["independent"]:
            continue
        if pd.notna(e.get("treated_HR")):
            et.append((e["treated_HR"], e["treated_p"], "treated"))
        if pd.notna(e.get("interaction_HR")):
            et.append((e["interaction_HR"], e["interaction_p"], "interaction"))
    vt, ve = _classify(tt), (_classify(et) if et else "no cohort")
    pos = ("predictive", "predictive (interaction only)", "survival-associated")
    if vt in pos and ve in pos:
        comb = "VALIDATED (TCGA + external)"
    elif ve in pos and vt in ("not testable", "trend", "not supported"):
        comb = "validated externally" if vt != "not supported" else "external only (TCGA not supported)"
    elif vt in pos and ve == "no cohort":
        comb = "TCGA held-out only (needs external)"
    elif vt in pos:
        comb = "TCGA only; NOT replicated externally"
    elif "contradicted" in (vt, ve):
        comb = "contradicted"
    elif "trend" in (vt, ve):
        comb = "trend only"
    elif vt == "not testable" and ve in ("no cohort", "not testable"):
        comb = "not testable"
    else:
        comb = "not supported"
    return vt, ve, comb


def short(p):
    return ", ".join(x.replace("REACTOME_", "").replace("GENE:", "").replace("SIG:", "signature: ").replace("_", " ").lower()[:45] for x in p) \
        if isinstance(p, list) else str(p)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    long = R.load_response_labels()
    long = long[long["therapy_type"].isin(THERAPY) & ~long["drug"].isin(SUPPORTIVE) & long["drug"].notna()]
    cdr = R.load_cdr()
    cancers = sorted({r[0] for r in ROWS})
    M, proj = RT.patient_matrix(cancers)
    Xc = T.center_within(M, proj)
    systemic = set(long["patient"])
    top = pd.concat([pd.read_csv(RESULTS / "atlas" / "cancer_specific" / "top10_by_cancer_drug.tsv", sep="\t"),
                     pd.read_csv(RESULTS / "atlas" / "cancer_specific" / "relaxed_top10.tsv", sep="\t")], ignore_index=True)
    loaders = dict(AE.cohorts())
    cache = {}
    for tag in {e[0] for r in ROWS for e in r[7]}:
        if tag in loaders:
            expr, clin = loaders[tag]()
            sc = AE.pathway_scores(tag.replace(":", "_"), expr)
            clin = clin.loc[[s for s in clin.index if s in sc.columns]]
            cache[tag] = (expr, clin, sc)
    wide, longrows = [], []
    for row in ROWS:
        cancer, drug, label, resp, res, evid, rule, _ = row
        if resp == "TOP10":
            g = top[(top["cancer"] == cancer) & (top["drug"] == drug)]
            rs, xs = list(g[g.direction == "response"].pathway), list(g[g.direction == "resistance"].pathway)
        else:
            rs, xs = resp, res
        r = {"cancer": cancer, "cancer_name": label, "drug": drug, "response_pathways": short(rs),
             "resistance_pathways": short(xs), "evidence_level": evid,
             "discovery": {"recist": "TCGA RECIST", "firstcourse": "TCGA first-course outcome", "gene": "external cohort",
                           "none": "external cohort" if (rs or xs) else "-"}[rule]}
        r.update(tcga_tests(row, long, cdr, Xc, proj, systemic, top))
        ext = external_tests(row, cache, top)
        for e in ext:
            longrows.append({"cancer": cancer, "drug": drug, **e})
        best = sorted([e for e in ext if e["independent"] and pd.notna(e.get("treated_HR"))], key=lambda e: e["treated_p"])
        if best:
            e = best[0]
            r.update({"ext_cohort": e["cohort"], "ext_endpoint": e["endpoint"], "ext_n": e["treated_n"], "ext_events": e["treated_events"],
                      "ext_HR": e["treated_HR"], "ext_p": e["treated_p"], "ext_interaction_HR": e.get("interaction_HR", np.nan),
                      "ext_interaction_p": e.get("interaction_p", np.nan), "ext_n_cohorts": len(ext)})
        inter = [e for e in ext if e["independent"] and pd.notna(e.get("interaction_HR"))]
        if inter:
            e = min(inter, key=lambda e: e["interaction_p"])
            r.update({"ext_best_interaction_cohort": e["cohort"], "ext_best_interaction_HR": e["interaction_HR"],
                      "ext_best_interaction_p": e["interaction_p"]})
        if rs or xs:
            r["verdict_tcga"], r["verdict_external"], r["verdict"] = verdict(r, ext)
        else:
            r["verdict_tcga"] = r["verdict_external"] = r["verdict"] = "no biomarker"
        wide.append(r)
        print(cancer, drug, r["verdict"], flush=True)
    W = pd.DataFrame(wide)
    W.to_csv(OUT / "FINAL_BIOMARKER_SURVIVAL.tsv", sep="\t", index=False)
    pd.DataFrame(longrows).to_csv(OUT / "final_external_long.tsv", sep="\t", index=False)
    try:
        with pd.ExcelWriter(OUT / "FINAL_BIOMARKER_SURVIVAL.xlsx") as xw:
            W.round(4).to_excel(xw, sheet_name="biomarkers", index=False)
            pd.DataFrame(longrows).round(4).to_excel(xw, sheet_name="external_cohorts", index=False)
    except Exception as exc:
        print("xlsx skipped", exc)
    write_md(W)


def f(x, d=2):
    return "–" if pd.isna(x) else (f"{x:.{d}g}" if abs(x) < 0.001 else f"{x:.{d}f}")


def write_md(W):
    L = ["# Final pathway biomarkers: survival validation (one table)", "",
         "Score = mean z(response pathways) - mean z(resistance pathways). HR per SD from Cox models adjusted for age "
         "and stage (TCGA) or age (external). **HR < 1 = higher score, longer survival.** Held-out = treated with the drug "
         "but not used for discovery. Interaction HR < 1 = the benefit is larger with the drug than without it "
         "(predictive). Log-rank = the original paper's median split. External = best independent cohort.", "",
         "| Cancer | Drug | Response ↑ | Resistance ↑ | Held-out n (ev) | OS HR [95% CI] p | log-rank p | PFI HR p | "
         "Untreated OS HR | Interaction HR p (OS / PFI) | External cohort: HR p; interaction | TCGA | External | Overall |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in W.itertuples(index=False):
        d = r._asdict()
        if d.get("verdict") == "no biomarker":
            L.append(f"| {d['cancer_name']} | {d['drug'].title()} | – | – | – | – | – | – | – | – | – | – | – | no biomarker |")
            continue
        os_ = f"{f(d.get('OS_heldout_HR'))} [{f(d.get('OS_heldout_lo'))}-{f(d.get('OS_heldout_hi'))}] p={f(d.get('OS_heldout_p'), 2)}" \
            if pd.notna(d.get("OS_heldout_HR")) else (f"all-treated* {f(d.get('OS_alltreated_HR'))} p={f(d.get('OS_alltreated_p'))}"
                                                     if pd.notna(d.get("OS_alltreated_HR")) else "–")
        pfi = f"{f(d.get('PFI_heldout_HR'))} p={f(d.get('PFI_heldout_p'))}" if pd.notna(d.get("PFI_heldout_HR")) else \
            (f"all-treated* {f(d.get('PFI_alltreated_HR'))} p={f(d.get('PFI_alltreated_p'))}" if pd.notna(d.get("PFI_alltreated_HR")) else "–")
        n = f"{int(d['tcga_n_heldout'])} ({int(d['OS_heldout_events'])})" if pd.notna(d.get("tcga_n_heldout")) and pd.notna(d.get("OS_heldout_events")) else "–"
        inter = f"{f(d.get('OS_interaction_HR'))} p={f(d.get('OS_interaction_p'))} / {f(d.get('PFI_interaction_HR'))} p={f(d.get('PFI_interaction_p'))}"
        ext = "–"
        if isinstance(d.get("ext_cohort"), str):
            ext = f"{d['ext_cohort']} ({d['ext_endpoint']}, n={int(d['ext_n'])}): {f(d['ext_HR'])} p={f(d['ext_p'])}"
            if pd.notna(d.get("ext_best_interaction_HR")):
                ext += f"; int {f(d['ext_best_interaction_HR'])} p={f(d['ext_best_interaction_p'])}"
        L.append(f"| {d['cancer_name']} | {d['drug'].title()} | {d['response_pathways']} | {d['resistance_pathways']} | {n} | {os_} | "
                 f"{f(d.get('OS_heldout_logrank_p'))} | {pfi} | {f(d.get('OS_untreated_HR'))} | {inter} | {ext} | {d['verdict_tcga']} | {d['verdict_external']} | **{d['verdict']}** |")
    L += ["", "\\* all treated patients (includes discovery patients; not independent) - shown only when the held-out set is too small.", "",
          "Verdict counts: " + ", ".join(f"{k} {v}" for k, v in W["verdict"].value_counts().items())]
    (OUT / "FINAL_BIOMARKER_SURVIVAL.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[-1:]))


if __name__ == "__main__":
    main()
