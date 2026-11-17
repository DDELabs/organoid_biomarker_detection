"""Extra paper panels for the GI sets (run after run_all.py and plots.py):

  Fig4b / SuppTable3   isogenic HCT-8 5-FU-sensitive vs -resistant cells (GSE81008/GSE81005; authors' NES + recomputed)
  SuppTable1           clinical features of predicted responders vs non-responders (chi-squared), every set x drug
  regimen split        colon 5-FU: 5-FU with oxaliplatin (FOLFOX-type) vs 5-FU without oxaliplatin
  Fig5a                cetuximab model vs BRAF V600E in TCGA-COAD (one-sided Mann-Whitney), organoids and PDX
  SuppFig15            5-FU models trained on GDSC1 cell lines (all lines; colorectal lines) applied to TCGA-COAD
  SuppFig1             PCA of lab models and patients before / after per-dataset scaling (pathway level)
  SuppFig10            pathway size vs single-pathway log-rank p
"""
import json
import pickle
import sys
import textwrap
import urllib.request
import zipfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2_contingency, mannwhitneyu  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import iso  # noqa: E402
from common import (C, DATA, PAPER_CUTOFF, RAW, REPO, ROOT, SETS, authors_inputs, cdr, gene_sets, proximity,  # noqa: E402
                    scores, scoring, set_response, treatment_table)
from engine import Combo, MLS, surv_stats  # noqa: E402
from plots import (BLUE, GREY, ORANGE, VERM, GREEN, COHORT_LABEL, fmt_hr, fmt_p, letter, load_all, p_km, rank1,  # noqa: E402
                   save, short)

T_ROOT = ROOT / "tables"
F_ROOT = ROOT / "figures"


def res_dir(r, kind):
    return ROOT / r["meta"]["folder"] / kind


# --------------------------------------------------------------------------- Fig 4b
def fig4b(allr):
    out_rows = []
    versions = {"authors": iso.authors_nes(), "current": iso.current_nes()}
    targets = []
    for tag, r in allr.items():
        if r["iso_tab"] is not None:
            targets.append((tag, r))
    fig, axs = plt.subplots(len(targets), 4, figsize=(18, 4.2 * len(targets)))
    axs = np.atleast_2d(axs)
    L = iter("abcdefghijklmnopqrstuvwxyzABCDEFGH")
    caps = []
    for i, (tag, r) in enumerate(targets):
        m = r["meta"]
        nes = versions["authors" if m["cohort"] == "COAD_authors" else "current"]
        tab, sens, res = iso.iso_table(nes)
        co = r["coefs"]["Ridge"].abs().sort_values(ascending=False)
        top3 = list(co.index[:3])
        bh3 = "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS"
        show = [bh3] + [p for p in top3 if p != bh3][:3]
        for j, pw in enumerate(show[:4]):
            ax = axs[i, j]
            if pw not in nes.index:
                ax.axis("off")
                continue
            a, b = nes.loc[pw, sens].values, nes.loc[pw, res].values
            ax.boxplot([a, b], widths=0.5, showfliers=False)
            ax.scatter(np.r_[np.ones(len(a)), 2 * np.ones(len(b))] + np.random.default_rng(0).uniform(-.08, .08, len(a) + len(b)),
                       np.r_[a, b], color=[BLUE] * len(a) + [VERM] * len(b), s=18, zorder=3)
            ax.set_xticks([1, 2])
            ax.set_xticklabels([f"5-FU-sensitive\nHCT-8 (n={len(a)})", f"5-FU-resistant\nHCT-8/5-FU (n={len(b)})"])
            ax.set_ylabel("ssGSEA score")
            rk = list(co.index).index(pw) + 1 if pw in co.index else None
            t = tab.loc[pw]
            ax.set_title(textwrap.fill(short(pw, 80), 40) + f"\n{'rank ' + str(rk) if rk else 'not in model'}; "
                         f"t-test p = {fmt_p(t.p)}", loc="left")
            letter(ax, next(L))
            coef = r["coefs"]["Ridge"].get(pw, np.nan)
            out_rows.append({"cohort": m["cohort"], "nes_version": "authors (MOESM3)" if m["cohort"] == "COAD_authors" else "recomputed (GSE81005)",
                             "pathway": pw, "ridge_rank": rk, "ridge_coef": coef, "proximity_z": r["z"].get(pw, np.nan),
                             "mean_sensitive": t.mean_sensitive, "mean_resistant": t.mean_resistant, "t": t.t, "p": t.p,
                             "direction_matches_model": (bool((t.mean_sensitive > t.mean_resistant) == (coef < 0))
                                                         if np.isfinite(coef) else None)})
        axs[i, 0].text(-0.35, 0.5, COHORT_LABEL.get(m["cohort"], m["cohort"]), transform=axs[i, 0].transAxes,
                       rotation=90, va="center", fontsize=10, fontweight="bold")
    fig.tight_layout()
    T = pd.DataFrame(out_rows)
    T.to_csv(ROOT / "COAD_vandewetering" / "tables" / "Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv", sep="\t", index=False)
    lines = "\n".join(f"- {r.cohort}: {short(r.pathway, 90)} (Ridge rank {r.ridge_rank}): sensitive mean {r.mean_sensitive:.4g}, "
                      f"resistant mean {r.mean_resistant:.4g}, p = {fmt_p(r.p)}; direction matches the organoid model: {r.direction_matches_model}"
                      for r in T.itertuples())
    save(fig, ROOT / "COAD_vandewetering" / "figures" / "Fig4b_SuppTable3_isogenic_HCT8_5FU", f"""# Fig. 4b (and Supp. Table 3): isogenic 5-FU-sensitive vs -resistant colorectal cells

HCT-8 colorectal cancer cells (5-FU-sensitive, n = 3) and HCT-8 cells made resistant by growing them in rising
5-FU doses (HCT-8/5-FU at 0, 24 and 48 h, n = 9), GEO GSE81008 (sub-series GSE81005), as in the paper. Each box
shows a pathway's ssGSEA score; dots are single arrays. p = two-sided two-sample Student's t-test (paper's test).
First column: the paper's biomarker (activation of BH3-only proteins). Other columns: the top Ridge pathways of the
model in that row. Row 1 uses the authors' own scores (Source Data sheet 'Figure 4b'); the other rows use scores
recomputed from GEO with this project's ssGSEA. "Direction matches" = the pathway is higher in sensitive cells
when a higher score predicts lower IC50 in organoids (negative Ridge coefficient).

{lines}
""")
    return T


# --------------------------------------------------------------------------- Supp Table 1
def clinical_tables(allr, cd):
    rows = []
    pat_auth = pd.read_csv(REPO / "data/TCGA_COAD/clinical_patient/4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt",
                           sep="\t").set_index("bcr_patient_barcode")
    for tag, r in allr.items():
        m = r["meta"]
        T, E, g, used = r["km"][("treated", 1)]
        pats = r["tp"]
        grp = pd.Series(np.where(g == 1, "non-responder", "responder"), index=pats)
        if m["cohort"] == "COAD_authors":
            age = pd.to_numeric(pat_auth.reindex(pats)["age_at_initial_pathologic_diagnosis"], errors="coerce")
            sex = pat_auth.reindex(pats)["gender"].str.upper()
            tst = pat_auth.reindex(pats)["ajcc_tumor_pathologic_pt"].astype(str).str[:2].map(
                {"T1": "T1-T2", "T2": "T1-T2", "T3": "T3-T4", "T4": "T3-T4"})
            nst = pat_auth.reindex(pats)["ajcc_nodes_pathologic_pn"].astype(str).str[:2].map({"N0": "N0", "N1": "N1", "N2": "N2"})
            feats = {"age": np.where(age < 60, "<60", ">=60"), "sex": sex, "T stage": tst, "N stage": nst}
        else:
            c = cd.reindex(pats)
            st = c["stage"].astype(str).str.upper().str.replace("STAGE ", "")
            stage = np.where(st.str.startswith("IV") | st.str.startswith("III"), "III-IV",
                             np.where(st.str.startswith("I"), "I-II", None))
            feats = {"age": np.where(c["age"] < 60, "<60", ">=60"), "sex": c["gender"].str.upper(), "stage": stage}
        for name, v in feats.items():
            v = pd.Series(v, index=pats).replace({"nan": np.nan, "None": np.nan})
            tab = pd.crosstab(v.dropna(), grp.loc[v.dropna().index])
            p = chi2_contingency(tab)[1] if tab.shape[0] > 1 and tab.shape[1] > 1 else np.nan
            for lev in tab.index:
                rows.append({"folder": m["folder"], "cohort": m["cohort"], "drug": m["model_drug"], "feature": name,
                             "level": lev, "responder": int(tab.loc[lev].get("responder", 0)),
                             "non_responder": int(tab.loc[lev].get("non-responder", 0)), "chi2_p": p})
    T = pd.DataFrame(rows)
    T.to_csv(T_ROOT / "SuppTable1_clinical_features_by_predicted_group.tsv", sep="\t", index=False)
    return T


# --------------------------------------------------------------------------- regimen split
def regimen_split(allr, cd):
    targets = [(t, r) for t, r in allr.items() if r["meta"]["patient_drug"] == "FLUOROURACIL"
               and r["meta"]["cohort"] in ("COAD_authors", "COAD", "READ", "COADREAD")]
    fig, axs = plt.subplots(len(targets), 3, figsize=(16, 4.3 * len(targets)))
    axs = np.atleast_2d(axs)
    L = iter("abcdefghijklmnopqrstuvwxyz")
    rows = []
    for i, (tag, r) in enumerate(targets):
        m = r["meta"]
        parts = ["COAD"] if m["cohort"] in ("COAD", "COAD_authors") else (["READ"] if m["cohort"] == "READ" else ["COAD", "READ"])
        tt = pd.concat([treatment_table(p) for p in parts])
        oxa = set(tt.loc[tt.drug == "OXALIPLATIN", "patient"])
        iri = set(tt.loc[tt.drug == "IRINOTECAN", "patient"])
        T, E, g, used = r["km"][("treated", 1)]
        pats = np.array(r["tp"])
        cls = np.array(["with oxaliplatin (FOLFOX-type)" if p in oxa else
                        ("with irinotecan, no oxaliplatin (FOLFIRI-type)" if p in iri else "no oxaliplatin or irinotecan record")
                        for p in pats])
        for j, c in enumerate(["with oxaliplatin (FOLFOX-type)", "no oxaliplatin or irinotecan record",
                               "with irinotecan, no oxaliplatin (FOLFIRI-type)"]):
            ax = axs[i, j]
            sel = cls == c
            if sel.sum() < 3 or g[sel].sum() in (0, sel.sum()):
                ax.axis("off")
                ax.text(0.5, 0.5, f"{c}: n = {sel.sum()} (too few to plot)", ha="center", transform=ax.transAxes)
                rows.append({"cohort": m["cohort"], "subgroup": c, "n": int(sel.sum())})
                continue
            st = p_km(ax, T[sel], E[sel], g[sel], f"{COHORT_LABEL.get(m['cohort'], m['cohort'])}\n5-FU {c}")
            letter(ax, next(L))
            rows.append({"cohort": m["cohort"], "subgroup": c, **st})
    fig.tight_layout()
    Tb = pd.DataFrame(rows)
    Tb.to_csv(ROOT / "COAD_vandewetering" / "tables" / "regimen_split_5FU_FOLFOX.tsv", sep="\t", index=False)
    lines = "\n".join(f"- {r.cohort}, {r.subgroup}: n = {r.n}" + (f", deaths = {int(r.events)}, p = {fmt_p(r.p)}, {fmt_hr(r.HR, r.HR_lo, r.HR_hi)}"
                                                                 if 'events' in Tb and np.isfinite(r.events) else "")
                      for r in Tb.itertuples())
    save(fig, ROOT / "COAD_vandewetering" / "figures" / "Fig2d_regimen_split_5FU_vs_FOLFOX", f"""# 5-FU-treated colorectal patients split by regimen

The paper's "5-FU-treated" patients are everyone with a 5-FU record. Most of them also received oxaliplatin
(FOLFOX). Here the same predicted groups (rank-1 Ridge pathway, median split over all 5-FU-treated patients) are
shown separately for patients who also had oxaliplatin, patients who had neither oxaliplatin nor irinotecan
(closest to "5-FU alone", often with leucovorin or capecitabine), and patients with irinotecan but no oxaliplatin.

{lines}
""")
    return Tb


# --------------------------------------------------------------------------- Fig 5a
def gdc_braf_v600e(project):
    f = DATA / f"braf_v600e_{project}.json"
    if f.exists():
        return json.loads(f.read_text())
    def post(endpoint, payload):
        req = urllib.request.Request(f"https://api.gdc.cancer.gov/{endpoint}", data=json.dumps(payload).encode(),
                                     headers={"Content-Type": "application/json"})
        return json.loads(urllib.request.urlopen(req, timeout=120).read())
    filt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "case.project.project_id", "value": [f"TCGA-{project}"]}},
        {"op": "in", "content": {"field": "ssm.consequence.transcript.gene.symbol", "value": ["BRAF"]}},
        {"op": "in", "content": {"field": "ssm.consequence.transcript.aa_change", "value": ["V600E"]}}]}
    d = post("ssm_occurrences", {"filters": filt, "fields": "case.submitter_id", "size": 2000})
    mut = sorted({h["case"]["submitter_id"] for h in d["data"]["hits"]})
    filt2 = {"op": "and", "content": [
        {"op": "in", "content": {"field": "project.project_id", "value": [f"TCGA-{project}"]}},
        {"op": "in", "content": {"field": "summary.experimental_strategies.experimental_strategy", "value": ["WXS"]}},
        {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Somatic Mutation"]}}]}
    d2 = post("cases", {"filters": filt2, "fields": "submitter_id", "size": 5000})
    allc = sorted({h["submitter_id"] for h in d2["data"]["hits"]})
    out = {"v600e": mut, "profiled": allc}
    f.write_text(json.dumps(out))
    return out


def fig5a(allr, cd):
    braf = gdc_braf_v600e("COAD")
    mutset, prof = set(braf["v600e"]), set(braf["profiled"])
    variants = []
    # organoids (paper: 3 hypermutated organoids removed before training), all organoids, and PDX
    for key, cohort, drop, label in (("vandewetering", "COAD_authors", ["P19A", "P19B", "P24A"], "organoids, authors' NES, 3 hypermutated removed (paper)"),
                                     ("vandewetering", "COAD_authors", [], "organoids, authors' NES, all 19"),
                                     ("vandewetering", "COAD", ["P19TA", "P19TB", "P24TA"], "organoids, current data, 3 hypermutated removed"),
                                     ("vandewetering", "COAD", [], "organoids, current data, all"),
                                     ("isella", "COAD", [], "PDX (Isella), current data")):
        z, _ = proximity("CETUXIMAB")
        if cohort == "COAD_authors":
            org, resp, pat, sur, treated, anydrug = authors_inputs("CETUXIMAB")
            osc, tsc, y = org, pat, resp
        else:
            import run_all as A
            osc, tsc, sur, tt, _ = A.cohort_inputs(key, "COAD", cd)
            y = set_response(key)["CETUXIMAB"]
        y = y.drop([d for d in drop if d in y.index])
        allpats = list(tsc.columns)
        surall = pd.DataFrame({"months": 1.0, "event": 0}, index=allpats)
        c = Combo(osc, y, tsc, surall, set(), set(), z)
        cf = c.coefs()["Ridge"]
        i = int(np.argmax(np.abs(cf)))
        pw = c.feats[i]
        score = cf[i] * c.P[pw]
        pats = [p for p in score.index if p in prof]
        mut = [p for p in pats if p in mutset]
        wt = [p for p in pats if p not in mutset]
        p1 = mannwhitneyu(score[mut], score[wt], alternative="greater").pvalue if mut else np.nan
        variants.append(dict(label=label, n_models=len(c.samples), pathway=pw, coef=cf[i], n_mut=len(mut), n_wt=len(wt),
                             p_one_sided=p1, mut=score[mut].values, wt=score[wt].values))
    fig, axs = plt.subplots(1, len(variants), figsize=(4.2 * len(variants), 5))
    for ax, v, L in zip(axs, variants, "abcde"):
        ax.boxplot([v["wt"], v["mut"]], widths=0.5, showfliers=False)
        rng = np.random.default_rng(0)
        ax.scatter(1 + rng.uniform(-.12, .12, len(v["wt"])), v["wt"], s=5, color=GREY, alpha=.6)
        ax.scatter(2 + rng.uniform(-.12, .12, len(v["mut"])), v["mut"], s=8, color=VERM)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"BRAF wild-type\n(n={v['n_wt']})", f"BRAF V600E\n(n={v['n_mut']})"])
        ax.set_ylabel("Predicted cetuximab resistance (coef x pathway score)")
        ax.set_title(textwrap.fill(v["label"], 34) + "\n" + textwrap.fill(short(v["pathway"], 70), 34) +
                     f"\none-sided p = {fmt_p(v['p_one_sided'])}", loc="left")
        letter(ax, L)
    fig.tight_layout()
    T = pd.DataFrame([{k: v for k, v in d.items() if k not in ("mut", "wt")} for d in variants])
    T.to_csv(ROOT / "COAD_vandewetering" / "tables" / "Fig5a_cetuximab_BRAF_V600E.tsv", sep="\t", index=False)
    lines = "\n".join(f"- {r.label}: {r.n_models} models; pathway {short(r.pathway, 80)}; V600E n = {r.n_mut}, wild-type n = {r.n_wt}; "
                      f"one-sided Mann-Whitney p = {fmt_p(r.p_one_sided)}" for r in T.itertuples())
    save(fig, ROOT / "COAD_vandewetering" / "figures" / "Fig5a_cetuximab_BRAF_V600E", f"""# Fig. 5a: cetuximab model vs BRAF V600E in TCGA colon tumours

BRAF V600E is a known marker of cetuximab resistance. A Ridge model was trained on cetuximab response of the lab
models using the cetuximab-proximal pathways; its rank-1 pathway gives each TCGA-COAD tumour a predicted-resistance
score (coefficient x scaled pathway score; higher = more resistant). BRAF V600E status comes from the GDC
open-access somatic mutation calls (cases with whole-exome mutation data). The test is the paper's one-sided
Mann-Whitney U (are V600E tumours predicted more resistant?). Paper: 47 V600E vs 328 wild-type, p = 0.037, using
organoids without the 3 hypermutated lines.

{lines}
""")
    return T


# --------------------------------------------------------------------------- Supp Fig 15
def gdsc_expression():
    f = DATA / "gdsc_5FU_expression.tsv.gz"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0)
    resp = pd.read_csv(RAW / "gdsc" / "GDSC1_5FU.tsv", sep="\t")
    want = set(resp["SANGER_MODEL_ID"])
    with zipfile.ZipFile(RAW / "gdsc" / "rnaseq_all_20220624.zip") as z:
        with z.open("rnaseq_tpm_20220624.csv") as fh:
            hdr = pd.read_csv(fh, nrows=0).columns
        keep = [0, 1] + [i for i, c in enumerate(hdr) if str(c).split(".")[0] in want]
        with z.open("rnaseq_tpm_20220624.csv") as fh:
            d = pd.read_csv(fh, usecols=keep, skiprows=[1, 2, 3, 4], low_memory=False)
    d = d.drop(columns=d.columns[0]).rename(columns={d.columns[1]: "gene"}).dropna(subset=["gene"])
    d = d.set_index("gene").apply(pd.to_numeric, errors="coerce")
    d.columns = [c.split(".")[0] for c in d.columns]
    d = d.T.groupby(level=0).mean().T                        # Sanger + Broad replicate columns -> mean
    d = np.log2(d.groupby(level=0).mean() + 1)
    d.to_csv(f, sep="\t")
    return d


def suppfig15(cd):
    resp = pd.read_csv(RAW / "gdsc" / "GDSC1_5FU.tsv", sep="\t")
    y_all = resp.groupby("SANGER_MODEL_ID")["LN_IC50"].median()
    tissue = resp.groupby("SANGER_MODEL_ID")["TCGA_DESC"].first()
    e = C.to_canonical(gdsc_expression())
    tag = "gdsc1_5FU_lines"
    texpr = C.tcga_star_fpkm_uq("COAD")
    genes = e.index.intersection(texpr.index)
    gs = {k: [g for g in v if g in set(genes)] for k, v in gene_sets(("REACTOME",)).items()}
    gs = {k: v for k, v in gs.items() if len(v) >= 5}
    lsc = scoring.cached(f"scores/{tag}_models", scoring.ssgsea, e.loc[genes], gs, threads=2)
    tsc = scoring.cached(f"scores/{tag}_COAD", scoring.ssgsea, texpr.loc[genes], gs, threads=2)
    import run_all as A
    _, _, sur, tt, untreated = A.cohort_inputs("vandewetering", "COAD", cd)
    treated = set(tt.loc[tt.drug == "FLUOROURACIL", "patient"])
    z, _ = proximity("FLUOROURACIL")
    org, resp_a, pat_a, sur_a, treated_a, anyd_a = authors_inputs("FLUOROURACIL")
    out, panels = [], []
    for label, ids in (("all GDSC1 lines", list(y_all.index)),
                       ("colorectal GDSC1 lines", [i for i in y_all.index if tissue.get(i) in ("COAD/READ", "COREAD")])):
        y = y_all.loc[[i for i in ids if i in lsc.columns]]
        for cohort, (tscx, surx, trx, unx) in (("COAD", (tsc, sur, treated, untreated)),
                                               ("COAD_authors", (pat_a, sur_a, treated_a, set(sur_a.index) - anyd_a))):
            c = Combo(lsc, y, tscx, surx, trx, unx, z)
            ev, cf = c.evaluate()
            for k in (1, 2, 3):
                # single pathway at rank k (paper's Supp Fig 15 shows the top 3 pathways one at a time)
                o = np.argsort(-np.abs(cf["Ridge"]), kind="stable")[k - 1]
                pw = c.feats[o]
                g = (cf["Ridge"][o] * c.P.loc[c.tp, pw].values > np.median(cf["Ridge"][o] * c.P.loc[c.tp, pw].values)).astype(int)
                st = surv_stats(c.T, c.E, g)
                out.append({"training": label, "n_lines": len(c.samples), "patients": cohort, "rank": k, "pathway": pw, **st})
                panels.append((label, cohort, k, pw, c.T, c.E, g, len(c.samples)))
    T = pd.DataFrame(out)
    T.to_csv(ROOT / "COAD_vandewetering" / "tables" / "SuppFig15_GDSC_cell_line_models_5FU.tsv", sep="\t", index=False)
    fig, axs = plt.subplots(4, 3, figsize=(16, 18))
    for ax, (label, cohort, k, pw, Tt, Et, g, nl), L in zip(axs.ravel(), panels, "abcdefghijkl"):
        p_km(ax, Tt, Et, g, f"{label} (n={nl}) -> {COHORT_LABEL[cohort]}\nrank {k}: " + textwrap.fill(short(pw, 70), 45))
        letter(ax, L)
    fig.tight_layout()
    lines = "\n".join(f"- {r.training} ({r.n_lines} lines) -> {r.patients}, rank {r.rank} {short(r.pathway, 80)}: n = {r.n}, "
                      f"deaths = {r.events}, p = {fmt_p(r.p)}, {fmt_hr(r.HR, r.HR_lo, r.HR_hi)}" for r in T.itertuples())
    save(fig, ROOT / "COAD_vandewetering" / "figures" / "SuppFig15_GDSC_cell_line_models_5FU", f"""# Supp. Fig. 15 (colon): 5-FU models trained on cancer cell lines instead of organoids

Same pipeline, but the Ridge model is trained on GDSC1 cell-line 5-FU response (ln IC50; GDSC release 8.5) with
cell-line RNA-seq (Cell Model Passports, log2 TPM+1), either all screened lines with RNA-seq or colorectal lines
only. The top three pathways are each used alone to split 5-FU-treated TCGA-COAD patients (current GDC data and the
authors' 2019 files). Paper: all 893 lines, top pathway p = 0.17; colorectal lines, top pathway p = 0.86.

{lines}
""")
    return T


# --------------------------------------------------------------------------- Supp Fig 1
def suppfig1():
    for key, (folder, name, projects, col, label) in SETS.items():
        sc = scores(key, projects)
        mods = sc["models"]
        pats = pd.concat([sc[p] for p in projects], axis=1)
        common = mods.index.intersection(pats.index)
        A = pd.concat([mods.loc[common], pats.loc[common]], axis=1).T
        lab = np.r_[np.ones(mods.shape[1]), np.zeros(pats.shape[1])]
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
        for ax, scaled, L in ((axs[0], False, "a"), (axs[1], True, "b")):
            X = A.values
            if scaled:
                X = np.vstack([StandardScaler().fit_transform(mods.loc[common].T.values),
                               StandardScaler().fit_transform(pats.loc[common].T.values)])
            pc = PCA(2).fit(X)
            Y = pc.transform(X)
            ax.scatter(Y[lab == 0, 0], Y[lab == 0, 1], s=6, color=GREY, alpha=.6, label=f"TCGA {'+'.join(projects)} (n={int((lab == 0).sum())})")
            ax.scatter(Y[lab == 1, 0], Y[lab == 1, 1], s=14, color=VERM, label=f"lab models (n={int(lab.sum())})")
            ev = pc.explained_variance_ratio_
            ax.set_xlabel(f"PC1 ({ev[0]:.0%})"); ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
            ax.set_title(("After" if scaled else "Before") + " per-dataset scaling (pathway scores)", loc="left")
            ax.legend(frameon=False, fontsize=8)
            letter(ax, L)
        fig.tight_layout()
        save(fig, ROOT / folder / "figures" / f"SuppFig1_PCA_models_vs_TCGA", f"""# Supp. Fig. 1: lab models vs patients, before and after scaling

Principal component analysis of the ssGSEA pathway scores ({len(common)} Reactome pathways) of the {int(lab.sum())}
lab models ({label}) and {int((lab == 0).sum())} TCGA tumours. (a) raw scores: models and tumours usually sit
apart. (b) after z-scoring every pathway within each dataset (what the method does before training and
prediction): the two groups overlap. This is the paper's batch-effect check.
""")


# --------------------------------------------------------------------------- Supp Fig 10
def suppfig10(allr):
    gsz = {k: len(v) for k, v in gene_sets(("REACTOME",)).items()}
    for tag, r in allr.items():
        m = r["meta"]
        co = r["coefs"]["Ridge"]
        sp = r["single_p"].loc[r["feats"]]
        p = np.where(co.loc[sp.index] > 0, sp.p_pos, sp.p_neg)
        size = np.array([gsz.get(f, np.nan) for f in sp.index])
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(size, -np.log10(p), color=BLUE, s=20)
        r2 = np.corrcoef(size, -np.log10(p))[0, 1] ** 2 if len(p) > 2 else np.nan
        ax.axhline(-np.log10(0.05), color="k", ls=":", lw=1)
        ax.set_xlabel("Pathway size (genes in MSigDB v6.1 Reactome set)")
        ax.set_ylabel("-log10 log-rank p, pathway used alone")
        ax.set_title(f"{m['model_drug'].lower()} | {COHORT_LABEL.get(m['cohort'], m['cohort'])}\nR2 = {r2:.3f}, "
                     f"{len(p)} proximal pathways", loc="left")
        fig.tight_layout()
        pd.DataFrame({"pathway": sp.index, "size": size, "p_alone": p}).to_csv(
            res_dir(r, "tables") / f"SuppFig10_size_vs_p_{m['cohort']}_{m['model_drug']}.tsv", sep="\t", index=False)
        save(fig, res_dir(r, "figures") / f"SuppFig10_size_vs_p_{m['cohort']}_{m['model_drug']}", f"""# Supp. Fig. 10: pathway size vs survival p ({m['model_drug'].lower()}, {COHORT_LABEL.get(m['cohort'], m['cohort'])})

Each dot is one proximal pathway used alone to split drug-treated patients at the median (coefficient sign from the
Ridge model). x = number of genes in the pathway; y = -log10 log-rank p. R2 = {r2:.3f} (squared Pearson
correlation). A value near 0 means pathway size does not drive the survival p-value (paper: R2 = 0.0026 for colon 5-FU).
""")


def main():
    allr = load_all()
    cd = cdr()
    what = sys.argv[1:] or ["fig4b", "clin", "regimen", "fig5a", "s1", "s10", "s15"]
    if "fig4b" in what:
        print(fig4b(allr).to_string())
    if "clin" in what:
        clinical_tables(allr, cd)
    if "regimen" in what:
        print(regimen_split(allr, cd).to_string())
    if "fig5a" in what:
        print(fig5a(allr, cd).to_string())
    if "s1" in what:
        suppfig1()
    if "s10" in what:
        suppfig10(allr)
    if "s15" in what:
        print(suppfig15(cd).to_string())


if __name__ == "__main__":
    main()
