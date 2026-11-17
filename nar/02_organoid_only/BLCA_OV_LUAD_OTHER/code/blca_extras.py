"""Bladder cisplatin: the paper's extra analyses, run after run_blca.py.

 1. the pathway the paper actually reported (amino acid synthesis and interconversion, rank 2 in the paper) - KM;
 2. paper's own bootstrap table (MOESM3 'Figure 4f') re-counted;
 3. isogenic T24 sensitive vs resistant cells (MOESM3 'Figure 4c' NES; Fig 4c + Supp Table 3);
 4. ERCC2 mutation concordance (Fig 5b; Supp Fig 13a-e comparison feature sets), TCGA open somatic mutations;
 5. comparison feature sets in cisplatin-treated patients (Fig 3b-d: whole transcriptome, whole pathway, 1st neighbours);
 6. drug-target genes alone (Supp Fig 8b); cross-test with the colorectal biomarker (Supp Fig 4b);
 7. clinical balance of the predicted groups (Supp Table 2);
 8. GDSC cell-line-trained models (Supp Fig 15c,d);
 9. our multi-pathway p-values next to the paper's source data (MOESM3 'FigureS5');
10. PCA of organoids vs patients before / after scaling (Supp Fig 1c,d,g,h).
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2_contingency, mannwhitneyu, ttest_ind  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import common as K  # noqa: E402
import engine as E  # noqa: E402
import run_blca as B  # noqa: E402

OUT = K.PARENT / "BLCA_lee2018"
FIG, TAB = OUT / "extras" / "figures", OUT / "extras" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
AA = "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION"
BH3 = "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS"
TARGETS = ["A2M", "ATOX1", "MPG", "TF"]
LOG = {}


def save(fig, name, caption):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    open(FIG / f"{name}.caption.md", "w").write(caption)


def zcoef(X, y):
    y = (y - y.mean()) / y.std()
    return RidgeCV(cv=3, alphas=np.arange(0.1, 1, 0.1)).fit(X, y).coef_


def setup(version):
    osc, tsc, surv = B.paper_version() if version == "paper" else B.current_version()
    y = B.response("CISPLATIN")
    z, _ = K.proximity("CISPLATIN")
    uni = sorted(set(osc.index) & set(tsc.index))
    F = [p for p in z.index[z <= -1.2816] if p in uni]
    S = [s for s in y.index if s in osc.columns]
    X = pd.DataFrame(StandardScaler().fit_transform(osc.loc[F, S].T), index=S, columns=F)
    c = pd.Series(zcoef(X.values, y.loc[S].values), index=F)
    pats = [p for p in tsc.columns if p in surv.index]
    P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[uni, pats].T), index=pats, columns=uni)
    tp = [p for p in pats if p in K.treated("BLCA", "CISPLATIN")]
    up = [p for p in pats if p not in K.any_drug("BLCA")]
    return dict(osc=osc, tsc=tsc, surv=surv, y=y.loc[S], F=F, X=X, c=c, P=P, tp=tp, up=up, uni=uni, S=S)


def TE(d, pats):
    return np.clip(d["surv"].loc[pats, "T"].values.astype(float), 0, None), d["surv"].loc[pats, "E"].values.astype(int)


# ---------------------------------------------------------------- 1 named biomarker
def named_biomarker(D):
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    rows = []
    for i, v in enumerate(("paper", "current")):
        d = D[v]
        rank = int((d["c"].abs() > abs(d["c"][AA])).sum() + 1)
        for j, (grp, pats) in enumerate((("treated", d["tp"]), ("untreated", d["up"]))):
            T, Ev = TE(d, pats)
            r, g = E.surv_test(d["c"][AA] * d["P"].loc[pats, AA].values, T, Ev)
            E.km_draw(axes[i, j], T, Ev, g, f"{v} data, {'cisplatin-treated' if grp == 'treated' else 'no drug record'}\n"
                      f"amino acid synthesis (our Ridge rank {rank})", r)
            E.panel_letter(axes[i, j], "abcd"[2 * i + j])
            rows.append(dict(version=v, group=grp, ridge_rank_of_AA=rank, coef=d["c"][AA], **r))
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "paper_biomarker_amino_acid_KM.tsv", sep="\t", index=False)
    fig.tight_layout()
    cap = ["# The pathway the paper reported for bladder cisplatin", "",
           "The paper reports 'amino acid synthesis and interconversion (transamination)' as the bladder cisplatin biomarker.",
           "In the paper's own coefficient plot (Fig. 2f) this pathway is rank 2, not rank 1, although the paper's code tests rank 1.",
           "Here the same pathway is used directly, with the sign of our organoid Ridge coefficient. Rows: paper source-data NES (top), current GDC data (bottom).", ""]
    for _, r in t.iterrows():
        cap.append(f"- {r.version}, {r.group}: our Ridge rank {r.ridge_rank_of_AA}; n = {r.n}, deaths = {r.events}, log-rank p = {r.p:.3g}, "
                   f"HR = {r.hr:.2f} ({r.hr_lo:.2f}-{r.hr_hi:.2f}); 5-year survival {r.surv5_resp:.2f} (pred. responders) vs {r.surv5_nonresp:.2f}.")
    cap.append("\nPaper: treated n = 77, p = 0.01; untreated n = 294, p = 0.066.")
    save(fig, "S_BLCA_paper_biomarker_KM", "\n".join(cap) + "\n")
    LOG["named_biomarker"] = t.to_dict("records")


# ---------------------------------------------------------------- 2 paper's bootstrap table
def paper_bootstrap():
    b = pd.read_csv(K.MOESM3 / "Figure_4f.tsv", sep="\t")
    tr = pd.to_numeric(b.top_rank, errors="coerce")
    out = {"iterations": len(b), "not_found": int(tr.isna().sum()), "rank1": int((tr == 1).sum()), "rank2": int((tr == 2).sum()),
           "empirical_p_rank<=1": float((tr <= 1).mean()), "empirical_p_rank<=2": float((tr <= 2).mean())}
    pd.DataFrame([out]).to_csv(TAB / "paper_bootstrap_Fig4f_recount.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    cnt = tr.dropna().astype(int).value_counts().sort_index()
    ax.bar(cnt.index, cnt.values, color=E.OK["sky"])
    ax.axvline(2, color=E.OK["verm"], lw=1.5, label="paper's biomarker: rank 2")
    ax.set_xlabel("Best rank meeting both criteria (random sets)")
    ax.set_ylabel("Iterations")
    ax.set_title(f"Paper's own bootstrap table (MOESM3 Fig 4f), recounted\nrank <= 2 in {out['rank1'] + out['rank2']} of {out['iterations']}: "
                 f"p = {out['empirical_p_rank<=2']:.4f}; none found in {out['not_found']}")
    ax.legend(frameon=False)
    save(fig, "S_BLCA_paper_bootstrap_recount",
         f"# Paper's bootstrap (Fig. 4f) recounted from its source data\n\nThe paper ran {out['iterations']} random pathway sets. "
         f"A set 'wins' if a pathway at a given coefficient rank predicts survival (p<0.05) and differs between resistant and sensitive "
         f"isogenic cells (p<0.05). The paper's biomarker reached rank 2. Random sets reaching rank 1: {out['rank1']}; rank 2: {out['rank2']}. "
         f"Empirical p (rank <= 2) = {out['empirical_p_rank<=2']:.4f} (paper: 0.014). No pathway met both criteria in {out['not_found']} sets.\n")
    LOG["paper_bootstrap"] = out


# ---------------------------------------------------------------- 3 isogenic
def isogenic(D):
    iso = B.isogenic()
    nes = iso["nes"]
    z, _ = K.proximity("CISPLATIN")
    zz = z[z.index.isin(nes.index)].sort_values()
    paths = [("paper biomarker", AA), ("our rank 1 (paper data)", D["paper"]["c"].abs().idxmax()),
             ("our rank 1 (current data)", D["current"]["c"].abs().idxmax())] + [(f"most proximal #{i + 1}", p) for i, p in enumerate(zz.index[:3])]
    rows = []
    fig, axes = plt.subplots(1, len(paths), figsize=(3.1 * len(paths), 4))
    for ax, (lab, p), L in zip(axes, paths, "abcdef"):
        a, b = nes.loc[p, iso["resistant"]].values, nes.loc[p, iso["sensitive"]].values
        t, pv = ttest_ind(b, a)
        rows.append(dict(label=lab, pathway=p, proximity_z=z.get(p, np.nan), mean_sensitive=b.mean(), mean_resistant=a.mean(),
                         t_sensitive_minus_resistant=t, p=pv, coef_paperdata=D["paper"]["c"].get(p, np.nan),
                         coef_currentdata=D["current"]["c"].get(p, np.nan)))
        ax.boxplot([b, a], widths=0.5)
        ax.scatter(np.r_[np.ones(3), 2 * np.ones(3)], np.r_[b, a], color=[E.OK["blue"]] * 3 + [E.OK["verm"]] * 3, zorder=3)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["T24 sensitive", "T24 resistant"])
        ax.set_title(f"{lab}\n{E.short(p, 30)}\nt-test p = {pv:.2g}", fontsize=8)
        ax.set_ylabel("ssGSEA NES")
        E.panel_letter(ax, L)
    fig.tight_layout()
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "isogenic_T24_tests.tsv", sep="\t", index=False)
    cap = ["# Isogenic cisplatin-sensitive vs -resistant bladder cells (paper Fig. 4c, Supp. Table 3)", "",
           "Pathway scores (ssGSEA NES) of 3 cisplatin-sensitive (T24S) and 3 resistant (T24R) samples, taken from the paper's source data",
           "(the raw data were a personal communication and are not public). Two-sided two-sample t-test, as in the paper.", ""]
    for _, r in t.iterrows():
        cap.append(f"- {r.label}: {E.short(r.pathway, 70)}; sensitive mean {r.mean_sensitive:.3f}, resistant {r.mean_resistant:.3f}, p = {r.p:.3g}. "
                   f"Organoid Ridge coefficient (paper data) {r.coef_paperdata:.3f}; a negative coefficient means higher activity = more sensitive.")
    cap.append("\nPaper: amino acid pathway lower in resistant cells, p = 0.00022.")
    save(fig, "S_BLCA_isogenic_T24", "\n".join(cap) + "\n")
    LOG["isogenic"] = t.to_dict("records")


# ---------------------------------------------------------------- gene-level data (current)
def gene_data():
    oexpr, _ = K.load_preclinical("bladder_lee2018")
    texpr = K.tcga_expr("BLCA")
    genes = oexpr.index.intersection(texpr.index)
    return oexpr.loc[genes], texpr.loc[genes]


def neighbours(k):
    g = K.paper_network()
    seen = set(t for t in TARGETS if t in g)
    front = set(seen)
    for _ in range(k):
        front = {n for f in front for n in g.neighbors(f)} - seen
        seen |= front
    return sorted(seen)


def feature_models(D, oexpr, texpr):
    """Return {name: (patients x 1 predicted resistance over ALL BLCA patients, n features)}."""
    y = D["current"]["y"]
    S = list(y.index)
    out = {}
    pats = list(texpr.columns)

    def ridge_pred(Xo, Xt):
        Xo = StandardScaler().fit_transform(Xo)
        Xt = StandardScaler().fit_transform(Xt)
        return Xt @ zcoef(Xo, y.values)
    var = oexpr.var(axis=1)
    g = oexpr.index[var > 0]
    out["whole transcriptome"] = (pd.Series(ridge_pred(oexpr.loc[g, S].T.values, texpr.loc[g, pats].T.values), index=pats), len(g))
    for v in ("paper", "current"):
        d = D[v]
        u = d["uni"]
        out[f"whole pathway ({v} data)"] = (pd.Series(ridge_pred(d["osc"].loc[u, S].T.values, d["tsc"].loc[u, d["tsc"].columns].T.values),
                                                      index=d["tsc"].columns), len(u))
    for k in (1, 2, 3):
        gg = [x for x in neighbours(k) if x in g]
        out[f"neighbours up to {k}"] = (pd.Series(ridge_pred(oexpr.loc[gg, S].T.values, texpr.loc[gg, pats].T.values), index=pats), len(gg))
    for v in ("paper", "current"):
        d = D[v]
        j = d["c"].abs().idxmax()
        full = pd.DataFrame(StandardScaler().fit_transform(d["tsc"].loc[[j, AA]].T), index=d["tsc"].columns, columns=[j, AA])
        out[f"network rank 1 ({v} data)"] = (d["c"][j] * full[j], 1)
        out[f"paper biomarker AA ({v} data)"] = (d["c"][AA] * full[AA], 1)
    return out


# ---------------------------------------------------------------- 4 ERCC2
def ercc2(models):
    m = pd.read_csv(K.RAW / "TCGA-BLCA.somaticmutation_wxs.tsv.gz", sep="\t", usecols=["sample", "gene", "effect"])
    m = m[m["sample"].str[13:15] == "01"]
    m["patient"] = m["sample"].str[:12]
    profiled = set(m.patient)
    silent = {"synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "upstream_gene_variant",
              "downstream_gene_variant", "intergenic_region", "non_coding_transcript_exon_variant"}
    mut = set(m.loc[(m.gene == "ERCC2") & ~m.effect.isin(silent), "patient"])
    rows = []
    names = list(models)
    fig, axes = plt.subplots(2, (len(names) + 1) // 2, figsize=(3.0 * ((len(names) + 1) // 2), 8.5))
    for ax, n, L in zip(axes.ravel(), names, "abcdefghijklmn"):
        s, nf = models[n]
        s = s[s.index.isin(profiled)]
        a, b = s[s.index.isin(mut)].values, s[~s.index.isin(mut)].values
        p = mannwhitneyu(a, b, alternative="less").pvalue
        rows.append(dict(feature_set=n, n_features=nf, n_mut=len(a), n_wt=len(b), median_mut=np.median(a), median_wt=np.median(b), p_one_sided=p))
        ax.boxplot([b, a], widths=0.5, showfliers=False)
        ax.scatter(1 + 0.08 * np.random.default_rng(0).standard_normal(len(b)), b, s=4, color=E.OK["grey"], alpha=0.5)
        ax.scatter(2 + 0.08 * np.random.default_rng(1).standard_normal(len(a)), a, s=8, color=E.OK["verm"])
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"WT\n(n={len(b)})", f"ERCC2 mut\n(n={len(a)})"])
        ax.set_title(f"{n}\n({nf} features)\none-sided MWU p = {p:.3g}", fontsize=8)
        ax.set_ylabel("Predicted cisplatin resistance")
        E.panel_letter(ax, L)
    for ax in axes.ravel()[len(names):]:
        ax.set_axis_off()
    fig.tight_layout()
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "ERCC2_mutation_concordance.tsv", sep="\t", index=False)
    cap = ["# ERCC2 mutation and predicted cisplatin resistance (paper Fig. 5b and Supp. Fig. 13)", "",
           "ERCC2 mutations make bladder tumours more sensitive to cisplatin in clinical studies. Each panel scores every TCGA-BLCA",
           "patient with mutation data (open GDC somatic mutations, Xena GDC hub; non-silent ERCC2 variants = mutant) using a model trained",
           "on the bladder organoids' cisplatin response, then asks whether mutant tumours get LOWER predicted resistance",
           "(one-sided Mann-Whitney U test, as in the paper).", ""]
    for _, r in t.iterrows():
        cap.append(f"- {r.feature_set} ({r.n_features} features): mutant n = {r.n_mut}, WT n = {r.n_wt}, p = {r.p_one_sided:.3g}.")
    cap.append("\nPaper: ERCC2 mut n = 37, WT n = 367; own model p = 0.002; whole transcriptome 0.04, whole pathway 0.037, "
               "neighbours 0.017 / 0.021 / 0.028.")
    save(fig, "S_BLCA_ERCC2_concordance", "\n".join(cap) + "\n")
    LOG["ercc2"] = t.to_dict("records")


# ---------------------------------------------------------------- 5 comparison feature sets in treated patients, 6 targets/cross
def comparisons(D, models, texpr):
    d = D["current"]
    rows = []
    sets = {k: v for k, v in models.items() if not k.startswith("paper biomarker")}
    # Supp Fig 4b cross-test: colorectal BH3 biomarker, sign from paper's COAD organoids
    o = pd.read_csv(K.MOESM3 / "Figure_2b.tsv", sep="\t", index_col=0)
    o.columns = [c.upper() for c in o.columns]
    ic = pd.read_csv(K.REPO / "data/organoid_COAD/drug_response/median_response.txt", sep="\t")
    ic = ic[ic["Drug Name"].str.upper() == "5-FLUOROURACIL"]
    ic = pd.Series(ic["median_IC50"].values, index=ic["Organoid"].str.upper().values)
    zc = K.paper_proximity()["FLUOROURACIL"]
    Fc = [p for p in zc.index[zc <= -1.2816] if p in o.index]
    Sc = [s for s in ic.index if s in o.columns]
    cc = pd.Series(zcoef(StandardScaler().fit_transform(o.loc[Fc, Sc].T), ic.loc[Sc].values), index=Fc)
    for v in ("paper", "current"):
        dd = D[v]
        full = pd.Series(StandardScaler().fit_transform(dd["tsc"].loc[[BH3]].T)[:, 0], index=dd["tsc"].columns)
        sets[f"cross-test: colorectal 5-FU biomarker BH3 ({v} data)"] = (cc[BH3] * full, 1)
    tex = texpr.loc[[t for t in TARGETS if t in texpr.index]]
    for t in tex.index:
        sets[f"target gene {t} (high = resistant)"] = (pd.Series(tex.loc[t].values, index=tex.columns), 1)
    n = len(sets)
    nc = 4
    fig, axes = plt.subplots(int(np.ceil(n / nc)), nc, figsize=(4.4 * nc, 4.3 * int(np.ceil(n / nc))))
    for ax, (name, (s, nf)), L in zip(axes.ravel(), sets.items(), "abcdefghijklmnopqrstuvwxyz"):
        ver = "paper" if "paper data" in name else "current"
        surv = D[ver]["surv"]
        pats = [p for p in s.index if p in surv.index and p in K.treated("BLCA", "CISPLATIN")]
        T = np.clip(surv.loc[pats, "T"].values.astype(float), 0, None)
        Ev = surv.loc[pats, "E"].values.astype(int)
        r, g = E.surv_test(s.loc[pats].values, T, Ev)
        E.km_draw(ax, T, Ev, g, f"{name} ({nf} feat.)", r)
        ax.legend(fontsize=6, frameon=False, loc="lower left")
        ax.title.set_fontsize(7.5)
        E.panel_letter(ax, L)
        rows.append(dict(model=name, n_features=nf, survival=("paper-era biotab" if ver == "paper" else "current GDC"), **r))
    for ax in axes.ravel()[n:]:
        ax.set_axis_off()
    fig.tight_layout()
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "comparison_feature_sets_treated_KM.tsv", sep="\t", index=False)
    cap = ["# Other ways of choosing features, cisplatin-treated TCGA-BLCA (paper Fig. 3b-d, Supp. Fig. 4b, 8b, 9)", "",
           "Each panel trains Ridge on the 9 bladder organoids with a different feature set (or uses one gene or the colorectal",
           "biomarker), scores the cisplatin-treated patients and splits them at the median.", ""]
    for _, r in t.iterrows():
        cap.append(f"- {r.model}: n = {r.n}, deaths = {r.events}, p = {r.p:.3g}, HR = {r.hr:.2f} ({r.hr_lo:.2f}-{r.hr_hi:.2f}).")
    cap.append("\nPaper (77 patients): whole transcriptome p = 0.79, whole pathway 0.82, first neighbours 0.37, 2nd/3rd neighbours 0.79/0.79, "
               "BH3 cross-test 0.86; target genes A2M 0.92, ATOX1 0.18, MPG 0.83, TF 0.48.")
    save(fig, "S_BLCA_comparison_feature_sets", "\n".join(cap) + "\n")
    LOG["comparisons"] = t.to_dict("records")


# ---------------------------------------------------------------- 7 clinical balance
def clinical(D):
    cl = pd.read_csv(K.RAW / "nationwidechildrens.org_clinical_patient_blca.txt", sep="\t", dtype=str, skiprows=[1, 2]).set_index("bcr_patient_barcode")
    rows = []
    for v in ("paper", "current"):
        d = D[v]
        for lab, j in (("rank 1", d["c"].abs().idxmax()), ("paper biomarker AA", AA)):
            s = d["c"][j] * d["P"].loc[d["tp"], j]
            g = (s > s.median()).map({True: "non-responder", False: "responder"})
            c = cl.reindex(g.index)
            feats = {"age": pd.to_numeric(c["age_at_diagnosis"], errors="coerce").map(lambda a: np.nan if pd.isna(a) else ("<60" if a < 60 else ">=60")),
                     "sex": c["gender"].str.upper().where(c["gender"].str.upper().isin(["MALE", "FEMALE"])),
                     "T stage": c["ajcc_tumor_pathologic_pt"].str.extract(r"(T\d)")[0].map({"T1": "T1-T2", "T2": "T1-T2", "T3": "T3-T4", "T4": "T3-T4"}),
                     "N stage": c["ajcc_nodes_pathologic_pn"].str.extract(r"(N\d)")[0]}
            for f, x in feats.items():
                tab = pd.crosstab(x, g)
                p = chi2_contingency(tab)[1] if tab.shape[0] > 1 else np.nan
                for lvl, rr in tab.iterrows():
                    rows.append(dict(version=v, biomarker=lab, feature=f, level=lvl, responders=int(rr.get("responder", 0)),
                                     non_responders=int(rr.get("non-responder", 0)), chi2_p=p))
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "SuppTable2_clinical_balance.tsv", sep="\t", index=False)
    LOG["clinical"] = t.groupby(["version", "biomarker", "feature"]).chi2_p.first().reset_index().to_dict("records")


# ---------------------------------------------------------------- 8 GDSC cell lines
def gdsc(D):
    s = pd.read_csv(K.PRECLIN / "ctrpv2_ccle" / "samples.tsv", sep="\t")
    expr, _ = K.load_preclinical("ctrpv2_ccle")
    cos2ach = {int(c): a for c, a in zip(s.COSMICID, s["sample"]) if pd.notna(c)}
    g = pd.read_csv(K.RAW / "GDSC2_fitted.tsv.gz", sep="\t")
    g = g[g.DRUG_NAME == "Cisplatin"].assign(ach=lambda x: x.COSMIC_ID.map(cos2ach)).dropna(subset=["ach"])
    g = g[g.ach.isin(expr.columns)]
    texpr = K.tcga_expr("BLCA")
    osc, tsc = K.paired_scores("BLCA_gdsc_ccle_all", expr[sorted(set(g.ach))], texpr)
    z, _ = K.proximity("CISPLATIN")
    surv = K.survival_current("BLCA")
    uni = sorted(set(osc.index) & set(tsc.index))
    F = [p for p in z.index[z <= -1.2816] if p in uni]
    pats = [p for p in tsc.columns if p in surv.index]
    P = pd.DataFrame(StandardScaler().fit_transform(tsc.loc[F, pats].T), index=pats, columns=F)
    tp = [p for p in pats if p in K.treated("BLCA", "CISPLATIN")]
    T, Ev = np.clip(surv.loc[tp, "T"].values, 0, None), surv.loc[tp, "E"].values.astype(int)
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for i, (lab, sub) in enumerate((("all GDSC lines", g), ("bladder GDSC lines", g[g.TCGA_DESC == "BLCA"]))):
        y = sub.groupby("ach")["LN_IC50"].median()
        X = StandardScaler().fit_transform(osc.loc[F, y.index].T)
        c = pd.Series(zcoef(X, y.values), index=F)
        for k, j in enumerate(c.abs().sort_values(ascending=False).index[:3]):
            r, gg = E.surv_test(c[j] * P.loc[tp, j].values, T, Ev)
            E.km_draw(axes[i, k], T, Ev, gg, f"{lab} (n={len(y)}), rank {k + 1}\n{E.short(j, 40)}", r)
            axes[i, k].legend(fontsize=6.5, frameon=False, loc="lower left")
            E.panel_letter(axes[i, k], "abcdef"[3 * i + k])
            rows.append(dict(training=lab, n_lines=len(y), rank=k + 1, pathway=j, coef=c[j], **r))
    fig.tight_layout()
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "SuppFig15_GDSC_trained.tsv", sep="\t", index=False)
    cap = ["# Models trained on cancer cell lines instead of organoids (paper Supp. Fig. 15c,d)", "",
           "Same procedure, but the training data are GDSC2 cisplatin ln IC50 values of cell lines that also have CCLE (DepMap) RNA-seq.",
           "Top: all such lines; bottom: bladder lines only. Each panel splits cisplatin-treated TCGA-BLCA patients (current GDC data) by one of the top-3 pathways.", ""]
    for _, r in t.iterrows():
        cap.append(f"- {r.training}, rank {r['rank']}: {E.short(r.pathway, 70)}; n = {r.n}, p = {r.p:.3g}, HR = {r.hr:.2f} ({r.hr_lo:.2f}-{r.hr_hi:.2f}).")
    cap.append("\nPaper: 831 lines -> p = 0.102, 0.0092, 0.84; 14 bladder lines -> p = 0.13, 0.23, 0.0092 (77 patients).")
    save(fig, "S_BLCA_GDSC_trained", "\n".join(cap) + "\n")
    LOG["gdsc"] = t.to_dict("records")


# ---------------------------------------------------------------- 9 multi-pathway vs paper
def multik_vs_paper():
    pap = pd.read_csv(K.MOESM3 / "FigureS5.tsv", sep="\t")
    pap = pap[pap.cancer_type == "blca"][["ML", "num_pathways", "pvalue"]].rename(columns={"ML": "model", "num_pathways": "k", "pvalue": "p_paper"})
    rows = []
    for v in ("paper", "current"):
        f = OUT / v / "tables" / "BLCA_CISPLATIN_multipathway_k.tsv"
        m = pd.read_csv(f, sep="\t")
        m = m[m.group == "treated"][["model", "k", "p"]].rename(columns={"p": f"p_ours_{v}_data"})
        pap = pap.merge(m, on=["model", "k"], how="left")
    pap.to_csv(TAB / "multipathway_k_vs_paper_sourcedata.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, ml, L in zip(axes, E.MODELS, "abc"):
        d = pap[pap.model == ml]
        for col, c, lab in (("p_paper", E.OK["black"], "paper source data (IC50)"), ("p_ours_paper_data", E.OK["blue"], "ours, paper NES + AUC"),
                            ("p_ours_current_data", E.OK["verm"], "ours, current GDC")):
            ax.plot(d.k, -np.log10(d[col]), "o-", color=c, label=lab)
        ax.axhline(-np.log10(0.05), ls=":", color=E.OK["grey"])
        ax.set_title(ml)
        ax.set_xlabel("k (top pathways)")
        ax.set_ylabel("-log10 log-rank p")
        E.panel_letter(ax, L)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "S_BLCA_multipathway_vs_paper",
         "# Multi-pathway models (k = 2..10) next to the paper's source data (MOESM3 'FigureS5', Supp. Fig. 5b)\n\n"
         "Black: the paper's own p-values (organoid IC50, 2019 TCGA files). Blue: our run on the paper's NES with CoderData AUC. "
         "Red: current GDC data. Numbers in `extras/tables/multipathway_k_vs_paper_sourcedata.tsv`. Paper Ridge k = 2: p = 0.022.\n")
    LOG["multik_vs_paper"] = pap.to_dict("records")


# ---------------------------------------------------------------- 10 PCA
def pca(D, oexpr, texpr):
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    rows = []
    for i, (lab, o, t) in enumerate((("genes (current data)", oexpr, texpr),
                                     ("pathways (paper NES)", D["paper"]["osc"].loc[D["paper"]["uni"]], D["paper"]["tsc"].loc[D["paper"]["uni"]]))):
        for j, scaled in enumerate((False, True)):
            A, Bm = o.T.values, t.T.values
            if scaled:
                A, Bm = StandardScaler().fit_transform(A), StandardScaler().fit_transform(Bm)
            M = np.vstack([Bm, A])
            M = M[:, np.std(M, axis=0) > 0]
            pc = PCA(2).fit(M)
            Z = pc.transform(M)
            ax = axes[i, j]
            ax.scatter(Z[: len(Bm), 0], Z[: len(Bm), 1], s=5, color=E.OK["grey"], alpha=0.5, label=f"TCGA-BLCA (n={len(Bm)})")
            ax.scatter(Z[len(Bm):, 0], Z[len(Bm):, 1], s=25, color=E.OK["verm"], label=f"organoids (n={len(A)})")
            ev = pc.explained_variance_ratio_ * 100
            ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
            ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
            ax.set_title(f"{lab}, {'z-scored within each dataset' if scaled else 'not scaled'}")
            ax.legend(frameon=False, fontsize=7)
            E.panel_letter(ax, "abcd"[2 * i + j])
            rows.append(dict(level=lab, scaled=scaled, pc1=ev[0], pc2=ev[1]))
    fig.tight_layout()
    pd.DataFrame(rows).to_csv(TAB / "SuppFig1_PCA_variance.tsv", sep="\t", index=False)
    save(fig, "S_BLCA_PCA_batch",
         "# Organoids vs patients before and after scaling (paper Supp. Fig. 1c,d,g,h)\n\nRed: bladder organoids; grey: TCGA-BLCA tumours. "
         "Left: raw values; right: each gene/pathway z-scored within its own dataset (the paper's correction). "
         "Explained variance in `extras/tables/SuppFig1_PCA_variance.tsv`. Paper: genes 21.6%/6.37% -> 9.50%/6.87%; pathways 27.0%/17.8% -> 20.7%/18.8%.\n")
    LOG["pca"] = rows


def main():
    D = {v: setup(v) for v in ("paper", "current")}
    named_biomarker(D)
    paper_bootstrap()
    isogenic(D)
    oexpr, texpr = gene_data()
    models = feature_models(D, oexpr, texpr)
    ercc2(models)
    comparisons(D, models, texpr)
    clinical(D)
    multik_vs_paper()
    pca(D, oexpr, texpr)
    gdsc(D)
    json.dump(LOG, open(TAB / "extras_log.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
