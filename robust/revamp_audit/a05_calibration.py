"""Audit 5: effect-size calibration.

(a) Literature anchors for transcriptomic / tissue predictive and prognostic markers (hand-curated,
    PubMed-checked where a DOI is given) -> out/literature_effect_sizes.tsv
(b) Retrospective design analysis (Gelman & Carlin 2014) of the nominally significant TCGA held-out
    results in FINAL_BIOMARKER_SURVIVAL: given the events available and a plausible true per-SD HR
    (0.8 or 0.7), what is the power, the type-S (wrong sign) rate and the expected exaggeration
    of a significant estimate -> out/retrodesign_final_rows.tsv
"""
import numpy as np
import pandas as pd

from common import OUT, REPO

LIT = [
    # marker, setting, design, effect, per-SD-ish translation, source DOI
    ("21-gene Recurrence Score (Oncotype DX)", "ER+ N0 breast, tamoxifen +/- chemo (NSABP B-20, n=651)", "RCT, treatment x marker",
     "chemo RR 0.26 for RS>=31 vs 1.31 for RS<18; interaction p=0.038", "strongest classic transcriptomic predictive marker reached only p=0.038 with 651 pts", "10.1200/JCO.2005.04.7985"),
    ("21-gene RS, TAILORx", "ER+ HER2- N0 breast, RS 11-25 (n~6,700 randomised)", "RCT",
     "endocrine alone non-inferior to chemo-endocrine overall (iDFS HR ~1.08); benefit limited to age<=50, RS 16-25",
     "intermediate-range predictive effect near null despite huge n", "10.1056/NEJMoa1804710"),
    ("70-gene signature (MammaPrint), MINDACT", "early breast, clinical-high/genomic-low (n=1,550)", "RCT, treatment by genomic risk",
     "5-y DMFS 94.7% without vs ~96.2% with chemo (1.5-point absolute difference)", "mainly prognostic; predictive component small", "10.1056/NEJMoa1602253"),
    ("Stromal TILs (continuous)", "early TNBC, anthracycline+/-taxane, pooled IPD (n=2,148; 736 iDFS events)", "prognostic, adjuvant",
     "iDFS HR 0.87 per 10% sTILs (SD 20%)", "~HR 0.76 per SD (prognostic, needed >700 events for a precise estimate)", "10.1200/JCO.18.01010"),
    ("MGMT promoter methylation", "GBM, RT vs RT+TMZ (EORTC/NCIC, n=206)", "RCT subgroup",
     "prognostic HR 0.45 irrespective of treatment; benefit of TMZ mainly in methylated", "binary marker with a very large prognostic effect; interaction modest", "10.1056/NEJMoa043331"),
    ("T-cell-inflamed GEP + TMB", "pembrolizumab, 22 tumour types (>300 pts, KEYNOTE)", "single-arm response",
     "both independently predictive of response; reported AUROCs generally ~0.6-0.75 by tumour", "best-in-class transcriptomic response marker for IO ~ AUROC 0.65-0.75", "10.1126/science.aar3593"),
    ("ERCC1 (IHC) for adjuvant cisplatin", "NSCLC (IALT-Bio discovery; JBR.10 + CALGB 9633 validation, n=494)", "RCT interaction",
     "original interaction did not validate; antibody did not detect the functional isoform", "cautionary: biologically plausible predictive marker failed replication", "10.1056/NEJMoa1214271"),
    ("Random signatures", "breast cancer outcome", "prognostic",
     ">90% of random signatures >100 genes significantly prognostic; driven by proliferation (meta-PCNA)", "prognostic association is the default for broad expression scores", "10.1371/journal.pcbi.1002240"),
]


def retro(se, true_log_hr, alpha=0.05, n_sim=200000, seed=0):
    """Gelman & Carlin retrodesign: power, type-S, exaggeration ratio (|est|/|true| given significance)."""
    from scipy import stats
    z = stats.norm.ppf(1 - alpha / 2)
    d = abs(true_log_hr)
    p_hi = 1 - stats.norm.cdf(z - d / se)
    p_lo = stats.norm.cdf(-z - d / se)
    power = p_hi + p_lo
    typeS = p_lo / power
    rng = np.random.default_rng(seed)
    est = d + se * rng.standard_normal(n_sim)
    sig = np.abs(est) > se * z
    exag = float(np.mean(np.abs(est[sig])) / d)
    return power, typeS, exag


def main():
    L = pd.DataFrame(LIT, columns=["marker", "setting", "design", "reported_effect", "calibration_note", "doi"])
    L.to_csv(OUT / "literature_effect_sizes.tsv", sep="\t", index=False)
    W = pd.read_csv(REPO / "robust/results/survival_validation/FINAL_BIOMARKER_SURVIVAL.tsv", sep="\t")
    rows = []
    for _, r in W.iterrows():
        for ep in ("OS", "PFI"):
            p, hr, ev = r.get(f"{ep}_heldout_p"), r.get(f"{ep}_heldout_HR"), r.get(f"{ep}_heldout_events")
            if pd.notna(p) and p < 0.05 and pd.notna(ev) and ev > 0:
                se = 1 / np.sqrt(ev)  # Schoenfeld SE for a per-SD covariate
                for true_hr in (0.8, 0.7):
                    pw, ts, ex = retro(se, np.log(true_hr))
                    rows.append({"cancer": r["cancer"], "drug": r["drug"], "endpoint": ep, "events": ev, "observed_HR": hr, "p": p,
                                 "assumed_true_HR": true_hr, "power": pw, "type_S": ts, "exaggeration_ratio": ex,
                                 "expected_sig_HR_if_true": float(np.exp(-abs(np.log(true_hr)) * ex))})
    Rt = pd.DataFrame(rows)
    Rt.to_csv(OUT / "retrodesign_final_rows.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(Rt.round(3).to_string())


if __name__ == "__main__":
    main()
