"""Is the TRIAD external signal drug-specific or generic chemosensitivity?

For each breast neoadjuvant cohort/arm where TRIAD validated, every drug's patient model
(fitted on all TCGA patients for that drug) is applied, and pCR is modelled as
  logit(pCR) ~ score + HR status + HER2 status + proliferation (Hallmark E2F/G2M)
Drug-specific signal = the matching drug's model beats the other drugs' models and keeps
an independent effect after adjustment for proliferation and receptor status.

Usage: python robust/run_triad_specificity.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import run_triad_external as EXT  # noqa: E402
from obd import RESULTS, nextgen  # noqa: E402
from obd import network as N  # noqa: E402
from obd import trial_cohorts as TC  # noqa: E402

DRUGS = ["PACLITAXEL", "DOXORUBICIN", "DOCETAXEL", "CISPLATIN", "FLUOROURACIL", "CARBOPLATIN", "GEMCITABINE"]
ARMS = [("GSE164458_brightness", "paclitaxel arm", lambda c: c["arm"].str.startswith("paclitaxel")),
        ("GSE164458_brightness", "all arms", lambda c: pd.Series(True, index=c.index)),
        ("GSE194040_ispy2", "control arm", lambda c: c["arm"] == "Ctr"),
        ("GSE194040_ispy2", "all arms", lambda c: pd.Series(True, index=c.index)),
        ("GSE25066_hatzis", "all", lambda c: pd.Series(True, index=c.index))]


def receptor_covariates(clin):
    out = pd.DataFrame(index=clin.index)
    for col, name in (("HR", "HR_pos"), ("er", "HR_pos"), ("HER2", "HER2_pos"), ("her2", "HER2_pos")):
        if col in clin and name not in out:
            v = clin[col].astype(str).str.upper().str.strip()
            out[name] = v.isin(["1", "1.0", "POS", "P", "POSITIVE", "+"]).astype(float)
            out.loc[clin[col].isna(), name] = np.nan
    return out.loc[:, out.nunique() > 1]


def main():
    graph = {}

    def g():
        if "g" not in graph:
            graph["g"] = N.load_string()
        return graph["g"]

    models = {d: EXT.fitted_models(d, g)["patient_only"] for d in DRUGS}
    rows = []
    for cid, label, sel in ARMS:
        sc, clin = EXT.cohort_scores(cid)
        expr, _ = TC.load_trial(cid, canonical=True)
        prolif = nextgen.proliferation_score(expr)
        c = clin[sel(clin)].dropna(subset=["responder"])
        cov = receptor_covariates(c)
        for d, coef in models.items():
            s = EXT.apply(coef, sc, list(c.index))
            s = (s - s.mean()) / s.std()
            X = pd.concat([s.rename("score"), prolif.reindex(c.index).rename("proliferation"), cov], axis=1).dropna()
            yv = c.loc[X.index, "responder"].astype(float)
            fit = sm.Logit(yv, sm.add_constant(X)).fit(disp=0)
            rows.append({"cohort": cid, "arm": label, "model_drug": d, "n": len(X), "pCR": int(yv.sum()),
                         "AUC": roc_auc_score(yv, s.loc[X.index]),
                         "AUC_proliferation_alone": roc_auc_score(yv, X["proliferation"]),
                         "adj_OR_per_SD": float(np.exp(fit.params["score"])), "adj_p": float(fit.pvalues["score"]),
                         "adjusted_for": ",".join(c2 for c2 in X.columns if c2 != "score")})
            print(rows[-1], flush=True)
    out = RESULTS / "triad" / "external"
    pd.DataFrame(rows).to_csv(out / "drug_specificity.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
