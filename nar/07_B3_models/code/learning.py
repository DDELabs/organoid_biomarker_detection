"""Learning curves vs number of training lab models (organoids / PDX / cell lines standing in for organoids).

For the larger sets: 5 outer folds (20 % of the lab models held out); inside the other 80 %, train on 25/50/75/100 %
(1 random draw per fold). Metrics: Spearman on the held-out lab models, and C-index in the treated TCGA patients
of the model trained on that subsample. Same subsamples for every model.
"""
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

import b3common as B
from bench import Ctx, fit_predict

PAIRS = [("COAD_isella", "CETUXIMAB"), ("LIHC_ji", "SORAFENIB"), ("PAAD_shi", "GEMCITABINE"),
         ("LUAD_cells", "PEMETREXED"), ("PDX_CRC", "FLUOROURACIL"), ("COAD_paper", "FLUOROURACIL")]
LC_MODELS = ["Paper: single network pathway", "Ridge: all network pathways", "All pathways (ridge)", "Elastic net (all pathways)",
             "Random forest (all pathways)", "XGBoost (all pathways)", "Fingerprint MLP: fine-tuned on lab models",
             "GCN on STRING (genes)"]


def main():
    rows = []
    rng = np.random.default_rng(B.SEED)
    for b, d in PAIRS:
        C = Ctx(b, d)
        n = len(C.S)
        for fo, (tr0, te) in enumerate(KFold(5, shuffle=True, random_state=0).split(np.arange(n))):
            for frac in (0.25, 0.5, 0.75, 1.0):
                for rep in range(1):
                    m = max(6, int(round(frac * len(tr0))))
                    tr = np.sort(rng.choice(tr0, m, replace=False)) if frac < 1 else tr0
                    for name in LC_MODELS:
                        try:
                            po, pp, _ = fit_predict(name, C, tr, te, full=True)
                        except Exception as e:  # noqa: BLE001
                            print("ERR", b, d, name, repr(e), flush=True)
                            continue
                        r = B.spearman_ci(po, C.y[te])["rho"]
                        c = B.cindex(C.T[C.ti], C.E[C.ti], pp[C.ti]) if np.std(pp[C.ti]) > 0 else np.nan
                        rows.append(dict(bundle=b, drug=d, fold=fo, frac=frac, rep=rep, n_train=len(tr), model=name,
                                         rho_heldout=r, C_index=c))
            print(b, d, "fold", fo, flush=True)
        pd.DataFrame(rows).to_csv(B.TABLES / "learning_curve_lab_models.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
