"""NIT: Network-Informed, deconfounded Transfer model (next-generation framework).

Keeps the foundations of the Kong et al. framework (drug-target network proximity,
pathway-level features, organoid-trained linear model, patient survival readout) and
fixes the failure modes found in the robust analysis:

1. Response deconfounding. In viability assays, fast-growing organoids look sensitive
   to almost every drug, and general sensitivity dominates single-drug responses
   (LICOB sorafenib AUC: rho = -0.35 with proliferation, 0.40 with PC1 of all drugs).
   In patients, proliferation is prognostic, so a model of raw AUC turns into
   "slow tumours do better" and the hazard ratio flips. We learn from the residual
   of the drug response after regressing out proliferation and leave-drug-out
   general sensitivity: the drug-specific part of the response.
2. Soft network prior. The hard proximity cut-off (z <= -1.28) is replaced by a
   smooth weight per pathway (logistic in z); distant pathways are shrunk, not dropped.
3. Transfer prior. When cell-line data for the drug exist, a pan-cancer cell-line
   model provides prior coefficients; organoid fits shrink toward it, not toward zero.
4. Patient side. Survival is re-based at treatment start (landmark) when start days
   are known, removing immortal-time bias; patient proliferation enters the Cox model,
   so the signature has to predict beyond general prognosis.
"""
import numpy as np
import pandas as pd

from . import EXTERNAL, models as M
from .cohorts import fetch
from .reference import gene_sets
from .scoring import rank_score

HALLMARK_URL = "https://raw.githubusercontent.com/saezlab/VisiumMS/main/config/h.all.v2023.1.Hs.symbols.gmt"
PROLIFERATION_SETS = ("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")


def proliferation_score(expr):
    """Per-sample proliferation (mean rank score of Hallmark E2F + G2M), z-scored."""
    fetch(HALLMARK_URL, EXTERNAL / "h.all.symbols.gmt")
    gene_sets.cache_clear()
    gs = gene_sets(("HALLMARK",))
    sc = rank_score(expr, {k: gs[k] for k in PROLIFERATION_SETS}, min_size=10).mean()
    return (sc - sc.mean()) / sc.std()


def general_sensitivity(response, drug, min_coverage=0.8):
    """First principal component of all OTHER drugs' responses (leave-drug-out)."""
    other = response.drop(columns=[drug])
    other = other.loc[:, other.notna().mean() >= min_coverage]
    if other.shape[1] < 5:
        return None
    z = ((other - other.mean()) / other.std()).fillna(0.0)
    u, s, _ = np.linalg.svd(z.values, full_matrices=False)
    pc1 = pd.Series(u[:, 0] * s[0], index=z.index)
    if np.corrcoef(pc1, z.mean(axis=1))[0, 1] < 0:
        pc1 = -pc1  # orient: high = generally resistant
    return (pc1 - pc1.mean()) / pc1.std()


def deconfound(response, drug, prolif=None, general=None):
    """Drug-specific response: residual of the response on proliferation and general sensitivity."""
    y = response[drug].dropna()
    cols = {}
    if prolif is not None:
        cols["proliferation"] = prolif
    if general is not None:
        cols["general"] = general
    if not cols:
        return y, {}
    Z = pd.DataFrame(cols).reindex(y.index).dropna()
    y = y.loc[Z.index]
    A = np.column_stack([np.ones(len(Z)), Z.values])
    beta, *_ = np.linalg.lstsq(A, y.values, rcond=None)
    resid = pd.Series(y.values - A @ beta, index=y.index)
    r2 = 1 - resid.var() / y.var()
    return resid, {"confounder_R2": float(r2), **{f"beta_{k}": float(b) for k, b in zip(Z.columns, beta[1:])}}


def network_weights(z, center=-1.2816, scale=0.5, floor=0.05):
    """Soft network prior: logistic weight in the proximity z-score (closer = higher)."""
    w = 1.0 / (1.0 + np.exp((z - center) / scale))
    return w[w >= floor]


def transfer_prior(cell_scores, cell_response, features, alpha=10.0):
    """Ridge model on cell lines (pathway scores -> response), standardised coefficients."""
    from sklearn.linear_model import Ridge

    y = cell_response.dropna()
    lines = [l for l in y.index if l in cell_scores.columns]
    feats = [f for f in features if f in cell_scores.index]
    if len(lines) < 30 or not feats:
        return None
    X = M.zscore(cell_scores.loc[feats, lines].T.values)
    yz = (y[lines].values - y[lines].mean()) / y[lines].std()
    coef = Ridge(alpha=alpha).fit(X, yz).coef_
    return pd.Series(coef, index=feats)
