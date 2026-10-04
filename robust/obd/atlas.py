"""ATLAS: pan-cancer, multi-drug pathway atlas of drug response.

One model for all drugs and cancers. Each patient x drug record (RECIST response) is
explained by pathway activity through three additive effect layers:

    logit P(response) = a_drug + covariates + x . (beta_shared + beta_class + beta_drug)

  beta_shared  pathways that govern response to systemic therapy in general
               ("overall response of the patient")
  beta_class   pathways specific to a mechanism class (platinum, taxane/microtubule, ...)
  beta_drug    drug-specific deviations

Implemented with the feature-augmentation trick for multi-task learning: each record's
pathway vector is copied into a shared block, its class block and its drug block, and a
single ridge-logistic model is fitted on the sparse augmented matrix. Block scales set
the relative penalties, so rare drugs borrow strength from their class and from the
shared layer. Pre-clinical (GDSC) and network priors enter the drug blocks as prior
means and per-pathway penalty weights, exactly as in TRIAD.

Features are Reactome rank scores centred within cancer type. Treatment setting and a
proliferation score are fitted as covariates but are not part of any pathway effect, so
pathway effects are "beyond proliferation and setting".
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import triad as T

DRUG_CLASS = {
    "CISPLATIN": "platinum", "CARBOPLATIN": "platinum", "OXALIPLATIN": "platinum",
    "FLUOROURACIL": "fluoropyrimidine", "CAPECITABINE": "fluoropyrimidine",
    "GEMCITABINE": "antimetabolite", "PEMETREXED": "antimetabolite",
    "PACLITAXEL": "microtubule", "DOCETAXEL": "microtubule", "VINORELBINE": "microtubule",
    "DOXORUBICIN": "topoisomerase_II", "EPIRUBICIN": "topoisomerase_II", "ETOPOSIDE": "topoisomerase_II",
    "IRINOTECAN": "topoisomerase_I",
    "CYCLOPHOSPHAMIDE": "alkylating", "TEMOZOLOMIDE": "alkylating", "DACARBAZINE": "alkylating",
    "IFOSFAMIDE": "alkylating",
    "BEVACIZUMAB": "antiangiogenic", "TAMOXIFEN": "hormonal", "BLEOMYCIN": "other_dna",
}


class AtlasDesign:
    """Builds the sparse augmented design for a set of records."""

    def __init__(self, pathways, drugs, scale_shared=1.0, scale_class=0.7, scale_drug=0.5):
        self.pathways = list(pathways)
        self.drugs = list(drugs)
        self.classes = sorted({DRUG_CLASS[d] for d in self.drugs})
        self.p = len(self.pathways)
        self.blocks = ["shared"] + [f"class:{c}" for c in self.classes] + [f"drug:{d}" for d in self.drugs]
        self.offset = {b: i * self.p for i, b in enumerate(self.blocks)}
        self.scale = {"shared": scale_shared, **{f"class:{c}": scale_class for c in self.classes},
                      **{f"drug:{d}": scale_drug for d in self.drugs}}
        self.n_path_cols = len(self.blocks) * self.p

    def matrix(self, X, drug, covariates):
        """X: records x pathways (array), drug: per-record drug names, covariates: records x k."""
        n = X.shape[0]
        rows, cols, vals = [], [], []
        for b_of in (lambda d: "shared", lambda d: f"class:{DRUG_CLASS[d]}", lambda d: f"drug:{d}"):
            for i in range(n):
                b = b_of(drug[i])
                rows.append(np.full(self.p, i))
                cols.append(self.offset[b] + np.arange(self.p))
                vals.append(X[i] * self.scale[b])
        A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, self.n_path_cols))
        # drug intercepts (one-hot, lightly penalised) + covariates
        D = sp.csr_matrix((np.ones(n), (np.arange(n), [self.drugs.index(d) for d in drug])), shape=(n, len(self.drugs)))
        C = sp.csr_matrix(np.asarray(covariates, float))
        return sp.hstack([A, D * 3.0, C * 3.0]).tocsr()

    def effects(self, b):
        """Split the fitted coefficient vector into per-block pathway effects (on feature scale)."""
        out = {}
        for blk in self.blocks:
            o = self.offset[blk]
            out[blk] = pd.Series(b[o:o + self.p] * self.scale[blk], index=self.pathways)
        return pd.DataFrame(out)

    def total_effect(self, eff, drug):
        return eff["shared"] + eff[f"class:{DRUG_CLASS[drug]}"] + eff[f"drug:{drug}"]


def fit(design, X, drug, covariates, y, lam=10.0, sample_weight=None, prior_m=None, prior_w=None):
    """Fit the multi-task model; prior_m / prior_w are dicts {drug: Series over pathways}."""
    M = design.matrix(X, drug, covariates)
    k = M.shape[1]
    m = np.zeros(k)
    w = np.ones(k)
    w[design.n_path_cols:] = 10.0  # intercepts and covariates nearly unpenalised
    for d in design.drugs:
        o = design.offset[f"drug:{d}"]
        sc = design.scale[f"drug:{d}"]
        if prior_m and d in prior_m and prior_m[d] is not None:
            m[o:o + design.p] = prior_m[d].reindex(design.pathways).fillna(0).values / sc
        if prior_w and d in prior_w and prior_w[d] is not None:
            w[o:o + design.p] = prior_w[d].reindex(design.pathways).fillna(0.25).values
    b0, b = T.fit_prior_logistic(M, np.asarray(y, float), lam, w, m, sample_weight)
    return b0, b


def predict(design, b, X, drug, covariates, include_covariates=False):
    M = design.matrix(X, drug, covariates)
    if not include_covariates:
        M = M[:, :design.n_path_cols]
        b = b[:design.n_path_cols]
    return M @ b
