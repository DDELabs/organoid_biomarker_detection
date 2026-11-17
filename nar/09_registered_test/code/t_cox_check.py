"""Check own Cox against statsmodels PHReg (entry + strata) on simulated data."""
import numpy as np
from statsmodels.duration.hazard_regression import PHReg
from common import cox_fit
rng = np.random.default_rng(1)
n = 3000
X = rng.standard_normal((n, 3))
T = rng.exponential(np.exp(-X @ [0.3, -0.2, 0.1]))
C = rng.uniform(0, 2, n)
t = np.round(np.minimum(T, C), 2) + 0.01
e = (T <= C).astype(float)
entry = np.where(rng.random(n) < 0.3, t * rng.random(n) * 0.9, 0)
st = rng.integers(0, 4, n)
b, se, _ = cox_fit(t, e, X, entry=entry, strata=st)
r = PHReg(t, X, status=e, entry=entry, strata=st, ties="breslow").fit()
print(b, r.params); print(se, r.bse)
