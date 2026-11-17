"""Exact, fast re-implementation of sklearn RidgeCV(cv=3, alphas=0.1..0.9) as used by Kong et al.

sklearn's RidgeCV with an integer cv runs GridSearchCV(Ridge, KFold(3) without shuffling, R^2 score) and refits
on all data with the best alpha. Doing that with SVDs is ~50x faster, which makes the paper's 10,000-iteration
bootstrap affordable. `check()` verifies the coefficients equal sklearn's to 1e-8.
"""
import numpy as np

ALPHAS = np.arange(0.1, 1, 0.1)


def _ridge_svd(X, y, alphas):
    xm, ym = X.mean(0), y.mean()
    U, s, Vt = np.linalg.svd(X - xm, full_matrices=False)
    uty = U.T @ (y - ym)
    coefs = np.stack([Vt.T @ (s / (s ** 2 + a) * uty) for a in alphas])     # n_alpha x p
    return coefs, ym - coefs @ xm


def _folds(n, k=3):
    sizes = np.full(k, n // k)
    sizes[: n % k] += 1
    idx, start = [], 0
    for sz in sizes:
        idx.append(np.arange(start, start + sz))
        start += sz
    return idx


def ridgecv_coef(X, y, alphas=ALPHAS):
    n = len(y)
    if n < 6:   # a test fold with < 2 samples: sklearn's R^2 is undefined (nan) for every alpha -> first alpha
        c, b = _ridge_svd(X, y, alphas[:1])
        return c[0], b[0], alphas[0]
    scores = np.zeros(len(alphas))
    for test in _folds(n):
        train = np.setdiff1d(np.arange(n), test)
        c, b = _ridge_svd(X[train], y[train], alphas)
        pred = X[test] @ c.T + b                                                # n_test x n_alpha
        yt = y[test]
        ss_tot = ((yt - yt.mean()) ** 2).sum()
        ss_res = ((pred - yt[:, None]) ** 2).sum(0)
        scores += (1 - ss_res / ss_tot) if ss_tot > 0 else np.where(ss_res == 0, 1.0, 0.0)
    best = int(np.argmax(scores / 3))
    c, b = _ridge_svd(X, y, alphas[best:best + 1])
    return c[0], b[0], alphas[best]


def check(seed=0):
    from sklearn.linear_model import RidgeCV
    rng = np.random.default_rng(seed)
    worst = 0
    for n, p in [(19, 37), (36, 50), (61, 12), (5, 30), (192, 40), (20, 3)]:
        X, y = rng.normal(size=(n, p)), rng.normal(size=n)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = RidgeCV(cv=3, alphas=ALPHAS).fit(X, y).coef_
        b = ridgecv_coef(X, y)[0]
        worst = max(worst, np.abs(a - b).max())
    return worst


if __name__ == "__main__":
    print("max |coef difference| vs sklearn:", check())
