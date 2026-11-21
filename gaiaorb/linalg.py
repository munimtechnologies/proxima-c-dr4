"""Small dense linear algebra in pure Python.

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.
Vectors are lists of floats; a design is a list of column vectors.
"""
import math
from operator import mul


def dot(a, b):
    return sum(map(mul, a, b))


def cholesky(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if not s > 0.0:
                    raise ValueError("matrix is not positive definite")
                L[i][i] = math.sqrt(s)
            else:
                L[i][j] = s / L[j][j]
    return L


def _forward(L, b):
    n = len(L)
    x = [0.0] * n
    for i in range(n):
        x[i] = (b[i] - sum(L[i][k] * x[k] for k in range(i))) / L[i][i]
    return x


def _backward_t(L, y):
    """Solve L^T x = y."""
    n = len(L)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(L[k][i] * x[k] for k in range(i + 1, n))) / L[i][i]
    return x


def chol_solve(L, b):
    return _backward_t(L, _forward(L, b))


def chol_inverse(L):
    n = len(L)
    cols = []
    for j in range(n):
        e = [0.0] * n
        e[j] = 1.0
        cols.append(chol_solve(L, e))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def wls(cols, y, sig):
    """Weighted least squares y ~ sum_a beta_a cols[a], weights 1/sig^2.
    Returns dict(beta, cov, chi2, resid, model)."""
    p = len(cols)
    N = len(y)
    w = [1.0 / (s * s) for s in sig]
    wc = [[c[i] * w[i] for i in range(N)] for c in cols]
    M = [[dot(wc[a], cols[b]) for b in range(p)] for a in range(p)]
    rhs = [dot(wc[a], y) for a in range(p)]
    L = cholesky(M)
    beta = chol_solve(L, rhs)
    cov = chol_inverse(L)
    model = [0.0] * N
    for a in range(p):
        ba = beta[a]
        ca = cols[a]
        for i in range(N):
            model[i] += ba * ca[i]
    resid = [y[i] - model[i] for i in range(N)]
    chi2 = sum(resid[i] * resid[i] * w[i] for i in range(N))
    return dict(beta=beta, cov=cov, chi2=chi2, resid=resid, model=model)


def extend_basis(base, cols, tol=1e-9):
    """Orthonormalise `cols` against the orthonormal vectors in `base` and
    against each other (modified Gram-Schmidt, two passes). Returns only the
    new unit vectors; columns that are (numerically) dependent are dropped."""
    new = []
    for c in cols:
        v = list(c)
        n0 = math.sqrt(dot(v, v))
        if n0 == 0.0:
            continue
        for _ in range(2):
            for b in base:
                pr = dot(b, v)
                v = [vi - pr * bi for vi, bi in zip(v, b)]
            for b in new:
                pr = dot(b, v)
                v = [vi - pr * bi for vi, bi in zip(v, b)]
        n = math.sqrt(dot(v, v))
        if n <= tol * n0:
            continue
        inv = 1.0 / n
        new.append([vi * inv for vi in v])
    return new
