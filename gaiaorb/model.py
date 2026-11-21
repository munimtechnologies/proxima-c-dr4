"""Along-scan models, jitter, grid statistic and orbit fitting.

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.

Single-star model with perspective acceleration, as in ESA's gaiasupdate
(GAIA-C3-TN-LU-LL-061, appendix E):
  w = L / D,  L = da* s + dd c + plx f + mua* t s + mud t c,
  D = 1 + (mu_r t + plx fr) k,  k = 1 mas in rad,  s = sin psi, c = cos psi.
Multiplying by D gives an exactly linear model
  w (1 + mu_r t k) = da* s + dd c + plx (f - w fr k) + mua* t s + mud t c,
with the radial proper motion mu_r = v_r plx / 4.740470446 (mas/yr) fixed
from the spectroscopic radial velocity. An orbit adds
  (B X + G Y) s + (A X + F Y) c   (relative size of the D factor < 1e-4).
"""
import math

from . import kepler
from .linalg import dot, wls, extend_basis

MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)
AU_KMS_YR = 4.740470446
PARAMS_5P = ("ra_offset_mas", "dec_offset_mas", "parallax_mas", "pmra_masyr", "pmdec_masyr")
PARAMS_TI = ("A", "B", "F", "G")


def mu_r_from_rv(rv_kms, plx_mas):
    return rv_kms * plx_mas / AU_KMS_YR


def five_cols(points, mu_r=0.0):
    k = MAS_TO_RAD
    y = [p["w"] * (1.0 + mu_r * p["t"] * k) for p in points]
    cols = [[p["s"] for p in points],
            [p["c"] for p in points],
            [p["f"] - p["w"] * p["fr"] * k for p in points],
            [p["t"] * p["s"] for p in points],
            [p["t"] * p["c"] for p in points]]
    return cols, y


def orbit_cols(points, P, e, Tp):
    """Columns for (A, B, F, G); t in days from J2017.5."""
    cA, cB, cF, cG = [], [], [], []
    for p in points:
        X, Y = kepler.xy(p["t"] * kepler.DAYS_PER_YEAR, P, e, Tp)
        cA.append(X * p["c"])
        cB.append(X * p["s"])
        cF.append(Y * p["c"])
        cG.append(Y * p["s"])
    return [cA, cB, cF, cG]


def orbit_signal(points, ti, P, e, Tp):
    out = []
    for p in points:
        X, Y = kepler.xy(p["t"] * kepler.DAYS_PER_YEAR, P, e, Tp)
        out.append(kepler.along_scan(ti, X, Y, p["s"], p["c"]))
    return out


def with_signal(points, sig_add):
    return [dict(p, w=p["w"] + s) for p, s in zip(points, sig_add)]


def _sig(base, eps):
    return [math.sqrt(s * s + eps * eps) for s in base]


def fit_jitter(cols, y, sig0, n_par, iters=60):
    """Smallest eps >= 0 with chi2 / (N - n_par) <= 1 (bisection)."""
    N = len(y)
    target = N - n_par
    if target <= 0:
        return 0.0

    def chi2(eps):
        return wls(cols, y, _sig(sig0, eps))["chi2"]

    if chi2(0.0) <= target:
        return 0.0
    srt = sorted(sig0)
    hi = max(1e-4, srt[len(srt) // 2])
    n = 0
    while chi2(hi) > target:
        hi *= 2.0
        n += 1
        if n > 80:
            raise ValueError("jitter search diverged")
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if chi2(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def clean_5p(points, mu_r=0.0, clip=5.0, jitter=True):
    """5-parameter fit with jitter and one pass of |r|/sigma > clip rejection."""
    pts = list(points)
    n_clipped = 0
    for rnd in range(2):
        cols, y = five_cols(pts, mu_r)
        base = [p["sig"] for p in pts]
        eps = fit_jitter(cols, y, base, 5) if jitter else 0.0
        sig = _sig(base, eps)
        fit = wls(cols, y, sig)
        if rnd == 1 or clip is None:
            break
        keep = [i for i in range(len(pts)) if abs(fit["resid"][i]) <= clip * sig[i]]
        if len(keep) == len(pts):
            break
        n_clipped = len(pts) - len(keep)
        pts = [pts[i] for i in keep]
    params = {name: (fit["beta"][a], math.sqrt(fit["cov"][a][a])) for a, name in enumerate(PARAMS_5P)}
    return dict(points=pts, sig=sig, eps=eps, fit=fit, params=params, n_clipped=n_clipped,
                chi2=fit["chi2"], n=len(pts), mu_r=mu_r)


class GridStat(object):
    """Delta chi^2 of (5p + orbit) vs 5p at every (P, e, Tp) grid point,
    computed as |U_g^T r|^2 with U_g an orthonormal basis of the whitened
    orbit columns after removing the whitened 5p column space."""

    def __init__(self, points, sig, grid, mu_r=0.0):
        self.points = points
        self.N = len(points)
        self.wt = [1.0 / s for s in sig]
        cols, _ = five_cols(points, mu_r)
        self.V = extend_basis([], [[c[i] * self.wt[i] for i in range(self.N)] for c in cols])
        self.grid = list(grid)
        self.U = []
        for (P, e, Tp) in self.grid:
            oc = orbit_cols(points, P, e, Tp)
            self.U.append(extend_basis(self.V, [[c[i] * self.wt[i] for i in range(self.N)]
                                                for c in oc]))

    def whiten(self, vec):
        return [v * w for v, w in zip(vec, self.wt)]

    def values(self, wvec):
        out = []
        for U in self.U:
            s = 0.0
            for u in U:
                d = dot(u, wvec)
                s += d * d
            out.append(s)
        return out

    def scan(self, wvec):
        best, bi = -1.0, -1
        for g, U in enumerate(self.U):
            s = 0.0
            for u in U:
                d = dot(u, wvec)
                s += d * d
            if s > best:
                best, bi = s, g
        return best, bi


def delta_chi2_fn(points, sig, mu_r, resid):
    """Returns f(P, e, Tp) -> delta chi^2 (orbit vs 5p) for fixed weights."""
    N = len(points)
    wt = [1.0 / s for s in sig]
    cols, _ = five_cols(points, mu_r)
    V = extend_basis([], [[c[i] * wt[i] for i in range(N)] for c in cols])
    rw = [r * w for r, w in zip(resid, wt)]

    def f(P, e, Tp):
        oc = orbit_cols(points, P, e, Tp)
        U = extend_basis(V, [[c[i] * wt[i] for i in range(N)] for c in oc])
        return sum(dot(u, rw) ** 2 for u in U)

    return f


def nelder_mead(f, x0, steps, max_iter=400, tol=1e-9):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        x = list(x0)
        x[i] += steps[i]
        pts.append(x)
    vals = [f(p) for p in pts]
    for _ in range(max_iter):
        order = sorted(range(n + 1), key=lambda k: vals[k])
        pts = [pts[k] for k in order]
        vals = [vals[k] for k in order]
        if abs(vals[-1] - vals[0]) <= tol * (abs(vals[0]) + tol):
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        worst = pts[-1]
        xr = [cen[j] + (cen[j] - worst[j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2.0 * (cen[j] - worst[j]) for j in range(n)]
            fe = f(xe)
            if fe < fr:
                pts[-1], vals[-1] = xe, fe
            else:
                pts[-1], vals[-1] = xr, fr
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            if fr < vals[-1]:
                xc = [cen[j] + 0.5 * (xr[j] - cen[j]) for j in range(n)]
            else:
                xc = [cen[j] + 0.5 * (worst[j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < min(fr, vals[-1]):
                pts[-1], vals[-1] = xc, fc
            else:
                for k in range(1, n + 1):
                    pts[k] = [pts[0][j] + 0.5 * (pts[k][j] - pts[0][j]) for j in range(n)]
                    vals[k] = f(pts[k])
    k = min(range(n + 1), key=lambda q: vals[q])
    return pts[k], vals[k]


def final_orbit_fit(points, sig, mu_r, P, e, Tp, plx_for_mass=None, m1_msun=None):
    """Full linear fit (5p + Thiele-Innes) at fixed (P, e, Tp)."""
    c5, y = five_cols(points, mu_r)
    oc = orbit_cols(points, P, e, Tp)
    fit = wls(c5 + oc, y, sig)
    b = fit["beta"]
    ti = tuple(b[5:9])
    a0, inc, om, Om = kepler.ti_to_campbell(*ti)
    # numerical propagation of the TI covariance to a0 and cos i
    covti = [row[5:9] for row in fit["cov"][5:9]]

    def g(t):
        a, i_, _, _ = kepler.ti_to_campbell(*t)
        return a, math.cos(math.radians(i_)) if a > 0 else 0.0

    base = g(ti)
    J = []
    for j in range(4):
        h = 1e-6 * max(1e-3, abs(ti[j]) + a0)
        tt = list(ti)
        tt[j] += h
        gj = g(tuple(tt))
        J.append(((gj[0] - base[0]) / h, (gj[1] - base[1]) / h))
    var_a = sum(J[p][0] * covti[p][q] * J[q][0] for p in range(4) for q in range(4))
    var_ci = sum(J[p][1] * covti[p][q] * J[q][1] for p in range(4) for q in range(4))
    out = dict(P=P, e=e, Tp=Tp, chi2=fit["chi2"], n=len(points),
               params={name: (b[a], math.sqrt(fit["cov"][a][a])) for a, name in enumerate(PARAMS_5P)},
               ti=dict(zip(PARAMS_TI, ti)),
               ti_err=dict(zip(PARAMS_TI, [math.sqrt(fit["cov"][5 + j][5 + j]) for j in range(4)])),
               a0_mas=a0, a0_err_mas=math.sqrt(max(0.0, var_a)), inc_deg=inc,
               cos_inc_err=math.sqrt(max(0.0, var_ci)), omega_deg=om, Omega_deg=Om,
               resid=fit["resid"])
    if plx_for_mass and m1_msun:
        out["mass_function_msun"] = kepler.mass_function(a0, plx_for_mass, P)
        out["companion_mass_msun"] = kepler.companion_mass(a0, plx_for_mass, P, m1_msun)
    return out


def fit_orbit(points, sig, mu_r, p_min, p_max, e_grid=(0.0, 0.25, 0.5, 0.75), n_phase=8,
              df_frac=0.2, n_starts=5, max_iter=300, e_max=0.95, P_bounds=None):
    """Global orbit search: frequency x e x Tp grid, then Nelder-Mead from the
    best distinct peaks. Returns final_orbit_fit(...) of the best solution
    plus 'dchi2' and the grid peaks."""
    c5, y = five_cols(points, mu_r)
    f5 = wls(c5, y, sig)
    dfun = delta_chi2_fn(points, sig, mu_r, f5["resid"])
    tdays = [p["t"] * kepler.DAYS_PER_YEAR for p in points]
    span = max(tdays) - min(tdays)
    df = df_frac / span
    fmin, fmax = 1.0 / p_max, 1.0 / p_min
    lo_b, hi_b = P_bounds if P_bounds else (0.8 * p_min, 1.25 * p_max)
    evals = []
    k = 0
    while True:
        fr = fmin + k * df
        if fr > fmax + 1e-15:
            break
        P = 1.0 / fr
        for e in e_grid:
            phases = [0.0] if e == 0.0 else [j / float(n_phase) for j in range(n_phase)]
            for ph in phases:
                evals.append((dfun(P, e, ph * P), P, e, ph))
        k += 1
    evals.sort(key=lambda r: -r[0])
    starts = []
    for v, P, e, ph in evals:
        if all(abs(1.0 / P - 1.0 / s[1]) > 2.0 * df for s in starts):
            starts.append((v, P, e, ph))
        if len(starts) >= n_starts:
            break

    def unpack(x):
        fr = x[0] * df
        if fr <= 0.0:
            return None
        P = 1.0 / fr
        if not (lo_b <= P <= hi_b):
            return None
        e = min(abs(x[1]), e_max)
        return P, e, (x[2] % 1.0) * P

    def obj(x):
        u = unpack(x)
        if u is None:
            return 1e300
        return -dfun(*u)

    best = None
    for v, P, e, ph in starts:
        x0 = [(1.0 / P) / df, e if e > 0 else 0.05, ph]
        x, val = nelder_mead(obj, x0, [0.3, 0.1, 0.05], max_iter=max_iter)
        if best is None or val < best[1]:
            best = (x, val)
    P, e, Tp = unpack(best[0])
    out = final_orbit_fit(points, sig, mu_r, P, e, Tp)
    out["dchi2"] = -best[1]
    out["chi2_5p"] = f5["chi2"]
    out["grid_peaks"] = [dict(dchi2=v, P=P_, e=e_, phase=ph_) for v, P_, e_, ph_ in starts]
    out["n_grid"] = len(evals)
    return out
