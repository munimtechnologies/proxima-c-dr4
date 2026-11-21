"""Kepler's equation, Thiele-Innes constants, Campbell elements and masses.

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.

Conventions follow Gaia DR3/DR4 (Halbwachs et al. 2023; the ESA DR4
pre-release Gaia BH3 notebook):
  X = cos E - e,  Y = sqrt(1 - e^2) sin E,  E - e sin E = 2 pi (t - Tp) / P
  d(alpha*) = B X + G Y,   d(delta) = A X + F Y
  along-scan signal = (B X + G Y) sin(psi) + (A X + F Y) cos(psi)
  A =  a (cos w cos W - sin w sin W cos i)
  B =  a (cos w sin W + sin w cos W cos i)
  F = -a (sin w cos W + cos w sin W cos i)
  G = -a (sin w sin W - cos w cos W cos i)
with w = argument of periastron, W = longitude of ascending node.
"""
import math

TWO_PI = 2.0 * math.pi
DAYS_PER_YEAR = 365.25
M_SUN_IN_EARTH = 332946.0487


def ecc_anomaly(M, e):
    """Solve E - e sin E = M (radians) for 0 <= e < 1."""
    M = (M + math.pi) % TWO_PI - math.pi
    if e == 0.0:
        return M
    E = M + (0.85 * e if math.sin(M) >= 0.0 else -0.85 * e)
    for _ in range(60):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = f / fp
        E -= dE
        if abs(dE) < 1e-14:
            break
    return E


def xy(t_days, P_days, e, Tp_days):
    E = ecc_anomaly(TWO_PI * (t_days - Tp_days) / P_days, e)
    return math.cos(E) - e, math.sqrt(1.0 - e * e) * math.sin(E)


def campbell_to_ti(a, inc_deg, omega_deg, Omega_deg):
    i = math.radians(inc_deg)
    w = math.radians(omega_deg)
    W = math.radians(Omega_deg)
    ci = math.cos(i)
    cw, sw, cW, sW = math.cos(w), math.sin(w), math.cos(W), math.sin(W)
    A = a * (cw * cW - sw * sW * ci)
    B = a * (cw * sW + sw * cW * ci)
    F = -a * (sw * cW + cw * sW * ci)
    G = -a * (sw * sW - cw * cW * ci)
    return A, B, F, G


def ti_to_campbell(A, B, F, G):
    """Returns (a, inc_deg, omega_deg, Omega_deg) with Omega in [0, 180)."""
    u = 0.5 * (A * A + B * B + F * F + G * G)
    v = A * G - B * F
    a2 = u + math.sqrt(max(0.0, (u + v) * (u - v)))
    a = math.sqrt(max(0.0, a2))
    if a == 0.0:
        return 0.0, float("nan"), float("nan"), float("nan")
    ci = max(-1.0, min(1.0, v / a2))
    inc = math.degrees(math.acos(ci))
    s = math.atan2(B - F, A + G)          # w + W
    d = math.atan2(-(B + F), A - G)       # w - W
    w = 0.5 * (s + d)
    W = 0.5 * (s - d)
    W = W % TWO_PI
    if W >= math.pi:
        W -= math.pi
        w -= math.pi
    w = w % TWO_PI
    return a, inc, math.degrees(w), math.degrees(W)


def mass_function(a0_mas, plx_mas, P_days):
    """(a0/plx)^3 / P^2 in solar masses (au, yr)."""
    a_au = a0_mas / plx_mas
    P_yr = P_days / DAYS_PER_YEAR
    return a_au ** 3 / P_yr ** 2


def companion_mass(a0_mas, plx_mas, P_days, m1_msun):
    """Companion mass (Msun) for a dark companion: m^3/(m1+m)^2 = f."""
    f = mass_function(a0_mas, plx_mas, P_days)
    if f <= 0.0:
        return 0.0

    def h(m):
        return m ** 3 / (m1_msun + m) ** 2 - f

    lo, hi = 0.0, 1.0
    while h(hi) < 0.0:
        hi *= 2.0
        if hi > 1e9:
            raise ValueError("no companion mass found")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reflex_mas(m_earth, P_days, mstar_msun, plx_mas):
    """Angular semi-major axis (mas) of the star's reflex orbit."""
    m = m_earth / M_SUN_IN_EARTH
    P_yr = P_days / DAYS_PER_YEAR
    a_rel = ((mstar_msun + m) * P_yr * P_yr) ** (1.0 / 3.0)
    return a_rel * m / (mstar_msun + m) * plx_mas


def along_scan(ti, X, Y, s, c):
    A, B, F, G = ti
    return (B * X + G * Y) * s + (A * X + F * Y) * c


def major_axis_time(ti, P_days, Tp_days):
    """For a circular orbit (e = 0), a time (days) at which the projected
    position is at an end of the projected ellipse's major axis (= the line
    of nodes). Defined modulo P/2; returned in [Tp, Tp + P/2)."""
    A, B, F, G = ti
    c1 = (B, A)   # coefficient of X = cos M
    c2 = (G, F)   # coefficient of Y = sin M
    n1 = c1[0] ** 2 + c1[1] ** 2
    n2 = c2[0] ** 2 + c2[1] ** 2
    x12 = c1[0] * c2[0] + c1[1] * c2[1]
    twoM = math.atan2(2.0 * x12, n1 - n2)
    M = (0.5 * twoM) % math.pi
    return Tp_days + M / TWO_PI * P_days


def phase_offset_days(t_a, t_b, half_period):
    """Signed difference t_a - t_b wrapped to [-half_period/2, half_period/2)."""
    d = (t_a - t_b) % half_period
    if d >= 0.5 * half_period:
        d -= half_period
    return d
