"""Validation gates on the Gaia DR4 pre-release epoch astrometry
(rules: PREREGISTRATION.md, section 3).

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.

References:
- Gaia BH3: ESA's own fit of the same pre-release epochs (esa/gaia-bhthree,
  branch gaia-dr4-prerelease, notebook Gaia-DR4-prerelease_BH3_fit_astrometric_orbit):
  5p fit ra 0.71001, dec 0.26549, plx 0.98329, pmra -30.94602, pmdec -148.54974;
  orbit P 4183.7 +- 72.1 d, e 0.72782 +- 0.00353, a0 27.110 +- 0.367 mas,
  i 110.5456 +- 0.0715, omega 78.219 +- 0.110, Omega 136.3603 +- 0.0925,
  plx 1.65980 +- 0.00652, M2 = 34.68 Msun for M1 = 0.76.
- HD 114762 and Gaia-4: Gaia DR3 nss_two_body_orbit (34-month data).
"""
import math
import random

from . import kepler, model

BH3 = 4318465066420528000
HD114762 = 3937211745905473024
GAIA4 = 1457486023639239296
HD183633 = 4181040337841125632
QSOS = (2237987199365376, 60730287810150016, 10973744521070720)

ESA_BH3_5P = dict(ra_offset_mas=0.71001, dec_offset_mas=0.26549, parallax_mas=0.98329,
                  pmra_masyr=-30.94602, pmdec_masyr=-148.54974)
ESA_BH3_ORBIT = dict(P=4183.7, P_err=72.1, e=0.72782, a0=27.110, inc=110.5456, omega=78.219,
                     Omega=136.3603, plx=1.65980, m2=34.68, m1=0.76)
DR3_HD114762 = dict(P=83.73818467439406, e=0.323893426184342,
                    ti=(1.7332884654363947, 0.13130528033229202, -0.03868510638699714, 1.7603556459114833))
DR3_GAIA4 = dict(P=564.0242807973492, P_err=11.195968,
                 ti=(0.2863899608855841, -0.12249067844532667, -0.07725115765112353, -0.11373315816651848))


def _ang(a, b):
    d = (a - b) % 360.0
    return min(d, 360.0 - d)


def gate(name, value, ref, tol, ok):
    return dict(gate=name, value=value, reference=ref, tolerance=tol, passed=bool(ok))


def esa_mode(points):
    """Formal errors, no jitter, no clipping (as in ESA's notebook)."""
    return [p["sig"] for p in points]


def check_bh3_5p(points):
    sig = esa_mode(points)
    c5, y = model.five_cols(points, 0.0)
    fit = model.wls(c5, y, sig)
    out = []
    for a, name in enumerate(model.PARAMS_5P):
        v = fit["beta"][a]
        out.append(gate("bh3_5p_" + name, v, ESA_BH3_5P[name], 0.05, abs(v - ESA_BH3_5P[name]) <= 0.05))
    return out


def check_bh3_orbit(points):
    sig = esa_mode(points)
    fit = model.fit_orbit(points, sig, 0.0, 200.0, 10000.0)
    r = ESA_BH3_ORBIT
    m2 = kepler.companion_mass(fit["a0_mas"], fit["params"]["parallax_mas"][0], fit["P"], r["m1"])
    g = [gate("bh3_P", fit["P"], r["P"], 3 * r["P_err"], abs(fit["P"] - r["P"]) <= 3 * r["P_err"]),
         gate("bh3_e", fit["e"], r["e"], 0.02, abs(fit["e"] - r["e"]) <= 0.02),
         gate("bh3_a0", fit["a0_mas"], r["a0"], 0.8, abs(fit["a0_mas"] - r["a0"]) <= 0.8),
         gate("bh3_inc", fit["inc_deg"], r["inc"], 1.0, abs(fit["inc_deg"] - r["inc"]) <= 1.0),
         gate("bh3_omega", fit["omega_deg"], r["omega"], 2.0, _ang(fit["omega_deg"], r["omega"]) <= 2.0),
         gate("bh3_Omega", fit["Omega_deg"], r["Omega"], 2.0, _ang(fit["Omega_deg"], r["Omega"]) <= 2.0),
         gate("bh3_plx", fit["params"]["parallax_mas"][0], r["plx"], 0.03,
              abs(fit["params"]["parallax_mas"][0] - r["plx"]) <= 0.03),
         gate("bh3_m2", m2, r["m2"], 2.0, abs(m2 - r["m2"]) <= 2.0)]
    fit.pop("resid", None)
    return g, fit


def check_dr3_orbit(points, ref, label, p_tol, a0_lo, a0_hi, e_tol=None, p_min=20.0, p_max=5000.0):
    sig = esa_mode(points)
    fit = model.fit_orbit(points, sig, 0.0, p_min, p_max)
    a0_ref = kepler.ti_to_campbell(*ref["ti"])[0]
    g = [gate(label + "_P", fit["P"], ref["P"], p_tol, abs(fit["P"] - ref["P"]) <= p_tol),
         gate(label + "_a0", fit["a0_mas"], a0_ref, [a0_lo, a0_hi], a0_lo <= fit["a0_mas"] <= a0_hi)]
    if e_tol is not None:
        g.append(gate(label + "_e", fit["e"], ref["e"], e_tol, abs(fit["e"] - ref["e"]) <= e_tol))
    fit.pop("resid", None)
    return g, fit


def circular_grid(p_min, p_max, span_days, df_frac=0.2):
    df = df_frac / span_days
    grid = []
    f = 1.0 / p_max
    while f <= 1.0 / p_min + 1e-15:
        grid.append((1.0 / f, 0.0, 0.0))
        f += df
    return grid


def check_quasar(points, label, n_sim=300, seed=7):
    cl = model.clean_5p(points, 0.0, clip=5.0)
    plx, eplx = cl["params"]["parallax_mas"]
    pts, sig = cl["points"], cl["sig"]
    span = (pts[-1]["t"] - pts[0]["t"]) * kepler.DAYS_PER_YEAR
    gs = model.GridStat(pts, sig, circular_grid(100.0, 4000.0, span), 0.0)
    rw = [r / s for r, s in zip(cl["fit"]["resid"], sig)]
    S = gs.scan(rw)[0]
    rng = random.Random(seed)
    n = len(pts)
    sims = [gs.scan([rng.gauss(0.0, 1.0) for _ in range(n)])[0] for _ in range(n_sim)]
    p = (1 + sum(1 for v in sims if v >= S)) / float(1 + n_sim)
    g = [gate(label + "_plx_over_err", plx / eplx, 0.0, 4.0, abs(plx / eplx) < 4.0),
         gate(label + "_orbit_search_p", p, None, 0.01, p > 0.01)]
    return g, dict(plx=plx, plx_err=eplx, S=S, p=p, n=n, jitter=cl["eps"])


def check_injection_recovery(points, n=10, seed=11, a0=1.0, P=700.0, e=0.3):
    rng = random.Random(seed)
    cl = model.clean_5p(points, 0.0, clip=5.0)
    base = cl["points"]
    ok = 0
    rows = []
    for _ in range(n):
        inc = math.degrees(math.acos(rng.uniform(-1.0, 1.0)))
        om, Om = rng.uniform(0.0, 360.0), rng.uniform(0.0, 180.0)
        Tp = rng.uniform(0.0, P)
        ti = kepler.campbell_to_ti(a0, inc, om, Om)
        pts = model.with_signal(base, model.orbit_signal(base, ti, P, e, Tp))
        # noise model of the un-injected data (formal errors + its 5p jitter)
        fit = model.fit_orbit(pts, cl["sig"], 0.0, 50.0, 5000.0)
        good = (abs(fit["P"] - P) / P < 0.03 and abs(fit["a0_mas"] - a0) / a0 < 0.15
                and abs(fit["e"] - e) < 0.1)
        ok += good
        rows.append(dict(inc=inc, omega=om, Omega=Om, Tp=Tp, P_fit=fit["P"], e_fit=fit["e"],
                         a0_fit=fit["a0_mas"], recovered=bool(good)))
    return [gate("hd183633_injection_recovery", ok, n, ">= 0.8 n", ok >= 0.8 * n)], rows
