"""The pre-registered Proxima c test (rules: PREREGISTRATION.md).

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.

Statistic S = max over a (P, e, Tp) grid inside the radial-velocity period
box of delta chi^2 (5p + orbit vs 5p), with per-transit errors plus one
jitter term fitted to the 5p model. S is calibrated by two nulls (Gaussian
noise and a permutation of the whitened 5p residuals across transits) and
by injections of Proxima c-like orbits into the real epochs.
"""
import math
import random

from . import kepler, model

J2017_5_JD = 2457936.875

PROXIMA = dict(
    name="Proxima Centauri",
    dr3_source_id=5853498713190525696,
    plx_mas=768.0665,          # Gaia DR3
    rv_kms=-22.204,            # Kervella, Thevenin & Lovis 2017
    mass_msun=0.1221,          # Kervella, Thevenin & Lovis 2017
    P_days=1900.0, P_plus=96.0, P_minus=82.0,       # Damasso et al. 2020
    tconj_bjd=2455892.0, tconj_sigma=101.5,          # Damasso et al. 2020
    msini_earth=5.8, msini_sigma=1.9,                # Damasso et al. 2020
)

DEFAULT_CFG = dict(
    P_min=1600.0, P_max=2200.0, P_step=10.0, e_grid=(0.0, 0.1, 0.2, 0.3), n_phase=12,
    n_null=2000, fap=0.001,
    masses=(1.0, 2.5, 4.0, 5.8, 8.0, 12.0, 18.0, 25.0), n_inj=200,
    families=("iso", "edge", "k20"), k20_inc_deg=152.0,
    n_jitter_check=50, clip=5.0, min_transits=30, seed=20261202,
)


def narrow_grid(cfg):
    grid = []
    n = int(round((cfg["P_max"] - cfg["P_min"]) / cfg["P_step"]))
    for k in range(n + 1):
        P = cfg["P_min"] + k * cfg["P_step"]
        for e in cfg["e_grid"]:
            if e == 0.0:
                grid.append((P, 0.0, 0.0))
            else:
                for j in range(cfg["n_phase"]):
                    grid.append((P, e, j * P / cfg["n_phase"]))
    return grid


def quantile(vals, q):
    s = sorted(vals)
    k = max(0, min(len(s) - 1, int(math.ceil(q * len(s))) - 1))
    return s[k]


def crossing(masses, fracs, level):
    """Smallest mass at which the fraction first reaches `level` (linear
    interpolation). Returns (value, flag) with flag 'interp', 'below_grid'
    (already reached at the smallest mass) or 'above_grid' (never)."""
    for k, (m, f) in enumerate(zip(masses, fracs)):
        if f >= level:
            if k == 0:
                return masses[0], "below_grid"
            m0, f0 = masses[k - 1], fracs[k - 1]
            return m0 + (level - f0) * (m - m0) / (f - f0), "interp"
    return None, "above_grid"


def draw_period(rng, star, cfg):
    z = rng.gauss(0.0, 1.0)
    P = star["P_days"] + z * (star["P_plus"] if z > 0 else star["P_minus"])
    return min(cfg["P_max"], max(cfg["P_min"], P))


def draw_orbit(rng, family, mass_earth, star, cfg):
    """Circular orbit, period and conjunction time from the RV prior. The
    star is at argument of latitude 90 deg (projected minor axis) at the
    planet's conjunction; the overall sign of the orbit does not change the
    detection statistics."""
    P = draw_period(rng, star, cfg)
    tconj = star["tconj_bjd"] - J2017_5_JD + rng.gauss(0.0, 1.0) * star["tconj_sigma"]
    a0 = kepler.reflex_mas(mass_earth, P, star["mass_msun"], star["plx_mas"])
    if family == "iso":
        inc = math.degrees(math.acos(rng.uniform(-1.0, 1.0)))
    elif family == "edge":
        inc = 90.0
    elif family == "k20":
        inc = cfg["k20_inc_deg"]
    else:
        raise ValueError(family)
    Om = rng.uniform(0.0, 360.0)
    return dict(P=P, e=0.0, Tp=tconj, a0=a0, inc=inc, omega=90.0, Omega=Om)


def orbit_ti(orb):
    return kepler.campbell_to_ti(orb["a0"], orb["inc"], orb["omega"], orb["Omega"])


def injected_signal(points, orb):
    return model.orbit_signal(points, orbit_ti(orb), orb["P"], orb["e"], orb["Tp"])


def detection_details(pts, sig, mu_r, star, cfg, g_best, grid):
    """Refine the best grid orbit inside the period box and describe it."""
    P0, e0, Tp0 = grid[g_best]
    c5, y = model.five_cols(pts, mu_r)
    f5 = model.wls(c5, y, sig)
    dfun = model.delta_chi2_fn(pts, sig, mu_r, f5["resid"])

    def obj(x):
        P = x[0]
        if not (cfg["P_min"] <= P <= cfg["P_max"]):
            return 1e300
        e = min(abs(x[1]), 0.5)
        return -dfun(P, e, (x[2] % 1.0) * P)

    x, val = model.nelder_mead(obj, [P0, e0 if e0 > 0 else 0.05, Tp0 / P0], [20.0, 0.05, 0.05])
    P, e, Tp = x[0], min(abs(x[1]), 0.5), (x[2] % 1.0) * x[0]
    fit = model.final_orbit_fit(pts, sig, mu_r, P, e, Tp)
    fit.pop("resid", None)
    m_earth = kepler.companion_mass(fit["a0_mas"], star["plx_mas"], P, star["mass_msun"]) * kepler.M_SUN_IN_EARTH
    sin_i = math.sin(math.radians(fit["inc_deg"])) if fit["a0_mas"] > 0 else 0.0
    fit.update(dchi2=-val, mass_earth=m_earth,
               mass_err_earth=m_earth * fit["a0_err_mas"] / fit["a0_mas"] if fit["a0_mas"] > 0 else None,
               msini_earth=m_earth * sin_i)
    # circular refit at the same period: node-time consistency with the RV phase
    circ = model.final_orbit_fit(pts, sig, mu_r, P, 0.0, 0.0)
    ti = tuple(circ["ti"][k] for k in model.PARAMS_TI)
    t_major = kepler.major_axis_time(ti, P, 0.0)
    t_pred = star["tconj_bjd"] - J2017_5_JD + 0.25 * P
    d = kepler.phase_offset_days(t_major, t_pred, 0.5 * P)
    fit["rv_phase_check"] = dict(t_major_days=t_major, t_pred_days=t_pred, offset_days=d,
                                 p_random_phase=abs(d) / (0.25 * P),
                                 note="offset is mod P/2; under a random phase |offset|/(P/4) is uniform on [0,1]")
    return fit


def run_test(points, star=PROXIMA, cfg=None, log=print, use_rv=True):
    cfg = dict(DEFAULT_CFG, **(cfg or {}))
    rng = random.Random(cfg["seed"])
    res = dict(star=dict(star), cfg={k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()})
    if len(points) < cfg["min_transits"]:
        res.update(status="too_few_transits", n_transits=len(points))
        return res
    mu_r = model.mu_r_from_rv(star["rv_kms"], star["plx_mas"]) if use_rv else 0.0
    clean = model.clean_5p(points, mu_r, clip=cfg["clip"])
    pts, sig = clean["points"], clean["sig"]
    log("5p: N=%d clipped=%d jitter=%.4f mas plx=%.4f+-%.4f" % (
        clean["n"], clean["n_clipped"], clean["eps"], clean["params"]["parallax_mas"][0],
        clean["params"]["parallax_mas"][1]))
    # free perspective-acceleration check (diagnostic): fit mu_r as a 6th parameter
    k = model.MAS_TO_RAD
    c5, y = model.five_cols(pts, mu_r)
    c6 = c5 + [[-p["w"] * p["t"] * k for p in pts]]
    f6 = model.wls(c6, y, sig)
    mu_r_free = mu_r + f6["beta"][5]
    mu_r_free_err = math.sqrt(f6["cov"][5][5])

    grid = narrow_grid(cfg)
    gs = model.GridStat(pts, sig, grid, mu_r)
    rw = [r / s for r, s in zip(clean["fit"]["resid"], sig)]
    vals = gs.values(rw)
    g_obs = max(range(len(vals)), key=lambda g: vals[g])
    S_obs = vals[g_obs]
    per_P = {}
    for (P, e, Tp), v in zip(grid, vals):
        per_P[P] = max(per_P.get(P, 0.0), v)
    log("S_obs=%.3f at P=%.0f e=%.2f" % (S_obs, grid[g_obs][0], grid[g_obs][1]))

    N = len(pts)
    null_g = []
    for _ in range(cfg["n_null"]):
        null_g.append(gs.scan([rng.gauss(0.0, 1.0) for _ in range(N)])[0])
    null_p = []
    for _ in range(cfg["n_null"]):
        perm = list(rw)
        rng.shuffle(perm)
        null_p.append(gs.scan(perm)[0])
    q = 1.0 - cfg["fap"]
    thr_g, thr_p = quantile(null_g, q), quantile(null_p, q)
    thr = max(thr_g, thr_p)
    p_g = (1 + sum(1 for v in null_g if v >= S_obs)) / float(1 + len(null_g))
    p_p = (1 + sum(1 for v in null_p if v >= S_obs)) / float(1 + len(null_p))
    log("thresholds gauss=%.3f perm=%.3f  p=%.4g/%.4g" % (thr_g, thr_p, p_g, p_p))

    inj_rows = []
    table = {}
    for fam in cfg["families"]:
        for m in cfg["masses"]:
            n_det = n_exc = 0
            for _ in range(cfg["n_inj"]):
                orb = draw_orbit(rng, fam, m, star, cfg)
                s = injected_signal(pts, orb)
                S = gs.scan([r + si / sg for r, si, sg in zip(rw, s, sig)])[0]
                n_det += S > thr
                n_exc += S > S_obs
            table[(fam, m)] = (n_det / float(cfg["n_inj"]), n_exc / float(cfg["n_inj"]))
            inj_rows.append(dict(family=fam, mass_earth=m,
                                 a0_mas_at_P0=kepler.reflex_mas(m, star["P_days"], star["mass_msun"], star["plx_mas"]),
                                 n=cfg["n_inj"], frac_detected=table[(fam, m)][0],
                                 frac_above_observed=table[(fam, m)][1]))
            log("inj %s m=%.1f det=%.3f exc=%.3f" % (fam, m, table[(fam, m)][0], table[(fam, m)][1]))
    masses = list(cfg["masses"])
    limits = {}
    for fam in cfg["families"]:
        det = [table[(fam, m)][0] for m in masses]
        exc = [table[(fam, m)][1] for m in masses]
        for name, fr, lev in (("m50_detect", det, 0.5), ("m90_detect", det, 0.9), ("m95_exclusion", exc, 0.95)):
            v, flag = crossing(masses, fr, lev)
            limits["%s_%s" % (name, fam)] = dict(mass_earth=v, flag=flag)

    # jitter re-estimation check at the isotropic 95% exclusion mass
    jcheck = None
    m95 = limits.get("m95_exclusion_iso", {}).get("mass_earth")
    if m95 is not None and cfg["n_jitter_check"] > 0:
        n_exc = 0
        for _ in range(cfg["n_jitter_check"]):
            orb = draw_orbit(rng, "iso", m95, star, cfg)
            s = injected_signal(pts, orb)
            cl = model.clean_5p(model.with_signal(pts, s), mu_r, clip=cfg["clip"])
            g2 = model.GridStat(cl["points"], cl["sig"], grid, mu_r)
            r2 = [r / sg for r, sg in zip(cl["fit"]["resid"], cl["sig"])]
            n_exc += g2.scan(r2)[0] > S_obs
        jcheck = dict(mass_earth=m95, n=cfg["n_jitter_check"],
                      frac_above_observed=n_exc / float(cfg["n_jitter_check"]))

    best = detection_details(pts, sig, mu_r, star, cfg, g_obs, grid)
    res.update(
        status="ok", n_transits=N, n_clipped=clean["n_clipped"], jitter_mas=clean["eps"],
        median_sigma_mas=sorted(p["sig"] for p in pts)[N // 2],
        span_years=pts[-1]["t"] - pts[0]["t"],
        params_5p=clean["params"], chi2_5p=clean["chi2"],
        mu_r_fixed=mu_r, mu_r_free=mu_r_free, mu_r_free_err=mu_r_free_err,
        S_obs=S_obs, grid_best=dict(P=grid[g_obs][0], e=grid[g_obs][1], Tp=grid[g_obs][2]),
        threshold=thr, threshold_gauss=thr_g, threshold_perm=thr_p, p_gauss=p_g, p_perm=p_p,
        null_quantiles=dict(gauss={str(qq): quantile(null_g, qq) for qq in (0.5, 0.9, 0.99, 0.999)},
                            perm={str(qq): quantile(null_p, qq) for qq in (0.5, 0.9, 0.99, 0.999)}),
        injections=inj_rows, limits=limits, jitter_check=jcheck, best_fit=best,
        periodogram=[dict(P=P, S=per_P[P]) for P in sorted(per_P)])
    return res


def verdict(res, validation_ok, star=PROXIMA):
    """Pre-registered verdict (PREREGISTRATION.md, section 5)."""
    if not validation_ok or res.get("status") != "ok":
        return "INCONCLUSIVE"
    if res["S_obs"] > res["threshold"]:
        lo = star["msini_earth"] - 3.0 * star["msini_sigma"]
        hi = star["msini_earth"] + 3.0 * star["msini_sigma"]
        msini = res["best_fit"]["msini_earth"]
        return "DETECTED" if lo <= msini <= hi else "SIGNAL_NOT_MATCHING_RV"
    for row in res["injections"]:
        if row["family"] == "edge" and abs(row["mass_earth"] - star["msini_earth"]) < 1e-9:
            if row["frac_above_observed"] >= 0.95:
                return "EXCLUDED"
    return "UPPER_LIMIT"
