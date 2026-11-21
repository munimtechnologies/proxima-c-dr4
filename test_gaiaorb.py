"""Tests for gaiaorb (proxima-c-dr4). Standard library only.

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Run: python -m unittest -v test_gaiaorb
"""
import csv
import math
import os
import random
import shutil
import tempfile
import unittest

from gaiaorb import epochs, kepler, linalg, model, prereg, validate

K = model.MAS_TO_RAD


def synth_points(n=120, seed=1, plx=50.0, pm=(100.0, -50.0), off=(1.0, -2.0), mu_r=0.0,
                 sigma=0.05, orbit=None, span=5.4, noise=True):
    """Synthetic per-transit points generated with the exact w = L / D model."""
    rng = random.Random(seed)
    pts = []
    ti = None
    if orbit:
        ti = kepler.campbell_to_ti(orbit["a0"], orbit["inc"], orbit["omega"], orbit["Omega"])
    for k in range(n):
        t = -0.5 * span + span * (k + 0.5 * rng.random()) / n
        psi = rng.uniform(-math.pi, math.pi)
        s, c = math.sin(psi), math.cos(psi)
        f, fr = rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)
        L = off[0] * s + off[1] * c + plx * f + pm[0] * t * s + pm[1] * t * c
        w = L / (1.0 + (mu_r * t + plx * fr) * K)
        if orbit:
            X, Y = kepler.xy(t * kepler.DAYS_PER_YEAR, orbit["P"], orbit["e"], orbit["Tp"])
            w += kepler.along_scan(ti, X, Y, s, c)
        if noise:
            w += rng.gauss(0.0, sigma)
        pts.append(dict(t=t, w=w, sig=sigma, s=s, c=c, f=f, fr=fr, n=9, infl=1.0))
    return pts


class TestLinalg(unittest.TestCase):
    def test_wls_exact(self):
        rng = random.Random(3)
        x = [rng.uniform(-1, 1) for _ in range(30)]
        cols = [[1.0] * 30, x, [v * v for v in x]]
        y = [2.0 - 3.0 * v + 0.5 * v * v for v in x]
        fit = linalg.wls(cols, y, [0.1] * 30)
        for got, want in zip(fit["beta"], (2.0, -3.0, 0.5)):
            self.assertAlmostEqual(got, want, places=9)
        self.assertLess(fit["chi2"], 1e-12)

    def test_extend_basis_orthonormal(self):
        rng = random.Random(4)
        cols = [[rng.gauss(0, 1) for _ in range(20)] for _ in range(4)]
        V = linalg.extend_basis([], cols[:2])
        U = linalg.extend_basis(V, cols[2:] + [cols[0]])   # cols[0] is dependent: dropped
        self.assertEqual(len(U), 2)
        for a in V + U:
            self.assertAlmostEqual(linalg.dot(a, a), 1.0, places=12)
        for i, a in enumerate(V + U):
            for b in (V + U)[i + 1:]:
                self.assertAlmostEqual(linalg.dot(a, b), 0.0, places=12)


class TestKepler(unittest.TestCase):
    def test_kepler_equation(self):
        for e in (0.0, 0.3, 0.7, 0.95):
            for k in range(-20, 21):
                M = 0.37 * k
                E = kepler.ecc_anomaly(M, e)
                self.assertLess(abs(math.remainder(E - e * math.sin(E) - M, 2 * math.pi)), 1e-11)

    def test_ti_round_trip(self):
        rng = random.Random(5)
        for _ in range(200):
            a, inc = rng.uniform(0.1, 10.0), rng.uniform(5.0, 175.0)
            om, Om = rng.uniform(0.0, 360.0), rng.uniform(0.0, 179.9)
            a2, i2, om2, Om2 = kepler.ti_to_campbell(*kepler.campbell_to_ti(a, inc, om, Om))
            self.assertAlmostEqual(a2, a, places=9)
            self.assertAlmostEqual(i2, inc, places=6)
            self.assertLess(validate._ang(om2, om), 1e-6)
            self.assertAlmostEqual(Om2, Om, places=6)

    def test_esa_bh3_thiele_innes(self):
        # ESA DR4 pre-release notebook: A, B, F, G -> a0, i, omega, Omega
        a, inc, om, Om = kepler.ti_to_campbell(2.4222, 10.560, 20.547, -16.909)
        self.assertAlmostEqual(a, 27.110, delta=0.01)
        self.assertAlmostEqual(inc, 110.5456, delta=0.05)
        self.assertAlmostEqual(om, 78.219, delta=0.1)
        self.assertAlmostEqual(Om, 136.3603, delta=0.1)

    def test_esa_bh3_masses(self):
        self.assertAlmostEqual(kepler.mass_function(27.110, 1.65980, 4183.7), 33.21, delta=0.05)
        self.assertAlmostEqual(kepler.companion_mass(27.110, 1.65980, 4183.7, 0.76), 34.68, delta=0.1)

    def test_dr3_reference_a0(self):
        self.assertAlmostEqual(kepler.ti_to_campbell(*validate.DR3_HD114762["ti"])[0], 1.79714, delta=0.002)
        self.assertAlmostEqual(kepler.ti_to_campbell(*validate.DR3_GAIA4["ti"])[0], 0.31285, delta=0.001)

    def test_proxima_reflex(self):
        a0 = kepler.reflex_mas(5.8, 1900.0, 0.1221, 768.0665)
        self.assertAlmostEqual(a0, 0.1632, delta=0.001)
        self.assertAlmostEqual(kepler.reflex_mas(12.0, 1900.0, 0.1221, 768.0665) / a0, 12.0 / 5.8, delta=1e-3)

    def test_major_axis_time_circular(self):
        P, Tp = 1900.0, -2044.875
        for inc, Om in ((60.0, 35.0), (90.0, 200.0), (130.0, 300.0)):
            ti = kepler.campbell_to_ti(1.0, inc, 90.0, Om)
            tm = kepler.major_axis_time(ti, P, Tp)
            self.assertLess(abs(kepler.phase_offset_days(tm, Tp + 0.25 * P, 0.5 * P)), 1e-6)

    def test_phase_offset(self):
        self.assertAlmostEqual(kepler.phase_offset_days(100.0, 0.0, 950.0), 100.0)
        self.assertAlmostEqual(kepler.phase_offset_days(900.0, 0.0, 950.0), -50.0)


def _transit_row(T0, bary, centroids, errs, used, psi=30.0, f=0.5, step=4.86e9):
    n = len(centroids)
    return dict(source_id=1, obs_time_tcb=[int(T0 + j * step) for j in range(n)], obs_time_bary_corr=bary,
                centroid_pos_al=centroids, centroid_pos_error_al=errs, scan_pos_angle=[psi] * n,
                parallax_factor_al=f, used_by_agis_al=used)


class TestEpochs(unittest.TestCase):
    T0 = int(7.5 * epochs.NS_PER_JYEAR)

    def rows(self):
        cen = [99.0, 10.0, 10.2, 9.8, 10.0, 10.0, 10.0, 10.0, 10.0, None]
        err = [0.1] + [0.2] * 8 + [None]
        used = [False] + [True] * 8 + [False]
        r1 = _transit_row(self.T0, 0.0, cen, err, used)
        r2 = _transit_row(self.T0 + 1e15, 0.0, [1.0, 2.0, 3.0], [0.2] * 3, [False, True, True])
        r3 = _transit_row(self.T0 + 2e15, 499.004783836156e9, [5.0, 5.0, 5.0], [0.3] * 3, [True] * 3)
        return [r1, r2, r3]

    def test_time_zero_at_j2017_5(self):
        self.assertAlmostEqual(epochs.rel_time_years(self.T0, 0.0), 0.0, places=9)
        self.assertAlmostEqual(epochs.rel_time_years(0, 0.0), -7.5, places=12)

    def test_normal_points(self):
        pts, counts = epochs.normal_points(self.rows())
        self.assertEqual(counts["transits"], 3)
        self.assertEqual(counts["kept"], 2)
        self.assertEqual(counts["too_few_ccd"], 1)
        p = pts[0]
        self.assertAlmostEqual(p["w"], 10.0, places=12)
        self.assertAlmostEqual(p["sig"], 0.2 / math.sqrt(8.0), places=12)
        self.assertEqual(p["n"], 8)
        self.assertAlmostEqual(p["t"], 4.5 * 4.86e9 / epochs.NS_PER_JYEAR, delta=1e-9)
        self.assertAlmostEqual(p["s"], 0.5, places=12)
        self.assertAlmostEqual(p["c"], math.sqrt(3.0) / 2.0, places=12)
        self.assertAlmostEqual(p["f"], 0.5, places=12)
        self.assertAlmostEqual(p["fr"], 0.0, places=12)
        self.assertAlmostEqual(pts[1]["fr"], -1.0, places=9)

    def test_inflation(self):
        r = _transit_row(self.T0, 0.0, [9.0, 11.0, 9.0, 11.0], [0.1] * 4, [True] * 4)
        p = epochs.normal_points([r])[0][0]
        # chi2_within = 4 * 1/0.01 = 400 over 3 dof
        self.assertAlmostEqual(p["sig"], math.sqrt(400.0 / 3.0) * 0.1 / 2.0, places=12)
        q = epochs.normal_points([r], inflate=False)[0][0]
        self.assertAlmostEqual(q["sig"], 0.05, places=12)

    def test_csv_round_trip(self):
        d = tempfile.mkdtemp()
        try:
            path = os.path.join(d, "ea.csv")
            keys = ["source_id", "obs_time_bary_corr", "centroid_pos_al", "centroid_pos_error_al",
                    "scan_pos_angle", "parallax_factor_al", "obs_time_tcb", "used_by_agis_al"]

            def cell(v):
                if isinstance(v, list):
                    return "(" + ", ".join("NaN" if x is None else ("true" if x is True else
                                           "false" if x is False else repr(x)) for x in v) + ")"
                return repr(v)

            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(keys)
                for r in self.rows():
                    w.writerow([cell(r[k]) for k in keys])
            rows = epochs.read_csv_rows(path)
            a = epochs.normal_points(rows)[0]
            b = epochs.normal_points(self.rows())[0]
            self.assertEqual(len(a), len(b))
            for x, y in zip(a, b):
                for k in ("t", "w", "sig", "s", "c", "f", "fr"):
                    self.assertAlmostEqual(x[k], y[k], places=12)
        finally:
            shutil.rmtree(d)

    def test_exploded_rows(self):
        rows = []
        for j, v in enumerate((1.0, 2.0, 3.0)):
            rows.append(dict(source_id=5, transit_id=77, obs_time_tcb=self.T0 + j, obs_time_bary_corr=0.0,
                             centroid_pos_al=v, centroid_pos_error_al=0.1, scan_pos_angle=0.0,
                             parallax_factor_al=0.1, used_by_agis_al=True))
        pts, counts = epochs.normal_points(rows)
        self.assertEqual(counts["kept"], 1)
        self.assertAlmostEqual(pts[0]["w"], 2.0, places=12)

    def test_parse_cell(self):
        self.assertEqual(epochs.parse_cell("(1, 2.5, NaN, true)"), [1, 2.5, None, True])
        self.assertEqual(epochs.parse_cell("0.68"), 0.68)
        self.assertIsNone(epochs.parse_cell(""))


class TestModel(unittest.TestCase):
    def test_perspective_acceleration_exact(self):
        mu_r = model.mu_r_from_rv(-22.204, 768.0665)
        self.assertAlmostEqual(mu_r, -3597.6, delta=0.5)
        pts = synth_points(plx=768.0665, pm=(-3781.7, 769.5), mu_r=mu_r, noise=False)
        cols, y = model.five_cols(pts, mu_r)
        fit = linalg.wls(cols, y, [p["sig"] for p in pts])
        for got, want in zip(fit["beta"], (1.0, -2.0, 768.0665, -3781.7, 769.5)):
            self.assertAlmostEqual(got, want, delta=1e-6)
        self.assertLess(max(abs(r) for r in fit["resid"]), 1e-6)
        # ignoring the radial term leaves a clearly visible residual
        cols0, y0 = model.five_cols(pts, 0.0)
        fit0 = linalg.wls(cols0, y0, [p["sig"] for p in pts])
        self.assertGreater(max(abs(r) for r in fit0["resid"]), 0.05)

    def test_gridstat_equals_chi2_difference(self):
        pts = synth_points(n=60, seed=9, orbit=dict(P=900.0, e=0.2, Tp=50.0, a0=0.2, inc=40.0,
                                                    omega=10.0, Omega=80.0))
        sig = [p["sig"] for p in pts]
        grid = [(900.0, 0.2, 50.0), (400.0, 0.0, 0.0), (1500.0, 0.5, 300.0)]
        gs = model.GridStat(pts, sig, grid, 0.0)
        c5, y = model.five_cols(pts, 0.0)
        f5 = linalg.wls(c5, y, sig)
        vals = gs.values([r / s for r, s in zip(f5["resid"], sig)])
        for (P, e, Tp), v in zip(grid, vals):
            f9 = linalg.wls(c5 + model.orbit_cols(pts, P, e, Tp), y, sig)
            self.assertAlmostEqual(v, f5["chi2"] - f9["chi2"], delta=1e-6 * max(1.0, v))
        self.assertEqual(max(range(3), key=lambda g: vals[g]), 0)

    def test_jitter(self):
        pts = synth_points(n=200, seed=2, sigma=0.05)
        for p in pts:
            p["sig"] = 0.025    # true noise 0.05 -> jitter ~ 0.043
        cl = model.clean_5p(pts, 0.0, clip=None)
        self.assertAlmostEqual(cl["chi2"], 200 - 5, delta=0.01)
        self.assertTrue(0.03 < cl["eps"] < 0.06)

    def test_fit_orbit_recovers(self):
        orb = dict(P=700.0, e=0.3, Tp=123.0, a0=1.0, inc=50.0, omega=40.0, Omega=110.0)
        pts = synth_points(n=120, seed=21, orbit=orb)
        fit = model.fit_orbit(pts, [p["sig"] for p in pts], 0.0, 100.0, 4000.0)
        self.assertLess(abs(fit["P"] - 700.0) / 700.0, 0.02)
        self.assertAlmostEqual(fit["a0_mas"], 1.0, delta=0.05)
        self.assertAlmostEqual(fit["e"], 0.3, delta=0.05)
        self.assertAlmostEqual(fit["inc_deg"], 50.0, delta=3.0)
        self.assertLess(validate._ang(fit["omega_deg"], 40.0), 5.0)
        self.assertLess(validate._ang(fit["Omega_deg"], 110.0), 5.0)


SMALL = dict(P_step=100.0, e_grid=(0.0, 0.2), n_phase=4, n_null=100, masses=(2.0, 60.0), n_inj=10,
             families=("iso", "edge"), n_jitter_check=3)


class TestPrereg(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(len(prereg.narrow_grid(prereg.DEFAULT_CFG)), 61 * 37)
        self.assertEqual(prereg.quantile(list(range(1, 1001)), 0.999), 999)
        self.assertAlmostEqual(prereg.crossing([1.0, 2.0, 4.0], [0.1, 0.5, 1.0], 0.9)[0], 3.6)
        self.assertEqual(prereg.crossing([1.0, 2.0], [0.96, 1.0], 0.95), (1.0, "below_grid"))
        self.assertEqual(prereg.crossing([1.0, 2.0], [0.1, 0.2], 0.95), (None, "above_grid"))

    def test_draw_orbit(self):
        rng = random.Random(1)
        for fam in ("iso", "edge", "k20"):
            o = prereg.draw_orbit(rng, fam, 5.8, prereg.PROXIMA, prereg.DEFAULT_CFG)
            self.assertTrue(1600.0 <= o["P"] <= 2200.0)
            self.assertAlmostEqual(o["a0"], kepler.reflex_mas(5.8, o["P"], 0.1221, 768.0665), places=12)
        self.assertEqual(prereg.draw_orbit(rng, "edge", 1.0, prereg.PROXIMA, prereg.DEFAULT_CFG)["inc"], 90.0)

    def _proxima_like(self, orbit=None, seed=31):
        mu_r = model.mu_r_from_rv(-22.204, 768.0665)
        return synth_points(n=100, seed=seed, plx=768.0665, pm=(-3781.7, 769.5), mu_r=mu_r,
                            sigma=0.05, orbit=orbit)

    def test_run_test_no_planet(self):
        res = prereg.run_test(self._proxima_like(), prereg.PROXIMA, SMALL, log=lambda m: None)
        self.assertEqual(res["status"], "ok")
        self.assertGreaterEqual(res["S_obs"], 0.0)
        self.assertGreater(res["threshold"], 0.0)
        big = [r for r in res["injections"] if r["mass_earth"] == 60.0]
        self.assertEqual(len(big), 2)
        self.assertTrue(all(r["frac_detected"] == 1.0 for r in big))
        self.assertLess(abs(res["mu_r_free"] / res["mu_r_fixed"] - 1.0), 0.2)
        self.assertAlmostEqual(res["params_5p"]["parallax_mas"][0], 768.0665, delta=0.05)
        self.assertIn(prereg.verdict(res, True), ("UPPER_LIMIT", "DETECTED", "SIGNAL_NOT_MATCHING_RV"))
        self.assertEqual(prereg.verdict(res, False), "INCONCLUSIVE")

    def test_run_test_detects_rv_like_planet(self):
        a0 = kepler.reflex_mas(6.0, 1900.0, 0.1221, 768.0665)
        orb = dict(P=1900.0, e=0.0, Tp=-2044.875, a0=a0, inc=90.0, omega=90.0, Omega=40.0)
        res = prereg.run_test(self._proxima_like(orb), prereg.PROXIMA, SMALL, log=lambda m: None)
        self.assertGreater(res["S_obs"], res["threshold"])
        self.assertEqual(prereg.verdict(res, True), "DETECTED")
        self.assertAlmostEqual(res["best_fit"]["mass_earth"], 6.0, delta=3.0)
        self.assertLess(res["best_fit"]["rv_phase_check"]["p_random_phase"], 0.5)

    def test_run_test_flags_planet_too_massive_for_rv(self):
        a0 = kepler.reflex_mas(30.0, 1900.0, 0.1221, 768.0665)
        orb = dict(P=1900.0, e=0.0, Tp=-2044.875, a0=a0, inc=60.0, omega=90.0, Omega=40.0)
        res = prereg.run_test(self._proxima_like(orb), prereg.PROXIMA, SMALL, log=lambda m: None)
        self.assertEqual(prereg.verdict(res, True), "SIGNAL_NOT_MATCHING_RV")

    def test_too_few_transits(self):
        res = prereg.run_test(self._proxima_like()[:10], prereg.PROXIMA, SMALL, log=lambda m: None)
        self.assertEqual(res["status"], "too_few_transits")
        self.assertEqual(prereg.verdict(res, True), "INCONCLUSIVE")


if __name__ == "__main__":
    unittest.main()
