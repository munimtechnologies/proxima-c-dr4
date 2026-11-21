# Pre-registration: does Proxima c show up in Gaia DR4 epoch astrometry?

MV2, Munim, Inc. (mv2@munimtech.com). Fixed on 2026-11-21, before Gaia DR4 (due 2026-12-02) and before
any Proxima Centauri epoch astrometry exists in public. Made with AI agents working in MV2's name; no
human expert has checked it yet. Code: `gaiaorb/` (v0.1.0) in this repository. Licence: Apache-2.0
(code) and CC BY 4.0 (text).

## 1. Question and why it is open

Proxima c is a planet candidate found in radial velocities (RV) by Damasso et al. 2020 (Sci. Adv. 6,
eaax7467). They report P = 1900 (+96/-82) d, K = 1.2 +- 0.4 m/s, m sin i = 5.8 +- 1.9 Earth masses, and a
time of inferior conjunction of BJD 2455892 (+101/-102), on a circular orbit. These numbers are quoted from
Table 3 of arXiv:2003.09305 (MNRAS 494, 1387), which lists the Damasso et al. 2020 values.

- Proxima's absolute radial velocity, -22.204 +- 0.032 km/s (corrected for gravitational redshift), and its
  mass, 0.1221 +- 0.0022 Msun, come from Kervella, Thevenin & Lovis 2017 (A&A 598, L7).

- Kervella, Arenou & Schneider 2020 (A&A 635, L14) used the Gaia DR2 proper-motion anomaly. They found
  i = 152 +- 14 deg and a mass of about 12 Earth masses.
- Benedict & McArthur 2020 used HST-FGS and found 18 +- 5 Earth masses.
- Kervella et al. 2022 (A&A 657, A7) found that the EDR3 anomaly had dropped to S/N 0.9. They wrote that
  the signal "will likely be detectable only in the epoch astrometry of Gaia".
- Later RV work has not confirmed the planet. ARES I (Libralato, Bedin & Burgasser 2025, arXiv:2512.08533)
  gives m_c = 3.4 (+5.2/-3.4) Earth masses from HST + Gaia.

Gaia DR4 (66 months, about one full orbit of c) is the first release with epoch astrometry for every
source. This test asks one question: **is the RV planet's reflex orbit in Proxima's DR4 along-scan
measurements, and if not, which masses are excluded?**

Expected size of the star's reflex orbit at P = 1900 d (M* = 0.1221 Msun, parallax 768.0665 mas):

| mass (Earth masses) | 2.5 | 5.8 | 12 | 18 |
|---|---|---|---|---|
| reflex semi-major axis (mas) | 0.070 | 0.163 | 0.338 | 0.507 |

## 2. Data

- **Target.** Gaia DR4 epoch astrometry of Proxima Centauri (DR3 source_id 5853498713190525696). The DR4
  source is found by a cone search with parallax > 700 mas.
- **Measurements.** Only CCD measurements with `used_by_agis_al == True` and finite values are used. Each
  transit's AF CCDs are collapsed to one weighted-mean along-scan point. Transits with fewer than 3 such
  CCDs are dropped.
- **Point error.** The error of the weighted mean, multiplied by sqrt(chi2_within / (n - 1)) when that
  factor is above 1.
- **Time.** As in ESA's gaiasupdate package: t = 2010.0 + (obs_time_tcb + obs_time_bary_corr) / (Julian
  year in ns) - 2017.5.
- **Allowed changes after release.** Only the I/O may be adapted to the format actually released: the
  source lookup, the download URLs and the column parsing.

## 3. Validation gates (on the June 2026 DR4 pre-release, run before 2026-12-02)

The fitter must pass all of these on ESA's public pre-release epoch astrometry. If any gate fails, the
Proxima verdict is INCONCLUSIVE until the cause is fixed and documented.

- **Gaia BH3, 5-parameter fit.** Compared with ESA's own fit of the same epochs (esa/gaia-bhthree,
  branch gaia-dr4-prerelease). All five parameters must be within 0.05 mas (or mas/yr).
- **Gaia BH3 orbit,** against ESA's notebook values:
  - P: within 3 sigma (216 d) of 4183.7 d
  - e: within 0.02 of 0.7278
  - a0: within 0.8 mas of 27.110 mas
  - i: within 1 deg of 110.55 deg
  - omega and Omega: within 2 deg of 78.22 and 136.36 deg
  - parallax: within 0.03 mas of 1.6598 mas
  - M2: within 2 Msun of 34.68 Msun (for M1 = 0.76 Msun)
- **HD 114762 orbit,** against Gaia DR3: P within 1 d of 83.74 d, e within 0.13 of 0.32, a0 in
  [1.35, 2.25] mas (DR3: 1.80).
- **Gaia-4 orbit,** against Gaia DR3: P within 36 d of 564 d, a0 in [0.16, 0.47] mas (DR3: 0.31).
- **Three CRF3 quasars.** For each:
  - |parallax / error| < 4
  - a circular orbit search over 100-4000 d is not significant (Gaussian-null p > 0.01, 300 simulations)
- **Injection recovery into HD 183633 (G ~ 9).** Inject 10 orbits with a0 = 1 mas, P = 700 d, e = 0.3 and
  random orientations. At least 8 must be recovered with P within 3%, a0 within 15% and e within 0.1.

The orbit gates use formal errors with no jitter and no clipping, as ESA's notebook does. The global orbit
search is a frequency x e x Tp grid followed by Nelder-Mead from the 5 best distinct peaks.

## 4. The test

- **Model.** Five astrometric parameters plus perspective acceleration in ESA's form (w = L / D). The
  radial proper motion is fixed from the RV: -22.204 km/s (Kervella, Thevenin & Lovis 2017) times parallax
  768.0665 mas, which gives -3597.6 mas/yr. The orbit adds Thiele-Innes terms (B X + G Y) sin psi +
  (A X + F Y) cos psi.
- **Noise.** One jitter term is added in quadrature. It is set so that the 5-parameter fit has reduced
  chi2 = 1. Transits with |residual| > 5 sigma are dropped once, and the jitter is refit.
- **Statistic.** S is the maximum, over a grid, of delta chi2 (orbit model minus 5-parameter model). The
  grid covers:
  - P = 1600-2200 d in 10 d steps (Damasso et al.'s narrow prior box)
  - e = 0, 0.1, 0.2, 0.3
  - 12 periastron phases for each e > 0
  - 2,257 grid points in all
- **Nulls.** There are 2,000 draws of each:
  - Gaussian noise with the adopted errors
  - random permutations of the whitened 5-parameter residuals across transits
- **Detection threshold.** The larger of the two nulls' 99.9% quantiles (false-alarm probability 0.001).
  p-values are reported for both.
- **Injections.** Circular orbits are injected into the real epochs:
  - P drawn from 1900 (+96/-82) d, clipped to the box
  - conjunction time drawn from BJD 2455892 +- 101.5
  - reflex size from M* and the parallax
  - masses 1, 2.5, 4, 5.8, 8, 12, 18 and 25 Earth masses, 200 draws each
- **Orientation families.** "iso" (isotropic inclination, random node), "edge" (i = 90 deg, so true mass =
  m sin i) and "k20" (i = 152 deg, Kervella et al. 2020).
- **Reported per mass:** the fraction above the detection threshold, and the fraction above the observed S
  (the exclusion fraction).
- **Jitter check.** At the isotropic 95% exclusion mass, 50 more injections are run with the jitter and
  clipping redone from scratch. They are reported, not used for the verdict.
- **Diagnostics,** reported but not part of the verdict:
  - a fit with the radial proper motion left free (it should be near -3598 mas/yr)
  - the best-fit orbit (P, e, a0, i, Omega, mass, m sin i)
  - for a circular refit, the offset of the projected major-axis time from the RV-predicted node time
    (conjunction + P/4, modulo P/2)

## 5. Verdicts (fixed now)

- **INCONCLUSIVE:** a validation gate failed, or fewer than 30 transits survive.
- **DETECTED:** S > threshold, and the fitted m sin i lies in [0.1, 11.5] Earth masses (the RV value +- 3
  sigma). The report gives the mass, the inclination and the RV-phase offset.
- **SIGNAL_NOT_MATCHING_RV:** S > threshold, but m sin i is outside that range.
- **EXCLUDED:** no detection, and at least 95% of the edge-on 5.8 Earth-mass injections give S above the
  observed S. The published planet would then be ruled out at 95% for every inclination, given a circular
  orbit with P inside the box. Edge-on is the hardest case at fixed m sin i.
- **UPPER_LIMIT:** otherwise. The report gives the 95% exclusion masses for the iso, edge and k20 families.

## 6. Rehearsal (run now)

Before the release, the whole test in section 4 runs on the pre-release star HD 183633 (G ~ 9, like
Proxima). Its own epochs serve as the data, with the radial term switched off (its parallax is about
1 mas), and Proxima c-like orbits are injected using Proxima's parallax and mass. This forecasts how well
DR4 can detect Proxima c at G ~ 9. Proxima's real transit count and noise will differ, so the forecast is
indicative only.

## 7. What this does not claim

- A non-detection rules out only the circular, P = 1600-2200 d planet as modelled. It does not rule out
  RV signals at other periods, or strongly eccentric orbits.
- Gaia DR4 systematics at G ~ 9 (calibration, attitude, the gating of bright stars) are treated as white
  noise per transit plus one jitter term.
- DPAC's own DR4 exoplanet pipeline may report on Proxima. This test is independent of it and fixed in
  advance.
