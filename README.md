# proxima-c-dr4

**Is Proxima c in Gaia DR4's epoch astrometry?** This is a pre-registered test, fixed before the data exist.

Proxima c is a planet candidate. It was found in radial velocities in 2020 (P of about 5.2 yr, m sin i of
5.8 +- 1.9 Earth masses) but has never been confirmed. Kervella et al. (2022) concluded that its signal "will
likely be detectable only in the epoch astrometry of Gaia". Gaia DR4, due on 2 December 2026, is the first
release of that epoch astrometry: about 5.5 years of individual along-scan measurements, roughly one full
orbit of c. Depending on its mass, the planet should move Proxima by 0.07-0.5 mas.

This repository contains:

- `gaiaorb/`, a pure-Python (standard library) fitter for Gaia DR4 along-scan epoch astrometry. It covers:
  - the five astrometric parameters plus perspective acceleration in ESA's form (this term is about 0.5 mas
    for Proxima, larger than the planet's signal)
  - Thiele-Innes orbits, with a global orbit search
  - a grid statistic calibrated by Gaussian and permutation nulls
  - injection-recovery of RV-prior orbits
  - the fixed verdict rules
- `PREREGISTRATION.md`, the rules, fixed on 2026-11-21: the validation gates, the statistic, the
  calibration and the verdicts (DETECTED / SIGNAL_NOT_MATCHING_RV / EXCLUDED / UPPER_LIMIT / INCONCLUSIVE).
- `validate_job.py`, which runs on ESA's public DR4 pre-release (12 sources, June 2026). It does two things:
  - **Validation.** The fitter must reproduce ESA's own Gaia BH3 fit of the same epochs, the Gaia DR3 orbits
    of HD 114762 and Gaia-4, null results on three quasars, and injected orbits in a G ~ 9 star.
  - **Rehearsal.** The full Proxima c test is run on that G ~ 9 star, as a sensitivity forecast.
- `proxima_job.py`, the test itself, to run on Proxima's DR4 epochs once they are public.
- `test_gaiaorb.py`, the tests (`python -m unittest -v test_gaiaorb`, standard library only).

Run:
- validation: `pip install -r requirements.txt && python validate_job.py --data-dir data --out-dir out`
- the test, after 2026-12-02: `python proxima_job.py --data-dir data --out-dir out --validated`

Credit: the Proxima c work this builds on is described in `PREREGISTRATION.md` and `notes/`:
- Damasso et al. 2020 (the RV detection)
- Kervella, Arenou & Schneider 2020 and Kervella et al. 2022 (the Gaia proper-motion anomaly)
- Benedict & McArthur 2020 (HST-FGS)
- Libralato, Bedin & Burgasser 2025 (ARES I)
- Kervella, Thevenin & Lovis 2017 (Proxima's radial velocity and mass)

The data model, time conventions and the BH3 reference fit come from ESA's Gaia DR4 pre-release and ESA's
open `gaiasupdate` and `gaia-bhthree` code. This work uses data from the European Space Agency mission Gaia,
processed by the Gaia DPAC.

By MV2, Munim, Inc. (mv2@munimtech.com). The work was made with AI agents working in MV2's name, and no
human expert has checked it yet.
Code: Apache-2.0. Text and result files: CC BY 4.0.
