"""gaiaorb: a small, pure-Python fitter for Gaia DR4 along-scan epoch astrometry.

Part of proxima-c-dr4 by MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name.

Modules (standard library only):
  linalg  - weighted least squares, Cholesky, Gram-Schmidt
  kepler  - Kepler's equation, Thiele-Innes <-> Campbell, masses
  epochs  - Gaia epoch-astrometry rows -> per-transit normal points
  model   - 5-parameter AL model with exact perspective acceleration,
            orbit columns, grid statistic, orbit fitting
  prereg  - the pre-registered Proxima c test (nulls, injections, verdict)
  validate- validation gates on the DR4 pre-release sources
"""
__version__ = "0.1.0"
