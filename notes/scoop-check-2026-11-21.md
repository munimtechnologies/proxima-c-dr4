# Scoop check, 2026-11-21 (before Gaia DR4)

Searched:
- the arXiv API (all: Proxima AND Gaia AND astrometry, newest first)
- web search for "Proxima c" "Gaia DR4" epoch astrometry 2026

Nothing analyses Proxima c with Gaia DR4 epoch astrometry. That is expected, because the data are not
public until 2026-12-02. The newest related work is:

- **ARES I** (Libralato, Bedin & Burgasser, arXiv:2512.08533, A&A). Its HST + Gaia DR3 proper-motion
  anomaly gives m_c = 3.4 (+5.2/-3.4) Earth masses, and it plans HST spatial scans.
- **Kervella et al. 2022** (A&A 657, A7). The EDR3 anomaly has S/N 0.9, and the signal is "likely
  detectable only in the epoch astrometry of Gaia".
- **Kervella, Arenou & Schneider 2020** (A&A 635, L14). From DR2: i = 152 +- 14 deg, about 12 Earth masses.
- **Benedict & McArthur 2020.** From HST-FGS: 18 +- 5 Earth masses.

The race: on release day, DPAC's own DR4 exoplanet products, the Kervella group, the Damasso/Faria/Suarez
Mascareno RV groups and the ARES team can all fit Proxima's epochs. This repository fixes the RV-prior test,
its calibration and its verdicts before any of those data exist.

The biggest scoop risk is DPAC itself. The DR4 product list (cosmos.esa.int/web/gaia/dr4, updated
2026-09-14) includes:
- `gaia_exoplanet_candidates`
- `rv_combined_orbital_solutions`
- `rv_constrained_single_systems` and `rv_constrained_multiple_systems`, from the paper "Gaia Data
  Release 4: The Exoplanet Candidate Catalogue"

DPAC may therefore publish its own RV-constrained solution for Proxima on release day. If it does, this test
becomes an independent check, fixed in advance, with a calibrated upper limit or detection, and it should
be described that way. Check these tables for Proxima first on 2026-12-02.

Epoch astrometry for most stars is in the DataLink-only `epoch_astrometry` table. The 49 of the 51
brightest stars are in `bright_source_epoch_astrometry`. Proxima (G ~ 9) should be in the main table.

Re-check right before publishing a result:
- arXiv listings astro-ph.EP and astro-ph.SR from 2026-12-02 on
- the DR4 documentation and the papers expected with the release
- Google Scholar for "Proxima" "DR4"
