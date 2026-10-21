"""
cortado: fit PAH features, atomic/molecular lines, a stellar+dust continuum
and mixed-geometry dust + H2O-ice attenuation to a JWST/MIRI LRS spectrum,
optionally joined with a Spitzer/IRS spectrum and broad-band photometry
(integrated through the real filter curves).

cortado is this file, measure_pahs_general.py, with its data files (see
README.md). It began as a generalized, cleaned-up version of the LRSz2 fitter
measure_pahs_lrs_and_irs.py, nicknamed "CAFE at home", hence the name. The
physical model is the same:

    F_nu(lam) = A(lam) * [ continuum + sum_k PAH_k + sum_j line_j ]
    A         = (1 - exp(-tau)) / tau                     (mixed geometry)
    tau       = tau97 * kvt(lam) + tau_ice305 * Drude(3.05 um, FWHM 0.4 um)

lam is the rest wavelength [um] and F_nu the observed flux density [mJy].
Continuum: SB99 10 Myr stars + power law ('stellar_pl'), or a polynomial
in lam ('poly', used for the AGN-dominated sources).
Blended PAH complexes (7.7, 11.3, 12.7 um; PAH_COMPLEXES) are also reported
as totals, with errors that include the covariance between their members.

Partly covered features
-----------------------
By default a PAH is fit only if the data reach half an (instrument-broadened)
FWHM on both sides of its centre. With partial_pahs=True (--partial-pahs) a
feature is also fit when at least min_coverage (default 0.3) of its profile
flux lands on the data, e.g. the blue side of the 7.7 um complex when the
spectrum stops at ~8 um rest. Its centre is then fixed, so the flux comes from
the covered part of the profile plus its fixed shape. The table flags these
(partial, coverage). Lines and the fully-covered criterion are unchanged.

Tested on LRSz2 galaxies whose LRS ends at 7.6-8.5 um rest, against their
LRS+IRS fits and with noiseless injections of those fits:
  * Report the complex total (PAH7.7), not its components: the blue side of
    the complex cannot separate 7.60 / 7.85 / 8.33 um (components off by up
    to x3 even without noise).
  * The total is good to ~2% when every member has coverage >~ 0.85 (data
    reaching ~8.5 um rest). At coverage ~0.5 (data stopping at ~7.9 um, the
    7.85 um peak) it carries a ~10% systematic, set by how the unfit 8.33 um
    wing is absorbed: +11% at min_coverage = 0.3, -10% at 0.1. This is NOT in
    the error bar; refit at both thresholds to bracket it.
  * Below 0.3 the 7.85 um member is not fit and no total is reported.
  * Without the option, a fully covered 7.60 um next to an unfit 7.85 um
    absorbs its wing (x8 too bright in the injection test).

Instrumental broadening
-----------------------
The old script fit the intrinsic (Smith+07) PAH profiles directly to the
data, i.e. it assumed infinite resolving power. R is dimensionless, so a
source at z is seen with rest-frame LSF FWHM = lam_rest / R(lam_rest (1+z)).
For the 3.3 um PAH (intrinsic FWHM 0.043 um) that is 0.059 um at z = 0.7
(R_LRS(5.6 um) = 56), 0.041 um at z = 1 and 0.021 um at z = 2; for 11.3 um
in IRS LL1 at z = 1 it is 0.17 um vs 0.135 um intrinsic for the 11.23 um
component. The model is now forward-modelled:

  * continuum, attenuation and PAHs are evaluated on a fine rest-frame grid
    and convolved with a Gaussian LSF whose width is set *per data pixel* by
    that pixel's instrument and R(lam_obs). The convolution is one sparse
    matrix, built once per spectrum.
  * lines are fit to the MIRI LRS pixels only (line_insts, --line-insts): the
    Spitzer/IRS spectra are only sensitive to the continuum and PAHs, so IRS
    pixels see no line emission and a line covered only by IRS is not fit.
  * lines are Gaussians of intrinsic FWHM `line_fwhm_kms` (default 300 km/s,
    i.e. unresolved), convolved analytically with the same per-pixel LSF.

Convolution conserves flux, so integrated fluxes come straight from the
intrinsic profile parameters.

R(lam_obs):
  LRS  P750L curve from the pandeia dispersion file in templates/
       (Kendrew+15: R ~ 40 at 5 um, ~160 at 10 um, ~250 at 14 um).
  IRS  Instrument Handbook ranges, linear in lam within each order:
       SL2 60-127 (5.2-7.7), SL1 61-120 (7.4-14.5), LL2 57-126 (14.0-21.3),
       LL1 58-112 (19.5-38.0 um); order switches at 7.5, 14.0 and 19.5 um.
Either can be replaced by a 4th column (R) in the input file or a callable,
and scaled with R_scale for systematics tests.

Photometry
----------
Photometric points are not treated like spectral pixels. A band given by a
filter name is the full model F_nu (attenuated continuum + PAHs + lines, at
intrinsic resolution) integrated through the filter, in the facility's
calibration convention:

    F_band = int F_nu(lam) T(lam) w(lam) dlam / int s(lam) T(lam) w(lam) dlam

with lam observed, T the SVO Filter Profile Service transmission, w = 1/lam
(photon counter) or 1/lam^2 (energy counter, SVO DetectorType), and s the
reference spectrum F_nu,ref(lam)/F_nu,ref(lam_ref): flat for JWST, a 10^4 K
blackbody for MIPS, nu F_nu = const for IRAC and PACS (FILTER_CALIBRATION).
F_band is therefore the quoted, not colour-corrected, catalogue value. Give
bands as (filter, flux, err) with a filter alias (FILTER_ALIASES, e.g. MIPS24,
MIPS70, PACS70, PACS100, F770W) or any SVO ID; curves are downloaded from SVO
on first use and cached in filters/svo/. The MIPS conventions reproduce the
MIPS Handbook colour-correction tables to <= 0.003; PACS matches Balog et al.
(2014) to <= 1.6%.

A band given as a number (e.g. 24, or the text '24') is NOT matched to any
filter: it is the observed wavelength in um of a delta-function measurement,
and the model is the value of the attenuated continuum + PAHs at exactly
lam_rest = lam_obs / (1+z), at intrinsic resolution (no filter, no LSF, no
colour correction). Emission lines do not contribute to such a point: at their
~300 km/s intrinsic width they would matter only within a hair of line centre.

Bugs fixed relative to measure_pahs_lrs_and_irs.py
--------------------------------------------------
* Luminosities: the old code multiplied by (1+z) where it should divide.
  L = 4 pi d_L^2 F_obs with F_obs = int F_nu dnu_obs = int F_nu dnu_rest/(1+z),
  so the old L were high by (1+z)^2 (0.4-1.1 dex across LRSz2).
* log-L errors were 0.434*log10(sigma_L)/log10(L); they are now
  0.434*sigma_L/L. Errors now come from the fit covariance (see
  Uncertainties) instead of the mean S/N near the feature, which was ~5x
  too small at 3.3 um.
* Photometry was evaluated monochromatically at the nominal wavelength; it is
  now integrated through the filter (see Photometry).
* The Chebyshev continuum was evaluated on raw lam (no domain mapping); it is
  now mapped to [-1, 1] over the fitted range (same cubic-in-lam model).
* (1 - e^-tau)/tau is finite at tau = 0 (tau97 = 0 is an allowed bound).
* Features need data on both sides (the old check only tested the red edge).
* No paths outside this directory; results never overwrite lrsz2_pah33.pkl.

Uncertainties
-------------
The stated pixel errors under-describe the misfit (chi2_nu ~ 3-20 for LRSz2),
so covariance errors are inflated per feature: by sqrt(max(1, chi2_nu)) within
+-1 observed FWHM of the feature, times sqrt(1 + 2 sum_k (1-k/4) rho_k) from the
lag-k autocorrelation of the normalized residuals within +-3 FWHM. Checked on
LRSz2 against a 4-pixel block residual bootstrap (100 refits per source): a
single global chi2_nu factor is right on average but ~1.4x too small at 3.3 um
and ~2x too large at 7.7 um. The local scheme removes those trends and sits
~1.15x above the bootstrap in the median. The correlation factor sets the
overall level but does not predict the per-feature bootstrap value (r ~ 0.2),
so expect ~0.25 dex feature-to-feature scatter in the error estimates
themselves. Choose with error_scaling = 'local' (default) / 'global' / 'none'.

Usage
-----
    python measure_pahs_general.py spec.txt --z 1.44
    python measure_pahs_general.py lrs_x1d.fits --z 1.44 --irs irs.txt
    python measure_pahs_general.py combined.dat --z 1.44 --seam 14.01 \\
        --phot MIPS24:0.45:0.02 --phot PACS100:12:2
    python measure_pahs_general.py --lrsz2             # the whole LRSz2 sample
    python measure_pahs_general.py --lrsz2 GN_IRS12    # selected IDs

Text inputs: observed um, F_nu, error [, R]; F_nu units set by --flux-unit
(mJy default). FITS inputs: JWST x1d (EXTRACT1D WAVELENGTH/FLUX/FLUX_ERROR).

From Python:
    from measure_pahs_general import build_spectrum, fit_spectrum
    spec = build_spectrum('combined.dat', seam=14.01)
    res = fit_spectrum(spec, z=1.44, phot=[('MIPS24', 0.45, 0.02)])
    res.table; res.summary; res.phot_table(); res.plot('fit.png')   # model_range='data' for the old style
"""
__version__ = '2026.10.20'

import os
import argparse
import warnings
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.cosmology import LambdaCDM
from scipy import sparse
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

warnings.filterwarnings('ignore', category=RuntimeWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
LRS_DISP = os.path.join(HERE, 'templates', 'jwst_miri_p750l_disp_20170404135013.fits')
SB99_FILE = os.path.join(HERE, 'sb99-10myr.txt')
SPECDIR = os.path.join(HERE, 'specsForCafe')
OUTDIR = os.path.join(HERE, 'pah_fits')

COSMO = LambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725, Ode0=0.7)  # as in the original
C_UM = 2.99792458e14    # um / s
C_KMS = 2.99792458e5
MJY = 1e-29             # W m^-2 Hz^-1
LSUN = 3.828e26         # W
FWHM2SIG = 1. / (2 * np.sqrt(2 * np.log(2)))
LRS_MIN_WAVE = 4.77     # obs um; LRS wavelength calibration unreliable below
CORR_BLOCK = 4          # pixels; residual-correlation length for the error scaling
LOCAL_MIN_PIX = 10      # minimum pixels in a feature's local chi2 window


# ---------------------------------------------------------------- resolution
_disp = fits.getdata(LRS_DISP, 1)
_LRS_W, _LRS_R = np.array(_disp['WAVELENGTH']), np.array(_disp['R'])


def R_lrs(lam_obs):
	"""MIRI LRS (P750L) resolving power at observed wavelength [um]."""
	return np.interp(lam_obs, _LRS_W, _LRS_R)


# IRS Instrument Handbook: (lam_lo, lam_hi, R_lo, R_hi), R linear in lam per order
IRS_ORDERS = {'SL2': (5.2, 7.7, 60., 127.), 'SL1': (7.4, 14.5, 61., 120.),
	'LL2': (14.0, 21.3, 57., 126.), 'LL1': (19.5, 38.0, 58., 112.)}
IRS_SWITCH = [('SL2', 7.5), ('SL1', 14.0), ('LL2', 19.5), ('LL1', np.inf)]


def R_irs(lam_obs):
	"""Spitzer IRS low-res resolving power at observed wavelength [um]."""
	lam = np.atleast_1d(np.asarray(lam_obs, dtype=float))
	R = np.empty_like(lam)
	lo = -np.inf
	for order, hi in IRS_SWITCH:
		l0, l1, r0, r1 = IRS_ORDERS[order]
		m = (lam >= lo) & (lam < hi)
		R[m] = r0 + (r1 - r0) * (lam[m] - l0) / (l1 - l0)
		lo = hi
	return R


# ---------------------------------------------------------------- attenuation
def _kvt(x):
	"""CAFE's KVT silicate profile (+18 um Drude, beta = 0.1 power law)."""
	kvt_wav = np.array([8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.7, 9.75, 9.8,
		10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2,
		12.4, 12.6, 12.7])
	kvt_int = np.array([.06, .09, .16, .275, .415, .575, .755, .895, .98, .99, 1.0, .99,
		.94, .83, .745, .655, .58, .525, .43, .35, .27, .20, .13, .09, .06, .045,
		.04314])
	# extend to short wavelengths, stitched to the first tabulated point
	xs = x[x < kvt_wav[0]]
	short = kvt_int.min() * np.exp(2.03 * (xs - kvt_wav[0]))
	short *= kvt_int[0] / short.max()
	sx, sy = np.r_[xs, kvt_wav], np.r_[short, kvt_int]
	xl = x[x >= kvt_wav[-1]]
	hw = 0.247 * 18.0 / 2
	ext = np.r_[np.interp(x[x < kvt_wav[-1]], sx, sy), 0.4 * hw**2 / ((xl - 18.0)**2 + hw**2)]
	beta = 0.1
	return (1 - beta) * ext + beta * (9.7 / x)**1.7


# Same grid step and smoothing (2 pixels) as the original, which tabulated the
# curve on np.linspace(1, 100, 1000); extended here to 0.5-200 um.
_KVT_DX = 99. / 999.
_KVT_X = 1. + _KVT_DX * np.arange(-5, 2010)
_KVT_Y = gaussian_filter(_kvt(_KVT_X), 2)


def tau_silicate(lam):
	"""Normalized silicate optical depth (~1 at 9.7 um)."""
	return np.interp(lam, _KVT_X, _KVT_Y)


def tau_ice(lam, tau305, mean=3.05, std=0.20):
	"""H2O ice optical depth, Drude profile peaking at tau305."""
	g = 2 * std / mean
	return tau305 * g**2 / ((lam / mean - mean / lam)**2 + g**2)


def mixed_attenuation(tau):
	"""(1 - e^-tau) / tau, finite at tau = 0."""
	tau = np.asarray(tau, dtype=float)
	out = np.ones_like(tau)
	m = np.abs(tau) > 1e-8
	out[m] = -np.expm1(-tau[m]) / tau[m]
	return out


def d_mixed_attenuation(tau):
	"""d/dtau of (1 - e^-tau) / tau."""
	tau = np.asarray(tau, dtype=float)
	out = np.asarray(-0.5 + tau / 3.)
	m = np.abs(tau) > 1e-4
	out[m] = (np.exp(-tau[m]) * (1 + tau[m]) - 1) / tau[m]**2
	return out


# ---------------------------------------------------------------- stellar continuum
_sb = ascii.read(SB99_FILE, comment=';')
_SB_W = np.array(_sb['col2']) / 1e4
_SB_F = 10.**np.array(_sb['col4']) * np.array(_sb['col2'])**2   # F_lambda -> F_nu shape
_m = _SB_W > 0.9
_SB_W, _SB_F = _SB_W[_m], _SB_F[_m]
_SB_F = _SB_F / _SB_F[np.argmin(np.abs(_SB_W - 3.005))]           # = 1 at 3 um
_o = np.argsort(_SB_W)
_SB_W, _SB_F = _SB_W[_o], _SB_F[_o]


def sb99(lam):
	return np.interp(lam, _SB_W, _SB_F)


# ---------------------------------------------------------------- features
# PAHs: rest um, intrinsic FWHM um (Smith+07 / CAFE)
PAH_FEATURES = {
	'PAH33': (3.29, 0.043), 'PAH34': (3.40, 0.031), 'PAH35': (3.47, 0.100),
	'PAH52': (5.24, 0.058), 'PAH53': (5.33, 0.043), 'PAH57': (5.70, 0.200),
	'PAH62': (6.22, 0.187), 'PAH69': (6.69, 0.468),
	'PAH7a': (7.42, 0.935), 'PAH7b': (7.60, 0.334), 'PAH7c': (7.85, 0.416),
	'PAH8a': (8.33, 0.417), 'PAH8b': (8.61, 0.336),
	'PAH11a': (11.23, 0.135), 'PAH11b': (11.33, 0.363),
	'PAH12a': (12.62, 0.530), 'PAH12b': (12.69, 0.165),
	'PAH19': (18.92, 0.359), 'PAH21': (21.65, 0.359),
}
# the set fit by the original script (no 6.9 or 7.42 um)
DEFAULT_PAHS = ['PAH33', 'PAH34', 'PAH35', 'PAH52', 'PAH53', 'PAH57', 'PAH62',
	'PAH7b', 'PAH7c', 'PAH8a', 'PAH8b', 'PAH11a', 'PAH11b', 'PAH12a', 'PAH12b',
	'PAH19', 'PAH21']

# Blended complexes, reported as the sum of whichever members were requested;
# errors include the (anti)correlation of the members from the fit covariance.
# The default PAH list has no 7.42 um component, so PAH7.7 = 7.60 + 7.85 um.
PAH_COMPLEXES = {'PAH7.7': ['PAH7a', 'PAH7b', 'PAH7c'],
	'PAH11.3': ['PAH11a', 'PAH11b'], 'PAH12.7': ['PAH12a', 'PAH12b']}

# Lines: rest vacuum um
LINES = {
	'Paa': 1.87561, 'Brb': 2.62587, 'Pfd': 3.29699, 'Pfg': 3.74056, 'Bra': 4.05226,
	'Pfb': 4.65378, 'Pfa': 7.45986, 'Hua': 12.3719,
	'H2S7': 5.51116, 'H2S6': 6.10856, 'H2S5': 6.90952, 'H2S4': 8.02505,
	'H2S3': 9.66491, 'H2S2': 12.27861, 'H2S1': 17.03484, 'H2S0': 28.21883,
	'FeII5': 5.34017, 'MgVII': 5.50300, 'MgV': 5.60985, 'ArII': 6.98527,
	'NeVI': 7.65240, 'ArIII': 8.99138, 'SIV': 10.51050, 'NeII': 12.81355,
	'NeV14': 14.32170, 'NeIII': 15.55510, 'SIII18': 18.71300, 'NeV24': 24.31750,
	'OIV': 25.89030, 'FeII26': 25.98829, 'SIII33': 33.48100, 'SiII': 34.81520,
}
# Lines that are not blended with a PAH or another line at R ~ 40-250. Left out:
# Pfd (on PAH33), FeII5 (PAH53), MgVII (H2S7), NeVI and Pfa (7.7 complex),
# Hua (H2S2/PAH12), FeII26 (OIV), and the weak Pfg/Pfb/MgV.
DEFAULT_LINES = ['Paa', 'Brb', 'Bra', 'H2S7', 'H2S6', 'H2S5', 'H2S4', 'H2S3', 'H2S2',
	'H2S1', 'H2S0', 'ArII', 'ArIII', 'SIV', 'NeII', 'NeV14', 'NeIII', 'SIII18',
	'NeV24', 'OIV', 'SIII33', 'SiII']


def lorentz(lam, amp, x0, fwhm):
	h = fwhm / 2
	return amp * h * h / ((lam - x0)**2 + h * h)


def drude(lam, amp, x0, fwhm):
	g = fwhm / x0
	return amp * g * g / ((lam / x0 - x0 / lam)**2 + g * g)


PROFILES = {'lorentz': lorentz, 'drude': drude}


# ---------------------------------------------------------------- photometry
# Filter curves come from the SVO Filter Profile Service and are cached as VOTables in
# filters/svo/ (downloaded on first use of a filter). A band measurement is modelled as
#     F_band = int F_nu(lam) T(lam) w(lam) dlam / int s(lam) T(lam) w(lam) dlam
# T: SVO transmission vs observed wavelength; w = 1/lam for photon counters and 1/lam^2
# for energy counters (SVO DetectorType 1 / 0); s: the facility's calibration reference
# spectrum in F_nu, normalized to 1 at its reference wavelength. F_band is the quoted,
# not colour-corrected, catalogue flux density, so catalogue values are fit as published.
#   JWST           flat F_nu (Gordon et al. 2022, AJ 163, 267, Sec. 4.2.3)
#   Spitzer/MIPS   10^4 K blackbody at 23.675, 71.44, 155.899 um (MIPS Instrument Handbook)
#   Spitzer/IRAC   nu F_nu = const at 3.550, 4.493, 5.731, 7.872 um
#   Herschel/PACS  nu F_nu = const at 70, 100, 160 um
# Checked: the MIPS 24/70 colour corrections this implies reproduce the MIPS Handbook
# tables (power laws and blackbodies >= 50 K) to <= 0.003; PACS 4000 K corrections match
# Balog et al. (2014) to 0.2/0.5/1.6% at 70/100/160 um.
FILTER_DIR = os.path.join(HERE, 'filters', 'svo')
SVO_FPS = 'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='
FILTER_ALIASES = {
	'MIPS24': 'Spitzer/MIPS.24mu', 'MIPS70': 'Spitzer/MIPS.70mu', 'MIPS160': 'Spitzer/MIPS.160mu',
	'PACS70': 'Herschel/Pacs.blue', 'PACS100': 'Herschel/Pacs.green', 'PACS160': 'Herschel/Pacs.red',
	'IRAC1': 'Spitzer/IRAC.I1', 'IRAC2': 'Spitzer/IRAC.I2', 'IRAC3': 'Spitzer/IRAC.I3', 'IRAC4': 'Spitzer/IRAC.I4',
}
FILTER_ALIASES.update({f: 'JWST/MIRI.' + f for f in
	('F560W', 'F770W', 'F1000W', 'F1130W', 'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W')})
# (reference spectrum, reference wavelength [um]); JWST filters default to flat F_nu
FILTER_CALIBRATION = {
	'Spitzer/MIPS.24mu': ('bb10000', 23.675), 'Spitzer/MIPS.70mu': ('bb10000', 71.44),
	'Spitzer/MIPS.160mu': ('bb10000', 155.899),
	'Herschel/Pacs.blue': ('nufnu', 70.), 'Herschel/Pacs.green': ('nufnu', 100.),
	'Herschel/Pacs.red': ('nufnu', 160.),
	'Spitzer/IRAC.I1': ('nufnu', 3.550), 'Spitzer/IRAC.I2': ('nufnu', 4.493),
	'Spitzer/IRAC.I3': ('nufnu', 5.731), 'Spitzer/IRAC.I4': ('nufnu', 7.872),
}


def _bnu(lam, T):
	"""Planck B_nu shape (arbitrary normalization), lam in um."""
	return lam**-3 / np.expm1(14387.769 / (lam * T))


@dataclass
class Filter:
	name: str             # SVO filter ID
	wave: np.ndarray      # observed um
	trans: np.ndarray
	photon: bool          # photon counter (else energy counter)
	reference: str        # 'flat' | 'nufnu' | 'bb10000'
	lam_ref: float        # um

	def ref_shape(self, lam):
		if self.reference == 'flat':
			return np.ones_like(lam)
		if self.reference == 'nufnu':
			return lam / self.lam_ref
		return _bnu(lam, 1e4) / _bnu(self.lam_ref, 1e4)

	def grid(self, dln=5e-4, cut=1e-3):
		"""Log-spaced observed wavelengths spanning T > cut * max(T)."""
		ok = self.trans > cut * self.trans.max()
		lo, hi = self.wave[ok].min(), self.wave[ok].max()
		return np.exp(np.arange(np.log(lo), np.log(hi) + dln, dln))

	def _density(self, lam):
		return np.interp(lam, self.wave, self.trans, left=0., right=0.) * lam**(-1. if self.photon else -2.)

	@property
	def norm(self):
		"""int s T w dlam on the band grid."""
		if not hasattr(self, '_norm'):
			g = self.grid()
			self._norm = np.sum(self._density(g) * self.ref_shape(g) * np.gradient(g))
		return self._norm

	def weights(self, lam):
		"""W such that F_band = sum W F_nu(lam), for ascending observed lam (um)."""
		return self._density(lam) * np.gradient(lam) / self.norm

	def density(self, lam):
		"""Weight per unit observed wavelength (for unresolved lines)."""
		return self._density(np.asarray(lam, float)) / self.norm

	def synphot(self, lam, fnu):
		"""Quoted band flux density of the spectrum fnu(lam) (observed um)."""
		return float(np.sum(self.weights(lam) * fnu))

	def halfmax(self):
		"""Observed wavelengths where T crosses half its maximum."""
		ok = self.wave[self.trans >= 0.5 * self.trans.max()]
		return ok.min(), ok.max()


_FILTERS = {}


def get_filter(name):
	"""Filter by alias (e.g. 'MIPS24', 'F770W'; case-insensitive) or SVO ID (e.g. 'Spitzer/MIPS.24mu')."""
	name = name.strip()
	aliases = {k.lower(): v for k, v in FILTER_ALIASES.items()}
	fid = aliases.get(name.lower(), name)
	if fid in _FILTERS:
		return _FILTERS[fid]
	path = os.path.join(FILTER_DIR, fid.replace('/', '_') + '.xml')
	if not os.path.isfile(path):
		from urllib.request import urlopen
		from urllib.parse import quote
		try:
			with urlopen(SVO_FPS + quote(fid, safe='/.'), timeout=60) as r:
				data = r.read()
		except OSError as err:
			raise ValueError('filter %r is not in %s and could not be downloaded from SVO (%s)'
				% (fid, FILTER_DIR, err)) from None
		if b'<TABLE' not in data:
			raise ValueError('unknown photometric band %r: give a number (observed wavelength in um), '
				'a filter alias (%s) or an SVO filter ID such as JWST/MIRI.F770W'
				% (name, ', '.join(sorted(FILTER_ALIASES))))
		os.makedirs(FILTER_DIR, exist_ok=True)
		with open(path, 'wb') as f:
			f.write(data)
	from astropy.io.votable import parse
	tab = parse(path).get_first_table()
	par = {p.name: p for p in tab.params}
	t = tab.to_table()
	wave = np.array(t['Wavelength'], float) * u.Unit(t['Wavelength'].unit or 'Angstrom').to(u.um)
	ref, lam_ref = FILTER_CALIBRATION.get(fid, (None, None))
	if ref is None:
		if not fid.startswith('JWST/'):
			warnings.warn('no calibration convention for %s; assuming flat F_nu' % fid)
		pv = par['WavelengthPivot']
		ref, lam_ref = 'flat', float(pv.value) * u.Unit(pv.unit or 'Angstrom').to(u.um)
	flt = Filter(fid, wave, np.array(t['Transmission'], float), int(par['DetectorType'].value) == 1,
		ref, lam_ref)
	_FILTERS[fid] = flt
	return flt


# ---------------------------------------------------------------- data
@dataclass
class Spectrum:
	wave: np.ndarray    # observed um, ascending
	flux: np.ndarray    # mJy
	err: np.ndarray     # mJy
	inst: np.ndarray    # 'LRS' / 'IRS' per pixel
	R: np.ndarray       # resolving power per pixel


def read_table(path, flux_unit='mJy'):
	"""(wave_obs [um], flux [mJy], err [mJy], R or None) from text or JWST x1d FITS."""
	if path.lower().endswith(('.fits', '.fit', '.fits.gz')):
		with fits.open(path) as h:
			hdu = h['EXTRACT1D']
			d = hdu.data
			unit = hdu.columns['FLUX'].unit or 'Jy'
			scale = u.Unit(unit).to(u.mJy)
			return (np.array(d['WAVELENGTH'], float), np.array(d['FLUX'], float) * scale,
				np.array(d['FLUX_ERROR'], float) * scale, None)
	s = np.loadtxt(path, ndmin=2)
	scale = u.Unit(flux_unit).to(u.mJy)
	R = s[:, 3] if s.shape[1] > 3 else None
	return s[:, 0], s[:, 1] * scale, s[:, 2] * scale, R


def build_spectrum(lrs, irs=None, seam=None, flux_unit='mJy', lrs_min_wave=LRS_MIN_WAVE,
		R_lrs=R_lrs, R_irs=R_irs, R_scale=1.):
	"""
	Combine an LRS spectrum with an optional IRS spectrum.

	lrs, irs : file path or (wave_obs, flux_mJy, err_mJy[, R]) tuple
	seam     : if set, `lrs` is a combined file; pixels at lam_obs > seam are IRS
	R_lrs, R_irs : callables R(lam_obs), used where no R column is given
	R_scale  : multiply every R by this (systematics tests)
	"""
	def load(x):
		t = read_table(x, flux_unit) if isinstance(x, str) else tuple(x) + (None,) * (4 - len(x))
		w, f, e, R = (None if a is None else np.asarray(a, float) for a in t[:4])
		return w, f, e, R

	parts = []
	w, f, e, R = load(lrs)
	if seam is not None:
		if irs is not None:
			raise ValueError('give either seam or irs, not both')
		m = w > seam
		parts.append(('IRS', w[m], f[m], e[m], None if R is None else R[m]))
		w, f, e, R = w[~m], f[~m], e[~m], None if R is None else R[~m]
	parts.append(('LRS', w, f, e, R))
	if irs is not None:
		parts.append(('IRS',) + load(irs))

	out = {k: [] for k in ('wave', 'flux', 'err', 'inst', 'R')}
	for inst, w, f, e, R in parts:
		good = np.isfinite(w) & np.isfinite(f) & np.isfinite(e) & (e > 0)
		if inst == 'LRS':
			good &= w >= lrs_min_wave
		if R is None:
			R = (R_lrs if inst == 'LRS' else R_irs)(w)
		out['wave'].append(w[good])
		out['flux'].append(f[good])
		out['err'].append(e[good])
		out['inst'].append(np.full(good.sum(), inst))
		out['R'].append(np.asarray(R, float)[good] * R_scale)
	cat = {k: np.concatenate(v) for k, v in out.items()}
	o = np.argsort(cat['wave'])
	return Spectrum(**{k: v[o] for k, v in cat.items()})


def segments(wave, inst, gap=4.):
	"""Index arrays of contiguous runs, split at instrument changes and gaps."""
	d = np.diff(wave)
	ref = np.minimum(np.r_[d[0], d[:-1]], np.r_[d[1:], d[-1]])
	brk = np.where((inst[1:] != inst[:-1]) | (d > gap * ref))[0] + 1
	return np.split(np.arange(wave.size), brk)


# ---------------------------------------------------------------- LSF operator
def fine_grid(lam, sig, min_hwhm_ln=np.inf, nsig=6.):
	"""Log-spaced rest grid resolving both the LSF and the narrowest PAH."""
	dln = np.clip(min(np.min(sig / lam), min_hwhm_ln) / 3., 1e-4, 2e-3)
	lo, hi = np.min(lam - nsig * sig), np.max(lam + nsig * sig)
	return np.exp(np.arange(np.log(lo), np.log(hi) + dln, dln))


def lsf_matrix(lam, sig, grid, nsig=5.):
	"""
	Sparse (n_pix x n_grid) operator: a model sampled on `grid` -> the same model
	seen through a Gaussian LSF of width sig[j] centred on pixel lam[j]. Rows sum
	to 1 (flux conserving). Where the LSF is narrower than the grid step the row
	is plain linear interpolation.
	"""
	dg = np.gradient(grid)
	rows, cols, vals = [], [], []
	for j, (l, s) in enumerate(zip(lam, sig)):
		i = np.clip(np.searchsorted(grid, l), 1, grid.size - 1)
		if s < dg[i]:
			t = (l - grid[i - 1]) / (grid[i] - grid[i - 1])
			idx, w = np.array([i - 1, i]), np.array([1 - t, t])
		else:
			a, b = np.searchsorted(grid, [l - nsig * s, l + nsig * s])
			idx = np.arange(a, b)
			w = np.exp(-0.5 * ((grid[idx] - l) / s)**2) * dg[idx]
			w /= w.sum()
		rows.append(np.full(idx.size, j))
		cols.append(idx)
		vals.append(w)
	return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
		shape=(lam.size, grid.size))


# ---------------------------------------------------------------- model
class MIRModel:
	"""
	Spectral model for one source. Parameters (free ones are fit):
	  cont_*            continuum ('stellar_pl': logAstar, logApl, alpha; 'poly': c0..cN)
	  tau97, tau_ice305 attenuation
	  <PAH>_amp, _x0    intrinsic peak F_nu [mJy] and centre [rest um]; FWHM fixed
	  <line>_flux, _x0  intrinsic int F_nu dlam_rest [mJy um] and centre
	"""

	def __init__(self, spec, z, phot=None, pahs=DEFAULT_PAHS, lines=DEFAULT_LINES,
			continuum='stellar_pl', poly_degree=3, profile='lorentz', fit_ice=True,
			line_fwhm_kms=300., pah_shift_tol=0.01, line_shift_tol=0.003,
			partial_pahs=False, min_coverage=0.3, line_insts=('LRS',)):
		if continuum not in ('stellar_pl', 'poly'):
			raise ValueError(continuum)
		self.spec, self.z = spec, z
		self.continuum, self.poly_degree = continuum, poly_degree
		self.profile = PROFILES[profile]
		self.profile_name = profile
		self.line_fwhm_kms = line_fwhm_kms

		self.lam = spec.wave / (1 + z)
		self.sig = self.lam / spec.R * FWHM2SIG        # rest-frame LSF sigma per pixel
		self._setup_photometry([] if phot is None else list(phot))   # list, tuple or (n, 3) array
		self.y = np.r_[spec.flux, self.phot[:, 1]]
		self.e = np.r_[spec.err, self.phot[:, 2]]
		self.segments = segments(spec.wave, spec.inst)
		lam_all = np.r_[self.lam, self.phot_grid]
		self.domain = (lam_all.min(), lam_all.max())

		self.requested = {'pah': list(pahs), 'line': list(lines)}
		# PAHs: fully covered ones always; with partial_pahs, also those with >= min_coverage
		# of their (broadened) profile flux on the data, e.g. the blue half of 7.7 um
		self.coverage = {p: self.coverage_fraction(*PAH_FEATURES[p]) for p in pahs}
		full = [p for p in pahs if self.covered(*PAH_FEATURES[p])]
		self.partial = [p for p in pahs if partial_pahs and p not in full
			and self.coverage[p] >= min_coverage]
		self.pahs = [p for p in pahs if p in full or p in self.partial]
		# atomic / H2 lines only where the spectrum is sensitive to them (line_insts):
		# fit only if covered by those pixels, and modelled only at those pixels
		self.line_insts = tuple(line_insts)
		self.line_pix = np.isin(spec.inst, self.line_insts)
		self.lines = [l for l in lines if self.covered(LINES[l], 0., self.line_pix)]
		# complexes with >= 2 requested members; reported only if all of them are fit
		self.requested['complex'] = []
		self.complexes = {}
		for c, mem in PAH_COMPLEXES.items():
			req = [x for x in mem if x in pahs]
			if len(req) >= 2:
				self.requested['complex'].append(c)
				if all(x in self.pahs for x in req):
					self.complexes[c] = req

		hwhm = [PAH_FEATURES[p][1] / 2 / PAH_FEATURES[p][0] for p in self.pahs]
		self.grid = fine_grid(self.lam, self.sig, min(hwhm, default=np.inf))
		self.K = lsf_matrix(self.lam, self.sig, self.grid)
		self._setup_params(fit_ice, pah_shift_tol, line_shift_tol)

	# ---- photometry
	def _setup_photometry(self, phot):
		"""
		phot: list / tuple / (n, 3) array of (band, flux_mJy, err_mJy). band is a filter alias
		or SVO ID, integrated through that filter; or a number (or numeric text) = observed um,
		a delta function: the model value at exactly that wavelength, never matched to a filter.
		Builds self.W (n_band x n_grid): F_band = W @ F_nu(phot_grid), phot_grid in rest um.
		"""
		z = self.z
		self.filters, self.phot_names, rows, grids, vals = [], [], [], [], []
		ref = []
		for i, (band, f, e) in enumerate(phot):
			if isinstance(band, str):
				try:
					band = float(band)      # a wavelength given as text, e.g. '24'
				except ValueError:
					pass
			if isinstance(band, str):
				flt = get_filter(band)
				g = flt.grid()
				w = flt.weights(g)
				ref.append(flt.lam_ref)
				self.phot_names.append(band)
			else:
				flt, g, w = None, np.array([float(band)]), np.array([1.])
				ref.append(float(band))
				self.phot_names.append('%g um' % float(band))
			self.filters.append(flt)
			rows.append(np.full(g.size, i))
			grids.append(g / (1 + z))
			vals.append(w)
		n = len(phot)
		self.phot = np.array([[r, f, e] for r, (_, f, e) in zip(ref, phot)], float).reshape(-1, 3)
		self.lam_phot = self.phot[:, 0] / (1 + z)          # reference wavelengths, rest um
		self.phot_grid = np.concatenate(grids) if n else np.zeros(0)
		self.W = sparse.csr_matrix((np.concatenate(vals) if n else [], (np.concatenate(rows) if n else [],
			np.arange(self.phot_grid.size))), shape=(n, self.phot_grid.size))

	def phot_line_unit(self, x0):
		"""Band response to a unit (1 mJy um, rest), unattenuated line at rest x0."""
		lo = x0 * (1 + self.z)
		return np.array([(1 + self.z) * float(flt.density(lo)) if flt else 0. for flt in self.filters])

	def photometry(self, p):
		"""Model band flux densities (quoted, not colour-corrected)."""
		m = self.W @ self.unconvolved(p, self.phot_grid)
		for n in self.lines:
			x0, flux = p[self.ix[n + '_x0']], p[self.ix[n + '_flux']]
			m = m + flux * mixed_attenuation(self.tau(p, x0)) * self.phot_line_unit(x0)
		return m

	# ---- coverage
	def lsf_fwhm(self, lam0):
		"""Rest-frame LSF FWHM at lam0, from the nearest data pixel."""
		return self.sig[np.argmin(np.abs(self.lam - lam0))] / FWHM2SIG

	def coverage_fraction(self, lam0, fwhm):
		"""Fraction of a Lorentzian feature's flux (FWHM and LSF added in quadrature) on the data."""
		w = np.hypot(fwhm, self.lsf_fwhm(lam0))
		cdf = lambda x: np.arctan(2 * (x - lam0) / w) / np.pi
		return float(sum(cdf(self.lam[s[-1]]) - cdf(self.lam[s[0]]) for s in self.segments))

	def covered(self, lam0, fwhm, pix=None):
		"""Data (pixels `pix`, default all) reach +-W/2 of lam0 (W = observed FWHM) with no
		gap wider than W."""
		lam = self.lam if pix is None else self.lam[pix]
		if lam.size == 0:
			return False
		w = np.hypot(fwhm, self.lsf_fwhm(lam0))
		lo, hi = lam0 - w / 2, lam0 + w / 2
		if lam[0] > lo or lam[-1] < hi:
			return False
		inside = lam[(lam >= lo) & (lam <= hi)]
		return inside.size > 0 and np.diff(np.r_[lo, inside, hi]).max() <= w

	# ---- parameters
	def _setup_params(self, fit_ice, pah_tol, line_tol):
		names, val, lo, hi, free = [], [], [], [], []

		def add(n, v, l=-np.inf, h=np.inf, f=True):
			names.append(n); val.append(v); lo.append(l); hi.append(h); free.append(f)

		norm = np.nanmedian(self.spec.flux)
		norm = norm if norm > 0 else np.nanmax(self.spec.flux)
		ln = np.log10(norm)
		if self.continuum == 'stellar_pl':
			add('cont_logAstar', ln - 0.5, ln - 4, ln + 2)
			add('cont_logApl', ln, ln - 4, ln + 4)
			add('cont_alpha', 2., 0., 5.)
		else:
			for k in range(self.poly_degree + 1):
				add('cont_c%d' % k, norm if k == 0 else 0.)
		add('tau97', 0.3, 0., 5.)
		add('tau_ice305', 0.5 if fit_ice else 0., 0., 10., fit_ice)
		for p in self.pahs:
			x0 = PAH_FEATURES[p][0]
			add(p + '_amp', 0.1 * norm, 0., np.inf)
			# a partly covered feature cannot constrain its centre: keep it fixed
			add(p + '_x0', x0, x0 * (1 - pah_tol), x0 * (1 + pah_tol), pah_tol > 0 and p not in self.partial)
		for l in self.lines:
			x0 = LINES[l]
			add(l + '_flux', 0.05 * norm * self.lsf_fwhm(x0), 0., np.inf)
			add(l + '_x0', x0, x0 * (1 - line_tol), x0 * (1 + line_tol), line_tol > 0)

		self.names = names
		self.ix = {n: i for i, n in enumerate(names)}
		self.p_init = np.array(val)
		self.free = np.array(free)
		self.lo, self.hi = np.array(lo)[self.free], np.array(hi)[self.free]
		span = np.where(np.isfinite(self.hi - self.lo), self.hi - self.lo, 1.)
		eps = 1e-6 * np.maximum(span, 1e-30)
		self.q_init = np.clip(self.p_init[self.free], self.lo + eps, self.hi - eps)

	def full(self, q):
		p = self.p_init.copy()
		p[self.free] = q
		return p

	# ---- evaluation
	def intrinsic_continuum(self, p, lam):
		if self.continuum == 'stellar_pl':
			return (10.**p[self.ix['cont_logAstar']] * sb99(lam)
				+ 10.**p[self.ix['cont_logApl']] * (lam / 10.)**p[self.ix['cont_alpha']])
		a, b = self.domain
		c = [p[self.ix['cont_c%d' % k]] for k in range(self.poly_degree + 1)]
		return np.polynomial.chebyshev.chebval((2 * lam - a - b) / (b - a), c)

	def tau(self, p, lam):
		return p[self.ix['tau97']] * tau_silicate(lam) + tau_ice(lam, p[self.ix['tau_ice305']])

	def pah(self, p, name, lam):
		return self.profile(lam, p[self.ix[name + '_amp']], p[self.ix[name + '_x0']],
			PAH_FEATURES[name][1])

	def line_profile(self, x0, lam, sig):
		"""Unit-flux Gaussian line of intrinsic FWHM line_fwhm_kms seen through LSF sigma."""
		s2 = sig**2 + (x0 * self.line_fwhm_kms / C_KMS * FWHM2SIG)**2
		return np.exp(-0.5 * (lam - x0)**2 / s2) / np.sqrt(2 * np.pi * s2)

	def line(self, p, name, lam, sig):
		"""Attenuated line seen through LSF sigma `sig` at pixels `lam`."""
		x0, flux = p[self.ix[name + '_x0']], p[self.ix[name + '_flux']]
		return flux * self.line_profile(x0, lam, sig) * mixed_attenuation(self.tau(p, x0))

	def unconvolved(self, p, lam):
		"""Attenuated continuum + PAHs at infinite resolution."""
		f = self.intrinsic_continuum(p, lam)
		for n in self.pahs:
			f = f + self.pah(p, n, lam)
		return mixed_attenuation(self.tau(p, lam)) * f

	def model(self, p):
		"""Model at the data: spectral pixels then photometry."""
		m = self.K @ self.unconvolved(p, self.grid)
		lp = self.line_pix
		for n in self.lines:
			m[lp] += self.line(p, n, self.lam[lp], self.sig[lp])
		return np.r_[m, self.photometry(p)]

	def _unconvolved_jac(self, p, lam):
		"""d unconvolved(p, lam) / dp, shape (lam.size, n_params)."""
		ix = self.ix
		J = np.zeros((lam.size, len(self.names)))
		tau = self.tau(p, lam)
		att = mixed_attenuation(tau)
		S = self.intrinsic_continuum(p, lam)
		for n in self.pahs:
			x0, amp, fw = p[ix[n + '_x0']], p[ix[n + '_amp']], PAH_FEATURES[n][1]
			shape = self.profile(lam, 1., x0, fw)
			S = S + amp * shape
			h = 1e-6 * x0
			J[:, ix[n + '_amp']] = att * shape
			J[:, ix[n + '_x0']] = att * amp * (self.profile(lam, 1., x0 + h, fw)
				- self.profile(lam, 1., x0 - h, fw)) / (2 * h)
		dS = d_mixed_attenuation(tau) * S
		J[:, ix['tau97']] = dS * tau_silicate(lam)
		J[:, ix['tau_ice305']] = dS * tau_ice(lam, 1.)
		if self.continuum == 'stellar_pl':
			pl = 10.**p[ix['cont_logApl']] * (lam / 10.)**p[ix['cont_alpha']]
			J[:, ix['cont_logAstar']] = att * np.log(10) * 10.**p[ix['cont_logAstar']] * sb99(lam)
			J[:, ix['cont_logApl']] = att * np.log(10) * pl
			J[:, ix['cont_alpha']] = att * pl * np.log(lam / 10.)
		else:
			a, b = self.domain
			x = (2 * lam - a - b) / (b - a)
			for k in range(self.poly_degree + 1):
				J[:, ix['cont_c%d' % k]] = att * np.polynomial.chebyshev.chebval(x, np.eye(k + 1)[k])
		return J

	def jacobian(self, p):
		"""d model(p) / dp for the free parameters (analytic; x0 terms by central difference)."""
		ix = self.ix
		J = self.K @ self._unconvolved_jac(p, self.grid)
		lp = self.line_pix
		lam, sig = self.lam[lp], self.sig[lp]
		for n in self.lines:
			x0, flux = p[ix[n + '_x0']], p[ix[n + '_flux']]
			g = self.line_profile(x0, lam, sig)
			t0 = self.tau(p, x0)
			J[lp, ix[n + '_flux']] = g * mixed_attenuation(t0)
			h = 1e-6 * x0
			pp, pm = p.copy(), p.copy()
			pp[ix[n + '_x0']] += h
			pm[ix[n + '_x0']] -= h
			J[lp, ix[n + '_x0']] = (self.line(pp, n, lam, sig) - self.line(pm, n, lam, sig)) / (2 * h)
			dg = flux * g * d_mixed_attenuation(t0)
			J[lp, ix['tau97']] += dg * tau_silicate(x0)
			J[lp, ix['tau_ice305']] += dg * tau_ice(x0, 1.)
		# photometry: filter-weighted sums of the same derivatives, plus the lines
		Jp = np.asarray(self.W @ self._unconvolved_jac(p, self.phot_grid)).reshape(self.W.shape[0], len(self.names))
		for n in self.lines:
			x0, flux = p[ix[n + '_x0']], p[ix[n + '_flux']]
			t0 = self.tau(p, x0)
			u0 = self.phot_line_unit(x0)
			Jp[:, ix[n + '_flux']] += mixed_attenuation(t0) * u0
			h = 1e-6 * x0
			Jp[:, ix[n + '_x0']] += flux * (mixed_attenuation(self.tau(p, x0 + h)) * self.phot_line_unit(x0 + h)
				- mixed_attenuation(self.tau(p, x0 - h)) * self.phot_line_unit(x0 - h)) / (2 * h)
			du = flux * d_mixed_attenuation(t0) * u0
			Jp[:, ix['tau97']] += du * tau_silicate(x0)
			Jp[:, ix['tau_ice305']] += du * tau_ice(x0, 1.)
		return np.vstack([J, Jp])[:, self.free]

	def components(self, p, lam, sig, K=None, grid=None):
		"""Each convolved component at pixels (lam, sig); sum = total model. `grid` (default:
		the fit grid) must extend ~6 sigma beyond lam."""
		grid = self.grid if grid is None else grid
		K = lsf_matrix(lam, sig, grid) if K is None else K
		att = mixed_attenuation(self.tau(p, grid))
		out = {'continuum': K @ (att * self.intrinsic_continuum(p, grid))}
		for n in self.pahs:
			out[n] = K @ (att * self.pah(p, n, grid))
		for n in self.lines:
			out[n] = self.line(p, n, lam, sig)
		return out

	def model_curve(self, p, cut=0.01):
		"""
		Model components on one continuous rest-frame grid spanning the spectra and the
		photometric filters (where T > cut * max T). The LSF width at each wavelength is
		that of the nearest data pixel (R interpolated in lam_obs, clamped beyond the
		data), so the curve matches the fitted model where there is data.
		Returns lam, sig, components, covered (True inside a spectral data segment).
		"""
		lo, hi = self.lam.min(), self.lam.max()
		for f in self.filters:
			if f is not None:
				w = f.wave[f.trans > cut * f.trans.max()] / (1 + self.z)
				lo, hi = min(lo, w.min()), max(hi, w.max())
		if self.lam_phot.size:      # monochromatic points
			lo, hi = min(lo, self.lam_phot.min()), max(hi, self.lam_phot.max())
		dln = np.clip(np.min(self.sig / self.lam) / 2, 2e-4, 1e-3)
		lam = np.exp(np.arange(np.log(lo), np.log(hi) + dln, dln))
		R = np.interp(lam * (1 + self.z), self.spec.wave, self.spec.R)
		sig = lam / R * FWHM2SIG
		hwhm = [PAH_FEATURES[n][1] / 2 / PAH_FEATURES[n][0] for n in self.pahs]
		grid = fine_grid(lam, sig, min(hwhm, default=np.inf))
		comp = self.components(p, lam, sig, grid=grid)
		covered = np.zeros(lam.size, bool)
		for seg in self.segments:
			covered |= (lam >= self.lam[seg[0]]) & (lam <= self.lam[seg[-1]])
		return lam, sig, comp, covered

	# ---- derived quantities [all fluxes W/m2, EW rest um]
	def derived(self, p):
		z = self.z
		to_si = MJY * C_UM / (1 + z)     # mJy um / um^2 (rest) -> W/m2 (observed)
		out = {}
		feats = [(n, 'pah') for n in self.pahs] + [(n, 'line') for n in self.lines]
		for n, kind in feats:
			x0 = p[self.ix[n + '_x0']]
			if kind == 'pah':
				I = np.pi * p[self.ix[n + '_amp']] * PAH_FEATURES[n][1] / 2   # mJy um
			else:
				I = p[self.ix[n + '_flux']]
			fc = I * to_si / x0**2
			out[n + ':lam'] = x0
			out[n + ':flux_corr'] = fc
			out[n + ':flux_obs'] = fc * mixed_attenuation(self.tau(p, x0))
			out[n + ':ew'] = I / self.intrinsic_continuum(p, x0)
		for c, mem in self.complexes.items():
			for k in ('flux_corr', 'flux_obs', 'ew'):
				out[c + ':' + k] = sum(out[n + ':' + k] for n in mem)
		ice = p[self.ix['tau_ice305']]
		out['tau_305'] = self.tau(p, 3.05)
		out['tau_305_ice'] = tau_ice(3.05, ice)
		out['tau_97'] = self.tau(p, 9.7)
		return out


# ---------------------------------------------------------------- fitting
@dataclass
class FitResult:
	model: MIRModel
	p: np.ndarray            # full parameter vector
	cov: np.ndarray          # covariance of the free parameters
	chi2: float
	dof: int
	success: bool
	message: str
	table: Table = None      # one row per requested feature
	summary: dict = None     # attenuation, chi2, ...
	flux_corr_rho: tuple = None   # (feature names, correlation matrix of flux_corr)

	def phot_table(self):
		"""Per-band data vs model. K = band model / monochromatic model at lam_ref (colour correction)."""
		m = self.model
		model = m.photometry(self.p)
		mono = m.unconvolved(self.p, m.lam_phot)
		return Table({'ID': [self.summary['ID']] * len(m.filters), 'band': m.phot_names,
			'filter': [f.name if f else '' for f in m.filters],
			'reference': [f.reference if f else 'monochromatic' for f in m.filters],
			'lam_ref': m.phot[:, 0] * u.um, 'flux': m.phot[:, 1] * u.mJy, 'err': m.phot[:, 2] * u.mJy,
			'model': model * u.mJy, 'pull': (m.phot[:, 1] - model) / m.phot[:, 2],
			'model_mono': mono * u.mJy, 'K': model / mono})

	def rho_table(self):
		"""Long table of flux_corr correlations: one row per pair of fitted features."""
		names, R = self.flux_corr_rho
		i, j = np.triu_indices(len(names), 1)
		return Table({'ID': [self.summary['ID']] * i.size, 'feature_1': [names[k] for k in i],
			'feature_2': [names[k] for k in j], 'rho': R[i, j]})

	def param(self, name):
		m = self.model
		v = self.p[m.ix[name]]
		if not m.free[m.ix[name]]:
			return v, 0.
		k = np.cumsum(m.free)[m.ix[name]] - 1
		return v, np.sqrt(self.cov[k, k])

	def plot(self, path=None, show=False, title=None, xlim=None, model_range='full'):
		return plot_fit(self, path=path, show=show, title=title, xlim=xlim, model_range=model_range)


def _propagate(fun, q, cov, return_cov=False):
	"""Values and 1-sigma errors of dict-valued fun(q), linear propagation
	(with return_cov, also (keys, covariance matrix) of all outputs)."""
	f0 = fun(q)
	keys = list(f0)
	sd = np.sqrt(np.clip(np.diag(cov), 0, None))
	D = np.zeros((len(keys), q.size))
	for i in range(q.size):
		h = 1e-3 * sd[i] if sd[i] > 0 else 1e-8 * max(1., abs(q[i]))
		dq = np.zeros_like(q)
		dq[i] = h
		fp, fm = fun(q + dq), fun(q - dq)
		D[:, i] = [(fp[k] - fm[k]) / (2 * h) for k in keys]
	C = D @ cov @ D.T
	err = np.sqrt(np.clip(np.diag(C), 0, None))
	out = {k: (f0[k], e) for k, e in zip(keys, err)}
	return (out, (keys, C)) if return_cov else out


def residual_correlation_factor(z, segs, block=CORR_BLOCK):
	"""
	sqrt(1 + 2 sum_{k<block} (1 - k/block) rho_k): variance inflation of a local
	sum of normalized residuals z with lag-k autocorrelation rho_k, pooled over
	the index arrays `segs` (never across a seam or gap). This is the factor a
	residual bootstrap in blocks of `block` pixels picks up. Floored at 1.
	"""
	den = sum(np.sum(z[s]**2) for s in segs)
	f2 = 1.
	for k in range(1, block):
		num = sum(np.sum(z[s][:-k] * z[s][k:]) for s in segs if s.size > k)
		f2 += 2 * (1 - k / block) * num / den
	return np.sqrt(max(f2, 1.))


def _window(lam, x0, width, n_min):
	"""Pixels within +-width of x0, or the n_min nearest if that is fewer."""
	d = np.abs(lam - x0)
	sel = np.where(d <= width)[0]
	return sel if sel.size >= n_min else np.sort(np.argsort(d)[:n_min])


def local_chi2nu(lam, z, x0, width, n_min=LOCAL_MIN_PIX):
	"""Mean z^2 over pixels within +-width of x0 (at least the n_min nearest)."""
	return float(np.mean(z[_window(lam, x0, width, n_min)]**2))


def local_correlation_factor(lam, z, segs, x0, width, n_min=2 * LOCAL_MIN_PIX):
	"""residual_correlation_factor over +-width of x0, split at segment edges; NaN if too short."""
	sel = _window(lam, x0, width, n_min)
	sub = [np.intersect1d(s, sel) for s in segs]
	sub = [s for s in sub if s.size > CORR_BLOCK]
	return residual_correlation_factor(z, sub) if sub else np.nan


def fit_spectrum(spec, z, phot=None, error_scaling='local', max_nfev=20000, name='', **model_kw):
	"""
	Fit one spectrum. model_kw go to MIRModel (pahs, lines, continuum, poly_degree,
	profile, fit_ice, line_fwhm_kms, pah_shift_tol, line_shift_tol, line_insts).

	error_scaling: how covariance errors are inflated for the misfit that the
	stated pixel errors do not describe.
	  'local'  (default) each feature's errors x sqrt(max(1, local chi2_nu)) x
	           a residual-correlation factor. local chi2_nu is over +-1 observed
	           FWHM of the feature (>= LOCAL_MIN_PIX pixels); the correlation
	           factor over +-3 FWHM (>= 2 LOCAL_MIN_PIX pixels), falling back to
	           the instrument-wide value. See the module docstring.
	  'global' every error x sqrt(max(1, chi2_nu))
	  'none'   unscaled covariance
	Attenuation parameters (tau*) and FitResult.cov always use the global scaling
	(unless 'none'): they are set by the whole spectrum, not one window.
	"""
	if error_scaling not in ('local', 'global', 'none'):
		raise ValueError(error_scaling)
	m = MIRModel(spec, z, phot=phot, **model_kw)
	resid = lambda q: (m.model(m.full(q)) - m.y) / m.e
	jac = lambda q: m.jacobian(m.full(q)) / m.e[:, None]
	r = least_squares(resid, m.q_init, jac=jac, bounds=(m.lo, m.hi), x_scale='jac', max_nfev=max_nfev)
	chi2 = float(np.sum(r.fun**2))
	dof = max(m.y.size - r.x.size, 1)
	gscale = 1. if error_scaling == 'none' else max(1., chi2 / dof)
	cov_raw = np.linalg.pinv(r.jac.T @ r.jac)
	res = FitResult(m, m.full(r.x), cov_raw * gscale, chi2, dof, bool(r.success), r.message)
	der, (keys, C) = _propagate(lambda q: m.derived(m.full(q)), r.x, cov_raw, return_cov=True)
	# correlation between the corrected fluxes of all fitted features/complexes; unchanged by
	# the per-feature error scaling, so for a ratio a/b use
	# sigma_log = 0.4343 sqrt(ea^2 + eb^2 - 2 rho ea eb) with ea, eb fractional table errors
	fk = [i for i, k in enumerate(keys) if k.endswith(':flux_corr')]
	sd = np.sqrt(np.clip(np.diag(C)[fk], 0, None))
	with np.errstate(invalid='ignore', divide='ignore'):
		rho = C[np.ix_(fk, fk)] / np.outer(sd, sd)
	res.flux_corr_rho = ([keys[i].split(':')[0] for i in fk], rho)

	# per-feature error scale
	zs = r.fun[:spec.wave.size]
	fcorr = {}
	for inst in ('LRS', 'IRS'):
		segs = [s for s in m.segments if spec.inst[s[0]] == inst]
		fcorr[inst] = residual_correlation_factor(zs, segs) if segs else np.nan
	def window(n):
		"""(lo, hi) rest um: +-1 observed FWHM around a fitted feature."""
		x0 = der[n + ':lam'][0]
		w = np.hypot(PAH_FEATURES[n][1] if n in PAH_FEATURES else 0., m.lsf_fwhm(x0))
		return x0 - w, x0 + w

	scales = {}
	wins = {n: window(n) for n in m.pahs + m.lines}
	for c, mem in m.complexes.items():     # union of the members' windows
		wins[c] = (min(wins[n][0] for n in mem), max(wins[n][1] for n in mem))
	for n, (lo, hi) in wins.items():
		x0, width = (lo + hi) / 2, (hi - lo) / 2
		c2 = local_chi2nu(m.lam, zs, x0, width)
		fc = local_correlation_factor(m.lam, zs, m.segments, x0, 3 * width)
		if not np.isfinite(fc):
			fc = fcorr[spec.inst[np.argmin(np.abs(m.lam - x0))]]
		s = {'local': np.sqrt(max(1., c2)) * fc, 'global': np.sqrt(gscale), 'none': 1.}[error_scaling]
		scales[n] = {'err_scale': s, 'chi2nu_local': c2, 'corr_factor': fc}
		for k in ('lam', 'flux_obs', 'flux_corr', 'ew'):
			if n + ':' + k in der:
				v, e = der[n + ':' + k]
				der[n + ':' + k] = (v, e * s)
	for k in ('tau_305', 'tau_305_ice', 'tau_97'):
		der[k] = (der[k][0], der[k][1] * np.sqrt(gscale))

	res.table = _feature_table(m, der, name, scales)
	res.summary = {'cortado_version': __version__, 'ID': name, 'z': z, 'continuum': m.continuum, 'profile': m.profile_name,
		'chi2': chi2, 'dof': dof, 'chi2_nu': chi2 / dof, 'success': res.success,
		'n_pix': spec.wave.size, 'n_phot': m.phot.shape[0],
		'has_irs': bool(np.any(spec.inst == 'IRS')), 'error_scaling': error_scaling,
		'corr_factor_LRS': fcorr['LRS'], 'corr_factor_IRS': fcorr['IRS'],
		'partial_pahs': ','.join(m.partial)}
	for k in ('tau97', 'tau_ice305'):
		res.summary[k], res.summary[k + '_err'] = res.param(k)
	for k in ('tau_305', 'tau_305_ice', 'tau_97'):
		res.summary[k], res.summary[k + '_err'] = der[k]
	return res


def _feature_table(m, der, name, scales):
	dl = COSMO.luminosity_distance(m.z).to(u.m).value
	lum = 4 * np.pi * dl**2 / LSUN       # W/m2 -> Lsun
	rows = []
	for kind, req in m.requested.items():
		fitted = {'pah': m.pahs, 'line': m.lines, 'complex': m.complexes}[kind]
		for n in req:
			if kind == 'complex':
				lam0, fw = float(n[3:]), np.nan
				members = '+'.join(m.complexes.get(n, []))
			else:
				lam0, fw = PAH_FEATURES[n] if kind == 'pah' else (LINES[n], 0.)
				members = ''
			row = {'ID': name, 'feature': n, 'kind': kind, 'members': members, 'fitted': n in fitted,
				'lam_rest': lam0, 'fwhm_int': fw, 'fwhm_lsf': m.lsf_fwhm(lam0)}
			for k in ('lam', 'flux_obs', 'flux_corr', 'ew'):
				v, e = der.get(n + ':' + k, (np.nan, np.nan))
				row[k], row[k + '_err'] = v, e
			for k in ('obs', 'corr'):
				# NaN when the amplitude is pinned at its lower bound (use logL_corr_ul3)
				f, e = row['flux_' + k], row['flux_%s_err' % k]
				det = f > 1e-6 * e
				row['logL_' + k] = np.log10(f * lum) if det else np.nan
				row['logL_%s_err' % k] = 0.4343 * e / f if det else np.nan
			row['logL_corr_ul3'] = np.log10(3 * row['flux_corr_err'] * lum) if row['fitted'] else np.nan
			# detection significance from the observed flux; the corrected-flux error also
			# carries the uncertainty of the attenuation correction (1.5x larger at 3.3 um)
			row['snr'] = row['flux_obs'] / row['flux_obs_err'] if row['fitted'] else np.nan
			row['snr_corr'] = row['flux_corr'] / row['flux_corr_err'] if row['fitted'] else np.nan
			row.update(scales.get(n, {'err_scale': np.nan, 'chi2nu_local': np.nan, 'corr_factor': np.nan}))
			if kind == 'pah':
				row['coverage'], row['partial'] = m.coverage[n], n in m.partial
			elif kind == 'complex':
				mem = [x for x in PAH_COMPLEXES[n] if x in m.requested['pah']]
				row['coverage'] = min(m.coverage[x] for x in mem)
				row['partial'] = any(x in m.partial for x in mem)
			else:
				row['coverage'], row['partial'] = np.nan, False
			rows.append(row)
	t = Table(rows=rows)
	t['lam'].unit = t['lam_err'].unit = t['lam_rest'].unit = u.um
	t['fwhm_int'].unit = t['fwhm_lsf'].unit = t['ew'].unit = t['ew_err'].unit = u.um
	for c in ('flux_obs', 'flux_obs_err', 'flux_corr', 'flux_corr_err'):
		t[c].unit = u.W / u.m**2
	t.meta['cortado_version'] = __version__
	t.meta['comments'] = [
		'flux_obs: observed (attenuated) line flux; flux_corr: corrected for the mixed attenuation.',
		'logL_*: log10 L/Lsun (H0=70, Om=0.3, OL=0.7); logL_corr_ul3: log10 of 3-sigma L_corr.',
		'ew: rest-frame equivalent width vs the intrinsic continuum at the feature centre.',
		'fwhm_lsf: rest-frame instrumental FWHM at lam_rest.',
		'snr: flux_obs / flux_obs_err (detection significance); snr_corr: flux_corr / flux_corr_err, '
		'which also includes the attenuation-correction uncertainty.',
		'kind=complex: sum of the listed members; errors include their covariance.',
		'coverage: fraction of the feature profile flux on the data (complex: lowest member); '
		'partial: fit from a partly covered profile (partial_pahs=True), flux relies on the fixed '
		'profile shape and centre.',
		'*_err: covariance errors x err_scale; for error_scaling=local, err_scale = '
		'sqrt(max(1, chi2nu_local)) x corr_factor (residual autocorrelation of the instrument).',
	]
	return t


# ---------------------------------------------------------------- plotting
def plot_fit(res, path=None, show=False, title=None, xlim=None, model_range='full'):
	"""
	model_range: 'full' (default) draws the model continuously from the bluest data / filter
	edge to the reddest, solid where there is spectral data and dashed elsewhere
	(extrapolated, smoothed at the nearest data's R); 'data' draws it only over the
	spectral segments.
	"""
	if model_range not in ('data', 'full'):
		raise ValueError(model_range)
	m, p, s = res.model, res.p, res.model.spec
	fig, ax = plt.subplots(figsize=(7, 3.6))
	ax2 = ax.twinx()
	mod = m.model(p)
	n_spec = s.wave.size
	for seg in m.segments:
		lam = m.lam[seg]
		c = 'k' if s.inst[seg[0]] == 'LRS' else '0.35'
		ax.step(lam, s.flux[seg], where='mid', c=c, lw=1.2)
		ax.fill_between(lam, s.flux[seg] - s.err[seg], s.flux[seg] + s.err[seg],
			step='mid', facecolor=c, alpha=.25, lw=0, zorder=-10)
		if model_range == 'full':
			continue
		# model on a 4x denser pixel grid with the segment's own LSF
		n = seg.size
		t = np.linspace(0, n - 1, 4 * n)
		lp = np.interp(t, np.arange(n), lam)
		sp = np.interp(t, np.arange(n), m.sig[seg])
		comp = m.components(p, lp, sp)
		cont = comp.pop('continuum')
		for k, v in comp.items():
			ax.plot(lp, cont + v, c='m' if k in m.pahs else 'g', lw=.5, ls='--' if k in m.partial else '-')
		ax.plot(lp, cont, c='b', lw=.7)
		ax.plot(lp, cont + sum(comp.values()), c='r', lw=.9)
		ax2.plot(lp, mixed_attenuation(m.tau(p, lp)), 'k:', lw=.7)
	ymax_model = 0.
	if model_range == 'full':
		lp, sp, comp, cov = m.model_curve(p)
		cont = comp.pop('continuum')
		for k, v in comp.items():
			ax.plot(lp, cont + v, c='m' if k in m.pahs else 'g', lw=.5, ls='--' if k in m.partial else '-')
		ax.plot(lp, cont, c='b', lw=.7)
		tot = cont + sum(comp.values())
		# solid on the data, dashed elsewhere; widen each mask by one point so the pieces join
		on = cov | np.r_[cov[1:], False] | np.r_[False, cov[:-1]]
		ax.plot(lp, np.where(on, tot, np.nan), c='r', lw=.9)
		ax.plot(lp, np.where(~cov | ~on, tot, np.nan), c='r', lw=.9, ls='--')
		ax2.plot(lp, mixed_attenuation(m.tau(p, lp)), 'k:', lw=.7)
		ymax_model = np.nanmax(tot)
		if xlim is None:
			xlim = (lp.min() / 1.03, lp.max() * 1.03)
	if m.phot.size:
		# filters: point at the calibration reference wavelength, bar = half-max width (rest)
		hm = np.array([np.array(f.halfmax()) / (1 + m.z) if f else (l, l) for f, l in zip(m.filters, m.lam_phot)])
		xerr = np.abs(hm.T - m.lam_phot)
		ax.errorbar(m.lam_phot, m.phot[:, 1], xerr=xerr, yerr=m.phot[:, 2], fmt='o', c='k', ms=4, elinewidth=.8)
		ax.plot(m.lam_phot, mod[n_spec:], 'rx', ms=6)

	y = np.r_[s.flux, m.phot[:, 1]]
	ax.set_yscale('log')
	ax.set_ylim(np.nanmin(y[y > 0]) * 0.4, max(np.nanmax(y), ymax_model) * 1.6)
	ax.set_xscale('log')
	if xlim:
		ax.set_xlim(*xlim)
	ax2.set_ylim(0, 1.05)
	ax.set_xlabel(r'Rest-frame wavelength [$\mu$m]')
	ax.set_ylabel(r'$F_\nu$ [mJy]')
	ax2.set_ylabel(r'$(1-e^{-\tau})/\tau$')
	ax.set_title(title or '%s  z=%.3f  $\\chi^2_\\nu$=%.2f' % (res.summary['ID'], m.z,
		res.summary['chi2_nu']), fontsize=10)
	for a in (ax, ax2):
		a.minorticks_on()
		a.tick_params(which='both', direction='in')
	fig.tight_layout()
	if path:
		fig.savefig(path, dpi=200)
	if show:
		plt.show()
	else:
		plt.close(fig)
	return fig


# ---------------------------------------------------------------- LRSz2 sample
LRSZ2_EXCLUDE = ['GN_IRS16', 'GN_IRS62', 'FLS_IRS19456000']   # not in the final 37
LRSZ2_POLY = ['GN_IRS12', 'GN_IRS19', 'FLS_IRS509', 'FLS_IRS8226', 'FLS_IRS16122']  # AGN: poly continuum
LRSZ2_SEAM = 14.01      # obs um; LRS <= seam < IRS in specsForCafe files


def lrsz2_photometry(src):
	"""Spitzer MIPS 24, MIPS 70 (else Herschel PACS 70) and PACS 100 as (filter, mJy, mJy)."""
	phot = []
	if src['S24'] > 0:
		phot.append(('MIPS24', src['S24'], src['S24_err']))
	if src['MIPS70'] > 0:
		phot.append(('MIPS70', src['MIPS70'] / 1e3, src['MIPS70_err'] / 1e3))   # stored in uJy
	elif src['PACS70'] > 0:
		phot.append(('PACS70', src['PACS70'], src['PACS70_err']))
	if src['PACS100'] > 0:
		phot.append(('PACS100', src['PACS100'], src['PACS100_err']))
	return phot


def run_lrsz2(ids=None, outdir=OUTDIR, use_irs=True, use_phot=True, show=False, model_range='full', **fit_kw):
	"""Fit the LRSz2 sample; writes per-source PNGs and two ECSV tables to outdir."""
	import pandas as pd
	import lrsz2_stacking as LS        # sample corrections (L_IR, z) shared with the stacks
	sample = LS.apply_overrides(pd.read_pickle(os.path.join(HERE, 'lrsz2_pah33.pkl')))
	sample = sample[~sample['ID'].isin(LRSZ2_EXCLUDE)]
	if ids:
		sample = sample[sample['ID'].isin(ids)]
	os.makedirs(outdir, exist_ok=True)
	tabs, summ, rhos, phots = [], [], [], []
	for _, src in sample.iterrows():
		path = os.path.join(SPECDIR, src['ID'] + '_SpecforCafe.dat')
		if not os.path.isfile(path):
			print(src['ID'], 'no spectrum, skipped')
			continue
		spec = build_spectrum(path, seam=LRSZ2_SEAM)
		if not use_irs:
			spec = Spectrum(*(a[spec.inst == 'LRS'] for a in
				(spec.wave, spec.flux, spec.err, spec.inst, spec.R)))
		kw = dict(continuum='poly' if src['ID'] in LRSZ2_POLY else 'stellar_pl')
		kw.update(fit_kw)
		res = fit_spectrum(spec, src['zspec'], phot=lrsz2_photometry(src) if use_phot else None,
			name=src['ID'], **kw)
		res.plot(os.path.join(outdir, src['ID'] + '.png'), show=show, model_range=model_range)
		t = res.table
		print('%-16s z=%.3f chi2_nu=%6.2f  logL3.3=%s' % (src['ID'], src['zspec'],
			res.summary['chi2_nu'], np.round(t['logL_corr'][t['feature'] == 'PAH33'][0], 3)))
		tabs.append(t)
		summ.append(res.summary)
		rhos.append(res.rho_table())
		if res.model.phot.size:
			phots.append(res.phot_table())
	feats = vstack(tabs)
	summary = Table(rows=summ)
	feats.write(os.path.join(outdir, 'lrsz2_features.ecsv'), overwrite=True)
	summary.write(os.path.join(outdir, 'lrsz2_fit_summary.ecsv'), overwrite=True)
	vstack(rhos).write(os.path.join(outdir, 'lrsz2_flux_correlations.ecsv'), overwrite=True)
	if phots:
		vstack(phots).write(os.path.join(outdir, 'lrsz2_photometry.ecsv'), overwrite=True)
	return feats, summary


# ---------------------------------------------------------------- CLI
def _cli():
	ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('--version', action='version', version='cortado %s' % __version__)
	ap.add_argument('spectrum', nargs='?', help='LRS (or combined LRS+IRS) spectrum: text or x1d FITS')
	ap.add_argument('--z', type=float, help='redshift')
	ap.add_argument('--irs', help='separate IRS spectrum')
	ap.add_argument('--seam', type=float, help='split a combined file: lam_obs > seam is IRS')
	ap.add_argument('--flux-unit', default='mJy', help='flux unit of text inputs (mJy, Jy, uJy)')
	ap.add_argument('--phot', action='append', default=[], metavar='BAND:F:ERR',
		help='photometry, same flux unit as the spectra; BAND = filter alias (%s), an SVO filter '
		'ID (e.g. JWST/MIRI.F770W), or a number = observed wavelength in um, fit as a delta '
		'function (model value at exactly that wavelength, no filter); '
		'repeatable' % ', '.join(sorted(FILTER_ALIASES)))
	ap.add_argument('--continuum', choices=['stellar_pl', 'poly'], default='stellar_pl')
	ap.add_argument('--poly-degree', type=int, default=3)
	ap.add_argument('--profile', choices=list(PROFILES), default='lorentz')
	ap.add_argument('--no-ice', action='store_true', help='fix the 3.05 um ice optical depth to 0')
	ap.add_argument('--pahs', nargs='+', default=DEFAULT_PAHS, choices=list(PAH_FEATURES))
	ap.add_argument('--lines', nargs='*', default=DEFAULT_LINES, choices=list(LINES))
	ap.add_argument('--line-fwhm', type=float, default=300., help='intrinsic line FWHM [km/s]')
	ap.add_argument('--line-insts', nargs='+', default=['LRS'], choices=['LRS', 'IRS'],
		help='instruments whose spectra include atomic/H2 lines (default LRS only)')
	ap.add_argument('--partial-pahs', action='store_true',
		help='also fit PAHs only partly covered by the data (e.g. the blue side of 7.7 um)')
	ap.add_argument('--min-coverage', type=float, default=0.3,
		help='with --partial-pahs: minimum fraction of a feature profile on the data')
	ap.add_argument('--R-scale', type=float, default=1., help='scale all resolving powers')
	ap.add_argument('--error-scaling', choices=['local', 'global', 'none'], default='local',
		help='per-feature local chi2 x residual correlation (default), one global chi2, or none')
	ap.add_argument('--name', default=None, help='source name for outputs')
	ap.add_argument('--outdir', default=OUTDIR)
	ap.add_argument('--show', action='store_true')
	ap.add_argument('--model-range', choices=['full', 'data'], default='full',
		help='plot the model over the full data + photometry range (default) or only on the spectra')
	ap.add_argument('--lrsz2', nargs='*', metavar='ID', help='run on the LRSz2 sample (optionally these IDs)')
	ap.add_argument('--no-irs', action='store_true', help='LRSz2 mode: LRS only')
	a = ap.parse_args()

	model_kw = dict(continuum=a.continuum, poly_degree=a.poly_degree, profile=a.profile,
		fit_ice=not a.no_ice, pahs=a.pahs, lines=a.lines, line_fwhm_kms=a.line_fwhm,
		partial_pahs=a.partial_pahs, min_coverage=a.min_coverage, line_insts=a.line_insts)
	if a.lrsz2 is not None:
		model_kw.pop('continuum')   # per-source choice (LRSZ2_POLY)
		run_lrsz2(a.lrsz2 or None, outdir=a.outdir, use_irs=not a.no_irs, show=a.show,
			model_range=a.model_range, error_scaling=a.error_scaling, **model_kw)
		return
	if a.spectrum is None or a.z is None:
		ap.error('need a spectrum and --z (or --lrsz2)')
	scale = u.Unit(a.flux_unit).to(u.mJy)
	phot = []
	for x in a.phot:
		band, f, e = x.rsplit(':', 2)
		try:
			band = float(band)
		except ValueError:
			pass
		phot.append((band, float(f) * scale, float(e) * scale))
	spec = build_spectrum(a.spectrum, irs=a.irs, seam=a.seam, flux_unit=a.flux_unit, R_scale=a.R_scale)
	name = a.name or os.path.splitext(os.path.basename(a.spectrum))[0]
	res = fit_spectrum(spec, a.z, phot=phot, name=name, error_scaling=a.error_scaling, **model_kw)
	os.makedirs(a.outdir, exist_ok=True)
	res.table.write(os.path.join(a.outdir, name + '_features.ecsv'), overwrite=True)
	res.plot(os.path.join(a.outdir, name + '.png'), show=a.show, model_range=a.model_range)
	cols = ['feature', 'lam', 'fwhm_lsf', 'flux_corr', 'flux_corr_err', 'logL_corr', 'logL_corr_err', 'ew', 'snr', 'err_scale']
	t = res.table[res.table['fitted']][cols]
	for c in cols[1:]:
		t[c].info.format = '.4g'
	t.pprint(max_lines=-1, max_width=-1)
	print({k: v for k, v in res.summary.items()})


if __name__ == '__main__':
	_cli()
