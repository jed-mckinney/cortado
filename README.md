# cortado

cortado fits PAH features, atomic and H2 emission lines, a stellar + dust
continuum and dust/ice attenuation to JWST/MIRI LRS spectra, optionally joined
with a Spitzer/IRS spectrum and broad-band photometry. It returns integrated
fluxes, luminosities, equivalent widths and uncertainties for every feature,
plus a figure of the fit. It was developed for, and is described in,
McKinney et al. 2025, "A JWST MIRI LRS Survey of 37 Massive Star-Forming
Galaxies and AGN at Cosmic Noon -- Overview and First Results", accepted to
ApJ ([arXiv:2510.07365](https://arxiv.org/abs/2510.07365)).

It grew out of the PAH fitter for the LRSz2 program, nicknamed "CAFE at home",
hence the name. It was shared earlier as the "LRSz2 PAH fitter"; apart from
the name and the version number, the code is the same.

Version: 2026.10.20 (`python measure_pahs_general.py --version`; changes since
the 2026-10-15 package are listed at the end)
Contact: Jed McKinney (jed.mckinney.astro@gmail.com)
License: MIT (see `LICENSE`)


## Contents

| Path | What it is |
|---|---|
| `measure_pahs_general.py` | cortado itself: command-line tool and Python module |
| `run_pah_fits.sh` | Runs cortado on every source in a CSV table |
| `README.md` | This file |
| `requirements.txt` | Python packages |
| `sb99-10myr.txt` | Starburst99 10 Myr stellar spectrum (continuum model) |
| `templates/jwst_miri_p750l_disp_*.fits` | MIRI LRS (P750L) resolving power vs wavelength |
| `filters/svo/*.xml` | Filter transmission curves from the SVO Filter Profile Service: Spitzer/MIPS 24/70/160, Spitzer/IRAC 5.8/8.0, Spitzer/IRS blue peak-up, Herschel/PACS 70/100/160, JWST/MIRI F560W–F2550W, JWST/NIRCam F444W |
| `example/` | Example table and two spectra from McKinney et al. 2025 (see Example data) |
| `LICENSE` | MIT license |

Keep this folder structure: cortado finds its data files relative to
`measure_pahs_general.py`. Every output table records the cortado version
that produced it (`cortado_version`); please quote it when reporting a
problem.


## Install

Python 3 with numpy, scipy, astropy and matplotlib:

```
pip install -r requirements.txt
```

The bash script needs bash 3.2 or later (the macOS default is fine).
Filters not in `filters/svo/` are downloaded from SVO the first time they are
used, which needs an internet connection.


## Quick start

From this folder, fit the two example galaxies:

```
./run_pah_fits.sh example/pah_fit_table_example.csv
```

This writes one figure (`<id>.png`) and one feature table
(`<id>_features.ecsv`) per galaxy to `./pah_fits_table/`, and logs to
`./pah_fits_table/logs/`. Each fit takes about 1–10 s. Add `-n` to print the
commands without running them.

One spectrum directly:

```
python measure_pahs_general.py example/spectra/GN_IRS12_SpecforCafe.dat --z 2.0091 \
    --seam 14.01 --continuum poly --phot MIPS24:1.224:0.00682 --phot PACS100:15.19:0.54 \
    --outdir my_fits
```


## Inputs

**Spectra.** Either a text file with columns: observed wavelength [µm],
F_ν, error, and optionally the resolving power R of each pixel; or a JWST
x1d FITS file (EXTRACT1D extension; flux units read from the header). Text
fluxes are mJy unless `--flux-unit Jy` or `--flux-unit uJy` is given.

- LRS only: give the spectrum.
- LRS + IRS in two files: `--irs irs.txt`.
- LRS + IRS in one file: `--seam 14.01` (pixels with λ_obs > seam are IRS).
- LRS pixels below 4.77 µm observed are dropped (unreliable wavelength
  calibration).

**Redshift:** `--z`.

**Photometry:** `--phot BAND:flux:error`, repeatable, in the same flux unit
as the spectra. BAND is either a filter or a wavelength, and the two are
treated differently:

- **A filter**, given as an alias (`MIPS24`, `MIPS70`, `MIPS160`,
  `PACS70`, `PACS100`, `PACS160`, `IRAC1`–`IRAC4`, `F560W`, `F770W`,
  `F1000W`, `F1130W`, `F1280W`, `F1500W`, `F1800W`, `F2100W`, `F2550W`;
  case-insensitive) or any SVO filter ID (e.g. `JWST/MIRI.F1800W`,
  downloaded on first use). The model (continuum, PAHs and lines) is
  integrated through the filter's transmission curve in the facility's
  calibration convention, so give catalogue flux densities as published
  (not colour-corrected).
- **A number**, e.g. `--phot 24:0.45:0.02`: the observed wavelength in µm of
  a delta-function measurement. It is never matched to a filter, even when it
  equals a filter's nominal wavelength. The model value used is the
  attenuated continuum + PAHs at exactly λ_rest = λ_obs / (1+z), at full
  intrinsic resolution: no filter, no instrumental smoothing, no colour
  correction. Emission lines do not contribute to such a point (at their
  ~300 km/s intrinsic width they would matter only within a hair of the line
  centre). Use this for a flux density that is already a monochromatic value;
  use the filter name for broad-band photometry.

So `--phot 24:…` and `--phot MIPS24:…` are different measurements: the
first is the model at exactly 24 µm, the second the model integrated over the
MIPS 24 µm band.


## Batch runs: `run_pah_fits.sh`

```
./run_pah_fits.sh [options] TABLE [-- options passed to every fit]
```

TABLE is a comma-separated file with one header row:

| Column | Meaning |
|---|---|
| `id` | source name, used for output file names (required) |
| `spec` | spectrum file (required) |
| `z` | redshift (required) |
| `irs` | separate IRS spectrum (optional) |
| `seam` | split a combined file at this observed wavelength (optional) |
| `args` | extra fitter options for this source only, e.g. `--continuum poly` (optional) |
| `BAND`, `BAND_err` | one pair of columns per photometric band, BAND as above |

Relative paths are taken from the table's folder. A band is skipped for a
source when its value or error is empty or `nan`, or the error is not
positive. Fields are split on commas without quoting, so values cannot
contain commas. Script options: `-n` dry run, `-o DIR` output folder,
`-m full|data` model range in the figures, `-p` Python interpreter. Run
`./run_pah_fits.sh --help` for details. The exit status is the number of
sources that failed.

Example, with partly covered PAH features allowed and fluxes in Jy:

```
./run_pah_fits.sh -o fits sources.csv -- --partial-pahs --flux-unit Jy
```


## Main options

| Option | Default | Meaning |
|---|---|---|
| `--continuum stellar_pl\|poly` | `stellar_pl` | stellar (SB99) + power law, or a cubic polynomial (for AGN-dominated spectra) |
| `--profile lorentz\|drude` | `lorentz` | PAH profile shape |
| `--no-ice` | off | fix the 3.05 µm ice absorption to zero |
| `--pahs ...` | 17 features | which PAH features to fit |
| `--lines ...` | 22 lines | which emission lines to fit (`--lines` with no names fits none) |
| `--line-insts LRS [IRS]` | `LRS` | which spectra contain emission lines (see below) |
| `--partial-pahs` | off | also fit PAH features only partly covered by the data (see below) |
| `--min-coverage` | 0.3 | with `--partial-pahs`: minimum fraction of a feature on the data |
| `--error-scaling local\|global\|none` | `local` | how uncertainties are scaled for fit quality (see below) |
| `--model-range full\|data` | `full` | model drawn over the full data + photometry range, or only on the spectra |
| `--R-scale` | 1 | multiply all resolving powers (for systematics tests) |
| `--outdir`, `--name` | | output folder and source name |

`python measure_pahs_general.py --help` lists everything. The `--lrsz2`
option only works inside the LRSz2 project and is not supported here.


## Outputs

**`<id>_features.ecsv`**, one row per requested feature (PAHs, lines, and
the 7.7, 11.3 and 12.7 µm complexes), including:

| Column | Meaning |
|---|---|
| `fitted` | False if the spectrum does not cover the feature |
| `flux_obs`, `flux_obs_err` | observed (attenuated) integrated flux [W m⁻²] |
| `flux_corr`, `flux_corr_err` | flux corrected for the fitted attenuation [W m⁻²] |
| `logL_corr`, `logL_corr_err` | log₁₀ of the corrected luminosity [L☉] (H₀ = 70, Ω_m = 0.3, Ω_Λ = 0.7) |
| `logL_corr_ul3` | log₁₀ of the 3σ luminosity upper limit |
| `ew`, `ew_err` | rest-frame equivalent width [µm] against the fitted continuum |
| `snr` | detection S/N, from the observed flux: a feature is detected if > 3 |
| `lam`, `fwhm_lsf` | fitted centre and instrumental FWHM [rest µm] |
| `coverage`, `partial` | fraction of the feature on the data; fit from a partial profile |
| `members` | for complexes, the components summed |

Use `snr`, not flux/error, to decide detection: the corrected-flux error also
includes the uncertainty of the attenuation correction. Use the complex rows
for 7.7, 11.3 and 12.7 µm totals rather than adding components yourself;
their errors include the strong anticorrelation between blended components.

**`<id>.png`:** data (LRS black, IRS grey) and photometry (at the filter's
calibration wavelength, with its half-maximum width). The model is red: solid
on the spectra, dashed in gaps and beyond (extrapolated). Red crosses are the
model integrated through each filter. Also shown: the continuum (blue), each
PAH (magenta, dashed if partial), each line (green) and the attenuation
factor (dotted, right axis). Photometry given as a wavelength is drawn at that
wavelength, without a filter-width bar.

**From Python** you get the same plus the fit itself:

```python
from measure_pahs_general import build_spectrum, fit_spectrum
spec = build_spectrum('spec.txt', seam=14.01)
res = fit_spectrum(spec, z=1.44, phot=[('MIPS24', 0.45, 0.02),   # filter-integrated
                                     (100.0, 12.0, 2.0)])       # delta function at 100 µm obs
res.table          # feature table
res.summary        # chi2, attenuation (tau_9.7, ice)
res.phot_table()   # photometry: data, model, colour correction per band
res.rho_table()    # correlations between feature fluxes, for ratio errors
res.plot('fit.png')
```


## What the model does

```
F_ν(λ) = A(λ) [ continuum + Σ PAH_k + Σ line_j ],   A = (1 − e^−τ)/τ
τ = τ_9.7 × silicate profile + τ_ice × 3.05 µm ice profile
```

- **Emission lines on LRS only (default):** atomic and H2 lines are fit to
  the MIRI LRS pixels only. The Spitzer/IRS spectra are treated as sensitive
  only to the continuum and PAHs, so the model has no line emission at IRS
  pixels, and a line that falls only in the IRS spectrum is not fit
  (`fitted` = False in the output). `--line-insts LRS IRS` fits lines in
  both spectra; `--line-insts IRS` in IRS only.
- **Instrumental resolution:** the model is smoothed pixel by pixel with each
  instrument's resolving power before it is compared with the data (MIRI LRS
  from the P750L curve; IRS from the Instrument Handbook values for each
  order). Smoothing conserves flux, so fluxes come from the intrinsic profiles.
- **Photometry:** the model is integrated through each filter curve, using
  the facility's calibration convention: a 10⁴ K blackbody for MIPS,
  νF_ν = constant for IRAC and PACS, flat F_ν for JWST.
- **Uncertainties:** from the fit covariance, scaled for each feature by the
  reduced χ² within one feature width and by the correlation of neighbouring
  residuals. In tests these are ~15–20 % larger than a residual bootstrap.
- **Partly covered features:** with `--partial-pahs`, features with at least
  `--min-coverage` of their profile on the data are fit with fixed centres.
  For a spectrum stopping at ~7.9 µm rest, report the 7.7 µm complex total,
  not its components, and allow a ~10 % systematic (compare
  `--min-coverage 0.3` and `0.1`).

The docstring at the top of `measure_pahs_general.py` gives the full
details and the validation tests.


## Caveats

- The LRS resolving-power curve is the pre-launch pandeia table (Kendrew et
  al. 2015). Changing R by 20–30 % moves PAH luminosities by ≤ 0.01 dex,
  except the 11.23 µm component on its own (~0.05 dex).
- For IRS, the order of each pixel is assumed from its wavelength (switching
  at 7.5, 14.0 and 19.5 µm). If your IRS data say otherwise, give R in a
  fourth column of the spectrum file.
- The 7.42 µm component of the 7.7 µm complex is not in the default PAH
  list, so the 7.7 µm complex total is 7.60 + 7.85 µm only. Add it with
  `--pahs` if you need the full PAHFIT-style complex.
- Some spectra are not described well by the continuum models (reduced χ²
  well above 10); check the figure before using the numbers.
- Extrapolated (dashed) parts of the model are not constrained by the
  spectra.
- Filters without a known calibration convention (here, the Spitzer/IRS blue
  peak-up filter, or any other SVO filter outside JWST, MIPS, IRAC and PACS)
  are treated as calibrated to a flat F_ν, with a warning. Check the
  facility's convention if the colour correction matters.


## Example data

`example/spectra/` holds two galaxies from the LRSz2 survey, published in
McKinney et al. 2025 (arXiv:2510.07365): GN_IRS12 (z = 2.009, AGN-dominated,
fit with `--continuum poly`) and FLS_IRS8040 (z = 0.76). Each file combines
the JWST/MIRI LRS spectrum (observed wavelengths ≤ 14.01 µm) with the
Spitzer/IRS spectrum above it; IRS pixels with S/N below ~1 were removed.
The MIRI LRS data are archived at MAST:
[doi:10.17909/18v0-8068](https://doi.org/10.17909/18v0-8068); the IRS spectra
are archival Spitzer data.
Columns: observed wavelength [µm], F_ν [mJy], uncertainty [mJy]; the header
lines (starting with `#`) give the source and citation. If you use these
spectra, please cite McKinney et al. 2025 and the MAST data DOI.

FLS_IRS8040 is one of the galaxies the continuum model fits poorly (reduced
χ² ≈ 76), so it also shows what such a fit looks like.


## Citing cortado

If cortado is useful for your work, please cite McKinney et al. 2025 (ApJ,
accepted; arXiv:2510.07365) and link to this repository. The paper's MIRI LRS spectra
are available from MAST at
[doi:10.17909/18v0-8068](https://doi.org/10.17909/18v0-8068).


## Credits

cortado bundles or relies on:

- filter transmission curves from the SVO Filter Profile Service
  (Rodrigo, Solano & Bayo 2012; Rodrigo & Solano 2020);
- the Starburst99 stellar population models (Leitherer et al. 1999);
- the MIRI LRS resolving power from the pandeia reference data
  (Pontoppidan et al. 2016), based on Kendrew et al. (2015);
- the Kemper, Vriend & Tielens (2004) silicate extinction profile as
  implemented in CAFE, and the PAH band parameters of Smith et al. (2007).


## Changes since the 2026-10-15 package

- Public release on GitHub under the MIT license. The examples are now two
  published spectra (GN_IRS1 removed), with citation headers.

- Renamed cortado; `--version` added, and the version is written to every
  output table and fit summary. File names are unchanged.
- Emission lines are fit to the LRS spectrum only by default (`--line-insts`).
- The fit figure no longer has a residual panel.
- Photometry input is more forgiving: aliases are case-insensitive; a
  wavelength may be given as text (`'24'`); photometry may be passed to
  `fit_spectrum` as an (n, 3) array; an unknown band name now gives an error
  listing the valid aliases.
- The documentation states explicitly how numeric photometry is treated (a
  delta function at that wavelength, never a filter). Its behaviour has not
  changed.
- Four more cached filter curves: IRAC 5.8 and 8.0 µm, IRS blue peak-up,
  NIRCam F444W.
