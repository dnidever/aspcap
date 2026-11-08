# aspcap

ASPCAP provides multi-order spectral data objects, external continuum
normalization, and element-abundance fitting with FERRE for HGE.

- [Spec1D objects](https://github.com/dnidever/aspcap/blob/main/docs/spectrum.rst): single or multiple orders, unequal chip
  lengths, masks, uncertainties, copying, and native-grid mapping.
- [External continuum normalization](https://github.com/dnidever/aspcap/blob/main/docs/continuum.rst): polynomial, spline,
  median, percentile, and model-ratio methods with diagnostics.
- [FERRE abundance fitting](https://github.com/dnidever/aspcap/blob/main/docs/quickstart.rst): fit an element using fixed
  global parameters and a supplied continuum.

## The ASPCAP workflow

ASPCAP (the APOGEE Stellar Parameters and Chemical Abundances Pipeline) compares observed spectra with synthetic spectral grids using FERRE. The workflow combines continuum normalization with two fitting stages: a global stellar-parameter fit and individual-element abundance fits.

1. **Prepare and normalize the spectrum.** Keep wavelengths, flux, uncertainties, and bad-pixel masks together, preserving separate orders or detector chips where appropriate. Estimate the continuum and divide both flux and uncertainty by it. Observations and models must use compatible normalization and wavelength conventions.

2. **Fit global stellar parameters.** Fit the usable spectrum jointly across a multidimensional synthetic grid. The planned HGE global fit has six free labels: effective temperature, surface gravity (log g), overall metallicity [M/H], alpha enhancement [alpha/M], carbon [C/M], and nitrogen [N/M]. Microturbulence is prescribed by a relation with log g rather than fitted as an independent parameter. These labels describe effects across the spectrum, including molecular absorption. The relation and its coefficients must be specified consistently with the synthetic grid; this overview does not prescribe a particular calibration.

3. **Fit individual elemental abundances.** Start from the global solution and hold all grid dimensions fixed except the abundance dimension, which is varied to fit the target element's abundance. Fit selected element-sensitive wavelength windows, with pixel weights emphasizing useful features and excluding unsuitable pixels. Alpha elements can use [alpha/M], carbon and nitrogen their respective coordinates, and other elements [M/H]. The coordinate is a proxy for the target element within its windows: changing an alpha or metallicity coordinate changes multiple elements in the synthetic grid. Convert the fitted coordinate to [X/H] according to the grid's abundance convention and the fixed global labels.

   We are also developing an alternative called mini-grids. Each mini-grid adds a dedicated seventh dimension to the six global stellar labels, varying the abundance of a single target element. During that element's abundance fit, the six global labels are held at their global-fit values and only the dedicated element-abundance dimension is varied. This provides an element-specific fitting coordinate instead of using the shared alpha or metallicity coordinate as a proxy.

See the [SDSS ASPCAP description](https://www.sdss4.org/dr17/irspec/aspcap/) for the survey workflow and its grid-dependent choices.

### What this package currently provides

This HGE package provides Spec1D objects, external continuum estimators, a Grid interface for FERRE interpolation and fits with multiple free labels, and fit_abundance plus a command-line interface for individual-element fits. A global fit can be performed with Grid.fit using a suitable supplied grid. An automatic end-to-end pipeline that selects grids, normalizes spectra, runs the global fit, and loops over elements is not yet provided.

The current fit_abundance interface accepts supplied global parameters and a fixed continuum. It uses unit weights inside the selected windows and zero outside, excludes bad pixels, and uses measurement uncertainties in the fit. Arbitrary element-specific weight curves are not yet exposed by this interface; Grid.fit accepts explicit pixel weights. External normalization is separate from FERRE's internal continuum controls.

Install the latest release from PyPI with:

```sh
python -m pip install --upgrade aspcap
```

The spectrum and continuum modules do not need a FERRE executable or Doppler
model files. Install [FERRE](https://github.com/dnidever/ferre) separately when
running spectral fits: `python -m pip install --upgrade ferre`.

## Try the worked example

Install Python 3.9 or later and [FERRE](https://github.com/dnidever/ferre).
The example was verified with FERRE 5.1.3. From the repository checkout, run:

```sh
python -m pip install -e .
aspcap abundances examples/one-star.json --output abundance.csv --pixels pixels.csv
```

If `ferre` is outside your PATH, add `--ferre /path/to/ferre`. The example
uses a synthetic spectrum and a toy grid. Expect a raw `ALPHA` coordinate of
about 0.20 and a reported `[Mg/H]` of about −0.05, using four pixels and fixed
`METALS = -0.25`. The pixel CSV shows the normalized flux, uncertainty, model,
and whether each pixel was used. Add `--keep-runs` to inspect FERRE's input
files after a successful fit. A failed fit keeps its files and reports their
location.

## Run your own fit

1. Copy [one-star.json](https://github.com/dnidever/aspcap/blob/main/examples/one-star.json) and point it to your inputs.
   Relative paths are resolved from the configuration file. Supply the global
   fit as a dictionary with one starting value per grid label.
2. Supply a spectrum CSV with `wavelength,flux,error,continuum` columns. The
   spectrum must already be on the native rest-frame grid in the library's
   wavelength units and air/vacuum convention. Give flux and one-sigma error
   in the same units as the positive global continuum. The wrapper divides
   both by that continuum and holds it fixed during the fit.
3. Copy [mg.json](https://github.com/dnidever/aspcap/blob/main/examples/mg.json) and set the species, the one free grid label,
   and one or more inclusive wavelength windows. Overlapping windows still
   give each pixel weight one. Nonfinite flux and nonpositive or nonfinite
   uncertainties get zero weight. Use `add_labels` to list any fixed logarithmic
   coordinates to add to the fitted coordinate to report `[X/H]`. Leave it empty
   if the fitted coordinate is already `[X/H]`. Verify this mapping and the
   wavelength windows against the scientific conventions of your grid.

This version supports single-segment, non-PCA FERRE libraries with a positive
spectral resolution, either formatted ASCII or matching `.hdr` and `.unf`
files. Point `grid` to the ASCII library or to the binary library's `.hdr`.
The fit uses linear interpolation. The reported error is conditional on the
fixed global parameters and continuum. The reduced χ² column comes from FERRE.

## Use from Python

```python
from aspcap import Grid, fit_abundance

grid = Grid.read("examples/toy-grid.dat")
result, pixels = fit_abundance(
    grid, wavelength, flux, error, continuum,
    {"ALPHA": 0.0, "METALS": -0.25},
    star_id="my-star", species="Mg", label="ALPHA",
    windows=[(15002, 15005)], add_labels=["METALS"],
)
```

Pass the four CSV columns as sequences. `result` is a dictionary and `pixels`
iterates over the diagnostic rows.

The FERRE protocol lives in `aspcap/ferre.py`, the abundance logic in
`aspcap/abundances.py`, and file input and output in `aspcap/cli.py`. Run
`python -m unittest discover -s tests -v` for the tests. Set
`ASPCAP_TEST_FERRE=/path/to/ferre` to include the integration test if FERRE
is outside your PATH.

## Documentation

The [user guide](https://github.com/dnidever/aspcap/blob/main/docs/index.rst) includes a quickstart, input reference,
[Spec1D guide](https://github.com/dnidever/aspcap/blob/main/docs/spectrum.rst), [continuum methods](https://github.com/dnidever/aspcap/blob/main/docs/continuum.rst), and
Python API. The repository includes Read the Docs build configuration. To
preview the guide locally after installing ASPCAP:

```sh
python -m pip install -r docs/requirements.txt
python -m sphinx -W -b html docs docs/_build/html
```


The full-wavelength APOGEE DR17 M-giant test grid and detector weights are
included in ASPCAP installations.

```sh
python -m pip install --upgrade aspcap
```

See the [Grid walkthrough](https://github.com/dnidever/aspcap/blob/main/docs/grid.rst#fit-the-full-apogee-dr17-wavelength-grid)
for an example using the installed files.


## Multi-order spectra and external continuum fitting

`Spec1D` follows Doppler's `(pixels, orders)` layout, accepts unequal-length
orders, preserves masks/metadata, and can map detector chips back to a native
wavelength grid. It imports no Doppler fitting code or model loaders. Existing
Doppler objects can be used directly with the continuum functions.

```python
from aspcap import Spec1D, continuum

spec = Spec1D(flux_orders, wave=wavelength_orders, err=error_orders)
result = continuum.fit(spec, method="chebyshev", degree=4)
normalized = result.apply(spec)  # independent copy
```

Choose the estimator with `continuum.fit(spec, method=...)`:

| Method | Useful options | Defaults |
|---|---|---|
| `chebyshev` / `polynomial` | `degree`, `weighted`, `clip` | Degree 4, unweighted, no clipping |
| `spline` | `degree`, `smoothing`, `weighted`, `clip` | Degree 3, noise-based smoothing |
| `median` | `window`, `mode` | 151 pixels, mirror edges |
| `percentile` | `percentile`, `bins` | 90th percentile, 20 bins |
| `model_ratio` | `model`, `window`, `ratio_bounds`, `mode` | Model required; 151 pixels; ratios 0.1–10 |

Methods operate independently per order. Use `result.fit_mask`, `result.flags`,
`result.parameters`, and `result.settings` to inspect the estimate. Invalid
orders raise an error by default; `on_failure="mask"` returns flagged, masked
results instead. `result.apply(spec)` returns a copy; `inplace=True` updates the
input explicitly. All normalization here is outside FERRE. Applying a
continuum rescales both flux and error; it does not include uncertainty in the
continuum estimate. See [the spectrum guide](https://github.com/dnidever/aspcap/blob/main/docs/spectrum.rst),
[the continuum guide](https://github.com/dnidever/aspcap/blob/main/docs/continuum.rst), and the runnable
[example](https://github.com/dnidever/aspcap/blob/main/examples/spectrum_continuum.py). The core design is adapted from
[Doppler Spec1D](https://github.com/dnidever/doppler/blob/6e165f87c6d766eff6eaa1fbee2c82be5bb30cc0/doppler/spec1d.py);
its MIT notice is retained in `LICENSE-DOPPLER`.
