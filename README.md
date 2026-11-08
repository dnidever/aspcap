# aspcap

APOGEE ASPCAP abundance software. This first workflow fits one element in one
prepared spectrum using the stellar parameters and continuum from a global fit.

## Try the worked example

Install Python 3.9 or later and [FERRE](https://github.com/callendeprieto/ferre).
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

1. Copy [one-star.json](examples/one-star.json) and point it to your inputs.
   Relative paths are resolved from the configuration file. Supply the global
   fit as a dictionary with one starting value per grid label.
2. Supply a spectrum CSV with `wavelength,flux,error,continuum` columns. The
   spectrum must already be on the native rest-frame grid in the library's
   wavelength units and air/vacuum convention. Give flux and one-sigma error
   in the same units as the positive global continuum. The wrapper divides
   both by that continuum and holds it fixed during the fit.
3. Copy [mg.json](examples/mg.json) and set the species, the one free grid label,
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

The [user guide](docs/index.rst) includes a quickstart, input reference, and
Python API. The repository includes Read the Docs build configuration. To
preview the guide locally after installing ASPCAP:

```sh
python -m pip install -r docs/requirements.txt
python -m sphinx -W -b html docs docs/_build/html
```


The full-wavelength APOGEE DR17 M-giant test grid and detector weights are
included when you install ASPCAP from this branch.

```sh
python -m pip install 'git+https://github.com/dnidever/aspcap.git@hge/first-abundance-fit'
```

See the [Grid walkthrough](docs/grid.rst#fit-the-full-apogee-dr17-wavelength-grid)
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

Polynomial, spline, median, binned-percentile and model-ratio methods operate
independently per order. All normalization here is outside FERRE. Applying a
continuum rescales both flux and error; it does not include uncertainty in the
continuum estimate. See [the spectrum guide](docs/spectrum.rst),
[the continuum guide](docs/continuum.rst), and the runnable
[example](examples/spectrum_continuum.py). The core design is adapted from
[Doppler Spec1D](https://github.com/dnidever/doppler/blob/6e165f87c6d766eff6eaa1fbee2c82be5bb30cc0/doppler/spec1d.py);
its MIT notice is retained in `LICENSE-DOPPLER`.
