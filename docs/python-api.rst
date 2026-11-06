Python API
==========

The Python interface accepts the same prepared inputs as the command-line
tool. The following example reads the bundled spectrum and returns one
result dictionary plus an iterator over the diagnostic pixels.

.. code-block:: python

   import csv

   from aspcap import Grid, fit_abundance

   grid = Grid.read("examples/toy-grid.dat")
   with open("examples/toy-spectrum.csv", newline="") as file:
       rows = list(csv.DictReader(file))

   wavelength, flux, error, continuum = [
       [float(row[name]) for row in rows]
       for name in ("wavelength", "flux", "error", "continuum")
   ]

   result, pixels = fit_abundance(
       grid, wavelength, flux, error, continuum,
       {"ALPHA": 0.0, "METALS": -0.25},
       star_id="synthetic-star",
       species="Mg",
       label="ALPHA",
       windows=[(15002, 15005)],
       add_labels=["METALS"],
   )

   print(result["abundance"])
   pixel_rows = list(pixels)

Pass ``executable="/path/to/ferre.x"`` if FERRE is outside your PATH.
Each pixel row contains wavelength, normalized flux, normalized uncertainty,
model flux, and weight.

Reference
---------

.. autoclass:: aspcap.Grid
   :members: read

.. autofunction:: aspcap.fit_abundance

Build this guide locally
------------------------

From the checkout, install the package and documentation requirements,
then build the HTML.

.. code-block:: sh

   python -m pip install -e .
   python -m pip install -r docs/requirements.txt
   python -m sphinx -W -b html docs docs/_build/html

Work directly with a grid
-------------------------

Construct a grid with its library path and, if needed, an executable path.
Use ``print(grid)`` to inspect the label order, ranges, and wavelength coverage.

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat", executable="/path/to/ferre.x")
   print(grid)
   spectrum = grid({"ALPHA": 0.20, "METALS": -0.25})
   batch = grid([(0.20, -0.25), (0.00, -0.25)])
   flux_at_two_pixels = grid.model([15002.5, 15004.5], 0.20, -0.25)

Calling the grid returns a dictionary with ``wave`` and ``flux`` tuples. For
a batch, ``flux`` contains one tuple per model. The optional ``wave`` argument
or constructor setting ``outwave`` linearly resamples the native model within
the grid's wavelength coverage. Flux stays on the library's scale.

Use ``grid.fit(flux, error)`` to fit a prepared spectrum on the native grid.
Supply per-pixel one-sigma errors. The fit starts at the grid midpoint and
varies all labels by default. To provide your own starting point or freeze
some labels, pass ``start`` and ``free_labels``. An optional ``weights`` vector
can give individual pixels zero weight. The result has ``parameters``,
``errors``, ``reduced_chi2``, and ``model`` attributes.

.. automethod:: aspcap.Grid.fit
