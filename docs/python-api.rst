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

Pass ``executable="/path/to/ferre"`` if FERRE is outside your PATH.
Each pixel row contains wavelength, normalized flux, normalized uncertainty,
model flux, and weight.

Reference
---------

See :doc:`grid` for model evaluation, fitting, and the Grid API reference.

.. autofunction:: aspcap.fit_abundance

Build this guide locally
------------------------

From the checkout, install the package and documentation requirements,
then build the HTML.

.. code-block:: sh

   python -m pip install -e .
   python -m pip install -r docs/requirements.txt
   python -m sphinx -W -b html docs docs/_build/html
