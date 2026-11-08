ASPCAP for HGE
==============

ASPCAP provides spectral data objects, external continuum normalization, and
FERRE-based element-abundance fitting for HGE.

Start with :doc:`spectrum` for Spec1D objects, multiple orders, masks, and
APOGEE chip mapping. See :doc:`continuum` for polynomial, spline, median,
percentile, and model-ratio estimators, their options, and worked examples.
These operations are outside FERRE and require no executable or Doppler models.

The abundance workflow runs FERRE on a prepared spectrum using parameters and
a fixed continuum from a global fit. The element configuration controls which
wavelengths contribute and which synthetic-grid label FERRE can vary.
See :doc:`quickstart` for that workflow.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   quickstart
   inputs
   grid
   spectrum
   continuum
   python-api
