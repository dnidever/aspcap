Observed Spectra
================

The supported input product for this ASPCAP workflow is **apStar**. You can fit either the **combined spectrum** or the
**individual visit spectra stored in the same file**. Both are already on the
common, stellar rest-frame wavelength grid.

Choosing a spectrum
-------------------

An apStar file contains multiple spectral rows. Select either a combined
spectrum or one of the individual-visit rows, along with its corresponding
uncertainties and pixel masks. The product can contain more than one combined
spectrum with different weighting schemes; use the release's data model to
identify the intended row rather than assuming every row is a visit.

The combined spectrum uses multiple visits to improve signal-to-noise.
Individual rows let you analyze visits separately while using the same
wavelength sampling. They remain individual measurements even though they
have already been resampled. Keep each visit's identity and relevant
instrumental information with its fit.

Preparing an apStar row for fitting
-----------------------------------

1. Select the combined or individual-visit row and its matching error and
   mask arrays. Retain the star and visit identifiers and inspect quality flags.
2. Confirm that the chosen synthetic library has compatible wavelength units,
   vacuum/air convention, coverage, and :doc:`lsf`. A combined spectrum can
   have an effective instrumental profile different from a particular visit.
3. Match the selected row to the library's native pixels. A library can contain
   selected detector segments rather than the full apStar array. Do not
   interpolate merely because a different row was selected: all these rows
   already share the apStar wavelength scale.
4. Exclude bad pixels and detector gaps, then apply :doc:`continuum`
   consistently with the synthetic library. An apStar rest-frame spectrum is
   not necessarily continuum normalized for ASPCAP fitting.
5. Supply :doc:`initial-estimates` and run :doc:`global-fitting`, followed
   by :doc:`element-fitting`, for the selected spectrum.

For the bundled full-range APOGEE fixture, the library uses the apStar grid
with zero-weight detector gaps. For other libraries, use their actual
wavelength arrays and modeled-pixel selection. The current fitting interfaces
require native library pixels in the correct order. If a library uses a genuinely
different sampling, an explicit preparation step is still needed; see
:doc:`grid` and :doc:`ferre`.

Context from earlier implementations
------------------------------------

DR17 ASPCAP's inspected IDL reader can select the combined spectrum and
individual-visit rows from an apStar file. The driver prepares the selected
rows for each synthetic-library segment. This supports the same useful
distinction here: choose which spectrum to fit without repeating the visit
combination pipeline's velocity correction and common-grid resampling.

The inspected Astra ASPCAP entry point targets coadded spectra in apStar,
and its FERRE preprocessing selects the modeled detector pixels from the
prepared arrays. That inspected main path does not establish an individual-visit
analysis interface.

Sources: `DR17 ASPCAP reader
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/aspcap_rdobsnew.pro>`_,
`DR17 ASPCAP driver
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro>`_,
`Astra ASPCAP entry point
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/__init__.py>`_, and
`Astra FERRE preprocessing
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/ferre/processing.py>`_.

Current Python interface
------------------------

The apStar-only scope above defines the input workflow; an automatic FITS
reader and row-selection interface are not yet implemented in this package.
Currently, extract the selected arrays before constructing a :doc:`spectrum`
object or calling the fitting interfaces. Spec1D supports orders, errors,
masks, metadata, and mapping of matching detector pixels onto a native grid;
it does not automatically select FITS rows or interpret velocity headers.

Preserve the original pixel flags, propagate uncertainties for any additional
resampling, and avoid interpolating across detector gaps. Previously resampled
visit rows can already contain neighboring-pixel covariance; the current
fitting interfaces use per-pixel errors rather than a full covariance matrix.
