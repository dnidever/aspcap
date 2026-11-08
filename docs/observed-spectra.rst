Observed Spectra
================

ASPCAP-style analysis can use an individual APOGEE visit spectrum or the
combined spectrum for a star. Both must ultimately be compared with synthetic
spectra on a compatible rest-frame wavelength scale. They differ in which
preparation steps have already been performed.

apVisit: an individual visit
----------------------------

An ``apVisit`` file (``asVisit`` for APOGEE-South) contains a visit's flux,
uncertainties, masks, and wavelength arrays for three detector chips. The
original visit wavelengths retain the observed Doppler shift. Prepare these
spectra by removing the measured shift and, for the current ASPCAP wrapper,
resampling onto the synthetic library's native pixels.

Velocity information is available in the processed visit headers and/or the
corresponding survey summary data, but keyword names and conventions depend
on the release. Historical readers used ``VRAD`` for the observed Doppler
shift. DR17 reports measured visit velocities as ``VREL``, the barycentric
correction as ``BC``, and the corrected velocity as ``VHELIO``. Despite its
name, ``VHELIO`` is barycentric. Do not apply it directly to observer-frame
wavelengths without accounting for the correction already included in it.
Check the release's velocity convention, available header fields, and RV flags.

For a receding velocity v that describes the shift of the actual wavelength
array, the survey combination prescription is:

.. math::

   \lambda_{\mathrm{rest}} = \frac{\lambda_{\mathrm{obs}}}{1 + v/c}.

Use consistent velocity units and conventions. First correcting the wavelength
array does not interpolate the measured flux; interpolation happens when the
spectrum is placed on a different set of pixels. Preserve the three chip
boundaries, masks, and uncertainties during that operation.

See the SDSS `visit RV guidance
<https://www.sdss4.org/dr17/irspec/use-radial-velocities/>`_ and
`RV pipeline description
<https://www.sdss4.org/dr17/irspec/radialvelocities/>`_.

apStar: a combined, rest-frame spectrum
---------------------------------------

The ``apStar``/``asStar`` product contains combined spectra and individual-visit
rows that have already been velocity-corrected and resampled. Select the intended
row: the combined product is not interchangeable with an original ``apVisit``
array. Do not apply the visit velocity correction a second time to those rows.

APOGEE combination removes each visit's measured Doppler shift, resamples with
sinc interpolation onto a common vacuum, logarithmic wavelength grid, and
combines the visits with weighting and bad-pixel treatment. The standard apStar
grid has 8,575 pixels with:

.. math::

   \log_{10}(\lambda_i / \mathrm{\mathring{A}}) = 4.179 + 6\times10^{-6}i,
   \qquad i=0,\ldots,8574.

A combined rest-frame spectrum still needs continuum normalization compatible
with the library. The combination process restores an average continuum shape;
``apStar`` does not mean an already ASPCAP-normalized spectrum. Inspect RV and
spectral quality flags as well as the pixel masks.

See `APOGEE visit combination
<https://www.sdss4.org/dr17/irspec/spectral_combination/>`_ and
`Using APOGEE Spectra <https://www.sdss4.org/dr17/irspec/spectra/>`_.

Which wavelength grid should be used for fitting?
-------------------------------------------------

The apStar grid is a useful standard, but the fitting target is the grid the
analysis actually uses. Some libraries contain only the modeled chip segments,
a subset of the apStar wavelengths, or a different sampling altogether. Correcting
a velocity, resampling a spectrum, and selecting library pixels are distinct
operations. Avoid repeated interpolation when direct preparation onto the
required fitting pixels is possible.

For the current wrapper, ``Grid.fit`` and ``fit_abundance`` require native
library pixels in the library's order. Our bundled full-range APOGEE fixture
uses the apStar wavelength grid, with zero-weight detector gaps, so apStar
sampling is the appropriate target for that fixture. This is not a requirement
for every possible FERRE library.

FERRE itself offers another route: supply the rest-frame observed wavelengths
and use ``WINTER=2`` to interpolate trial models onto them. Then the observed
flux need not be resampled. The current ASPCAP wrapper does not expose that
internal resampling configuration; see :doc:`ferre`.

What the earlier implementations did
------------------------------------

QASPCAP
~~~~~~~

The inspected QASPCAP reads either apVisit or apStar data and obtains ``VRAD``
from the header. It applies the first-order correction
``x2 = x * (1 - VRAD/c)``. Its usual ``WINTER=0`` path interpolates flux and
errors onto wavelengths read from the synthetic library, rather than hard-coding
the full apStar grid. Its ``lsf2d`` path sets ``WINTER=2``, retains selected
corrected observed wavelengths, and writes them for FERRE's model resampling.
This is a historical implementation description, not a recommended method for
interpolating modern bitmasks or propagating resampling uncertainties.

Sources: `qaspcap.pro
<https://github.com/callendeprieto/qaspcap/blob/636af1374158acce27cb02dfcf3966efd08eb1f1/qaspcap.pro#L213-L270>`_ and
`rd_apvisit.pro
<https://github.com/callendeprieto/qaspcap/blob/636af1374158acce27cb02dfcf3966efd08eb1f1/rd_apvisit.pro>`_.

DR17 ASPCAP
~~~~~~~~~~~

The inspected IDL reader accepts apVisit data and apStar spectra, including
optional visit rows already stored in apStar. It sorts the original apVisit
chips into increasing wavelength. The driver applies the header ``VRAD``
correction, handles air/vacuum compatibility, and interpolates each chip onto
wavelengths constructed from that library segment's header, in log wavelength.
Thus the target is the library's sampling, which can correspond to apStar chip
pixels without being the entire 8,575-pixel array. Header corrections must be
interpreted in the context of the product and release, avoiding a second
correction to already rest-frame data.

Sources: `aspcap_rdobsnew.pro
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/aspcap_rdobsnew.pro>`_ and
`doaspcap.pro
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro#L245-L285>`_.

Astra ASPCAP
~~~~~~~~~~~~

The inspected main ASPCAP entry point targets ``ApogeeCoaddedSpectrumInApStar``
objects. FERRE preprocessing copies those flux/error arrays and selects modeled
APOGEE pixels through its chip mask; wavelength interpolation defaults to zero.
This path uses already prepared coadded spectra. It should not be described as
a general original-apVisit Doppler-correction and resampling routine.

Sources: `ASPCAP entry point
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/__init__.py>`_ and
`FERRE preprocessing
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/ferre/processing.py>`_.

Current support and preparation responsibilities
------------------------------------------------

This package currently provides :doc:`spectrum` containers and explicit mapping
of already aligned chips onto native grid pixels. It does not yet provide
an apVisit/apStar FITS reader, automatic header-velocity correction, or a general
observed-spectrum resampling pipeline. ``Spec1D.from_native`` and ``to_native``
map matching pixels; they do not perform Doppler shifts or interpolate arbitrary
visit sampling. ``Grid`` output resampling resamples models, not observations.

For either product, preserve provenance, per-chip :doc:`lsf` information, and
bad-pixel flags. Propagate uncertainties using the chosen interpolation weights;
resampling can introduce neighboring-pixel covariance. Treat bitmasks as flags,
not continuous flux values, and do not interpolate across detector gaps. Then
apply :doc:`continuum` and proceed to :doc:`initial-estimates`,
:doc:`global-fitting`, and :doc:`element-fitting`.
