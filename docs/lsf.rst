Line-spread function (LSF)
==========================

The line-spread function describes how the instrument spreads an intrinsically
narrow spectral feature over detector pixels. Its width and shape affect line
profiles and therefore the comparison between observations and synthetic spectra.
An LSF can vary with wavelength, detector chip, and observing configuration.

Where it enters the workflow
----------------------------

Before :doc:`global-fitting` or :doc:`element-fitting`, the model spectra must
have an instrumental profile compatible with the observation. Continuum
normalization does not correct an LSF mismatch. Wavelength resampling alone
also does not perform instrumental broadening.

For a prebroadened library, choose a grid with the appropriate instrumental
profile. If starting with higher-resolution models, apply the appropriate
broadening before comparing them with the data. Avoid broadening an already
matched library a second time. A scalar resolving power is useful metadata,
but does not describe every wavelength-dependent or non-Gaussian LSF.

Two ways to account for the LSF
-------------------------------

Pre-convolved spectral grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Convolve every synthetic spectrum with a representative APOGEE LSF when
building the library, then interpolate and fit that pre-convolved library.
Separate versions of a grid can represent different instrumental profiles.

Historically, APOGEE used five LSF versions of each grid: one for an average
LSF and four for different fiber ranges. Thus the five versions were not five
distinct fiber bins. The appropriate version approximates the instrumental
profile of the observation. For combined spectra, the survey used mean fiber
information to select among the fiber-dependent grids.

This approach moves convolution out of the fitting loop, making repeated fits
less expensive, at the cost of storing several libraries and approximating
the individual observation's LSF by a representative profile. Keep the LSF
group, telescope, and library provenance with the results.

See the historical discussion in `Holtzman et al. (2018)
<https://research.monash.edu/files/255726790/255350468_oa.pdf>`_ and the
`DR17 spectral-library documentation
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_. The latter describes
four fiber regions; the exact library set and selection rules depend on the
survey release.

On-the-fly LSF convolution
~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, keep a sufficiently high-resolution synthetic library,
interpolate a trial spectrum at the requested stellar labels, and convolve
that model with the observation's LSF before evaluating its residuals.
Apply the appropriate wavelength-dependent profile separately across orders
or detector segments, with sampling handled consistently.

This approach can use an observation-specific LSF without storing a separate
grid for every profile, but adds work to model evaluation. Kernels or convolution
operators can be prepared once and reused when the LSF is fixed during a fit.
Using a known LSF on the fly does not imply fitting the LSF parameters themselves.

Both approaches must account for the model's existing broadening and sampling.
Avoid convolving an already matched library a second time. Broadening cannot
recover resolution absent from the input models. The same LSF treatment should
be used consistently in global and individual-element fits.

Current package support
-----------------------

``Spec1D`` accepts an optional ``lsf`` attribute and preserves it as opaque
metadata. ASPCAP does not interpret that object, construct an LSF, convolve
spectra, or automatically select a matching library. Supplying ``lsf`` therefore
does not change a model evaluation or fit.

.. code-block:: python

   from aspcap import Spec1D

   spec = Spec1D(
       flux, wave=wavelength, err=error,
       lsf={"description": "User-supplied instrumental profile metadata"},
   )

This is a storage example, not an LSF model API. A Doppler LSF object can also
be retained as metadata without importing Doppler's fitting or model loaders.

``Grid.resolution`` records the library header's positive resolving power.
``Grid.__call__``, ``Grid.model``, and ``Grid.fit`` use the library as supplied;
they do not match it to ``spec.lsf``. The bundled APOGEE example already comes
from a specific SDSS LSF-a library; see :doc:`grid` for its provenance. It is
an example fixture rather than an automatic LSF match for every APOGEE spectrum.

Preparing an observation
------------------------

Record the provenance of the model broadening and the observation's instrumental
profile. Check compatibility separately for each order or chip, retain masks
for detector gaps and unusable pixels, and compare model and observed line
profiles before interpreting an abundance residual as a composition change.
Automatic LSF matching and convolution are future extensions to this package.
