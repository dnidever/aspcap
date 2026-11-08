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
