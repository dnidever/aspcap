Prepare your inputs
===================

Start with the example configuration files in the repository's
``examples/`` directory. Paths inside the run configuration are relative to
that configuration file.

Run configuration
-----------------

.. literalinclude:: ../examples/one-star.json
   :language: json

* ``star_id`` identifies the star in the output.
* ``grid`` points to a FERRE ASCII library or the ``.hdr`` of a matching
  ``.hdr`` and ``.unf`` pair. PCA-compressed, multi-segment, transposed, and
  sparse libraries are also accepted; see :doc:`grid`. The spectrum CSV must
  contain the reconstructed native pixels, concatenated in segment order,
  rather than PCA coefficients. The fit uses linear interpolation.
* ``spectrum`` points to the CSV described below.
* ``species_file`` points to the element configuration.
* ``wavelength_frame`` declares the air or vacuum convention and must match
  the grid header. Use ``air`` or ``vacuum``.
* ``parameters`` supplies the global fit as one value per grid label.
  Label names must match the header and values must lie within the grid.
* The optional ``ferre`` field specifies the executable. A path containing a
  slash is relative to the configuration file. ``--ferre`` overrides it.

Spectrum CSV
------------

Provide one row per native grid pixel in wavelength order, with columns
named ``wavelength,flux,error,continuum``. Wavelengths must already be in the
stellar rest frame and use the same units and air or vacuum convention as
the library.

Give flux and its one-sigma uncertainty in the same units as the global
continuum. The continuum must be finite and positive at every pixel. ASPCAP
divides both flux and uncertainty by it, then keeps the continuum fixed
throughout the fit. Nonfinite normalized flux and nonfinite or nonpositive
normalized uncertainties receive zero weight.

Element configuration
---------------------

The recommended element prescription is **one wavelength-weight file per
element**: positive weights define its windows, and zero excludes pixels.
No separate window file is needed. Load it with
``aspcap.windows.load_weights`` and pass the curve as ``weights=`` to ``Grid.fit_abundance``.
See :doc:`element-windows` for the formats and legacy readers.

.. literalinclude:: ../examples/mg.json
   :language: json

``species`` names the measurement and ``label`` names the only grid axis
FERRE can vary. Each pair in ``windows`` gives inclusive wavelength bounds.
Usable pixels in any window receive weight one, including pixels in
overlapping windows.

``add_labels`` lists fixed logarithmic coordinates to add to the fitted
value when reporting ``[X/H]``. In the example, the Mg-sensitive windows
measure a shared ``ALPHA`` coordinate relative to overall metallicity.
Adding the fixed ``METALS`` value converts it to ``[Mg/H]``. Omit
``add_labels`` when the fitted coordinate already represents ``[X/H]``.
Check that the chosen windows isolate the species and that this conversion
matches the conventions of your library.
