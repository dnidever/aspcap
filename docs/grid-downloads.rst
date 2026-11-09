Browse and Download Grids
=========================

Use the grid catalog to see available libraries, narrow the selection, inspect
file sizes, and download a complete FERRE bundle. Browsing the installed catalog
works offline. Downloading requires internet access; FERRE is not needed until
you interpolate or fit spectra.

The initial catalog contains 162 verified DR17 bundles: Synspec NLTE (the
primary DR17 analysis), Synspec LTE, Turbospectrum 20, Turbospectrum 20
plane-parallel comparisons, and ASSeT hot-star grids. It includes APOGEE-N
and APOGEE-S LSF variants. Entries come from the pinned DR17 pipeline grid
configurations and were checked against the public archive. This is a catalog
of production FERRE libraries, not every raw synthesis or every survey release.
See the `SDSS library description
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_.

Command line
------------

List the catalog or filter to one family and instrument:

.. code-block:: console

   aspcap grids list
   aspcap grids list --family M-giant --synthesis synspec-nlte --instrument apogee-n
   aspcap grids list --family M-giant --synthesis synspec-nlte --instrument apogee-n --lsf a

The table shows the selection ID, family, dimension count, and total size.
``--json`` returns complete metadata. Family filters are M-giant, GK-giant,
M-dwarf, GK-dwarf, F-dwarf, and BA. Synthesis filters are synspec-nlte,
synspec-lte, turbo20, turbo20-pp, and asset. Filters ignore case.

Inspect a selection and check its current remote availability before downloading:

.. code-block:: console

   aspcap grids info dr17-synspec-nlte-apogee-n-mg-a
   aspcap grids info dr17-synspec-nlte-apogee-n-mg-a --online
   aspcap grids download dr17-synspec-nlte-apogee-n-mg-a --dry-run

``info`` includes published temperature coverage, source URLs, and file sizes.
The catalog's sizes are snapshots; ``--online`` and ``--dry-run`` check the
server with HEAD requests without transferring the grid payload. Unknown
remote sizes are reported explicitly. Read the downloaded header for exact
axis coverage, labels, and conventions.

Download the chosen grid:

.. code-block:: console

   aspcap grids download dr17-synspec-nlte-apogee-n-mg-a --directory ./grids

The command prints the local header path when complete. Multiple IDs can be
supplied to download several selections serially. Relative directories resolve
from the working directory. Without ``--directory``, the parent directory is
``$XDG_CACHE_HOME/aspcap/grids`` or ``~/.cache/aspcap/grids``.
Each selection is stored in its own subdirectory.

Python
------

.. code-block:: python

   from aspcap import Grid, grids

   choices = grids.available(
       family="M-giant", synthesis="synspec-nlte", instrument="apogee-n", lsf="a"
   )
   selection = choices[0]
   print(selection["id"], selection["size"], selection["teff"])
   details = grids.info(selection["id"], online=True)
   header = grids.download(selection["id"], directory="grids")
   grid = Grid(header)

``available`` returns independent dictionaries and never downloads files.
``info`` describes a selection; ``download`` returns a pathlib.Path suitable
for Grid, Pipeline, or a recipe's grid-path setting. The optional progress
callback receives (filename, bytes_present, total_bytes_or_None).
A grid download does not select the correct LSF for an observation.
See :doc:`lsf` and :doc:`grid` before fitting.

Storage, interruption, and integrity
------------------------------------

Production grids are large. For example, the catalog's APO M-giant NLTE LSF-a
selection is about 11.70 GiB, including its PCA header/basis and coefficient
file. Disk space is checked against current known sizes before transfer.
Downloading does not establish that enough memory is available for fitting.
For initial tests, the :doc:`grid` page describes the small bundled fixtures;
they require no production download.

The downloader preserves filenames and fetches the complete listed bundle.
For these DR17 libraries that means the .hdr file (including its PCA basis)
and the .unf coefficient data. Additional listed sidecars are also retained.
Files remain in an ID.partial directory until all files pass length/checksum
checks and Grid can read the header. Completion publishes the bundle directory
and a download.json receipt with URLs, sizes, ETags, and local SHA256 hashes.
No archive extraction occurs.

Run the same download command after interruption. Complete staged files are
reused after verification; partial files resume when the server's ETag and
Range support allow it. Otherwise that file restarts. A completed cached bundle
is verified and reused without network access. Publisher SHA256 values are
checked when supplied by a catalog; the shipped DR17 catalog has size/ETag
metadata, so its locally computed hashes detect later corruption rather than
independently authenticating the first download. Corrupt completed bundles
raise an error and are not silently overwritten.

A lock prevents two writers downloading the same selection at once. A hard
process kill can leave a stale ID.lock directory: remove it only after
confirming the download is no longer running, then retry.

Extending the catalog
---------------------

``aspcap/data/grid-catalog.json`` records provenance and explicit file sets.
To use another collection without modifying ASPCAP, supply ``catalog=path``
to available, info, and download. A version-1 catalog has an entries list;
each entry requires a unique safe id, header filename, and files list with
name/url pairs. Optional file size, etag, and publisher sha256 fields support
inspection and verification. Include family, synthesis, instrument, lsf,
release, dimensions, and description to make entries browsable. Every file
needed by the library must be listed, preserving relative sidecar names as
simple filenames. The current downloader handles flat bundles; nested layouts
must be prepared separately.

Reference
---------

.. autofunction:: aspcap.grids.available

.. autofunction:: aspcap.grids.info

.. autofunction:: aspcap.grids.download
