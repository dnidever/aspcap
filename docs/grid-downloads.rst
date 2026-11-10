Browse and Download Grids
=========================

Use the grid catalog to see available libraries, narrow the selection, inspect
file sizes, and download a complete FERRE bundle. Browsing the installed catalog
works offline. Downloading requires internet access; FERRE is not needed until
you interpolate or fit spectra.

The catalog includes the initial derived HGE giant test grid,
``hge-giants-test-v1``, and 162 verified DR17 bundles: Synspec NLTE (the
primary DR17 analysis), Synspec LTE, Turbospectrum 20, Turbospectrum 20
plane-parallel comparisons, and ASSeT hot-star grids. It includes APOGEE-N
and APOGEE-S LSF variants. Entries come from the pinned DR17 pipeline grid
configurations and were checked against the public archive. This is a catalog
of selected FERRE libraries, not every raw synthesis or every survey release.
See the `SDSS library description
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_.

HGE giant test grid
-------------------

This new interface is available in the current GitHub version; see
:doc:`installation` for development-version installation.

The versioned selection ``hge-giants-test-v1`` is the derived six-label grid:
3000–5500 K, raw log g up to 3, approved excluded corners, prescribed
microturbulence, and shared PCA105. Its sources are DR17 Synspec NLTE M/GK
libraries preconvolved with **APOGEE-N LSF-a**. It is an initial test version,
not an official DR17 production grid. See :doc:`hge-giant-grid` for accuracy
limits and provenance.

.. code-block:: python

   from aspcap import Grid

   grid = Grid.from_download("hge-giants-test-v1")
   spec = grid({"C": 0, "N": 0, "O Mg Si S Ca Ti": 0,
                "METALS": -0.5, "LOGG": 1.5, "TEFF": 4500})

The bundle is about 3.20 GiB including provenance and training files. The first
call downloads from public Google Drive, checks publisher SHA256 hashes, and
opens the grid using direct disk access (``f_access=1``). Subsequent calls
verify and reuse it offline. No Google sign-in, API key, connected account,
or browser cookies are used. ``directory="grids"`` chooses a different cache
parent. Progress is printed by default; ``progress=False`` is quiet and
``progress=callback`` supplies custom reporting.

The versioned ID pins exact file content. Future test versions receive new IDs;
there is no ambiguous ``hge-giants`` or automatic latest-version substitution.
If Drive sharing or traffic limits block a download, an explicit error directs
you to the source link; HTML error pages are never accepted as grid data.

Command line
------------

List the catalog or filter to one family and instrument:

.. code-block:: console

   aspcap grids list
   aspcap grids info hge-giants-test-v1
   aspcap grids download hge-giants-test-v1
   aspcap grids list --family M-giant --synthesis synspec-nlte --instrument apogee-n
   aspcap grids list --family M-giant --synthesis synspec-nlte --instrument apogee-n --lsf a

The table shows the selection ID, family, dimension count, and total size.
``--json`` returns complete metadata. Family filters are M-giant, GK-giant,
M-dwarf, GK-dwarf, F-dwarf, BA, and HGE-giant. Synthesis filters are synspec-nlte,
synspec-lte, turbo20, turbo20-pp, and asset. Filters ignore case.

Inspect a selection and check its current remote availability before downloading:

.. code-block:: console

   aspcap grids info dr17-synspec-nlte-apogee-n-mg-a
   aspcap grids info dr17-synspec-nlte-apogee-n-mg-a --online
   aspcap grids download dr17-synspec-nlte-apogee-n-mg-a --dry-run

``info`` includes published temperature coverage, source URLs, and file sizes.
The catalog's sizes are snapshots; ``--online`` and ``--dry-run`` check the
server with HEAD requests for ordinary HTTP files and one-byte GET probes
for Drive files. Unknown
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
   grid = Grid(header, f_access=1)

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
Range support allow it. Drive files can also resume without an ETag when
their pinned publisher checksum identity matches. A server ignoring Range
causes a safe restart. Otherwise that file restarts. A completed cached bundle
is verified and reused without network access. Publisher SHA256 values are
checked when supplied by a catalog; the historical DR17 entries have size/ETag
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
simple filenames. Drive files additionally require ``provider="google-drive"``,
``drive_file_id``, a matching ``https://drive.google.com/file/d/ID/view`` URL,
a positive byte size, and a publisher SHA256. Prefer versioned IDs and never
replace a published release's checksum identity. The downloader handles flat bundles; nested layouts
must be prepared separately.

Reference
---------

.. autofunction:: aspcap.grids.available

.. autofunction:: aspcap.grids.info

.. autofunction:: aspcap.grids.download

.. automethod:: aspcap.Grid.from_download
