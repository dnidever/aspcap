"""Build the full-wavelength 5 x 5 x 5 test grid from the SDSS DR17 archive."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np


SOURCE = (
    "https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/"
    "synspec/marcs/giantisotopes/sgM_200921nlte_lsfa/"
    "p_apssgM_200921nlte_lsfa_012_075"
)
STEM = "apogee-dr17-mgiant-full"
SOURCE_SHAPE = (5, 11, 6, 8, 15, 8, 11)
COMPONENTS = 75
BLOCK_LENGTHS = (627,) * 11 + (617,)
RECORD_BYTES = len(BLOCK_LENGTHS) * COMPONENTS * 4
SOURCE_BYTES = 12_545_280_000
# Zero-based offsets on the APOGEE apStar wavelength grid.
SEGMENTS = ((246, 3028), (3585, 2495), (6344, 1991))
NPIX = 8575

HEADER = """&SYNTH
 ID = 'APOGEE DR17 M giant 5x5x5'
 N_OF_DIM = 3
 N_P = 5 5 5
 LABEL(1) = 'METALS'
 LABEL(2) = 'LOGG'
 LABEL(3) = 'TEFF'
 LLIMITS = -1.0 0.0 3400
 STEPS = 0.25 0.5 100
 NPIX = 8575
 WAVE = 4.179000 0.000006
 LOGW = 1
 VACUUM = 1
 RESOLUTION = 22500
/
"""


def read_coefficients(indices):
    """Read one little-endian PCA record via an HTTP range request."""
    # Fix LOG10VDOP, C, N, and the shared alpha axis at the chosen source nodes.
    record = np.ravel_multi_index((3, 6, 1, 3, *indices), SOURCE_SHAPE)
    first = record * RECORD_BYTES
    last = first + RECORD_BYTES - 1
    request = Request(SOURCE + ".unf", headers={"Range": f"bytes={first}-{last}"})
    with urlopen(request, timeout=60) as response:
        expected = f"bytes {first}-{last}/{SOURCE_BYTES}"
        if response.status != 206 or response.headers.get("Content-Range") != expected:
            raise RuntimeError(f"SDSS did not return the expected range for record {record}")
        data = response.read(RECORD_BYTES + 1)
    if len(data) != RECORD_BYTES:
        raise RuntimeError(f"Incomplete SDSS PCA record {record}")
    return np.frombuffer(data, dtype="<f4").reshape(len(BLOCK_LENGTHS), COMPONENTS)


def build_grid(output):
    """Reconstruct the retained models and place them on the apStar pixel grid."""
    with urlopen(SOURCE + ".hdr", timeout=60) as response:
        header = response.read().decode("ascii")
    # Four namelists precede the mean, eigenvalues, and basis vectors.
    parts = header.split("\n /\n", 4)
    if len(parts) != 5:
        raise RuntimeError("Unexpected SDSS PCA header structure")
    wave_count = sum(BLOCK_LENGTHS)
    pca = np.fromstring(parts[-1], sep=" ", dtype=np.float32)
    if pca.size != (COMPONENTS + 2) * wave_count:
        raise RuntimeError("Unexpected number of values in the SDSS PCA basis")
    pca = pca.reshape(-1, wave_count)
    mean, basis = pca[0], pca[2:]
    # FERRE stores the last axis (TEFF) fastest. Indices select five
    # adjacent source nodes on each axis.
    nodes = product(range(6, 11), range(1, 6), range(4, 9))
    with ThreadPoolExecutor(max_workers=8) as pool:
        coefficients = np.stack(list(pool.map(read_coefficients, nodes)))

    blocks = []
    offset = 0
    for i, count in enumerate(BLOCK_LENGTHS):
        section = slice(offset, offset + count)
        blocks.append(mean[section] + coefficients[:, i, :] @ basis[:, section])
        offset += count
    source_flux = np.concatenate(blocks, axis=1)
    if not np.isfinite(source_flux).all():
        raise RuntimeError("SDSS PCA reconstruction produced nonfinite flux")

    # Keep every detector sample. Unity is a placeholder where there is no
    # model coverage. The accompanying weights exclude those pixels in fits.
    flux = np.ones((len(coefficients), NPIX), dtype="<f4")
    weights = [0] * NPIX
    offset = 0
    for start, count in SEGMENTS:
        flux[:, start:start + count] = source_flux[:, offset:offset + count]
        weights[start:start + count] = [1] * count
        offset += count
    output.mkdir(parents=True, exist_ok=True)
    (output / f"{STEM}.unf").write_bytes(flux.tobytes())
    (output / f"{STEM}.hdr").write_text(HEADER, encoding="ascii")
    (output / f"{STEM}-weights.json").write_text(json.dumps(weights) + "\n")
    print(f"Wrote {len(coefficients)} spectra and detector weights to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="directory for the grid and weights")
    build_grid(parser.parse_args().output)
