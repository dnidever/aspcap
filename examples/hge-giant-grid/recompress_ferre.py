"""Train/recompress spectrum-major float32 FERRE PCA libraries in bounded memory.

This reader supports the specific DR17/HGE layout: .hdr with ASCII mean,
variance and basis rows, adjacent .unf with little-endian float32 coefficient
records. It does not read transposed, big-endian, formatted or sequential
Fortran-record payloads. Use PCABasis directly for other input formats.
"""

from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import f90nml
from aspcap.pca import PCABasis
from aspcap._grid_headers import _header, read_library


def read_source(header):
    """Return metadata, spectral basis and memory-mapped source coefficients."""
    header = Path(header).resolve()
    metadata = read_library(header)
    with header.open("rb") as stream:
        root = _header(stream)
        segments = [_header(stream) for _ in range(int(root.get("multi", 0)))]
        lengths = np.atleast_1d(root.get("npca", [])).astype(int)
        stored = int(root["npix"])
        if not len(lengths) or np.any(lengths <= 0) or stored % len(lengths):
            raise ValueError(
                "Expected a PCA library with equal component counts per block"
            )
        rank = stored // len(lengths)
        values = np.fromstring(stream.read().decode("ascii"), sep=" ")
    if values.size != (rank + 2) * sum(lengths):
        raise ValueError("Unexpected PCA basis layout")
    values = values.reshape(rank + 2, sum(lengths))
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite source PCA values")
    # FERRE physical reconstruction uses mean and basis; row 1 is VPCA metadata.
    payload = header.with_suffix(".unf")
    if payload.stat().st_size != metadata["nspectra"] * stored * 4:
        raise ValueError("Payload size does not match spectrum-major float32 records")
    coefficients = np.memmap(
        payload, dtype="<f4", mode="r", shape=(metadata["nspectra"], len(lengths), rank)
    )
    return metadata, root, segments, values[0], values[2:], lengths, coefficients


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def train(header, destination, components, rows=None, batch_size=1024):
    """Train only on selected zero-based record indices (all records by default)."""
    _, _, _, mean, vectors, lengths, coefficients = read_source(header)
    if Path(destination).exists():
        raise FileExistsError(destination)
    indices = (
        np.arange(len(coefficients))
        if rows is None
        else np.load(rows, allow_pickle=False)
    )
    if (
        indices.ndim != 1
        or indices.dtype.kind not in "iu"
        or not len(indices)
        or np.any(indices < 0)
        or np.any(indices >= len(coefficients))
        or len(np.unique(indices)) != len(indices)
    ):
        raise ValueError("Training rows must be unique valid integer record indices")
    batches = (
        coefficients[indices[i : i + batch_size]]
        for i in range(0, len(indices), batch_size)
    )
    basis = PCABasis.from_batches(
        batches, components, lengths, source_mean=mean, source_components=vectors
    )
    basis.save(destination)
    Path(str(destination) + ".json").write_text(
        json.dumps(
            dict(
                source_header=str(Path(header).resolve()),
                components=components,
                training_rows=len(indices),
                row_indices_sha256=hashlib.sha256(
                    indices.astype("<i8").tobytes()
                ).hexdigest(),
                note="Uniform physical-flux PCA; no region selection beyond the supplied row list",
            ),
            indent=2,
        )
        + "\n"
    )


def compress(header, basis_file, directory, batch_size=1024):
    """Write a new bundle, preserving coordinate order and wavelength metadata.

    The destination must not exist. The source is never changed or deleted.
    Inspect compression.json and validate with FERRE before replacing a grid.
    """
    meta, root, segments, mean, vectors, lengths, coefficients = read_source(header)
    basis = PCABasis.load(basis_file)
    if not np.array_equal(lengths, basis.block_lengths):
        raise ValueError("Source and target wavelength blocks must match")
    folder = Path(directory).resolve()
    if folder.exists():
        raise FileExistsError(folder)
    folder.mkdir(parents=True)
    stem = Path(header).stem
    root = dict(root)
    root["npix"] = basis.nblocks * basis.ncomponents
    if root.get("nodefile"):
        source = Path(root["nodefile"])
        if not source.is_absolute():
            source = Path(header).resolve().parent / source
        root["nodefile"] = stem + ".nodes"
        shutil.copy2(source, folder / root["nodefile"])
    with (folder / (stem + ".hdr")).open("w") as stream:
        stream.write(str(f90nml.Namelist(synth=root)) + "\n")
        for segment in segments:
            stream.write(str(f90nml.Namelist(synth=segment)) + "\n")
        np.savetxt(stream, basis.mean[None], fmt="%.10e")
        np.savetxt(stream, np.ones((1, len(mean))), fmt="%.10e")
        np.savetxt(stream, basis.components, fmt="%.10e")
    # Cache affine transforms instead of expanding flux or recomputing matrix products per batch.
    transforms = [vectors[:, sl] @ basis.components[:, sl].T for sl in basis.slices]
    offsets = [
        (mean[sl] - basis.mean[sl]) @ basis.components[:, sl].T for sl in basis.slices
    ]
    with (folder / (stem + ".unf")).open("wb") as stream:
        for first in range(0, len(coefficients), batch_size):
            chunk = np.asarray(coefficients[first : first + batch_size], dtype=float)
            if not np.isfinite(chunk).all():
                raise ValueError("Nonfinite source coefficients")
            result = np.stack(
                [
                    chunk[:, k] @ transform + offsets[k]
                    for k, transform in enumerate(transforms)
                ],
                axis=1,
            ).astype("<f4")
            if not np.isfinite(result).all():
                raise ValueError("Nonfinite compressed coefficients")
            stream.write(result.tobytes())
            if first % (batch_size * 64) == 0:
                print(
                    f"Compressed {first:,} / {len(coefficients):,} records", flush=True
                )
    basis.save(folder / "pca-basis.npz")
    files = {
        path.name: dict(bytes=path.stat().st_size, sha256=sha256(path))
        for path in [folder / (stem + ".hdr"), folder / (stem + ".unf")]
        + list(folder.glob("*.nodes"))
    }
    (folder / "SHA256SUMS").write_text(
        "".join(v["sha256"] + "  " + k + "\n" for k, v in files.items())
    )
    (folder / "compression.json").write_text(
        json.dumps(
            dict(
                spectra=meta["nspectra"],
                components_per_block=basis.ncomponents,
                wavelength_pixels=len(mean),
                block_lengths=lengths.tolist(),
                basis_sha256=sha256(basis_file),
                source_header=str(Path(header).resolve()),
                source_header_sha256=sha256(header),
                files=files,
                status="written_pending_FERRE_validation",
                layout="spectrum-major little-endian float32; original coordinate order retained",
            ),
            indent=2,
        )
        + "\n"
    )
    print("Written:", folder, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("train")
    p.add_argument("header", type=Path)
    p.add_argument("basis", type=Path)
    p.add_argument("--components", type=int, required=True)
    p.add_argument("--rows", type=Path)
    p.add_argument("--batch-size", type=int, default=1024)
    p = sub.add_parser("compress")
    p.add_argument("header", type=Path)
    p.add_argument("basis", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--batch-size", type=int, default=1024)
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("batch-size must be positive")
    if args.command == "train":
        train(args.header, args.basis, args.components, args.rows, args.batch_size)
    else:
        compress(args.header, args.basis, args.output, args.batch_size)
