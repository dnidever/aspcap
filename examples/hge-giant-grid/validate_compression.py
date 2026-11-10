"""Validate a recompressed sparse grid against its original PCA representation."""

from pathlib import Path
from unittest.mock import patch
import argparse, json, subprocess
import numpy as np
import f90nml
from aspcap import Grid, Spec1D, PCABasis
from aspcap.windows import load_dr17
from recompress_ferre import read_source, sha256


def validate(source_header, compressed_header, executable="ferre"):
    meta, root, segments, mean, vectors, lengths, coefficients = read_source(
        source_header
    )
    folder = Path(compressed_header).resolve().parent
    basis = PCABasis.load(folder / "pca-basis.npz")
    grid = Grid(compressed_header, executable=executable, f_access=1)
    if not root.get("nodefile"):
        raise ValueError("This validation example requires NODEFILE coordinates")
    nodes = Path(root["nodefile"])
    if not nodes.is_absolute():
        nodes = Path(source_header).resolve().parent / nodes
    coords = np.loadtxt(nodes)
    if sha256(nodes) != sha256(folder / (Path(compressed_header).stem + ".nodes")):
        raise ValueError("Coordinate order changed")
    rng = np.random.default_rng(901)
    sample = list(map(int, rng.choice(len(coefficients), 32, replace=False)))
    truth = {
        "C": 0.0,
        "N": 0.0,
        "O Mg Si S Ca Ti": 0.0,
        "METALS": -0.5,
        "LOGG": 1.5,
        "TEFF": 4500.0,
    }
    center = int(
        np.flatnonzero(
            np.all(
                np.isclose(coords, np.array([truth[x] for x in grid.labels])), axis=1
            )
        )[0]
    )
    sample.append(center)
    expected = basis.reconstruct(
        basis.project_coefficients(coefficients[sample], mean, vectors)
    )
    models = grid([dict(zip(grid.labels, coords[i])) for i in sample])
    error = max(
        float(np.max(np.abs(model.flux - reference)))
        for model, reference in zip(models, expected)
    )
    if error > 3e-5:
        raise ValueError(f"FERRE/Python reconstruction mismatch: {error}")
    # Independent original representation supplies the synthetic observation.
    source = []
    start = 0
    for k, n in enumerate(lengths):
        source.append(
            mean[start : start + n]
            + coefficients[center, k] @ vectors[:, start : start + n]
        )
        start += n
    observation = Spec1D(
        np.concatenate(source), wave=grid.wavelengths, err=np.full(grid.npix, 0.01)
    )
    initial = dict(truth, C=0.02, N=0.04, METALS=-0.47, LOGG=1.46, TEFF=4480.0)
    controls = dict(stopcr=1e-9, simp=1e-9, scope=0.05)
    original = subprocess.run

    def precise(*args, **kwargs):
        p = Path(kwargs["cwd"]) / "input.nml"
        n = f90nml.read(p)
        n["lista"].update(controls)
        n.write(p, force=True)
        return original(*args, **kwargs)

    with patch("aspcap.ferre.subprocess.run", side_effect=precise):
        fit = grid.fit(observation, start=initial)
    errors = {key: fit.parameters[key] - value for key, value in truth.items()}
    if abs(errors["TEFF"]) > 10 or any(
        abs(value) > 0.01 for key, value in errors.items() if key != "TEFF"
    ):
        raise ValueError("Synthetic global recovery failed")
    _, weights = load_dr17("Mg", wavelength=grid.wavelengths)
    initial = dict(truth)
    initial["O Mg Si S Ca Ti"] = 0.03
    with patch("aspcap.ferre.subprocess.run", side_effect=precise):
        element = grid.fit(
            observation,
            start=initial,
            free_labels=("O Mg Si S Ca Ti",),
            weights=weights,
        )
    if abs(element.parameters["O Mg Si S Ca Ti"]) > 0.01:
        raise ValueError("Weighted Mg recovery failed")
    coverage = []
    for t, g, want in [
        (3300, 0, True),
        (3300, 3, False),
        (5000, 0, False),
        (5000, 2, True),
        (5500, 3.5, False),
    ]:
        actual = grid.contains(dict(truth, TEFF=t, LOGG=g))
        coverage.append(dict(teff=t, logg=g, covered=actual))
        if actual != want:
            raise ValueError("Sparse coverage changed")
    manifest = json.loads((folder / "compression.json").read_text())
    for name, entry in manifest["files"].items():
        if sha256(folder / name) != entry["sha256"]:
            raise ValueError("Checksum mismatch")
    result = dict(
        status="validated",
        spectra=grid.nspectra,
        components_per_block=basis.ncomponents,
        wavelength_pixels=grid.npix,
        reconstruction_samples=len(sample),
        max_FERRE_Python_difference=error,
        global_fit=dict(truth=truth, recovered=fit.parameters, signed_errors=errors),
        Mg_fit=dict(truth=0, recovered=element.parameters["O Mg Si S Ca Ti"]),
        coverage=coverage,
        validation_controls=controls,
        note="Compression loss is characterized separately by held-out spectra and element-window surveys; this check verifies the written library and representative FERRE fits.",
    )
    (folder / "compression-validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    manifest["status"] = "validated"
    manifest["validation_file"] = "compression-validation.json"
    (folder / "compression.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("compressed", type=Path)
    p.add_argument("--ferre", default="ferre")
    a = p.parse_args()
    validate(a.source, a.compressed, a.ferre)
