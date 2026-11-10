"""Numerical checks against independent SVD, not implementation-shaped fixtures."""

import numpy as np
import pytest
from aspcap.pca import PCABasis


def test_streamed_physical_pca_matches_svd():
    rng = np.random.default_rng(18)
    flux = 1e6 + rng.normal(size=(73, 15)) * np.arange(1, 16)
    model = PCABasis.fit(flux, 3, [8, 7], batch_size=11)
    reconstruction = model.reconstruct(model.project(flux))
    expected = []
    for sl in (slice(0, 8), slice(8, 15)):
        values = flux[:, sl]
        mean = values.mean(axis=0)
        _, _, vectors = np.linalg.svd(values - mean, full_matrices=False)
        expected.append(mean + (values - mean) @ vectors[:3].T @ vectors[:3])
    np.testing.assert_allclose(
        reconstruction, np.concatenate(expected, axis=1), atol=1e-7, rtol=0
    )


def test_nonorthogonal_source_recompression_matches_physical_training():
    rng = np.random.default_rng(7)
    B = rng.normal(size=(5, 19))
    B[4] += 3 * B[0]
    mean = rng.normal(size=19)
    c = rng.normal(size=(90, 2, 5))
    flux = np.concatenate(
        [mean[:10] + c[:, 0] @ B[:, :10], mean[10:] + c[:, 1] @ B[:, 10:]], axis=1
    )
    direct = PCABasis.fit(flux, 3, [10, 9], batch_size=17)
    compressed = PCABasis.fit_coefficients(c, mean, B, 3, [10, 9], batch_size=13)
    np.testing.assert_allclose(
        compressed.reconstruct(compressed.project_coefficients(c, mean, B)),
        direct.reconstruct(direct.project(flux)),
        atol=1e-11,
    )
    np.testing.assert_allclose(
        compressed.project_coefficients(c, mean, B),
        compressed.project(flux),
        atol=1e-12,
    )


def test_roundtrip_and_single_spectrum(tmp_path):
    rng = np.random.default_rng(3)
    flux = rng.normal(size=(30, 12))
    basis = PCABasis.fit(flux, 6, [6, 6])
    p = tmp_path / "basis.npz"
    basis.save(p)
    loaded = PCABasis.load(p)
    np.testing.assert_allclose(
        loaded.reconstruct(loaded.project(flux[0])), flux[0], atol=1e-12
    )
    assert loaded.project(flux[0]).shape == (2, 6)


@pytest.mark.parametrize(
    "data,rank,lengths",
    [
        (np.ones((1, 8)), 1, [8]),
        (np.ones((4, 8)), 5, [4, 4]),
        (np.ones((4, 8)), 1, [3, 4]),
        (np.full((4, 8), np.nan), 1, [8]),
    ],
)
def test_rejects_invalid_training(data, rank, lengths):
    with pytest.raises(ValueError):
        PCABasis.fit(data, rank, lengths)


def test_ferre_recompression_preserves_physical_models(tmp_path):
    import importlib.util
    from pathlib import Path
    import f90nml

    script = (
        Path(__file__).resolve().parents[1]
        / "examples/hge-giant-grid/recompress_ferre.py"
    )
    spec = importlib.util.spec_from_file_location("compress_example", script)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    rng = np.random.default_rng(14)
    source_basis = PCABasis.fit(rng.normal(size=(40, 8)), 3, [4, 4])
    old_c = rng.normal(size=(20, 2, 3)).astype("<f4")
    physical = source_basis.reconstruct(old_c)
    new_basis = PCABasis.fit(physical, 2, [4, 4])
    saved = tmp_path / "basis.npz"
    new_basis.save(saved)
    header = tmp_path / "source.hdr"
    root = dict(
        n_of_dim=1,
        label=["TEFF"],
        n_p=[20],
        llimits=[3000],
        steps=[100],
        npix=6,
        npca=[4, 4],
        wave=[15000, 1],
        logw=0,
        resolution=22500,
    )
    with header.open("w") as f:
        f.write(str(f90nml.Namelist(synth=root)) + "\n")
        np.savetxt(
            f, np.vstack([source_basis.mean, np.ones(8), source_basis.components])
        )
    header.with_suffix(".unf").write_bytes(old_c.tobytes())
    destination = tmp_path / "compressed"
    example.compress(header, saved, destination, batch_size=7)
    meta, _, _, mean, vectors, _, new_c = example.read_source(
        destination / "source.hdr"
    )
    actual = np.concatenate(
        [
            mean[:4] + new_c[:, 0] @ vectors[:, :4],
            mean[4:] + new_c[:, 1] @ vectors[:, 4:],
        ],
        axis=1,
    )
    np.testing.assert_allclose(
        actual, new_basis.reconstruct(new_basis.project(physical)), atol=2e-7
    )
    assert meta["nspectra"] == 20
    assert meta["axes"][0] == tuple(range(3000, 5000, 100))
    with pytest.raises(FileExistsError):
        example.compress(header, saved, destination)
