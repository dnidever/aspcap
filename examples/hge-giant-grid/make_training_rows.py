"""Reproduce the adopted HGE abundance-held-out training selection."""

from pathlib import Path
import argparse, json
import numpy as np


def make_rows(nodefile, destination, proposal=None):
    """Save zero-based rows for 6,322 complete training abundance combinations.

    Input columns are C, N, alpha, METALS, LOGG, TEFF in spectrum-record order.
    The input must be the approved HGE grid, not an arbitrary source library.
    """
    coordinates = np.loadtxt(nodefile, ndmin=2)
    if coordinates.shape != (673200, 6):
        raise ValueError("Expected the full approved HGE coordinate table")
    proposal = (
        Path(proposal)
        if proposal
        else Path(__file__).with_name("hge-giant-grid-proposal.json")
    )
    allowed = {
        tuple(pair) for pair in json.loads(proposal.read_text())["teff_logg_nodes"]
    }
    if {(row[5], row[4]) for row in coordinates} != allowed:
        raise ValueError("Temperature/gravity region differs from the approved grid")
    combos, inverse, counts = np.unique(
        coordinates[:, :4], axis=0, return_inverse=True, return_counts=True
    )
    if (
        len(combos) != 7920
        or np.any(counts != 85)
        or len(np.unique(coordinates, axis=0)) != len(coordinates)
    ):
        raise ValueError(
            "Abundance combinations must each contain the approved 85 nodes"
        )
    corners = np.flatnonzero(
        np.all((combos == combos.min(axis=0)) | (combos == combos.max(axis=0)), axis=1)
    )
    validation = np.unique(
        np.r_[
            np.random.default_rng(20261109).choice(7920, 1584, replace=False), corners
        ]
    )
    rows = np.flatnonzero(~np.isin(inverse, validation))
    if Path(destination).exists():
        raise FileExistsError(destination)
    with Path(destination).open("wb") as f:
        np.save(f, rows.astype("<i8"), allow_pickle=False)
    print(
        f"Saved {len(rows):,} training records; {len(coordinates)-len(rows):,} held out."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("nodefile", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--proposal", type=Path)
    a = p.parse_args()
    make_rows(a.nodefile, a.output, a.proposal)
