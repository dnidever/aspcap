"""Run the two fitting stages on a synthetic observation from the toy grid."""

from pathlib import Path
import sys

import numpy as np
from aspcap import Grid, Pipeline


def main():
    grid = Grid(Path(__file__).with_name("toy-grid.dat"))
    spectrum = grid({"ALPHA": 0.2, "METALS": -0.25})
    # Give the synthetic observation a known continuum and one-sigma errors.
    spectrum.flux *= 2
    spectrum.err = np.full(grid.npix, 0.02)
    pipeline = Pipeline(
        grid,
        elements=[
            dict(
                species="Mg",
                label="ALPHA",
                windows=[(15002, 15005)],
                add_labels=["METALS"],
            )
        ],
        continuum=np.full(grid.npix, 2.0),
    )
    result = pipeline.fit(
        spectrum, initial={"ALPHA": 0.0, "METALS": -0.25}, star_id="synthetic-star"
    )
    if not result.success:
        print(result.errors, file=sys.stderr)
        return 1
    print("Global parameters:", result.parameters)
    print("Abundances [X/H]:", result.abundances)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
