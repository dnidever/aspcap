"""Command-line interface for a single HGE abundance fit."""

import argparse
import csv
import json
from pathlib import Path

from . import Grid


def _run(args):
    path = args.config.resolve()
    config = json.loads(path.read_text())
    species = json.loads((path.parent / config["species_file"]).read_text())
    grid = Grid.read(path.parent / config["grid"])
    if grid.vacuum is None:
        raise ValueError(
            "The CSV workflow requires segments with a common wavelength frame"
        )
    expected_frame = "vacuum" if grid.vacuum else "air"
    if config["wavelength_frame"] != expected_frame:
        raise ValueError(f"Grid uses {expected_frame} wavelengths")
    with (path.parent / config["spectrum"]).open(newline="") as file:
        rows = list(csv.DictReader(file))
    columns = [
        [float(row[name]) for row in rows]
        for name in ("wavelength", "flux", "error", "continuum")
    ]
    executable = args.ferre or config.get("ferre", "ferre")
    if not args.ferre and "/" in executable:
        executable = path.parent / executable
    result, pixels = grid.fit_abundance(
        *columns,
        config["parameters"],
        star_id=config["star_id"],
        species=species["species"],
        label=species["label"],
        windows=species["windows"],
        add_labels=species.get("add_labels", ()),
        executable=executable,
        keep_runs=args.keep_runs,
        run_root=args.run_root,
    )
    if args.pixels:
        with args.pixels.open("w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                ("wavelength", "normalized_flux", "normalized_error", "model", "weight")
            )
            writer.writerows(pixels)
    with args.output.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=result)
        writer.writeheader()
        writer.writerow(result)
    print(f"Wrote {args.output}")


def _recipe_run(args):
    """Execute a text recipe on a prepared, flat CSV spectrum."""
    from .recipe import PipelineRecipe
    from .spectrum import Spec1D

    recipe = PipelineRecipe.read(args.config, executable=args.ferre)
    with args.spectrum.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    spectrum = Spec1D(
        [float(row["flux"]) for row in rows],
        wave=[float(row["wavelength"]) for row in rows],
        err=[float(row["error"]) for row in rows],
        mask=[row.get("mask", "0").lower() in ("1", "true") for row in rows],
        wavevac=args.wavelength_frame == "vacuum",
    )
    result = recipe.fit(spectrum, star_id=args.star_id)

    def serial(value):
        if hasattr(value, "tolist"):
            return value.tolist()
        raise TypeError(f"Cannot serialize {type(value)}")

    args.output.write_text(
        json.dumps(
            dict(
                success=result.success,
                parameters=result.parameters,
                abundances=result.abundances,
                errors=result.errors,
                history=result.history,
                config=result.config,
            ),
            indent=2,
            default=serial,
            allow_nan=False,
        )
        + "\n"
    )
    if not result.success:
        raise RuntimeError(
            f"Recipe failed; partial results saved in {args.output}: {result.errors}"
        )
    print(f"Wrote {args.output}")


def main():
    parser = argparse.ArgumentParser(prog="aspcap", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    abundances = commands.add_parser("abundances", help="fit an individual species")
    abundances.add_argument("config", type=Path)
    abundances.add_argument("--ferre", help="FERRE executable name or path")
    abundances.add_argument("--output", type=Path, default=Path("abundance.csv"))
    abundances.add_argument(
        "--pixels", type=Path, help="save normalized data, model, and weights"
    )
    abundances.add_argument("--keep-runs", action="store_true")
    abundances.add_argument(
        "--run-root", type=Path, help="existing directory for FERRE run files"
    )
    run = commands.add_parser("run", help="run an ordered JSON pipeline recipe")
    run.add_argument("config", type=Path)
    run.add_argument(
        "spectrum", type=Path, help="prepared CSV: wavelength,flux,error[,mask]"
    )
    run.add_argument("--output", type=Path, default=Path("pipeline-result.json"))
    run.add_argument("--ferre", help="FERRE executable")
    run.add_argument("--star-id", default="star")
    run.add_argument("--wavelength-frame", choices=("air", "vacuum"), default="vacuum")
    args = parser.parse_args()
    try:
        if args.command == "run":
            if args.output.resolve() in (
                args.config.resolve(),
                args.spectrum.resolve(),
            ):
                raise ValueError(
                    "Output must differ from the recipe and input spectrum"
                )
            _recipe_run(args)
            return
        if args.pixels and args.pixels.resolve() == args.output.resolve():
            raise ValueError("--pixels and --output must use different files")
        _run(args)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        parser.exit(1, f"aspcap: {error}\n")


if __name__ == "__main__":
    main()
