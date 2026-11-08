"""Command-line interface for a single HGE abundance fit."""

import argparse
import csv
import json
from pathlib import Path

from . import Grid, fit_abundance


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
    result, pixels = fit_abundance(
        grid,
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
    args = parser.parse_args()
    try:
        if args.pixels and args.pixels.resolve() == args.output.resolve():
            raise ValueError("--pixels and --output must use different files")
        _run(args)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        parser.exit(1, f"aspcap: {error}\n")


if __name__ == "__main__":
    main()
