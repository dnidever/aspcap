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
    """Run one input or an ordered batch, preserving read/fit failures."""
    from .recipe import PipelineRecipe
    from .readers import read_spectrum

    recipe = PipelineRecipe.read(args.config, executable=args.ferre)
    if args.input_list and args.spectrum:
        raise ValueError("Use a spectrum argument or --input-list, not both")
    list_file = args.input_list
    if args.spectrum and args.spectrum.startswith("@"):
        list_file = Path(args.spectrum[1:])
    if list_file is not None:
        list_file = list_file.resolve()
        inputs = [
            list_file.parent / line.strip()
            for line in list_file.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        batch = True
    elif args.spectrum:
        inputs = [
            Path(name.strip()) for name in args.spectrum.split(",") if name.strip()
        ]
        batch = len(inputs) > 1
    else:
        raise ValueError("Supply FITS/CSV paths or --input-list")
    if not inputs:
        raise ValueError("Input list is empty")
    protected = [args.config.resolve(), *(p.resolve() for p in inputs)]
    if list_file is not None:
        protected.append(list_file)
    if args.output.resolve() in protected:
        raise ValueError("Output must differ from recipe, list, and input spectra")
    if args.stop_after is not None and not 1 <= args.stop_after <= len(
        recipe.config["steps"]
    ):
        raise ValueError("--stop-after must identify an existing one-based step")
    grid = recipe._grid(recipe.config["steps"][0]["grid"])
    results = []
    for index, path in enumerate(inputs):
        entry = dict(index=index, source=str(path), success=False)
        try:
            spectrum = read_spectrum(
                path,
                grid=grid,
                row=args.fits_row,
                bad_bits=args.bad_bits,
                wavevac=args.wavelength_frame == "vacuum",
            )
            identity = args.star_id or spectrum.metadata.get("star_id") or path.stem
            result = recipe.fit(spectrum, star_id=identity, stop_after=args.stop_after)
            entry.update(
                star_id=str(identity),
                success=result.success,
                parameters=result.parameters,
                abundances=result.abundances,
                errors=result.errors,
                history=result.history,
                config=result.config,
            )
        except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
            entry["errors"] = {"input": str(error)}
        results.append(entry)
    success = all(item["success"] for item in results)
    payload = (
        dict(success=success, count=len(results), results=results)
        if batch
        else results[0]
    )

    def serial(value):
        if hasattr(value, "tolist"):
            return value.tolist()
        raise TypeError(f"Cannot serialize {type(value)}")

    args.output.write_text(
        json.dumps(payload, indent=2, default=serial, allow_nan=False) + "\n"
    )
    if not success:
        raise RuntimeError(
            f"{sum(not item['success'] for item in results)} of {len(results)} inputs failed; results saved in {args.output}"
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
        "spectrum", nargs="?", help="FITS/CSV path, comma-separated paths, or @list.txt"
    )
    run.add_argument("--output", type=Path, default=Path("pipeline-result.json"))
    run.add_argument("--ferre", help="FERRE executable")
    run.add_argument("--star-id", help="override input identities")
    run.add_argument(
        "--input-list", type=Path, help="one path per line, relative to the list file"
    )
    run.add_argument(
        "--fits-row",
        type=int,
        default=0,
        help="zero-based apStar row, default combined row 0",
    )
    run.add_argument(
        "--bad-bits",
        type=lambda value: int(value, 0),
        default=255,
        help="FITS bad pixel bits, default 255",
    )
    run.add_argument(
        "--stop-after", type=int, help="run through this one-based recipe step"
    )
    run.add_argument("--wavelength-frame", choices=("air", "vacuum"), default="vacuum")
    args = parser.parse_args()
    try:
        if args.command == "run":
            _recipe_run(args)
            return
        if args.pixels and args.pixels.resolve() == args.output.resolve():
            raise ValueError("--pixels and --output must use different files")
        _run(args)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        parser.exit(1, f"aspcap: {error}\n")


if __name__ == "__main__":
    main()
