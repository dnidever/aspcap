"""Validate checkpoint integrity, reconstruction, interpolation and a FERRE fit."""

from pathlib import Path
import importlib.util, json, hashlib, time
from unittest.mock import patch
import f90nml
import subprocess
import numpy as np
from aspcap import Grid, Spec1D

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("builder", HERE / "build_grid.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)
EXE = "ferre"
PROPOSAL = HERE / "hge-giant-grid-proposal.json"
if not PROPOSAL.exists():
    PROPOSAL = HERE.parent / "giant-grid-3000-5500/proposal.json"


def continuous(source, comp, temperature, gravity):
    """Independent source-grid linear interpolation at analytic vmicro(g)."""

    def bracket(axis, value):
        right = int(np.searchsorted(axis, value, side="right"))
        left = max(0, min(right - 1, len(axis) - 2))
        f = (value - axis[left]) / (axis[left + 1] - axis[left])
        return [(left, 1 - f), (left + 1, f)]

    micro = np.log10(builder.A + builder.B * gravity)
    weights = {}
    for vi, vw in bracket(source.axes[0], micro):
        for gi, gw in bracket(source.axes[-2], gravity):
            for ti, tw in bracket(source.axes[-1], temperature):
                record = (
                    (vi * builder.NCOMP + comp) * source.ng * source.nt
                    + gi * source.nt
                    + ti
                )
                weights[record] = vw * gw * tw
    records = source.fetch(weights)
    coefficient = sum(
        weight * records[record].astype(float) for record, weight in weights.items()
    ).reshape(12, 75)
    return source.flux(coefficient)


def run(mini=False):
    m, g = builder.sources()
    folder = builder.WORK / "mini" if mini else HERE
    grid = Grid(folder / (builder.STEM + ".hdr"), executable=EXE, f_access=1)
    proposal = json.loads(PROPOSAL.read_text())
    tg = [tuple(pair) for pair in proposal["teff_logg_nodes"]]
    if not mini:
        array = np.memmap(
            folder / (builder.STEM + ".unf"),
            dtype="<f4",
            mode="r",
            shape=(builder.NCOMP, len(tg), 12, 151),
        )
        markers = list((folder / "completed-chunks").glob("*.json"))
        if len(markers) != 990:
            raise ValueError("Incomplete build checkpoints")
        for path in markers:
            record = json.loads(path.read_text())
            start, count = record["start"], record["count"]
            if not np.isfinite(array[start : start + count]).all():
                raise ValueError("Nonfinite coefficients in checkpoint: " + str(path))
            if (
                hashlib.sha256(array[start : start + count].tobytes()).hexdigest()
                != record["sha256"]
            ):
                raise ValueError("Checkpoint checksum mismatch: " + str(path))
        del array
    compositions = (
        [int(np.ravel_multi_index((6, 1, 3, 8), (11, 6, 8, 15)))]
        if mini
        else [0, 1234, 4567, 7919]
    )
    requests = []
    expected = []
    descriptions = []
    for comp in compositions:
        combination = [
            m.axes[i][j]
            for i, j in zip(range(1, 5), np.unravel_index(comp, (11, 6, 8, 15)))
        ]
        for source, point in [
            (m, (3300, 0)),
            (m, (4000, 2.5)),
            (g, (4500, 1)),
            (g, (5500, 3)),
        ]:
            t, gravity = point
            coefficient = source.coefficients([comp], [point])[0, 0]
            requests.append(dict(zip(grid.labels, (*combination, gravity, t))))
            expected.append(source.flux(coefficient))
            descriptions.append(
                dict(source=source.name, composition=comp, teff=t, logg=gravity)
            )
    models = grid(requests)
    exact = []
    for model, reference, description in zip(models, expected, descriptions):
        difference = model.flux - reference
        exact.append(
            dict(
                description,
                max_absolute=float(np.max(np.abs(difference))),
                rms=float(np.sqrt(np.mean(difference**2))),
            )
        )
    if max(r["max_absolute"] for r in exact) > 3e-5:
        raise ValueError("Stored node reconstruction mismatch")
    comp = compositions[len(compositions) // 2]
    combination = [
        m.axes[i][j]
        for i, j in zip(range(1, 5), np.unravel_index(comp, (11, 6, 8, 15)))
    ]
    points = [(m, 3550, 0.75), (m, 3850, 1.25), (g, 4625, 1.75), (g, 5375, 2.75)]
    requests = [
        dict(zip(grid.labels, (*combination, gravity, t)))
        for source, t, gravity in points
    ]
    models = grid(requests)
    interpolation = []
    for model, (source, t, gravity) in zip(models, points):
        reference = continuous(source, comp, t, gravity)
        difference = model.flux - reference
        interpolation.append(
            dict(
                source=source.name,
                teff=t,
                logg=gravity,
                rms=float(np.sqrt(np.mean(difference**2))),
                max_absolute=float(np.max(np.abs(difference))),
            )
        )
    truth = dict(zip(grid.labels, (*combination, 1.25, 4500.0)))
    spectrum = grid(truth)
    spectrum.err = np.full(grid.npix, 0.01)
    start = dict(truth, LOGG=1.15, TEFF=4475.0)
    fit = grid.fit(spectrum, start=start, free_labels=("LOGG", "TEFF"))
    if (
        abs(fit.parameters["LOGG"] - truth["LOGG"]) > 0.02
        or abs(fit.parameters["TEFF"] - truth["TEFF"]) > 5
    ):
        raise ValueError("Synthetic fit recovery failed")
    coverage = [
        dict(
            teff=t,
            logg=gravity,
            covered=grid.contains(dict(zip(grid.labels, (*combination, gravity, t)))),
        )
        for t, gravity in [(3300, 0), (3300, 3), (5000, 0), (5000, 2), (5500, 3.5)]
    ]
    if [c["covered"] for c in coverage] != [True, False, False, True, False]:
        raise ValueError("Coverage/exclusion mismatch")
    global_test = None
    if not mini:
        global_truth = {
            "C": 0.0,
            "N": 0.0,
            "O Mg Si S Ca Ti": 0.0,
            "METALS": -0.5,
            "LOGG": 1.25,
            "TEFF": 4500.0,
        }
        coefficient_file = folder / (builder.STEM + ".unf")
        identity = dict(
            build_state_sha256=hashlib.sha256(
                (folder / "build-state.json").read_bytes()
            ).hexdigest(),
            coefficient_file_size=coefficient_file.stat().st_size,
            coefficient_file_mtime_ns=coefficient_file.stat().st_mtime_ns,
        )
        cache = folder / "six-label-validation.json"
        cached = json.loads(cache.read_text()) if cache.exists() else {}
        if all(cached.get(k) == v for k, v in identity.items()):
            global_test = cached
        else:
            comp = int(np.ravel_multi_index((6, 1, 3, 8), (11, 6, 8, 15)))
            coefficients = g.coefficients([comp], [(4500, 1), (4500, 1.5)])[0]
            flux = 0.5 * (g.flux(coefficients[0]) + g.flux(coefficients[1]))
            observation = Spec1D(
                flux, wave=grid.wavelengths, err=np.full(grid.npix, 0.01)
            )
            estimate = dict(
                global_truth, C=0.02, N=0.04, METALS=-0.47, LOGG=1.2, TEFF=4475.0
            )
            controls = dict(stopcr=1e-9, simp=1e-9, scope=0.05)
            original_run = subprocess.run

            def precise(*args, **kwargs):
                path = Path(kwargs["cwd"]) / "input.nml"
                namelist = f90nml.read(path)
                namelist["lista"].update(controls)
                namelist.write(path, force=True)
                return original_run(*args, **kwargs)

            with patch("aspcap.ferre.subprocess.run", side_effect=precise):
                global_fit = grid.fit(observation, start=estimate)
            global_test = dict(
                truth=global_truth,
                recovered=global_fit.parameters,
                absolute_errors={
                    label: abs(global_fit.parameters[label] - value)
                    for label, value in global_truth.items()
                },
                reduced_chi2=global_fit.reduced_chi2,
                validation_controls=controls,
                **identity,
            )
            cache.write_text(json.dumps(global_test, indent=2) + "\n")
        errors = global_test["absolute_errors"]
        if errors["TEFF"] > 10 or any(
            value > 0.03 for label, value in errors.items() if label != "TEFF"
        ):
            raise ValueError("Six-label synthetic recovery failed: " + str(errors))
    result = dict(
        status="validated",
        mini=mini,
        spectra=grid.nspectra,
        labels=grid.labels,
        wavelength_pixels=grid.npix,
        six_label_synthetic_fit=global_test,
        source_node_checks=exact,
        off_node_analytic_relation_comparison=interpolation,
        synthetic_fit=dict(
            truth=truth, recovered=fit.parameters, reduced_chi2=fit.reduced_chi2
        ),
        coverage_checks=coverage,
        comment="Node checks compare original source coefficients after micro interpolation. Off-node comparisons measure interpolation of the collapsed relation against original 7D linear interpolation with analytic vmicro(g).",
    )
    files = {}
    for extension in ("hdr", "nodes", "unf"):
        path = folder / (builder.STEM + "." + extension)
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        files[path.name] = dict(bytes=path.stat().st_size, sha256=digest.hexdigest())
    result["files"] = files
    (folder / "SHA256SUMS").write_text(
        "".join(value["sha256"] + "  " + name + "\n" for name, value in files.items())
    )
    metadata_path = HERE / "build-metadata.json"
    if not mini and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        metadata.update(
            status="validated", files=files, validation_file="validation.json"
        )
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    (HERE / ("mini-full-validation.json" if mini else "validation.json")).write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(
        json.dumps(
            dict(
                status="validated",
                mini=mini,
                spectra=grid.nspectra,
                max_node_error=max(r["max_absolute"] for r in exact),
                interpolation=interpolation,
                fit=fit.parameters,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("output", nargs="?", type=Path)
    p.add_argument("--cache", type=Path)
    p.add_argument("--proposal", type=Path)
    p.add_argument("--ferre", default="ferre")
    p.add_argument("--mini", action="store_true")
    args = p.parse_args()
    if args.output:
        HERE = args.output.resolve()
    if args.cache:
        builder.WORK = args.cache.resolve()
    if args.proposal:
        PROPOSAL = args.proposal.resolve()
    EXE = args.ferre
    run(args.mini)
