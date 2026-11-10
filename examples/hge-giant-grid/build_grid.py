"""Build the approved six-label, sparse DR17 giant library reproducibly.

Downloads only required PCA records using validated HTTP multipart byte ranges.
Retains both source bases exactly; only logarithmic-microturbulence interpolation
and float32 coefficient storage are introduced. Resume with the same command.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from email.parser import BytesParser
from email.policy import default
from itertools import product
from pathlib import Path
from urllib.request import Request, urlopen
import argparse, hashlib, json, time, os
import numpy as np
import f90nml

HERE = Path(__file__).resolve().parent
WORK = (
    Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    / "aspcap/hge-grid-source"
)
PROPOSAL = HERE / "hge-giant-grid-proposal.json"
if not PROPOSAL.exists():
    PROPOSAL = HERE.parent / "giant-grid-3000-5500/proposal.json"
EXE = "ferre"
BASE = "https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/synspec/marcs/giantisotopes/"
STEM = "hge-giants-dr17-3000-5500-apo-a"
A, B = 2.3379795498263993, -0.43433033061467907
NCOMP = 11 * 6 * 8 * 15


class Source:
    def __init__(self, name):
        self.name = name
        self.url = (
            BASE + f"sg{name}_200921nlte_lsfa/p_apssg{name}_200921nlte_lsfa_012_075"
        )
        path = WORK / (name + ".hdr")
        if not path.exists():
            with urlopen(self.url + ".hdr", timeout=120) as r:
                path.write_bytes(r.read())
        parts = path.read_text().split("\n /\n", 4)
        if len(parts) != 5:
            raise ValueError("Expected four source namelists")
        self.headers = [dict(f90nml.reads(p + "\n/")["synth"]) for p in parts[:4]]
        self.h = self.headers[0]
        self.shape = tuple(self.h["n_p"])
        self.ng = self.shape[-2]
        self.nt = self.shape[-1]
        self.nvar = self.h["npix"] // len(self.h["npca"])
        self.npix = sum(self.h["npca"])
        self.pca = np.fromstring(parts[4], sep=" ", dtype=np.float64).reshape(
            self.nvar + 2, self.npix
        )
        self.size = int(np.prod(self.shape)) * self.h["npix"] * 4
        self.recordbytes = self.h["npix"] * 4
        self.lower = np.array(self.h["llimits"])
        self.steps = np.array(self.h["steps"])
        self.axes = [
            lo + step * np.arange(n)
            for lo, step, n in zip(self.lower, self.steps, self.shape)
        ]

    def mapping(self, tg):
        records = []
        for t, g in tg:
            gi = int(round((g - self.lower[-2]) / self.steps[-2]))
            ti = int(round((t - self.lower[-1]) / self.steps[-1]))
            if not (0 <= gi < self.ng and 0 <= ti < self.nt):
                raise ValueError("Node outside source")
            z = np.log10(A + B * g)
            x = (z - self.lower[0]) / self.steps[0]
            low = int(np.floor(x))
            fraction = x - low
            if not 0 <= low < self.shape[0] - 1:
                raise ValueError("Microturbulence outside source")
            records.append((low, gi * self.nt + ti, float(fraction)))
        return records

    def fetch(self, record_numbers):
        unique = sorted(set(map(int, record_numbers)))
        # Merge contiguous source records, reducing HTTP range headers.
        intervals = []
        for record in unique:
            if intervals and record == intervals[-1][1] + 1:
                intervals[-1][1] = record
            else:
                intervals.append([record, record])
        result = {}
        for offset in range(0, len(intervals), 96):
            group = intervals[offset : offset + 96]
            ranges = [
                (a * self.recordbytes, (b + 1) * self.recordbytes - 1) for a, b in group
            ]
            request = Request(
                self.url + ".unf",
                headers={
                    "Range": "bytes=" + ",".join(f"{a}-{b}" for a, b in ranges),
                    "Accept-Encoding": "identity",
                },
            )
            for attempt in range(5):
                try:
                    with urlopen(request, timeout=120) as response:
                        if response.status != 206:
                            raise ValueError("Archive ignored byte ranges")
                        content_type = response.headers.get("Content-Type", "")
                        data = response.read(
                            sum(b - a + 1 for a, b in ranges)
                            + len(ranges) * 1024
                            + 4096
                        )
                        if content_type.startswith("multipart/"):
                            message = BytesParser(policy=default).parsebytes(
                                (
                                    "Content-Type: "
                                    + content_type
                                    + "\r\nMIME-Version: 1.0\r\n\r\n"
                                ).encode()
                                + data
                            )
                            blocks = [
                                (part["Content-Range"], part.get_payload(decode=True))
                                for part in message.iter_parts()
                            ]
                        else:
                            blocks = [(response.headers.get("Content-Range"), data)]
                    found = set()
                    for content_range, raw in blocks:
                        unit, span = content_range.split()
                        ab, total = span.split("/")
                        a, b = map(int, ab.split("-"))
                        if (
                            unit != "bytes"
                            or int(total) != self.size
                            or (a, b) not in ranges
                            or len(raw) != b - a + 1
                        ):
                            raise ValueError("Incorrect range response")
                        values = np.frombuffer(raw, dtype="<f4").reshape(
                            -1, self.h["npix"]
                        )
                        for j, row in enumerate(values):
                            result[a // self.recordbytes + j] = row.copy()
                        found.add((a, b))
                    if found != set(ranges):
                        raise ValueError("Incomplete multipart response")
                    break
                except Exception:
                    if attempt == 4:
                        raise
                    time.sleep(min(2**attempt, 8))
        if set(result) != set(unique):
            raise ValueError("Missing source records")
        return result

    def coefficients(self, compositions, tg):
        mapping = self.mapping(tg)
        records = []
        stride = NCOMP * self.ng * self.nt
        for comp in compositions:
            for low, gt, f in mapping:
                records.extend(
                    [
                        low * stride + comp * self.ng * self.nt + gt,
                        (low + 1) * stride + comp * self.ng * self.nt + gt,
                    ]
                )
        source = self.fetch(records)
        result = np.empty(
            (len(compositions), len(tg), len(self.h["npca"]), self.nvar), dtype="<f4"
        )
        for i, comp in enumerate(compositions):
            for j, (low, gt, f) in enumerate(mapping):
                first = low * stride + comp * self.ng * self.nt + gt
                result[i, j] = (
                    (1 - f) * source[first].astype(float)
                    + f * source[first + stride].astype(float)
                ).reshape(-1, self.nvar)
        return result

    def flux(self, coeff):
        result = []
        start = 0
        for b, n in enumerate(self.h["npca"]):
            result.append(
                self.pca[0, start : start + n]
                + coeff[b] @ self.pca[2:, start : start + n]
            )
            start += n
        return np.concatenate(result)


def sources():
    WORK.mkdir(parents=True, exist_ok=True)
    return Source("M"), Source("GK")


def check_overlap(m, g):
    # Exact common source nodes; identical C,N,alpha,MH and adopted micro relation.
    points = [(3500, 0), (3500, 1), (4000, 0), (4000, 1), (4000, 2)]
    compositions = [
        np.ravel_multi_index((6, 1, 3, index), (11, 6, 8, 15)) for index in (6, 10)
    ]
    cm = m.coefficients(compositions, points)
    cg = g.coefficients(compositions, points)
    rows = []
    for i, comp in enumerate(compositions):
        for j, (t, gravity) in enumerate(points):
            diff = m.flux(cm[i, j]) - g.flux(cg[i, j])
            rows.append(
                dict(
                    teff=t,
                    logg=gravity,
                    metals=float(m.axes[4][(6, 10)[i]]),
                    rms=float(np.sqrt(np.mean(diff**2))),
                    p99_absolute=float(np.percentile(np.abs(diff), 99)),
                    max_absolute=float(np.max(np.abs(diff))),
                )
            )
    (HERE / "overlap-check.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(
        "Overlap RMS range:",
        min(r["rms"] for r in rows),
        max(r["rms"] for r in rows),
        flush=True,
    )
    return rows


def header(m, g, labels, directory=None):
    directory = HERE if directory is None else directory
    h = dict(
        id=STEM,
        multi=3,
        npca=m.h["npca"],
        npix=12 * 151,
        n_of_dim=6,
        label=labels,
        nodefile=STEM + ".nodes",
    )
    with (directory / (STEM + ".hdr")).open("w") as stream:
        stream.write(str(f90nml.Namelist(synth=h)) + "\n")
        for segment in m.headers[1:]:
            selected = {
                k: v
                for k, v in segment.items()
                if k in ("npix", "wave", "logw", "vacuum", "resolution")
            }
            selected.update(n_of_dim=6, label=labels)
            stream.write(str(f90nml.Namelist(synth=selected)) + "\n")
        np.savetxt(stream, m.pca[0:1], fmt="%.10e")
        np.savetxt(stream, np.ones((1, m.npix)), fmt="%.10e")
        np.savetxt(
            stream,
            np.concatenate(
                [m.pca[2:], g.pca[2:], (g.pca[0] - m.pca[0])[None, :]], axis=0
            ),
            fmt="%.10e",
        )


def mini():
    m, g = sources()
    tg = [tuple(p) for p in json.loads(PROPOSAL.read_text())["teff_logg_nodes"]]
    mi = [i for i, (t, _) in enumerate(tg) if t <= 4000]
    gi = [i for i, (t, _) in enumerate(tg) if t > 4000]
    comp = int(np.ravel_multi_index((6, 1, 3, 8), (11, 6, 8, 15)))
    mc = m.coefficients([comp], [tg[i] for i in mi])
    gc = g.coefficients([comp], [tg[i] for i in gi])
    directory = WORK / "mini"
    directory.mkdir(exist_ok=True)
    header(m, g, m.h["label"][1:], directory)
    block = np.zeros((1, len(tg), 12, 151), dtype="<f4")
    block[:, mi, :, :75] = mc
    block[:, gi, :, 75:150] = gc
    block[:, gi, :, 150] = 1
    (directory / (STEM + ".unf")).write_bytes(block.tobytes())
    combination = [m.axes[i][j] for i, j in zip(range(1, 5), (6, 1, 3, 8))]
    (directory / (STEM + ".nodes")).write_text(
        "".join(
            " ".join(map(str, (*combination, gravity, t))) + "\n" for t, gravity in tg
        )
    )
    from aspcap import Grid

    grid = Grid(directory / (STEM + ".hdr"), executable=EXE, f_access=1)
    errors = []
    for source, coef, indices in [(m, mc, mi), (g, gc, gi)]:
        for j in [0, len(indices) // 2, len(indices) - 1]:
            t, gravity = tg[indices[j]]
            model = grid(dict(zip(grid.labels, (*combination, gravity, t))))
            expected = source.flux(coef[0, j])
            errors.append(float(np.max(np.abs(model.flux - expected))))
    if max(errors) > 3e-5:
        raise ValueError("Union PCA reconstruction failed: " + str(errors))
    parameters = dict(zip(grid.labels, (*combination, 1.0, 3900.0)))
    spec = grid(parameters)
    spec.err = np.full(grid.npix, 0.01)
    fit = grid.fit(
        spec, start=dict(parameters, LOGG=0.9, TEFF=3875), free_labels=("LOGG", "TEFF")
    )
    (HERE / "mini-validation.json").write_text(
        json.dumps(
            dict(
                max_flux_errors=errors, fit_parameters=fit.parameters, truth=parameters
            ),
            indent=2,
        )
        + "\n"
    )
    print("MINI VALIDATION", errors, fit.parameters, flush=True)


def build(workers=6):
    m, g = sources()
    if m.h["npca"] != g.h["npca"] or m.nvar != 75 or g.nvar != 75:
        raise ValueError("Incompatible source PCA geometry")
    for first, second in zip(m.headers[1:], g.headers[1:]):
        for key in ("npix", "wave", "logw", "vacuum", "resolution"):
            if first.get(key) != second.get(key):
                raise ValueError("Source wavelength segments or resolution differ")
    for i in range(1, 5):
        if not np.array_equal(m.axes[i], g.axes[i]):
            raise ValueError("Abundance axes differ")
    overlap = check_overlap(m, g)
    if max(r["rms"] for r in overlap) > 0.02:
        raise ValueError("Large M/GK overlap mismatch; review before joining")
    proposal = json.loads(PROPOSAL.read_text())
    tg = [tuple(p) for p in proposal["teff_logg_nodes"]]
    mi = [i for i, (t, _) in enumerate(tg) if t <= 4000]
    gi = [i for i, (t, _) in enumerate(tg) if t > 4000]
    signature = dict(
        proposal_sha256=hashlib.sha256(PROPOSAL.read_bytes()).hexdigest(),
        source_header_sha256={
            source.name: hashlib.sha256(
                (WORK / (source.name + ".hdr")).read_bytes()
            ).hexdigest()
            for source in (m, g)
        },
        intercept=A,
        slope=B,
        components_per_block=151,
    )
    statepath = HERE / "build-state.json"
    if statepath.exists():
        if json.loads(statepath.read_text()) != signature:
            raise ValueError("Build settings changed; use a new output directory")
    else:
        if (HERE / "completed-chunks").exists() and any(
            (HERE / "completed-chunks").iterdir()
        ):
            raise ValueError("Existing checkpoints have no configuration fingerprint")
        statepath.write_text(json.dumps(signature, indent=2) + "\n")
    labels = m.h["label"][1:]
    header(m, g, labels)
    combinations = list(product(*m.axes[1:5]))
    nodepath = HERE / (STEM + ".nodes")
    with nodepath.open("w") as stream:
        stream.write("# " + " | ".join(labels) + "\n")
        for combination in combinations:
            for t, gravity in tg:
                stream.write(
                    " ".join(f"{v:.9g}" for v in (*combination, gravity, t)) + "\n"
                )
    target = HERE / (STEM + ".unf")
    shape = (NCOMP, len(tg), 12, 151)
    if target.exists() and target.stat().st_size != int(np.prod(shape)) * 4:
        raise ValueError("Existing coefficient file has the wrong length")
    output = np.memmap(
        target, dtype="<f4", mode="r+" if target.exists() else "w+", shape=shape
    )
    checkpoints = HERE / "completed-chunks"
    checkpoints.mkdir(exist_ok=True)
    chunk = 8
    starts = list(range(0, NCOMP, chunk))

    def process(start):
        marker = checkpoints / f"{start:05d}.json"
        if marker.exists():
            stored = json.loads(marker.read_text())
            existing = np.asarray(output[start : start + stored["count"]])
            if hashlib.sha256(existing.tobytes()).hexdigest() != stored["sha256"]:
                raise ValueError("Checkpoint data failed checksum")
            return start, False
        comps = list(range(start, min(start + chunk, NCOMP)))
        mc = m.coefficients(comps, [tg[i] for i in mi])
        gc = g.coefficients(comps, [tg[i] for i in gi])
        if not np.isfinite(mc).all() or not np.isfinite(gc).all():
            raise ValueError("Nonfinite source coefficients")
        block = np.zeros((len(comps), len(tg), 12, 151), dtype="<f4")
        block[:, mi, :, :75] = mc
        block[:, gi, :, 75:150] = gc
        block[:, gi, :, 150] = 1
        output[start : start + len(comps)] = block
        output.flush()
        marker.write_text(
            json.dumps(
                dict(
                    start=start,
                    count=len(comps),
                    sha256=hashlib.sha256(block.tobytes()).hexdigest(),
                )
            )
            + "\n"
        )
        return start, True

    complete = 0
    began = time.time()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(process, start) for start in starts]):
            start, fresh = future.result()
            complete += 1
            if complete % 10 == 0 or complete == len(starts):
                print(
                    f"{complete}/{len(starts)} chunks complete; {time.time()-began:.0f}s elapsed",
                    flush=True,
                )
    output.flush()
    del output
    metadata = dict(
        status="built_pending_validation",
        source_headers={s.name: str(WORK / (s.name + ".hdr")) for s in (m, g)},
        source_urls={s.name: s.url for s in (m, g)},
        labels=labels,
        spectra=NCOMP * len(tg),
        coefficient_count=1812,
        teff_logg_nodes=tg,
        vmicro=dict(intercept=A, slope=B, gravity="raw LOGG_SPEC", units="km/s"),
        format="sparse NODEFILE, PCA union basis, PCACHI=0 and PCAPROJECT=0 required",
        basis_note="75 M basis vectors, 75 GK vectors, one GK-minus-M mean vector per block; no additional PCA truncation",
        source_micro_interpolation="linear in LOG10VDOP coefficients",
        source_switch="M <=4000 K; GK >=4250 K",
        overlap=overlap,
    )
    (HERE / "build-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("BUILD COMPLETE", target, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--proposal", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--ferre", default="ferre")
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--mini", action="store_true")
    args = parser.parse_args()
    if args.output:
        HERE = args.output.resolve()
        HERE.mkdir(parents=True, exist_ok=True)
    if args.proposal:
        PROPOSAL = args.proposal.resolve()
    if args.cache:
        WORK = args.cache.resolve()
    EXE = args.ferre
    if args.preflight:
        check_overlap(*sources())
    elif args.mini:
        mini()
    else:
        build(args.workers)
