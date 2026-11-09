"""Browse a verified SDSS grid catalog and download complete FERRE bundles."""

from copy import deepcopy
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import re
import shutil
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


def _entries(catalog=None):
    path = (
        files("aspcap").joinpath("data/grid-catalog.json")
        if catalog is None
        else Path(catalog)
    )
    with path.open() as stream:
        data = json.load(stream)
    if data.get("version") != 1:
        raise ValueError("Unsupported grid catalog version")
    entries = data["entries"]
    if len({e["id"] for e in entries}) != len(entries):
        raise ValueError("Catalog identifiers must be unique")
    for entry in entries:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", entry["id"]):
            raise ValueError("Unsafe catalog identifier")
        names = [f["name"] for f in entry["files"]]
        if not names or len(set(names)) != len(names) or entry["header"] not in names:
            raise ValueError("Catalog bundle must have unique files and a header")
        for file in entry["files"]:
            if (
                "/" in file["name"]
                or "\\" in file["name"]
                or file["name"] in ("", ".", "..")
            ):
                raise ValueError("Catalog file names must be simple local names")
            if urlsplit(file["url"]).scheme not in ("http", "https"):
                raise ValueError("Only HTTP(S) grid downloads are supported")
    return entries


def available(
    *,
    family=None,
    synthesis=None,
    instrument=None,
    lsf=None,
    release=None,
    catalog=None,
):
    """Return catalog entries matching case-insensitive exact filters.

    No network access occurs. Each entry includes id, description, family,
    dimensions, published temperature range, synthesis, instrument, LSF group,
    source URLs, and snapshot file/total sizes in bytes. This initial catalog
    covers available DR17 production FERRE bundles and comparison syntheses,
    not every raw high-resolution spectrum or every survey. Use info(online=True)
    to check current availability and sizes. Returned mappings are independent.
    catalog optionally names an alternative version-1 JSON catalog.
    """
    filters = dict(
        family=family,
        synthesis=synthesis,
        instrument=instrument,
        lsf=lsf,
        release=release,
    )
    return deepcopy(
        [
            e
            for e in _entries(catalog)
            if all(
                value is None or str(e[key]).lower() == str(value).lower()
                for key, value in filters.items()
            )
        ]
    )


def _remote(file, timeout):
    with urlopen(
        Request(file["url"], method="HEAD", headers={"Accept-Encoding": "identity"}),
        timeout=timeout,
    ) as response:
        length = response.headers.get("Content-Length")
        return dict(
            file,
            size=None if length is None else int(length),
            etag=response.headers.get("ETag"),
        )


def info(identifier, *, online=False, catalog=None, timeout=30):
    """Describe one selection; online=True checks HEAD metadata, not payloads.

    Published catalog sizes are snapshots. Online errors propagate explicitly.
    A None total means at least one remote file does not publish its length.
    """
    matches = [e for e in available(catalog=catalog) if e["id"] == identifier]
    if not matches:
        raise ValueError(f"Unknown grid {identifier!r}; use grids.available()")
    entry = matches[0]
    if online:
        entry["files"] = [_remote(f, timeout) for f in entry["files"]]
    sizes = [f.get("size") for f in entry["files"]]
    entry["size"] = None if None in sizes else sum(sizes)
    return entry


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _receipt(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def _verify(directory, entry):
    receipt = json.loads((directory / "download.json").read_text())
    if receipt.get("id") != entry["id"]:
        raise ValueError("Downloaded bundle identity does not match catalog")
    for file in entry["files"]:
        record = receipt["files"][file["name"]]
        path = directory / file["name"]
        if (
            record["url"] != file["url"]
            or path.stat().st_size != record["size"]
            or _sha256(path) != record["sha256"]
            or (file.get("sha256") is not None and record["sha256"] != file["sha256"])
        ):
            raise ValueError(f"Downloaded file failed integrity check: {path}")
    return directory / entry["header"]


def download(identifier, directory=None, *, progress=None, catalog=None, timeout=60):
    """Download/resume a complete catalog bundle and return its header Path.

    Parameters
    ----------
    identifier : str
        Exact id from available(). The selection is explicit; no other grids
        are downloaded automatically.
    directory : path-like, optional
        Parent directory. Defaults to $XDG_CACHE_HOME/aspcap/grids or
        ~/.cache/aspcap/grids. Each selection has its own subdirectory.
    progress : callable, optional
        Called as progress(filename, bytes_present, total_or_None) while each
        file is transferred. Exceptions interrupt the transfer safely.
    catalog : path-like, optional
        Alternative catalog, useful for additional hosted FERRE libraries.
    timeout : float, optional
        Timeout per network operation, in seconds.

    Notes
    -----
    Existing complete downloads are checked with local SHA256 receipts and
    reused without network access. New downloads check current remote sizes,
    use resumable .partial staging directories and Range/If-Range when an ETag
    is available, and only publish a complete bundle after Grid can read its
    header. A lock prevents concurrent writers to the same selection. Remote
    errors leave partial files for a retry. No archives are extracted.
    Checksums detect subsequent local corruption; without publisher checksums
    they do not independently authenticate the initial server payload.
    Large production libraries can require many GB of storage and RAM.
    """
    entry = info(identifier, catalog=catalog)
    if directory is None:
        directory = (
            Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
            / "aspcap/grids"
        )
    parent = Path(directory).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    final, stage, lock = (
        parent / identifier,
        parent / (identifier + ".partial"),
        parent / (identifier + ".lock"),
    )
    if final.exists():
        return _verify(final, entry)
    try:
        lock.mkdir()
    except FileExistsError:
        raise RuntimeError(
            f"Download already locked: {lock}; remove a stale lock only after confirming no download is running"
        ) from None
    try:
        if final.exists():
            return _verify(final, entry)
        entry = info(identifier, online=True, catalog=catalog, timeout=timeout)
        stage.mkdir(exist_ok=True)
        receipt_path = stage / "download.json"
        receipt = (
            json.loads(receipt_path.read_text())
            if receipt_path.exists()
            else dict(id=identifier, files={})
        )
        if receipt.get("id") != identifier:
            raise ValueError("Partial download identity mismatch")
        state_path = stage / "resume.json"
        state = json.loads(state_path.read_text()) if state_path.exists() else {}
        required = sum(
            (
                max(
                    0,
                    (f.get("size") or 0)
                    - (stage / (f["name"] + ".part")).stat().st_size,
                )
                if (stage / (f["name"] + ".part")).exists()
                else (f.get("size") or 0)
            )
            for f in entry["files"]
            if not (stage / f["name"]).exists()
        )
        if shutil.disk_usage(parent).free < required:
            raise OSError(
                f"Insufficient space: approximately {required} additional bytes required"
            )
        for file in entry["files"]:
            name, expected, etag = file["name"], file["size"], file.get("etag")
            target, part = stage / name, stage / (name + ".part")
            previous = receipt["files"].get(name)
            if (
                target.exists()
                and previous
                and previous["url"] == file["url"]
                and previous.get("etag") == etag
                and target.stat().st_size == expected
                and _sha256(target) == previous["sha256"]
            ):
                continue
            if target.exists():
                target.unlink()
            resume = dict(url=file["url"], etag=etag, size=expected)
            if part.exists() and (state.get(name) != resume or not etag):
                part.unlink()
            state[name] = resume
            _receipt(state_path, state)
            offset = part.stat().st_size if part.exists() else 0
            if expected is not None and offset > expected:
                part.unlink()
                offset = 0
            if offset != expected or not part.exists():
                headers = {"Accept-Encoding": "identity"}
                if offset:
                    headers.update(Range=f"bytes={offset}-", **{"If-Range": etag})
                with urlopen(
                    Request(file["url"], headers=headers), timeout=timeout
                ) as response:
                    status = response.status
                    if status == 206:
                        match = re.fullmatch(
                            r"bytes (\d+)-(\d+)/(\d+)",
                            response.headers.get("Content-Range", ""),
                        )
                        if (
                            not match
                            or int(match[1]) != offset
                            or (expected is not None and int(match[3]) != expected)
                        ):
                            raise ValueError(
                                "Server returned an inconsistent byte range"
                            )
                    elif status == 200:
                        offset = 0
                    else:
                        raise ValueError(f"Unexpected download status {status}")
                    if expected is None:
                        if status == 206:
                            expected = int(match[3])
                        elif response.headers.get("Content-Length") is not None:
                            expected = int(response.headers["Content-Length"])
                    with part.open("ab" if offset else "wb") as stream:
                        if progress:
                            progress(name, offset, expected)
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            stream.write(block)
                            offset += len(block)
                            if progress:
                                progress(name, offset, expected)
            actual = part.stat().st_size
            if expected is not None and actual != expected:
                raise ValueError(
                    f"Truncated download for {name}: {actual} != {expected}"
                )
            digest = _sha256(part)
            # Optional publisher SHA256 entries are checked when provided.
            if file.get("sha256") is not None and digest != file["sha256"]:
                raise ValueError(f"Publisher checksum mismatch for {name}")
            part.replace(target)
            receipt["files"][name] = dict(
                url=file["url"], size=actual, etag=etag, sha256=digest
            )
            _receipt(receipt_path, receipt)
        from .ferre import Grid

        Grid(stage / entry["header"])
        _verify(stage, entry)
        state_path.unlink(missing_ok=True)
        stage.rename(final)
        return final / entry["header"]
    finally:
        lock.rmdir()
