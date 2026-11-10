"""HTTP bundle integrity, resuming, catalog filters, and cache reuse."""

from contextlib import contextmanager
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
import shutil
from pathlib import Path
import struct
import tempfile
from threading import Thread
import unittest

from aspcap import Grid, grids

ROOT = Path(__file__).resolve().parents[1]


@contextmanager
def server():
    text = (ROOT / "examples/toy-grid.dat").read_text()
    header, body = text.split("/\n", 1)
    values = list(map(float, body.split()))
    content = {
        "grid.hdr": (header + "/\n").encode(),
        "grid.unf": struct.pack(f"={len(values)}f", *values),
    }
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_HEAD(self):
            data = content.get(self.path[1:])
            if data is None:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("ETag", '"test-v1"')
            self.end_headers()

        def do_GET(self):
            data = content[self.path[1:]]
            value = self.headers.get("Range")
            requests.append((self.path, value))
            first = int(value[6:].split("-")[0]) if value else 0
            end = value[6:].split("-")[1] if value else ""
            last = int(end) if end else len(data) - 1
            self.send_response(206 if value else 200)
            if value:
                self.send_header("Content-Range", f"bytes {first}-{last}/{len(data)}")
            self.send_header("Content-Length", str(last + 1 - first))
            self.end_headers()
            self.wfile.write(data[first : last + 1])

    http = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=http.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{http.server_port}/", content, requests
    finally:
        http.shutdown()
        http.server_close()
        thread.join()


class DownloadTests(unittest.TestCase):
    def catalog(self, root, url, content):
        entry = dict(
            id="test-grid",
            family="M-giant",
            dimensions=2,
            release="TEST",
            synthesis="test",
            instrument="test",
            lsf="a",
            header="grid.hdr",
            files=[
                dict(
                    name=name,
                    url=url + name,
                    size=len(data),
                    sha256=hashlib.sha256(data).hexdigest(),
                )
                for name, data in content.items()
            ],
        )
        path = root / "catalog.json"
        path.write_text(json.dumps(dict(version=1, entries=[entry])))
        return path

    def test_download_cache_and_corruption(self):
        with tempfile.TemporaryDirectory() as root, server() as (
            url,
            content,
            requests,
        ):
            root = Path(root)
            catalog = self.catalog(root, url, content)
            self.assertEqual(len(grids.available(family="m-GIANT", catalog=catalog)), 1)
            self.assertEqual(grids.available(family="GK-giant", catalog=catalog), [])
            self.assertEqual(
                grids.info("test-grid", online=True, catalog=catalog)["size"],
                sum(map(len, content.values())),
            )
            header = grids.download("test-grid", root / "cache", catalog=catalog)
            self.assertEqual(Grid(header).npix, 8)
            executable = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
            if executable:
                downloaded = Grid(header, executable=executable)(
                    {"ALPHA": 0.2, "METALS": -0.25}
                )
                original = Grid(ROOT / "examples/toy-grid.dat", executable=executable)(
                    {"ALPHA": 0.2, "METALS": -0.25}
                )
                for actual, expected in zip(downloaded.flux, original.flux):
                    self.assertAlmostEqual(actual, expected, places=6)
            before = len(requests)
            self.assertEqual(
                grids.download("test-grid", root / "cache", catalog=catalog), header
            )
            self.assertEqual(len(requests), before)
            header.with_suffix(".unf").write_bytes(b"corrupted")
            with self.assertRaises(ValueError):
                grids.download("test-grid", root / "cache", catalog=catalog)

    def test_resume_retains_binary_alignment(self):
        with tempfile.TemporaryDirectory() as root, server() as (
            url,
            content,
            requests,
        ):
            root = Path(root)
            catalog = self.catalog(root, url, content)
            stage = root / "cache/test-grid.partial"
            stage.mkdir(parents=True)
            (stage / "grid.unf.part").write_bytes(content["grid.unf"][:100])
            (stage / "resume.json").write_text(
                json.dumps(
                    {
                        "grid.unf": dict(
                            url=url + "grid.unf",
                            etag='"test-v1"',
                            size=len(content["grid.unf"]),
                        )
                    }
                )
            )
            header = grids.download("test-grid", root / "cache", catalog=catalog)
            self.assertEqual(
                header.with_suffix(".unf").read_bytes(), content["grid.unf"]
            )
            self.assertIn(("/grid.unf", "bytes=100-"), requests)
            self.assertFalse(stage.exists())

    def test_interrupted_transfer_and_checksum_failure(self):
        with tempfile.TemporaryDirectory() as root, server() as (
            url,
            content,
            requests,
        ):
            root = Path(root)
            catalog = self.catalog(root, url, content)

            def stop(name, done, total):
                if done:
                    raise RuntimeError("interrupted")

            with self.assertRaises(RuntimeError):
                grids.download(
                    "test-grid", root / "cache", catalog=catalog, progress=stop
                )
            self.assertFalse((root / "cache/test-grid").exists())
            self.assertFalse((root / "cache/test-grid.lock").exists())
            self.assertTrue(
                grids.download("test-grid", root / "cache", catalog=catalog).exists()
            )
            doc = json.loads(catalog.read_text())
            doc["entries"][0]["files"][1]["sha256"] = "0" * 64
            catalog.write_text(json.dumps(doc))
            with self.assertRaisesRegex(ValueError, "checksum"):
                grids.download("test-grid", root / "other", catalog=catalog)

    def test_catalog_and_selection_validation(self):
        entries = grids.available()
        self.assertGreaterEqual(len(entries), 40)
        self.assertEqual(len({e["id"] for e in entries}), len(entries))
        self.assertTrue(grids.available(family="M-giant", synthesis="synspec-nlte"))
        with self.assertRaises(ValueError):
            grids.info("not-a-grid")
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "bad.json"
            path.write_text(
                json.dumps(
                    dict(
                        version=1,
                        entries=[
                            dict(
                                id="test",
                                header="../evil",
                                files=[
                                    dict(name="../evil", url="https://example.org/file")
                                ],
                            )
                        ],
                    )
                )
            )
            with self.assertRaises(ValueError):
                grids.available(catalog=path)
