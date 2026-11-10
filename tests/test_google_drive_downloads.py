"""Public Drive confirmation, integrity-pinned resume and one-call grid loading."""

from io import BytesIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch, Mock
from urllib.request import Request, urlopen
import pytest
from aspcap import Grid, grids
from aspcap._google_drive import open_file
from test_grid_downloads import server
import test_grid_downloads as helpers


class Response(BytesIO):
    def __init__(self, data, headers=None, status=200):
        super().__init__(data)
        self.headers = headers or {}
        self.status = status

    def geturl(self):
        return "https://drive.usercontent.google.com/download?id=test"


def test_public_confirmation_preserves_range_without_browser_cookies():
    html = b"""<form id="download-form" action="https://drive.usercontent.google.com/download"><input type="hidden" name="id" value="test"><input type="hidden" name="confirm" value="approved-token"></form>"""
    transport = Mock()
    transport.open.side_effect = [
        Response(html, {"Content-Type": "text/html"}),
        Response(
            b"bytes",
            {
                "Content-Disposition": 'attachment; filename="grid.unf"',
                "Content-Type": "application/octet-stream",
            },
            206,
        ),
    ]
    file = dict(
        drive_file_id="test",
        name="grid.unf",
        url="https://drive.google.com/file/d/test/view",
    )
    with patch("aspcap._google_drive.build_opener", return_value=transport):
        with open_file(file, headers={"Range": "bytes=10-"}) as response:
            assert response.read() == b"bytes"
    assert len(transport.open.call_args_list) == 2
    request = transport.open.call_args_list[-1].args[0]
    assert request.get_header("Range") == "bytes=10-"
    assert "confirm=approved-token" in request.full_url


@pytest.mark.parametrize(
    "page",
    [
        b"<html>Too many users have viewed or downloaded this file.</html>",
        b'<form id="download-form" action="https://evil.example/download"><input type="hidden" name="id" value="test"></form>',
        b'<form id="download-form" action="https://drive.usercontent.google.com/download"><input type="hidden" name="id" value="different"></form>',
    ],
)
def test_quota_and_unexpected_forms_never_become_file_bytes(page):
    transport = Mock()
    transport.open.return_value = Response(page, {"Content-Type": "text/html"})
    with patch("aspcap._google_drive.build_opener", return_value=transport):
        with pytest.raises(RuntimeError):
            open_file(
                dict(
                    drive_file_id="test",
                    name="grid.unf",
                    url="https://drive.google.com/file/d/test/view",
                )
            )
    assert transport.open.call_count == 1


def test_from_download_is_direct_access_and_reuses_cache_offline():
    with TemporaryDirectory() as temporary, server() as (url, content, requests):
        root = Path(temporary)
        catalog = helpers.DownloadTests().catalog(root, url, content)
        grid = Grid.from_download(
            "test-grid", root / "cache", catalog=catalog, executable="custom-ferre"
        )
        assert (
            grid.npix == 8 and grid.f_access == 1 and grid.executable == "custom-ferre"
        )
        with patch(
            "aspcap.grids.urlopen", side_effect=AssertionError("No network allowed")
        ):
            same = Grid.from_download(
                "test-grid", root / "cache", catalog=catalog, f_access=0
            )
        assert same.path == grid.path and same.f_access == 0


@pytest.mark.parametrize("honor_range", [True, False])
def test_drive_pinned_resume_with_no_etag_or_safe_restart(honor_range):
    with TemporaryDirectory() as temporary, server() as (url, content, requests):
        root = Path(temporary)
        catalog = helpers.DownloadTests().catalog(root, url, content)
        doc = json.loads(catalog.read_text())
        files = doc["entries"][0]["files"]
        for f in files:
            f.update(provider="google-drive", drive_file_id=f["name"].replace(".", "_"))
            f["url"] = f"https://drive.google.com/file/d/{f['drive_file_id']}/view"
        catalog.write_text(json.dumps(doc))
        stage = root / "cache/test-grid.partial"
        stage.mkdir(parents=True)
        f = files[1]
        (stage / "grid.unf.part").write_bytes(content["grid.unf"][:17])
        (stage / "resume.json").write_text(
            json.dumps(
                {
                    "grid.unf": dict(
                        url=f["url"], etag=None, size=f["size"], sha256=f["sha256"]
                    )
                }
            )
        )

        def direct(file, headers, timeout):
            if not honor_range:
                headers = {k: v for k, v in headers.items() if k != "Range"}
            return urlopen(
                Request(url + file["name"], headers=headers), timeout=timeout
            )

        with patch("aspcap.grids._open_file", side_effect=direct):
            header = grids.download("test-grid", root / "cache", catalog=catalog)
        assert header.with_suffix(".unf").read_bytes() == content["grid.unf"]
        assert (
            ("/grid.unf", "bytes=17-") in requests
            if honor_range
            else all(value is None for _, value in requests)
        )


def test_catalog_pins_test_version_and_rejects_unpinned_drive_files(tmp_path):
    release = grids.info("hge-giants-test-v1")
    assert (
        release["status"] == "test"
        and release["lsf"] == "a"
        and release["instrument"] == "apogee-n"
    )
    assert release["pca_components_per_block"] == 105 and release["spectra"] == 673200
    assert len(release["files"]) == 15 and all(
        len(f["sha256"]) == 64 for f in release["files"]
    )
    release["files"][0].pop("sha256")
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(dict(version=1, entries=[release])))
    with pytest.raises(ValueError, match="pinned size and SHA256"):
        grids.available(catalog=path)


def test_loader_progress_control_is_explicit(capsys):
    with patch("aspcap.grids.download") as download, patch.object(
        Grid, "__init__", return_value=None
    ):
        Grid.from_download("test-grid")
        callback = download.call_args.kwargs["progress"]
        callback("grid.unf", 0, 1024**3)
        callback("grid.unf", 1024, 1024**3)
        callback("grid.unf", 1024**3, 1024**3)
        assert capsys.readouterr().err.count("grid.unf") == 2
        Grid.from_download("test-grid", progress=False)
        assert download.call_args.kwargs["progress"] is None
        with pytest.raises(TypeError):
            Grid.from_download("test-grid", progress="yes")
