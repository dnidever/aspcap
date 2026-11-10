"""Anonymous downloads of public Drive files; no browser or account cookies.

Drive may return an HTML download-confirmation form instead of file bytes.
Only that specific form, pointing back to a Google download host, is followed.
Login, quota and other HTML pages are errors, never grid payloads.
"""

from html.parser import HTMLParser
from http.cookiejar import CookieJar
from urllib.error import HTTPError
from urllib.parse import parse_qsl, urlencode, urljoin, urlsplit, urlunsplit
from urllib.request import HTTPCookieProcessor, Request, build_opener


class _Confirmation(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.action = None
        self.fields = {}
        self.inside = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "form" and attrs.get("id") == "download-form":
            self.inside = True
            self.action = attrs.get("action")
        elif (
            self.inside
            and tag == "input"
            and attrs.get("type") == "hidden"
            and attrs.get("name")
        ):
            self.fields[attrs["name"]] = attrs.get("value", "")

    def handle_endtag(self, tag):
        if tag == "form":
            self.inside = False

    def url(self, previous, file_id):
        if not self.action:
            raise ValueError("Drive did not return a download confirmation form")
        url = urljoin(previous, self.action)
        parts = urlsplit(url)
        if (
            parts.scheme != "https"
            or parts.username
            or parts.password
            or parts.port not in (None, 443)
            or parts.hostname
            not in (
                "drive.google.com",
                "docs.google.com",
                "drive.usercontent.google.com",
            )
            or parts.path not in ("/download", "/uc")
        ):
            raise ValueError("Unexpected Drive download confirmation destination")
        fields = dict(parse_qsl(parts.query))
        fields.update(self.fields)
        if fields.get("id") != file_id:
            raise ValueError("Drive confirmation refers to a different file")
        return urlunsplit(parts._replace(query=urlencode(fields), fragment=""))


def open_file(file, *, headers=None, timeout=60):
    """Return an open binary HTTP response for a pinned public Drive file.

    The caller must close the response. Range/If-Range headers survive the
    confirmation flow. HTML is bounded to 512 KiB and never written to disk.
    This handles a download acknowledgement, not CAPTCHAs or login screens.
    """
    file_id = file["drive_file_id"]
    url = "https://drive.usercontent.google.com/download?" + urlencode(
        dict(id=file_id, export="download", confirm="t")
    )
    opener = build_opener(HTTPCookieProcessor(CookieJar()))
    request_headers = {"Accept-Encoding": "identity", **(headers or {})}
    for _ in range(3):
        try:
            response = opener.open(
                Request(url, headers=request_headers), timeout=timeout
            )
        except HTTPError as error:
            error.close()
            raise RuntimeError(
                f"Google Drive download failed (HTTP {error.code}) for {file['name']}. The public file may be unavailable or temporarily rate-limited; retry later or use {file['url']}"
            ) from None
        content_type = (
            response.headers.get("Content-Type", "").split(";", 1)[0].strip().lower()
        )
        attachment = (
            response.headers.get("Content-Disposition", "")
            .lower()
            .startswith("attachment")
        )
        if attachment and content_type not in ("text/html", "application/xhtml+xml"):
            return response
        with response:
            if content_type not in ("text/html", "application/xhtml+xml"):
                raise RuntimeError(
                    f"Drive returned no file attachment for {file['name']}"
                )
            page = response.read(512 * 1024 + 1)
            if len(page) > 512 * 1024:
                raise RuntimeError("Drive confirmation page was unexpectedly large")
            parser = _Confirmation()
            parser.feed(page.decode("utf-8", errors="replace"))
            try:
                url = parser.url(response.geturl(), file_id)
            except ValueError as error:
                raise RuntimeError(
                    f"Cannot download {file['name']} from public Google Drive: {error}. Check public sharing or traffic limits at {file['url']}"
                ) from None
    raise RuntimeError(
        f"Drive repeatedly returned confirmation pages for {file['name']}"
    )
