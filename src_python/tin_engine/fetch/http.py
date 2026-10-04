"""Byte ranges over HTTP, stdlib only (increment 23a-2).

`docs/increments/23-basin-scale.md`, "Downloading". The one module in
`tin_engine` that imports `urllib.request` or `http.client` (F11), so network
code stays in `fetch/` and nothing on `rasputin mesh`'s path can reach it.

A range answer must be 206 with `Content-Range: bytes start-(stop-1)/total`
and exactly that many bytes; a 200 is refused without reading its body (it
would be the whole file). One range per request (decided 3). Timeouts,
connection errors, short bodies and 5xx are retried after each of `delays`;
4xx and a malformed answer never are.

Every request names its client, `User-Agent: rasputin/<version>`, instead of
urllib's default (increment 29, "Data use"). `query_url` builds the query
strings of `fetch/nve.py`, which may not import `urllib` itself.
"""

from __future__ import annotations

import http.client
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Mapping
from dataclasses import dataclass

from tin_engine import installed_version

_CONTENT_RANGE = re.compile(r"bytes (\d+)-(\d+)/(\d+)")


def user_agent() -> str:
    """`rasputin/<version>`: the client names itself to every server."""
    return f"rasputin/{installed_version()}"


def query_url(base: str, params: Mapping[str, str]) -> str:
    """`base?params`, each value percent-encoded."""
    return f"{base}?{urllib.parse.urlencode(params)}"


class FetchError(Exception):
    """A refusal of the fetch step, naming the URL and, for a range, the range."""


class _TransientError(Exception):
    """Worth another try: a timeout, a dropped connection, a short body, a 5xx."""


@dataclass(frozen=True, slots=True)
class RangeResponse:
    """A range's bytes, the file's total length and its `Last-Modified`."""

    data: bytes
    total: int
    last_modified: str


class RangeClient:
    """Blocking GETs; `fetch/run.py` runs them in threads. Holds no connection."""

    def __init__(self, *, delays: tuple[float, ...] = (1, 2, 4), timeout: float = 60.0) -> None:
        self.delays, self.timeout = delays, timeout

    def get(self, url: str, start: int, stop: int) -> RangeResponse:
        """Bytes `[start, stop)` of `url`. A range reaching past the file's end
        is answered up to the end (RFC 9110 §14.1.2), so a file shorter than a
        header prefix is still read whole."""
        return self._retry(f"{url} bytes {start}-{stop - 1}", lambda: self._get(url, start, stop))

    def get_text(self, url: str) -> str:
        """The whole of `url` as UTF-8 (GLO-30's tile list)."""

        def call() -> str:
            request = urllib.request.Request(url, headers={"User-Agent": user_agent()})
            with self._open(request, url) as response:
                return self._read(response, url).decode()

        return self._retry(url, call)

    def _get(self, url: str, start: int, stop: int) -> RangeResponse:
        where = f"{url} bytes {start}-{stop - 1}"
        request = urllib.request.Request(
            url, headers={"Range": f"bytes={start}-{stop - 1}", "User-Agent": user_agent()}
        )
        with self._open(request, where) as response:
            if response.status != 206:  # the body is left unread; closing drops it
                raise FetchError(f"{where}: HTTP {response.status}, not 206 Partial Content")
            got = _CONTENT_RANGE.fullmatch(response.headers.get("Content-Range", ""))
            first, last, total = (int(g) for g in got.groups()) if got else (-1, -1, -1)
            if first != start or not (last == stop - 1 or (last == total - 1 and stop > total)):
                shown = response.headers.get("Content-Range")
                raise FetchError(f"{where}: Content-Range {shown!r} is not the range asked")
            data = self._read(response, where)
            if len(data) != last + 1 - start:
                raise _TransientError(f"{len(data)} bytes, not {last + 1 - start}")
            return RangeResponse(data, total, response.headers.get("Last-Modified", ""))

    def _open(self, request: urllib.request.Request, where: str) -> http.client.HTTPResponse:
        try:
            response: http.client.HTTPResponse = urllib.request.urlopen(
                request, timeout=self.timeout
            )
        except urllib.error.HTTPError as exc:
            exc.close()
            if exc.code >= 500:
                raise _TransientError(f"HTTP {exc.code}") from exc
            raise FetchError(f"{where}: HTTP {exc.code} {exc.reason}") from exc
        except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
            raise _TransientError(str(exc)) from exc
        return response

    @staticmethod
    def _read(response: http.client.HTTPResponse, where: str) -> bytes:
        try:
            return response.read()
        except (OSError, http.client.HTTPException) as exc:
            raise _TransientError(f"the body broke off: {exc!r}") from exc

    def _retry[T](self, where: str, call: Callable[[], T]) -> T:
        for delay in (*self.delays, None):
            try:
                return call()
            except _TransientError as exc:
                if delay is None:
                    raise FetchError(f"{where}: {exc}, after {len(self.delays) + 1} tries") from exc
                time.sleep(delay)
        raise AssertionError("unreachable")


__all__ = ["FetchError", "RangeClient", "RangeResponse", "query_url", "user_agent"]
