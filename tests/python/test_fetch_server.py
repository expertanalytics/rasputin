"""F1: the range server can serve, and it can fail in every way F6 needs (23a-2).

`docs/increments/23-basin-scale.md`, "Tests @tester can write red", 23a-2 F1.
Each fault is shown here against a plain `http.client` caller, which accepts
what it is given: a fault the fetch step must refuse is therefore one the
fixture really provokes, not a test that passes because nothing went wrong.
This file is green before 23a-2 lands; it tests the fixture, not the code.
"""

from __future__ import annotations

import http.client
from collections.abc import Iterator
from urllib.parse import urlsplit

import pytest

from fetch_fixtures import CHANGED, LAST_MODIFIED, MIB, RangeServer, padded

DATA = padded(bytes(range(256)) * 64, 2 * MIB)


@pytest.fixture
def server() -> Iterator[RangeServer]:
    with RangeServer() as running:
        running.put("f.tif", DATA)
        yield running


def get(server: RangeServer, path: str, header: str | None) -> http.client.HTTPResponse:
    """One GET with an optional `Range` header; the caller reads or closes."""
    url = urlsplit(server.url(path))
    connection = http.client.HTTPConnection(url.hostname or "", url.port, timeout=10)
    connection.request("GET", url.path, headers={"Range": header} if header else {})
    return connection.getresponse()


class TestTheFixtureServes:
    def test_a_range_is_206_with_its_content_range_and_bytes(self, server: RangeServer) -> None:
        response = get(server, "f.tif", "bytes=100-199")
        assert response.status == 206
        assert response.headers["Content-Range"] == f"bytes 100-199/{len(DATA)}"
        assert response.headers["Last-Modified"] == LAST_MODIFIED
        assert response.read() == DATA[100:200]
        assert server.ranges("f.tif") == [(100, 200)]

    def test_a_range_past_the_end_is_clipped_to_the_file(self, server: RangeServer) -> None:
        response = get(server, "f.tif", f"bytes={len(DATA) - 10}-{len(DATA) + 99}")
        assert response.status == 206
        assert response.headers["Content-Range"].endswith(f"-{len(DATA) - 1}/{len(DATA)}")
        assert response.read() == DATA[-10:]

    def test_an_unknown_path_is_404(self, server: RangeServer) -> None:
        assert get(server, "absent.tif", "bytes=0-9").status == 404

    def test_the_log_keeps_every_request_in_order(self, server: RangeServer) -> None:
        for header in ("bytes=0-9", "bytes=10-19", "bytes=5-5"):
            get(server, "f.tif", header).read()
        assert server.ranges("f.tif") == [(0, 10), (10, 20), (5, 6)]
        assert server.block_ranges("f.tif") == [(10, 20), (5, 6)]


class TestTheFixtureFails:
    def test_fail_after_k_responses(self, server: RangeServer) -> None:
        server.faults.fail_after = 2
        statuses = [get(server, "f.tif", "bytes=0-9").status for _ in range(4)]
        assert statuses == [206, 206, 500, 500]

    @pytest.mark.parametrize("status", [500, 404])
    def test_a_fixed_status(self, server: RangeServer, status: int) -> None:
        server.faults.status = status
        assert get(server, "f.tif", "bytes=0-9").status == status

    def test_ignoring_range_sends_200_and_the_whole_body(self, server: RangeServer) -> None:
        """What a client that does not check the status would take as its range."""
        server.faults.ignore_range = True
        response = get(server, "f.tif", "bytes=0-9")
        assert response.status == 200
        assert "Content-Range" not in response.headers
        assert response.read() == DATA

    def test_a_client_that_hangs_up_on_a_200_is_sent_less_than_the_body(
        self, server: RangeServer
    ) -> None:
        big = padded(b"", 16 * MIB)
        server.put("big.tif", big)
        server.faults.ignore_range = True
        response = get(server, "big.tif", "bytes=0-9")
        assert response.status == 200
        response.close()
        server.wait_done()
        (entry,) = server.log
        assert entry.sent < len(big)

    def test_a_short_body(self, server: RangeServer) -> None:
        server.faults.short = True
        response = get(server, "f.tif", "bytes=100-199")
        with pytest.raises(http.client.IncompleteRead):
            response.read()

    def test_a_wrong_content_range(self, server: RangeServer) -> None:
        server.faults.wrong_content_range = True
        response = get(server, "f.tif", "bytes=100-199")
        assert response.headers["Content-Range"] == f"bytes 101-199/{len(DATA)}"

    def test_a_changed_total_or_date_on_block_responses_only(self, server: RangeServer) -> None:
        server.faults.block_total, server.faults.block_last_modified = 1, CHANGED
        prefix, block = get(server, "f.tif", "bytes=0-9"), get(server, "f.tif", "bytes=100-199")
        assert prefix.headers["Content-Range"] == f"bytes 0-9/{len(DATA)}"
        assert prefix.headers["Last-Modified"] == LAST_MODIFIED
        assert block.headers["Content-Range"] == f"bytes 100-199/{len(DATA) + 1}"
        assert block.headers["Last-Modified"] == CHANGED

    def test_reset_clears_faults_and_log(self, server: RangeServer) -> None:
        server.faults.status = 500
        get(server, "f.tif", "bytes=0-9").read()
        server.reset()
        assert server.log == []
        assert get(server, "f.tif", "bytes=0-9").status == 206
