"""Increment 69, T1 and T2: the bearer tokens, read from the environment and verified.

`docs/increments/69-mcp-http.md`, sections 6.1 and 9. `tokens_from_env` and
`StaticTokens` are tested as pure code; the 401s and the per-user tool list
over a real loopback server.
"""

from __future__ import annotations

import hashlib
import hmac
from collections.abc import Iterator
from typing import Any

import pytest

pytest.importorskip("mcp")

from mcp_driver import tool_names
from mcp_http_driver import TOKENS, HttpServer, client, post, rasputin_server
from tin_engine.cli import app

VARIABLE = "RASPUTIN_MCP_TOKENS"
#: Two tokens at the floor and above it: 32 and 43 characters.
OLA = "A" * 32
KARI = "kari-" + "b" * 38


def digest(token: str) -> bytes:
    return hashlib.sha256(token.encode()).digest()


def http_module() -> Any:
    from tin_engine.mcp import http

    return http


# ---------------------------------------------------------------- T1


class TestTokensFromTheEnvironment:
    def test_two_users_parse_to_the_sha256_of_each_token(self) -> None:
        tokens = http_module().tokens_from_env({VARIABLE: f"ola={OLA}\n  kari={KARI}\t"})
        assert tokens == {"ola": digest(OLA), "kari": digest(KARI)}

    def test_a_32_character_token_is_accepted(self) -> None:
        assert http_module().tokens_from_env({VARIABLE: f"ola={'x' * 32}"}) == {
            "ola": digest("x" * 32)
        }

    @pytest.mark.parametrize(
        "environ", [{}, {VARIABLE: ""}, {VARIABLE: "  \n"}], ids=["missing", "empty", "blank"]
    )
    def test_a_missing_or_empty_variable_is_refused_naming_it(
        self, environ: dict[str, str]
    ) -> None:
        with pytest.raises(ValueError, match=VARIABLE):
            http_module().tokens_from_env(environ)

    @pytest.mark.parametrize(
        ("value", "user", "secret"),
        [
            pytest.param(f"Ola={OLA}", "Ola", OLA, id="capital in the name"),
            pytest.param(f"../x={OLA}", "../x", OLA, id="dots and slash in the name"),
            pytest.param(f"ola={OLA} ola={KARI}", "ola", KARI, id="repeated name"),
            pytest.param(f"ola={OLA} kari={OLA}", "kari", OLA, id="repeated token"),
            pytest.param(f"ola={'s' * 31}", "ola", "s" * 31, id="31-character token"),
        ],
    )
    def test_a_bad_pair_is_refused_naming_variable_and_user_never_the_token(
        self, value: str, user: str, secret: str
    ) -> None:
        with pytest.raises(ValueError) as caught:
            http_module().tokens_from_env({VARIABLE: value})
        words = str(caught.value)
        assert VARIABLE in words, words
        assert user in words, words
        assert secret not in words, f"the refusal shows a token: {words}"


# ---------------------------------------------------------------- T2


class TestTheVerifier:
    @pytest.fixture
    def verifier(self) -> Any:
        return http_module().StaticTokens({"ola": digest(OLA), "kari": digest(KARI)})

    async def test_the_right_token_gives_its_user_and_holds_no_token(self, verifier: Any) -> None:
        for user, token in (("ola", OLA), ("kari", KARI)):
            found = await verifier.verify_token(token)
            assert found is not None
            assert found.client_id == user
            assert found.token == user, "the token itself is kept in the request's scope"

    @pytest.mark.parametrize(
        "token", ["c" * 40, "", OLA + "A", OLA[:-1]], ids=["wrong", "empty", "one more", "one less"]
    )
    async def test_any_other_token_gives_none(self, verifier: Any, token: str) -> None:
        assert await verifier.verify_token(token) is None

    async def test_every_stored_digest_is_compared_with_no_early_exit(
        self, verifier: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        calls: list[tuple[Any, Any]] = []
        original = hmac.compare_digest

        def counting(a: Any, b: Any) -> bool:
            calls.append((a, b))
            return original(a, b)

        monkeypatch.setattr(hmac, "compare_digest", counting)
        found = await verifier.verify_token(OLA)  # the first user: a match ends nothing
        assert found is not None and found.client_id == "ola"
        assert len(calls) == 2, calls


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    with rasputin_server(tmp_path_factory.mktemp("tokens")) as running:
        yield running


LIST = {"jsonrpc": "2.0", "id": 1, "method": "tools/list"}


class TestOverHttp:
    def test_no_header_is_401(self, server: HttpServer) -> None:
        assert post(server, LIST).status_code == 401

    def test_a_wrong_token_is_401(self, server: HttpServer) -> None:
        response = post(server, LIST, {"Authorization": "Bearer " + "c" * 40})
        assert response.status_code == 401

    async def test_both_users_see_the_same_tools_download_link_included(
        self, server: HttpServer
    ) -> None:
        seen = {}
        for user in TOKENS:
            async with client(server, user) as c:
                seen[user] = [t.name for t in (await c.list_tools()).tools]
        assert seen["ola"] == seen["kari"]
        assert set(seen["ola"]) == {*tool_names(app), "download_link"}
