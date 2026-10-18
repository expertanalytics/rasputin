"""Increment 69, T1 and T2: the bearer tokens, read from the environment and verified.

`docs/increments/69-mcp-http.md`, sections 6.1 and 9. `tokens_from_env` and
`StaticTokens` are tested as pure code; the 401s and the per-user tool list
over a real loopback server.
"""

from __future__ import annotations

import hashlib
import hmac
import re
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
#: A reversed pair's left side: 32 hex digits starting with a letter, as `openssl rand
#: -hex 16` gives about 6 times in 16, which passes the name check.
HEX = "deadbeef" * 4
#: A reversed pair's left side as `secrets.token_urlsafe` writes one: capitals, which
#: fail the name check.
URLSAFE = "Xq3vB9_kLm2-Rt7sWp4zYc8nHd5fGj1aEo6u"
NAME_RULE = (
    "the name is not a lowercase letter then up to 31 of a-z, 0-9, _ and - "
    "(a pair is user=token, the name first)"
)


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
        ("value", "refusal", "hidden"),
        [
            pytest.param(
                f"Ola={OLA}", f"pair 1: {NAME_RULE}", ["Ola", OLA], id="capital in the name"
            ),
            pytest.param(
                f"../x={OLA}",
                f"pair 1: {NAME_RULE}",
                ["../x", OLA],
                id="dots and slash in the name",
            ),
            pytest.param(
                f"ola={OLA} ola={KARI}",
                "pair 2 repeats the name or the token of an earlier pair",
                ["ola", KARI],
                id="repeated name",
            ),
            pytest.param(
                f"ola={OLA} kari={OLA}",
                "pair 2 repeats the name or the token of an earlier pair",
                ["kari", OLA],
                id="repeated token",
            ),
            pytest.param(
                f"ola={'s' * 31}",
                "pair 1: the token is under 32 characters",
                ["ola", "s" * 31],
                id="31-character token",
            ),
            pytest.param(
                f"ola={OLA} {KARI}", "pair 2 is not user=token", [KARI], id="no equals sign"
            ),
            pytest.param(
                f"{HEX}=ola",
                "pair 1: the token is under 32 characters",
                [HEX, "ola"],
                id="reversed, hex token",
            ),
            pytest.param(
                f"{URLSAFE}=ola",
                f"pair 1: {NAME_RULE}",
                [URLSAFE, "ola"],
                id="reversed, urlsafe token",
            ),
            pytest.param(
                f"ola={OLA} Kari={KARI}",
                f"pair 2: {NAME_RULE}",
                ["Kari", KARI],
                id="a bad pair after a good one",
            ),
        ],
    )
    def test_a_bad_pair_is_refused_naming_variable_and_position_never_either_side(
        self, value: str, refusal: str, hidden: list[str]
    ) -> None:
        with pytest.raises(ValueError) as caught:
            http_module().tokens_from_env({VARIABLE: value})
        words = str(caught.value)
        for side in hidden:
            assert not shows(side, words), f"the refusal shows {side!r}: {words}"
        assert words == f"{VARIABLE}: {refusal}", words


def shows(word: str, text: str) -> bool:
    """Whether ``word`` is in ``text`` as a whole word: delimited by the text's ends, space,
    or punctuation, so ``ola`` is not found inside an ordinary word."""
    return re.search(rf"(?:^|[\s('\"=])({re.escape(word)})(?=$|[\s,.;:)'\"=])", text) is not None


def test_shows_finds_whole_words_only() -> None:
    assert (
        shows("ola", "user ola repeats") and shows("ola", "x=ola") and shows("../x", "user ../x is")
    )
    assert not shows("ola", "a payolas list") and not shows("ola", "olaf")


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
