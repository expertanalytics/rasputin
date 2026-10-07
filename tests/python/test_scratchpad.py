"""`tools/scratchpad.py`: is a path under a session scratchpad (h16 §3, G3a).

`docs/increments/h16-harness-fixes.md` §2 G3a fixes the pattern,
`^/private/tmp/claude-\\d+/[^/]+/[^/]+/scratchpad(/|$)`, matched against the
path's real path, and §3 says `under` takes absolute paths only. Any session's
scratchpad counts, not only the current one. These paths need not exist:
`os.path.realpath` resolves what does and keeps the rest.
"""

from __future__ import annotations

import sys

import pytest

from harness_fixtures import Tool

scratchpad = Tool("scratchpad")

PAD = "/private/tmp/claude-501/-Users-x-project/0f1e2d3c-session/scratchpad"

UNDER = (
    PAD,
    f"{PAD}/",
    f"{PAD}/r/.git/config",
    f"{PAD}/copy/CLAUDE.md",
    "/private/tmp/claude-0/p/s/scratchpad/x",
)

NOT_UNDER = (
    f"{PAD}s/x",  # a sibling whose name only starts with scratchpad
    "/private/tmp/claude-501/p/scratchpad/x",  # one level short
    "/private/tmp/claude-501/p/s/t/scratchpad/x",  # one level deep
    "/private/tmp/claude-abc/p/s/scratchpad/x",  # not a uid
    "/private/tmp/x/claude-501/p/s/scratchpad/x",  # not at /private/tmp
    f"{PAD}/../../../../x",  # resolves to /private/tmp/x
    "/Users/x/scratchpad/y",
    "scratchpad/x",  # relative: the guard does not track cd
    "r/.git/config",
)


@pytest.mark.parametrize("path", UNDER)
def test_a_path_in_a_session_scratchpad_is_under_it(path: str) -> None:
    assert scratchpad.under(path) is True


@pytest.mark.parametrize("path", NOT_UNDER)
def test_a_path_elsewhere_is_not_under_a_scratchpad(path: str) -> None:
    assert scratchpad.under(path) is False


@pytest.mark.skipif(sys.platform != "darwin", reason="/tmp links to /private/tmp on macOS only")
def test_tmp_resolves_to_private_tmp() -> None:
    assert scratchpad.under("/tmp/claude-501/p/s/scratchpad/x") is True
