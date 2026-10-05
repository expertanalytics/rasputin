"""Is a path under a Claude Code session scratchpad? (docs/increments/h16-harness-fixes.md §3)

The scratchpad is a session's temporary directory, `<tmp>/claude-<uid>/<project>/
<session>/scratchpad`. Any session's counts, not only the current one: all are
temporary, and nothing the harness reads lives there. guard_governance.py imports
this (G3a), so it is governed. Standard library only.
"""

import os
import re

#: Matched against the real path, so /tmp (a link to /private/tmp on macOS) and a
#: symlink out of the scratchpad resolve before the match. The tests rewrite the
#: one literal below to point the pattern elsewhere, so it is stated once.
PATTERN = re.compile(r"^/private/tmp/claude-\d+/[^/]+/[^/]+/scratchpad(/|$)")


def under(path: str) -> bool:
    """True for an absolute path whose real path lies under a scratchpad; relative paths never."""
    return os.path.isabs(path) and PATTERN.match(os.path.realpath(path)) is not None
