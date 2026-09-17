#!/usr/bin/env python3
"""Flag the prose citations a branch has put at risk.

Governance and increment docs cite each other and the source by line number --
`01-predicates.md:151`, `detria.hpp:1240-1242`. Line numbers drift, and worse,
the *text* at a line can be rewritten to say the opposite while the number
still resolves. On 2026-09-17 a sweep resolved all 46 citations in
`kernel-sufficiency-audit.md` and was still wrong three times, because three of
them pointed at text the same branch had rewritten. `@reviewer` put it exactly:
resolving a line number is not the same as resolving the quotation.

So this script does not claim to find wrong citations. It computes the two
things that *are* mechanically decidable:

  broken   -- the file is missing, or the line is past its end. Cheap, certain.
  at-risk  -- the cited file is modified by this branch, so the number may
              still resolve while the quotation no longer holds. These must be
              re-read against the new text; no script can read them for you.

Exit status is 1 only for `broken`. `at-risk` is a worklist, not a failure:
a branch that edits a cited file and updates the citing prose correctly is
right, and the script cannot tell that from the wrong case.

Usage: python3 tools/check_citations.py [--base master] [--paths docs .claude]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# A citation is a path-like token with a line or line range, inside backticks:
# `01-predicates.md:151`, `include/terrain/core/segment.hpp:68-74`. Bare
# parenthesised forms like (detria.hpp:391) are also used, so backticks are not
# required -- but an extension is, to avoid matching ratios and timestamps.
CITATION = re.compile(
    r"(?<![\w/.-])([\w./-]+\.(?:md|py|hpp|cpp|h|yaml|yml|toml|txt|cmake)):(\d+)(?:-(\d+))?"
)

# Where a bare basename may live. More than one hit is reported as ambiguous
# rather than picked between.
SEARCH_ROOTS = ("docs", ".claude", "include", "src", "src_python", "tests", "tools", "lib")

SCAN_SUFFIXES = (".md",)


def changed_files(base: str) -> set[str] | None:
    """Repo-relative paths this branch modifies, or None if base is unknown."""
    for ref in (f"{base}...HEAD", base):
        result = subprocess.run(
            ["git", "diff", "--name-only", ref],
            cwd=REPO,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return {line for line in result.stdout.split("\n") if line}
    return None


def resolve(cited: str) -> Path | list[Path] | None:
    """A cited path or bare basename, as a file, an ambiguity, or nothing.

    An ambiguous basename is returned as the candidate list rather than
    silently resolved to the first hit: guessing which `README.md:40` meant
    would make this script commit the error it exists to catch.
    """
    direct = REPO / cited
    if direct.is_file():
        return direct
    name = Path(cited).name
    hits = [
        candidate
        for root in SEARCH_ROOTS
        for candidate in sorted((REPO / root).rglob(name))
        if candidate.is_file()
    ]
    if len(hits) == 1:
        return hits[0]
    return hits or None


def line_count(path: Path) -> int:
    return len(path.read_text(errors="replace").splitlines())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="master")
    parser.add_argument("--paths", nargs="*", default=["docs", ".claude", "CLAUDE.md"])
    args = parser.parse_args()

    maybe_touched = changed_files(args.base)
    diffable = maybe_touched is not None
    if maybe_touched is None:
        print(f"warning: cannot diff against '{args.base}'; at-risk check disabled")
        touched: set[str] = set()
    else:
        touched = maybe_touched

    sources: list[Path] = []
    for entry in args.paths:
        target = REPO / entry
        if target.is_file():
            sources.append(target)
        elif target.is_dir():
            sources.extend(
                p for p in sorted(target.rglob("*")) if p.suffix in SCAN_SUFFIXES
            )

    broken: list[str] = []
    at_risk: list[str] = []

    for source in sources:
        rel_source = source.relative_to(REPO) if source.is_relative_to(REPO) else source
        for number, line in enumerate(
            source.read_text(errors="replace").splitlines(), start=1
        ):
            for cited, start, end in CITATION.findall(line):
                where = f"{rel_source}:{number}"
                target = resolve(cited)
                if target is None:
                    broken.append(f"{where}: cites '{cited}' -- no such file")
                    continue
                if isinstance(target, list):
                    names = ", ".join(str(p.relative_to(REPO)) for p in target)
                    broken.append(
                        f"{where}: cites '{cited}' -- ambiguous ({names}); "
                        f"cite the path, not the basename"
                    )
                    continue
                last = max(int(start), int(end or 0))
                total = line_count(target)
                if last > total:
                    broken.append(
                        f"{where}: cites '{cited}:{start}' -- file has {total} lines"
                    )
                    continue
                rel_target = str(target.relative_to(REPO))
                if rel_target in touched:
                    at_risk.append(
                        f"{where}: cites '{cited}:{start}', and this branch edits "
                        f"{rel_target} -- re-read the quotation, not the number"
                    )

    if at_risk:
        print(f"== at risk ({len(at_risk)}) -- cited text changed on this branch ==")
        for item in at_risk:
            print(f"  {item}")
    if broken:
        print(f"\n== broken ({len(broken)}) ==")
        for item in broken:
            print(f"  {item}")
        return 1
    if not at_risk:
        print(
            "All citations resolve"
            + (
                ", and none point into a file this branch edits."
                if diffable
                else ". At-risk was not evaluated -- see the warning above."
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
