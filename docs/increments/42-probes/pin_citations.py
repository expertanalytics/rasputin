"""Pin every unpinned line citation into `src_python/` (increment 42, section 7).

A citation `<path>:<n>` whose target resolves under `src_python/` becomes
`<full path>@<commit>:<n>`, where `<commit>` is the commit that last changed
the citing line (`git blame`). That is the tree its author was reading. The
citation's own patterns and path resolution are `tools/check_citations.py`'s,
imported, so this script finds exactly what that gate checks.

    python3 docs/increments/42-probes/pin_citations.py          # list only
    python3 docs/increments/42-probes/pin_citations.py --write  # rewrite

A pin whose file is missing or too short at that commit is not written; it is
listed under "refused" for a person to judge. Rerun after `--write`: it then
finds nothing.
"""

from __future__ import annotations

import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "tools"))
import check_citations as cc  # noqa: E402


def blame(rel: str, number: int) -> str:
    out = subprocess.run(
        ["git", "blame", "--porcelain", "-L", f"{number},{number}", "--", rel],
        cwd=REPO, capture_output=True, text=True, check=True,
    ).stdout
    return out.split(" ", 1)[0][:8]


def main() -> int:
    write = "--write" in sys.argv[1:]
    edits: dict[Path, list[tuple[int, str, str]]] = defaultdict(list)
    refused: list[str] = []
    for source in cc.default_sources(None):
        rel_source = source.relative_to(REPO).as_posix()
        try:
            text_by_line = cc.scanned_text(source)
        except (SyntaxError, ValueError):
            continue
        for number, fragments in sorted(text_by_line.items()):
            for fragment in fragments:
                rest = fragment = cc.ABSOLUTE.sub("", fragment)
                for m in cc.PINNED.finditer(fragment):
                    rest = rest[: m.start()] + " " * len(m.group(0)) + rest[m.end() :]
                for m in cc.CITATION.finditer(rest):
                    cited, start, end = m.groups()
                    if cited.startswith(cc.LEGACY_PREFIX):
                        continue
                    target = cc.resolve(cited)
                    if not isinstance(target, Path):
                        continue
                    rel_target = target.relative_to(REPO).as_posix()
                    if not rel_target.startswith("src_python/"):
                        continue
                    sha = blame(rel_source, number)
                    last = max(int(start), int(end or 0))
                    old = m.group(0)
                    new = f"{rel_target}@{sha}:{start}" + (f"-{end}" if end else "")
                    if cc.check_pinned(rel_target, sha, last) is not None:
                        refused.append(f"{rel_source}:{number}: {old} -> {new}")
                        continue
                    edits[source].append((number, old, new))
    total = sum(len(v) for v in edits.values())
    for source, items in sorted(edits.items()):
        lines = source.read_text().split("\n")
        for number, old, new in items:
            print(f"{source.relative_to(REPO)}:{number}: {old} -> {new}")
            line = lines[number - 1]
            # The cited token, whole: not preceded or followed by a path character.
            at = next(
                i for i in range(len(line))
                if line.startswith(old, i)
                and (i == 0 or not (line[i - 1].isalnum() or line[i - 1] in "_./-"))
                and not line[i + len(old) : i + len(old) + 1].isdigit()
            )
            lines[number - 1] = line[:at] + new + line[at + len(old) :]
        if write:
            source.write_text("\n".join(lines))
    print(f"TOTAL {total} in {len(edits)} files; refused {len(refused)}")
    for item in refused:
        print(f"refused: {item}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
