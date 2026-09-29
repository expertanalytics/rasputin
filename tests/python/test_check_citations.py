"""`tools/check_citations.py` resolving `legacy/` citations through a tag.

Release hygiene (`docs/increments/release-hygiene.md`, section 3 item 3) takes
`legacy/` out of the working tree. The prose still cites it by line, 147
times, as the increment records' evidence, so instead of rewriting them the
tool resolves a citation whose path starts with `legacy/` against the
annotated tag `legacy-archive` (`git cat-file -p legacy-archive:<path>`), and
reports it broken if the tag or the file is missing or the line is past the
end.

The tool finds its repository from its own path (`REPO`), so each test copies
it into a scratch repository and runs it there as a subprocess: the tests see
only its exit status and its report, never its internals.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[2] / "tools" / "check_citations.py"

TAG = "legacy-archive"
ARCHIVED = "legacy/rasputin/reader.py"
ARCHIVED_LINES = 5


@dataclass(frozen=True)
class Run:
    code: int
    out: str


def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, text=True)


def _commit_all(repo: Path, message: str) -> None:
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", message)


def _cite(repo: Path, *citations: str) -> None:
    """Commit a doc citing each of `citations`, one per line, in backticks."""
    body = "".join(f"See `{c}`.\n" for c in citations)
    (repo / "docs" / "record.md").write_text(f"# Record\n\n{body}")
    _commit_all(repo, "cite")


def _check(repo: Path) -> Run:
    result = subprocess.run(
        [sys.executable, str(repo / "tools" / "check_citations.py"), "--paths", "docs"],
        cwd=repo,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return Run(result.returncode, result.stdout + result.stderr)


def _repo(tmp_path: Path, *, tag: bool) -> Path:
    """A repository on master whose `legacy/` exists only in history.

    The archived file has `ARCHIVED_LINES` lines. With `tag`, the annotated tag
    `legacy-archive` points at the commit that still has `legacy/`; the next
    commit deletes it, as step 3 of the plan will.
    """
    repo = tmp_path / "repo"
    (repo / "tools").mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "src_python").mkdir()
    _git(repo, "init", "-q", "-b", "master")
    shutil.copy(TOOL, repo / "tools" / "check_citations.py")
    (repo / "src_python" / "live.py").write_text("a = 1\nb = 2\nc = 3\n")
    archived = repo / ARCHIVED
    archived.parent.mkdir(parents=True)
    archived.write_text("".join(f"line_{i} = {i}\n" for i in range(1, ARCHIVED_LINES + 1)))
    _commit_all(repo, "with legacy")
    if tag:
        _git(repo, "-c", "user.name=t", "-c", "user.email=t@t", "tag", "-a", TAG, "-m", "archive")
    shutil.rmtree(repo / "legacy")
    _commit_all(repo, "legacy leaves the tree")
    return repo


@pytest.fixture
def archived_repo(tmp_path: Path) -> Path:
    return _repo(tmp_path, tag=True)


@pytest.fixture
def untagged_repo(tmp_path: Path) -> Path:
    return _repo(tmp_path, tag=False)


def _broken_section(out: str) -> str:
    """The report after its `broken` header, or '' when nothing is broken."""
    _, header, rest = out.partition("== broken")
    return header + rest


# -- resolved through the tag -------------------------------------------------


@pytest.mark.parametrize(
    "citation",
    [f"{ARCHIVED}:1", f"{ARCHIVED}:{ARCHIVED_LINES}", f"{ARCHIVED}:2-{ARCHIVED_LINES}"],
    ids=["first-line", "last-line", "range-ending-on-last-line"],
)
def test_a_legacy_citation_resolves_through_the_tag(archived_repo: Path, citation: str) -> None:
    _cite(archived_repo, citation)

    run = _check(archived_repo)

    assert run.code == 0, run.out
    assert "All citations resolve" in run.out


def test_the_tag_is_read_even_where_the_working_tree_still_has_legacy(
    archived_repo: Path,
) -> None:
    # The plan resolves `legacy/` against the tag, not the working tree: a
    # stray or regrown `legacy/` in the checkout must not make a citation past
    # the end of the archived file resolve.
    regrown = archived_repo / ARCHIVED
    regrown.parent.mkdir(parents=True)
    regrown.write_text("x = 0\n" * (ARCHIVED_LINES + 10))
    _commit_all(archived_repo, "legacy regrows")
    _cite(archived_repo, f"{ARCHIVED}:{ARCHIVED_LINES + 3}")

    run = _check(archived_repo)

    assert run.code == 1, run.out
    assert f"{ARCHIVED}:{ARCHIVED_LINES + 3}" in _broken_section(run.out)


def test_a_live_citation_still_resolves_in_the_working_tree(archived_repo: Path) -> None:
    _cite(archived_repo, "src_python/live.py:3", f"{ARCHIVED}:2")

    run = _check(archived_repo)

    assert run.code == 0, run.out


# -- broken through the tag ---------------------------------------------------


@pytest.mark.parametrize(
    "citation",
    [f"{ARCHIVED}:{ARCHIVED_LINES + 1}", f"{ARCHIVED}:3-{ARCHIVED_LINES + 1}"],
    ids=["line-past-the-end", "range-past-the-end"],
)
def test_a_legacy_citation_past_the_end_of_the_tagged_file_is_broken(
    archived_repo: Path, citation: str
) -> None:
    _cite(archived_repo, citation)

    run = _check(archived_repo)

    assert run.code == 1, run.out
    broken = _broken_section(run.out)
    assert citation.split("-")[0] in broken
    # The same wording as a live file past its end: the file was found, in the
    # tag, and its length is reported -- not "no such file".
    assert f"file has {ARCHIVED_LINES} lines" in broken


def test_a_legacy_file_missing_from_the_tag_is_broken(archived_repo: Path) -> None:
    _cite(archived_repo, "legacy/rasputin/never_existed.py:1")

    run = _check(archived_repo)

    assert run.code == 1, run.out
    assert "legacy/rasputin/never_existed.py" in _broken_section(run.out)


def test_with_the_tag_absent_a_legacy_citation_is_reported_not_passed(
    untagged_repo: Path,
) -> None:
    _cite(untagged_repo, f"{ARCHIVED}:1")

    run = _check(untagged_repo)

    assert run.code == 1, run.out
    broken = _broken_section(run.out)
    assert ARCHIVED in broken
    # A clone without the tag (CI's shallow checkout) must be told what is
    # missing, so the report names the tag rather than blaming the citation.
    assert TAG in broken
