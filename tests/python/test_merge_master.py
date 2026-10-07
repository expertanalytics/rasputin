"""Tests for `tools/merge_master.py` (h18, T2).

The spec is `docs/increments/h18-worktree-and-merge-tools.md` §4, and the
tests are its §4.9, numbered here as there. No network, no C++ build, no real
suite: `git` runs for real against a temporary repository whose `origin` is a
bare repository on disk, whose master gets first-parent merges with subjects
`Merge pull request #N from x/pr-N`; the gates, `cmake` and the venv's
extension suffix are answered by `MergeTable` (`worktree_fixtures.Recorder`).
The worktree is a linked worktree of the temporary main checkout, on branch
`worktree-x`, with a stub `.venv/bin/python` that is never run.

The current directory is the main checkout, not the worktree: the tool takes
the worktree from its argument only (§2), so a tool that read the current
directory would act on `master` and fail these tests.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import pytest

from worktree_fixtures import (
    TRAILER,
    Answer,
    Outcome,
    Recorder,
    Repos,
    commit,
    failing_fetch,
    git,
    git_dir,
    git_sub,
    head,
    invoke,
    is_venv_python,
    isolated_imports,  # noqa: F401  (autouse: undoes the tool copy's imports)
    load,
    make_repos,
    program,
    push_guard_findings,
    stub_venv,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

TOOL = "tools/merge_master.py"
#: The venv's extension suffix as the table answers it; the dummy `_core` carries it.
SUFFIX = ".cpython-3x-h18test.so"
FAST_GATES = (
    "prohibited_deps", "detria_boundary", "citations", "ruff_check", "ruff_format", "mypy",
)  # fmt: skip
CORE_BYTES = b"h18 dummy _core\n"


def gate_name(argv: list[str]) -> str | None:
    """Which gate a venv-python argv runs, or None."""
    if not is_venv_python(argv):
        return None
    rest = argv[1:]
    scripts = {
        "check_prohibited_deps.py": "prohibited_deps",
        "check_detria_boundary.py": "detria_boundary",
        "check_citations.py": "citations",
    }
    if rest and Path(rest[0]).name in scripts:
        return scripts[Path(rest[0]).name]
    if rest[:2] == ["-m", "ruff"] and "--version" not in rest:
        return "ruff_format" if "format" in rest else "ruff_check"
    if rest[:2] == ["-m", "mypy"] and "--version" not in rest:
        return "mypy"
    if rest[:2] == ["-m", "pytest"]:
        return "pytest"
    return None


def gate_output(name: str, red: bool) -> str:
    return f"h18 gate {name} says {'RED: 1 error' if red else 'ok'}\n"


@dataclass
class MergeTable:
    """Answers the gates, `cmake --build` and the extension suffix for merge_master.py."""

    red: set[str] = field(default_factory=set)
    build_rc: int = 0
    suffix: str = SUFFIX

    def __call__(self, argv: list[str], cwd: Path | None) -> Answer:
        gate = gate_name(argv)
        if gate is not None:
            red = gate in self.red
            return Answer(1 if red else 0, gate_output(gate, red))
        if is_venv_python(argv) and argv[1:2] == ["-c"] and "EXT_SUFFIX" in argv[2]:
            return Answer(0, self.suffix + "\n")
        if program(argv) == "cmake" and "--build" in argv:
            return Answer(self.build_rc, "h18 cmake build output\n")
        pytest.fail(f"merge_master ran {argv} (cwd {cwd}), which the table does not answer")


@dataclass
class Setup:
    repos: Repos
    wt: Path
    old_head: str
    capsys: pytest.CaptureFixture[str]

    @property
    def git_dir(self) -> Path:
        return git_dir(self.wt)

    @property
    def state(self) -> Path:
        return self.git_dir / "merge_master.json"

    @property
    def message_file(self) -> Path:
        return self.git_dir / "MERGE_MASTER_MSG"

    @property
    def merge_head(self) -> Path:
        return self.git_dir / "MERGE_HEAD"

    def run(
        self, args: list[str], table: MergeTable | None = None, *, fail_fetch: bool = False
    ) -> Outcome:
        tool: ModuleType = load(self.repos.main, TOOL)
        recorder = Recorder(table or MergeTable(), failing_fetch if fail_fetch else None)
        return invoke(tool, args, recorder, self.capsys)

    def start(self, *extra: str, table: MergeTable | None = None, **kw: bool) -> Outcome:
        return self.run([str(self.wt), "--persona", "developer", "--trailer", TRAILER, *extra],
                        table, **kw)  # fmt: skip

    def resume(self, *extra: str, table: MergeTable | None = None) -> Outcome:
        return self.run(
            [str(self.wt), "--continue", "--persona", "developer", "--trailer", TRAILER, *extra],
            table,
        )

    def origin_master(self) -> str:
        return git(self.repos.origin, "rev-parse", "master").strip()

    def body(self) -> list[str]:
        return git(self.wt, "log", "-1", "--format=%B").rstrip("\n").splitlines()

    def subject(self) -> str:
        return self.body()[0]


@pytest.fixture
def setup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> Setup:
    """A worktree `worktree-x` off the root commit, with one commit of its own and a stub venv."""
    repos = make_repos(tmp_path)
    wt = repos.worktree("x")
    git(repos.main, "worktree", "add", "-q", "-b", "worktree-x", str(wt), "origin/master")
    commit(wt, {"mine.py": "mine = 1\n"}, "the branch's own work")
    stub_venv(wt)
    (wt / "build-pyext").mkdir()
    monkeypatch.chdir(repos.main)
    return Setup(repos, wt, head(wt), capsys)


def bring_in(setup: Setup, *prs: tuple[int, dict[str, str | None]]) -> None:
    for number, changes in prs:
        setup.repos.pr(number, changes)
    setup.repos.publish()


def conflict_on_a_py_and_claude_md(setup: Setup) -> None:
    """The branch and PR #12 both change a.py's first line and CLAUDE.md's last; #12
    also changes b.py, which the branch does not touch."""
    commit(setup.wt, {"a.py": "x = 'branch'\ny = 2\n",
                      "CLAUDE.md": "# rules\n\nbranch line\n"}, "branch edits")  # fmt: skip
    setup.old_head = head(setup.wt)
    bring_in(setup, (12, {"a.py": "x = 'master'\ny = 2\n",
                          "CLAUDE.md": "# rules\n\nmaster line\n", "b.py": "b = 2\n"}))  # fmt: skip


def resolve_conflicts(setup: Setup) -> None:
    (setup.wt / "a.py").write_text("x = 'both'\ny = 2\n")
    (setup.wt / "CLAUDE.md").write_text("# rules\n\nboth lines\n")
    git(setup.wt, "add", "a.py", "CLAUDE.md")


def assert_no_merge_started(setup: Setup, outcome: Outcome) -> None:
    """The start form wrote nothing: no MERGE_HEAD, no state file, the head unchanged."""
    assert not setup.merge_head.exists()
    assert not setup.state.exists()
    assert head(setup.wt) == setup.old_head
    assert "merge" not in outcome.recorder.git_subs()


def assert_left_in_progress(setup: Setup) -> None:
    assert setup.merge_head.exists()
    assert setup.state.exists()
    assert head(setup.wt) == setup.old_head


def refusal(outcome: Outcome) -> str:
    """Exit 2 and the reason as the last line on stderr."""
    assert outcome.code == 2, outcome.text
    assert outcome.err_lines(), "a refusal writes its reason on stderr"
    return outcome.err_lines()[-1]


def gates_called(outcome: Outcome) -> list[str]:
    return [g for g in (gate_name(a) for a in outcome.recorder.argvs()) if g is not None]


def builds(outcome: Outcome) -> list[list[str]]:
    return [a for a in outcome.recorder.argvs() if program(a) == "cmake"]


def same(a: str | Path, b: Path) -> bool:
    return Path(a).resolve() == b.resolve()


# --- 1. refusals of the start form: none starts a merge -------------------------------------


def test_a_directory_inside_the_worktree_is_refused(setup: Setup) -> None:
    outcome = setup.run([str(setup.wt / "tools"), "--persona", "developer", "--trailer", TRAILER])
    line = refusal(outcome)
    assert line.startswith("merge_master: ")
    assert line.endswith(" is not the top of a checkout of this repository")
    assert_no_merge_started(setup, outcome)


def test_a_checkout_of_another_repository_is_refused(setup: Setup, tmp_path: Path) -> None:
    stranger = tmp_path / "link" / "stranger"
    stranger.mkdir()
    git(stranger, "init", "-q", "-b", "side")
    git(stranger, "commit", "-q", "--allow-empty", "-m", "elsewhere")
    stub_venv(stranger)
    outcome = setup.run([str(stranger), "--persona", "developer", "--trailer", TRAILER])
    line = refusal(outcome)
    assert line.endswith(" is not the top of a checkout of this repository")
    assert not (stranger / ".git" / "MERGE_HEAD").exists()
    assert not (stranger / ".git" / "merge_master.json").exists()
    assert "merge" not in outcome.recorder.git_subs()


def test_a_detached_head_is_refused(setup: Setup) -> None:
    git(setup.wt, "checkout", "-q", "--detach")
    outcome = setup.start()
    line = refusal(outcome)
    assert line.startswith("merge_master: ") and "detached HEAD" in line
    assert line.endswith("; run it on a worktree's own branch")
    assert_no_merge_started(setup, outcome)


def test_the_master_branch_is_refused(setup: Setup) -> None:
    stub_venv(setup.repos.main)
    bring_in(setup, (12, {"notes.txt": "pr 12\n"}))
    main_head = head(setup.repos.main)
    outcome = setup.run([str(setup.repos.main), "--persona", "developer", "--trailer", TRAILER])
    line = refusal(outcome)
    assert line.startswith("merge_master: ") and " is on master" in line
    assert line.endswith("; run it on a worktree's own branch")
    assert head(setup.repos.main) == main_head
    assert not (setup.repos.main / ".git" / "MERGE_HEAD").exists()
    assert not (setup.repos.main / ".git" / "merge_master.json").exists()


def test_a_worktree_without_a_venv_is_refused_naming_new_worktree(setup: Setup) -> None:
    (setup.wt / ".venv" / "bin" / "python").unlink()
    outcome = setup.start()
    line = refusal(outcome)
    assert line.startswith("merge_master: ")
    assert " has no .venv; run: python3 tools/new_worktree.py --existing " in line
    assert same(line.rsplit("--existing ", 1)[1], setup.wt)
    assert_no_merge_started(setup, outcome)


def test_uncommitted_changes_are_refused_naming_the_file_and_the_stash(setup: Setup) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    (setup.wt / "notes.txt").write_text("not committed\n")
    outcome = setup.start()
    line = refusal(outcome)
    assert line.startswith("merge_master: ")
    assert " has uncommitted changes (notes.txt); commit them first. " in line
    assert line.endswith("Do not use git stash: every worktree shares one stash")
    assert_no_merge_started(setup, outcome)
    assert (setup.wt / "notes.txt").read_text() == "not committed\n"


def test_a_merge_already_in_progress_is_refused(setup: Setup) -> None:
    # The branch and master add the same file, so the hand-started merge
    # changes nothing: the tree is clean and only MERGE_HEAD shows the merge.
    commit(setup.wt, {"same.txt": "same\n"}, "the branch adds same.txt")
    setup.old_head = head(setup.wt)
    bring_in(setup, (12, {"same.txt": "same\n"}))
    git(setup.wt, "fetch", "-q", "origin", "master")
    git(setup.wt, "merge", "-q", "--no-ff", "--no-commit", "origin/master")
    assert git(setup.wt, "status", "--porcelain") == ""
    before = setup.merge_head.read_text()
    outcome = setup.start()
    line = refusal(outcome)
    assert line.startswith("merge_master: a merge is already in progress in ")
    assert line.endswith("; finish it with --continue or abandon it with --abort")
    assert setup.merge_head.read_text() == before
    assert not setup.state.exists()
    assert head(setup.wt) == setup.old_head


@pytest.mark.parametrize("persona", ["reviewer", "nobody"])
def test_a_persona_that_is_not_a_writer_is_refused(setup: Setup, persona: str) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    outcome = setup.run([str(setup.wt), "--persona", persona, "--trailer", TRAILER])
    line = refusal(outcome)
    prefix = f"merge_master: --persona {persona}: give one of "
    assert line.startswith(prefix)
    named = {n.strip() for n in line.removeprefix(prefix).split(",")}
    assert named == {"developer", "tester", "architect", "perf", "orchestrator"}
    assert_no_merge_started(setup, outcome)


@pytest.mark.parametrize(
    "trailer",
    ["Co-Authored-By: nobody", "Signed-off-by: A <a@b.c>", f"{TRAILER} extra", ""],
)
def test_a_trailer_that_is_not_a_co_authored_by_line_is_refused(setup: Setup, trailer: str) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    outcome = setup.run([str(setup.wt), "--persona", "developer", "--trailer", trailer])
    assert refusal(outcome) == (
        "merge_master: --trailer must be the Co-Authored-By line from your system context"
    )
    assert_no_merge_started(setup, outcome)


@pytest.mark.parametrize("missing", ["--persona", "--trailer"])
def test_the_first_form_without_persona_or_trailer_is_refused_up_front(
    setup: Setup, missing: str
) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    given = {"--persona": "developer", "--trailer": TRAILER}
    del given[missing]
    outcome = setup.run([str(setup.wt), *[w for kv in given.items() for w in kv]])
    assert refusal(outcome) == missing_option_message(missing)
    assert_no_merge_started(setup, outcome)


def missing_option_message(missing: str) -> str:
    """§4.6's row for a first form or `--continue` without `--persona` or `--trailer`."""
    return (
        f"merge_master: {missing} is missing;"
        " starting or continuing a merge needs --persona and --trailer"
    )


def test_a_failed_fetch_is_refused_as_blocked_on_network(setup: Setup) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    outcome = setup.start(fail_fetch=True)
    line = refusal(outcome)
    assert line.startswith("merge_master: could not fetch origin/master (")
    assert "h18 test network down" in line
    assert line.endswith("); blocked on network: stop and hand back")
    assert_no_merge_started(setup, outcome)


def test_no_build_on_the_first_form_refuses_cpp_from_master_before_merging(setup: Setup) -> None:
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start("--no-build")
    assert refusal(outcome) == (
        "merge_master: origin/master changed C++ (include/core.h); the suite needs a rebuilt"
        " _core and this run may not build C++. No merge was started; hand back"
    )
    assert_no_merge_started(setup, outcome)
    assert builds(outcome) == []
    assert gates_called(outcome) == []


# --- 1. refusals of --continue -------------------------------------------------------------


def test_continue_without_a_merge_this_tool_started_is_refused(setup: Setup) -> None:
    outcome = setup.resume()
    line = refusal(outcome)
    assert line.startswith("merge_master: no merge started by this tool in ")
    assert line.endswith("; start one with: python3 tools/merge_master.py " + line.split()[-1])
    assert same(line.split()[-1], setup.wt)
    assert head(setup.wt) == setup.old_head


@pytest.mark.parametrize("missing", ["--persona", "--trailer"])
def test_continue_without_persona_or_trailer_is_refused_and_leaves_the_merge(
    setup: Setup, missing: str
) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    given = {"--persona": "developer", "--trailer": TRAILER}
    del given[missing]
    outcome = setup.run([str(setup.wt), "--continue", *[w for kv in given.items() for w in kv]])
    assert refusal(outcome) == missing_option_message(missing)
    assert_left_in_progress(setup)
    assert gates_called(outcome) == []
    assert "commit" not in outcome.recorder.git_subs()


def test_continue_with_unmerged_paths_is_refused_naming_them(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    outcome = setup.resume()
    line = refusal(outcome)
    assert line.startswith("merge_master: still unresolved: ")
    assert line.endswith("; resolve with Edit, then git add")
    assert "a.py" in line and "CLAUDE.md" in line
    assert_left_in_progress(setup)
    assert gates_called(outcome) == []


def test_continue_with_a_tracked_change_not_added_is_refused(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    (setup.wt / "notes.txt").write_text("changed, not added\n")
    outcome = setup.resume()
    assert refusal(outcome) == (
        "merge_master: notes.txt is changed but not added; git add it if it belongs to the merge"
    )
    assert_left_in_progress(setup)
    assert gates_called(outcome) == []


def test_continue_with_a_wrong_merge_head_is_refused(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    setup.merge_head.write_text(setup.old_head + "\n")
    outcome = setup.resume()
    line = refusal(outcome)
    assert line.startswith("merge_master: MERGE_HEAD is ")
    assert line.endswith("; stop and hand back")
    found, recorded = line.removeprefix("merge_master: MERGE_HEAD is ").split(
        ", but this tool merged "
    )
    assert len(found.strip()) >= 7 and setup.old_head.startswith(found.strip())
    assert setup.origin_master().startswith(recorded.removesuffix("; stop and hand back").strip())
    assert head(setup.wt) == setup.old_head
    assert gates_called(outcome) == []


# --- 2. up to date ---------------------------------------------------------------------------


def test_a_branch_that_already_contains_master_is_left_alone(setup: Setup) -> None:
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    assert "already contains origin/master" in outcome.out
    assert "nothing to merge" in outcome.out
    assert head(setup.wt) == setup.old_head
    assert not setup.state.exists() and not setup.merge_head.exists()
    assert gates_called(outcome) == []


# --- 3. a clean merge ------------------------------------------------------------------------


def test_a_clean_merge_runs_the_gates_in_order_in_the_worktree_venv_and_commits(
    setup: Setup,
) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}), (13, {"notes.txt": "pr 13\n"}))
    outcome = setup.start()
    assert outcome.code == 0, outcome.text

    assert gates_called(outcome) == [*FAST_GATES, "pytest"]
    py = setup.wt / ".venv" / "bin" / "python"
    for call in outcome.recorder.calls:
        if gate_name(call.argv) is not None:
            assert same(call.argv[0], py), call.argv
            assert call.cwd is not None and same(call.cwd, setup.wt), call
    citations = next(a for a in outcome.recorder.argvs() if gate_name(a) == "citations")
    assert citations[-2:] == ["--base", "origin/master"]
    assert builds(outcome) == []

    parents = git(setup.wt, "rev-list", "--parents", "-n", "1", "HEAD").split()[1:]
    assert parents == [setup.old_head, setup.origin_master()]
    subject = setup.subject()
    assert subject.startswith("Merge origin/master ") and " into worktree-x: brings in " in subject
    assert "#12" in subject and "#13" in subject
    assert subject.endswith("(@developer)")
    assert setup.body()[-1] == TRAILER
    assert not setup.state.exists() and not setup.message_file.exists()
    assert not setup.merge_head.exists()


def test_the_merge_list_is_git_log_first_parent_merges_since_the_base(setup: Setup) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}), (13, {"notes.txt": "pr 13\n"}))
    git(setup.wt, "fetch", "-q", "origin", "master")  # the test's own copy of the commits
    base = git(setup.wt, "merge-base", "HEAD", setup.origin_master()).strip()
    expected = git(setup.wt, "log", "--merges", "--first-parent", "--oneline",
                   f"{base}..{setup.origin_master()}").splitlines()  # fmt: skip
    assert len(expected) == 2
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    body = [line.strip() for line in setup.body()]
    for line in expected:
        assert line in body, f"{line!r} missing from the commit body"
        assert line in outcome.out, f"{line!r} was not printed"
    assert not any(line.startswith("#") for line in setup.body()[1:])
    assert "Conflicts: none" in body
    assert "_core rebuilt: no, no C++ came in." in body, body
    assert not any(line.startswith("Changed by hand beyond git's merge:") for line in body)


def test_more_than_six_prs_name_six_and_count_the_rest(setup: Setup) -> None:
    bring_in(setup, *[(n, {f"n{n}.txt": f"{n}\n"}) for n in range(20, 28)])
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    subject = setup.subject()
    assert len(re.findall(r"#\d+", subject)) == 6
    assert "and 2 more" in subject


def test_commits_without_pr_subjects_are_counted(setup: Setup) -> None:
    setup.repos.master_commit({"n1.txt": "1\n"}, "a plain commit")
    setup.repos.master_commit({"n2.txt": "2\n"}, "another plain commit")
    setup.repos.publish()
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    assert "brings in 2 commits" in setup.subject()


# --- 4. a conflict stops the tool ----------------------------------------------------------


def test_a_conflict_stops_with_exit_3_and_resolves_nothing(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    outcome = setup.start()
    assert outcome.code == 3, outcome.text
    assert head(setup.wt) == setup.old_head
    assert_left_in_progress(setup)
    assert "merge_master: stopped for conflict resolution in 2 files:" in outcome.text
    listed = {line.split()[0]: line for line in outcome.text.splitlines()
              if line.startswith("  ") and line.split()}  # fmt: skip
    assert "a.py" in listed and "CLAUDE.md" in listed
    assert "states rules" in listed["CLAUDE.md"]
    assert "states rules" not in listed["a.py"]
    assert "--continue" in outcome.text and "--abort" in outcome.text
    assert (setup.wt / "a.py").read_text().startswith("<<<<<<< ")
    assert "<<<<<<< " in (setup.wt / "CLAUDE.md").read_text()
    unmerged = git(setup.wt, "diff", "--name-only", "--diff-filter=U").split()
    assert sorted(unmerged) == ["CLAUDE.md", "a.py"]
    assert gates_called(outcome) == []
    assert "add" not in outcome.recorder.git_subs()


def test_a_conflict_against_a_master_newer_than_recorded_is_refused_not_stopped(
    setup: Setup,
) -> None:
    """Code review round 1, suggestion "merge-head-check-at-stop" (§10 item 8):
    another worktree's fetch moves the shared origin/master between the tool's
    `rev-parse` and its `git merge`, and the merge conflicts. The tool must
    refuse with §4.6's MERGE_HEAD message before inviting a resolution that
    --continue would then refuse; the merge is left in progress for --abort."""
    conflict_on_a_py_and_claude_md(setup)
    recorded = setup.origin_master()
    moved: list[str] = []

    def fetch_elsewhere(argv: list[str]) -> Answer | None:
        if git_sub(argv) == "merge" and not moved:
            setup.repos.master_commit({"notes.txt": "master moved again\n"}, "a later master")
            setup.repos.publish()
            git(setup.repos.main, "fetch", "-q", "origin", "master")
            moved.append(setup.origin_master())
        return None

    tool: ModuleType = load(setup.repos.main, TOOL)
    recorder = Recorder(MergeTable(), fetch_elsewhere)
    args = [str(setup.wt), "--persona", "developer", "--trailer", TRAILER]
    outcome = invoke(tool, args, recorder, setup.capsys)
    assert moved and moved[0] != recorded, "the probe never moved origin/master"
    assert setup.merge_head.read_text().strip() == moved[0], "git merged the moved master"
    line = refusal(outcome)
    assert line.startswith("merge_master: MERGE_HEAD is ")
    assert line.endswith("; stop and hand back")
    found, merged = line.removeprefix("merge_master: MERGE_HEAD is ").split(
        ", but this tool merged "
    )
    assert moved[0].startswith(found.strip())
    assert recorded.startswith(merged.removesuffix("; stop and hand back").strip())
    assert "stopped for conflict resolution" not in outcome.text
    assert setup.merge_head.exists(), "the merge is left in progress for --abort"
    assert head(setup.wt) == setup.old_head
    assert gates_called(outcome) == []


# --- 5. --continue after a conflict --------------------------------------------------------


def test_continue_refuses_a_conflict_marker_naming_file_and_line(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    (setup.wt / "CLAUDE.md").write_text("# rules\n\nboth lines\n")
    marked = (setup.wt / "a.py").read_text().splitlines()
    first = next(i for i, line in enumerate(marked, 1) if line.startswith("<<<<<<< "))
    git(setup.wt, "add", "a.py", "CLAUDE.md")
    outcome = setup.resume()
    assert refusal(outcome) == f"merge_master: a.py:{first} still holds a conflict marker"
    assert_left_in_progress(setup)
    assert gates_called(outcome) == []


def test_continue_commits_and_lists_conflicts_and_hand_changes(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    (setup.wt / "b.py").write_text("b = 3  # changed by hand after the merge\n")
    git(setup.wt, "add", "b.py")
    outcome = setup.resume()
    assert outcome.code == 0, outcome.text
    body = [line.strip() for line in setup.body()]
    conflicts = next(line for line in body if line.startswith("Conflicts resolved by hand: "))
    assert set(conflicts.removeprefix("Conflicts resolved by hand: ").split(", ")) == {
        "a.py",
        "CLAUDE.md",
    }
    assert "Changed by hand beyond git's merge: b.py" in body
    parents = git(setup.wt, "rev-list", "--parents", "-n", "1", "HEAD").split()[1:]
    assert parents == [setup.old_head, setup.origin_master()]
    assert (setup.wt / "a.py").read_text() == "x = 'both'\ny = 2\n"
    assert not setup.state.exists() and not setup.message_file.exists()


def test_continue_without_a_hand_change_has_no_hand_change_line(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    outcome = setup.resume()
    assert outcome.code == 0, outcome.text
    body = [line.strip() for line in setup.body()]
    assert any(line.startswith("Conflicts resolved by hand: ") for line in body)
    assert not any(line.startswith("Changed by hand beyond git's merge:") for line in body)


def test_a_modify_delete_conflict_resolved_by_git_rm_commits(setup: Setup) -> None:
    commit(setup.wt, {"gone.py": "g = 2  # the branch changed it\n"}, "branch edits gone.py")
    setup.old_head = head(setup.wt)
    bring_in(setup, (12, {"gone.py": None}))
    stopped = setup.start()
    assert stopped.code == 3, stopped.text
    listed = [line.split()[0] for line in stopped.text.splitlines()
              if line.startswith("  ") and line.split()]  # fmt: skip
    assert "gone.py" in listed
    git(setup.wt, "rm", "-q", "gone.py")
    outcome = setup.resume()
    assert outcome.code == 0, outcome.text
    body = [line.strip() for line in setup.body()]
    conflicts = next(line for line in body if line.startswith("Conflicts resolved by hand: "))
    assert "gone.py" in conflicts
    assert "gone.py" not in git(setup.wt, "ls-tree", "--name-only", "HEAD").split()


# --- 6. red gates ------------------------------------------------------------------------------


def test_a_red_fast_gate_runs_the_rest_skips_pytest_and_commits_nothing(setup: Setup) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    outcome = setup.start(table=MergeTable(red={"ruff_check"}))
    assert outcome.code == 4, outcome.text
    assert gates_called(outcome) == list(FAST_GATES)
    assert gate_output("ruff_check", True) in outcome.text, "the gate's output, unchanged"
    assert head(setup.wt) == setup.old_head
    assert setup.merge_head.exists()
    assert "commit" not in outcome.recorder.git_subs()


def test_a_red_pytest_commits_nothing(setup: Setup) -> None:
    bring_in(setup, (12, {"b.py": "b = 2\n"}))
    outcome = setup.start(table=MergeTable(red={"pytest"}))
    assert outcome.code == 4, outcome.text
    assert gates_called(outcome) == [*FAST_GATES, "pytest"]
    assert head(setup.wt) == setup.old_head
    assert setup.merge_head.exists()
    assert "commit" not in outcome.recorder.git_subs()


# --- 7. C++ from master ------------------------------------------------------------------------


def put_core(setup: Setup, *names: str) -> list[Path]:
    """Dummy `_core` files in build-pyext, dated 2001 so a copy's touch is visible."""
    made = []
    for name in names:
        path = setup.wt / "build-pyext" / name
        path.write_bytes(CORE_BYTES)
        os.utime(path, (1e9, 1e9))
        made.append(path)
    return made


def site_package(setup: Setup) -> Path:
    return setup.wt / ".venv" / "lib" / "python3.12" / "site-packages" / "tin_engine"


def test_cpp_from_master_rebuilds_core_between_the_fast_gates_and_pytest(setup: Setup) -> None:
    put_core(setup, f"_core{SUFFIX}")
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    order = [gate_name(a) or ("build" if program(a) == "cmake" else None)
             for a in outcome.recorder.argvs()]  # fmt: skip
    assert [o for o in order if o is not None] == [*FAST_GATES, "build", "pytest"]
    (build,) = builds(outcome)
    assert build[:2] == ["cmake", "--build"] and same(build[2], setup.wt / "build-pyext")
    assert build[3:] == ["-j", "--target", "_core"]
    copied = site_package(setup) / f"_core{SUFFIX}"
    assert copied.read_bytes() == CORE_BYTES
    assert copied.stat().st_mtime > 1e9 + 1, "the copy's time stamp is touched"
    body = [line.strip() for line in setup.body()]
    assert "_core rebuilt: yes." in body, body


def test_no_build_on_continue_refuses_cpp_and_leaves_the_merge(setup: Setup) -> None:
    put_core(setup, f"_core{SUFFIX}")
    conflict_on_a_py_and_claude_md(setup)
    bring_in(setup, (13, {"include/core.h": "#pragma once\nint f();\n"}))
    assert setup.start().code == 3
    resolve_conflicts(setup)
    outcome = setup.resume("--no-build")
    assert refusal(outcome) == (
        "merge_master: origin/master changed C++ (include/core.h); the suite needs a rebuilt"
        " _core and this run may not build C++. The merge is left in progress; hand back"
    )
    assert builds(outcome) == []
    assert gates_called(outcome) == []
    assert_left_in_progress(setup)


def test_a_failed_build_is_exit_4_without_pytest_or_commit(setup: Setup) -> None:
    put_core(setup, f"_core{SUFFIX}")
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start(table=MergeTable(build_rc=2))
    assert outcome.code == 4, outcome.text
    assert "pytest" not in gates_called(outcome)
    assert_left_in_progress(setup)
    assert "commit" not in outcome.recorder.git_subs()


@pytest.mark.parametrize(
    "names",
    [(f"_core{SUFFIX}", "_core.cpython-311-other.so"), ("_core.cpython-311-other.so",), ()],
    ids=["two", "other-suffix", "none"],
)
def test_build_pyext_must_hold_one_core_with_the_venv_suffix(
    setup: Setup, names: tuple[str, ...]
) -> None:
    put_core(setup, *names)
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start()
    line = refusal(outcome)
    assert line.startswith("merge_master: ") and "build-pyext holds " in line
    assert f", not one _core{SUFFIX}; remove build-pyext and run: " in line
    assert line.endswith("The merge is left in progress")
    assert list(site_package(setup).iterdir()) == [], "nothing copied"
    assert "pytest" not in gates_called(outcome)
    assert_left_in_progress(setup)


def assert_no_build_pyext_refusal(setup: Setup, line: str, ending: str) -> None:
    """§4.6's row "C++ came in and no `build-pyext`", the worktree named twice."""
    prefix = "merge_master: origin/master changed C++ (include/core.h) and "
    middle = (
        " has no build-pyext to rebuild _core in; run: python3 tools/new_worktree.py --existing "
    )
    assert line.startswith(prefix) and middle in line, line
    assert line.endswith(f". {ending}"), line
    named, again = line.removeprefix(prefix).removesuffix(f". {ending}").split(middle)
    assert same(named, setup.wt) and same(again, setup.wt), line


def test_cpp_from_master_without_build_pyext_is_refused_naming_new_worktree(
    setup: Setup,
) -> None:
    (setup.wt / "build-pyext").rmdir()
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start()
    assert_no_build_pyext_refusal(setup, refusal(outcome), "No merge was started")
    assert_no_merge_started(setup, outcome)
    assert builds(outcome) == []
    assert gates_called(outcome) == []


def test_continue_without_build_pyext_refuses_cpp_and_leaves_the_merge(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    bring_in(setup, (13, {"include/core.h": "#pragma once\nint f();\n"}))
    assert setup.start().code == 3
    resolve_conflicts(setup)
    (setup.wt / "build-pyext").rmdir()
    outcome = setup.resume()
    assert_no_build_pyext_refusal(setup, refusal(outcome), "The merge is left in progress")
    assert builds(outcome) == []
    assert gates_called(outcome) == []
    assert_left_in_progress(setup)


def test_no_build_and_no_build_pyext_on_continue_is_the_no_build_refusal(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    bring_in(setup, (13, {"include/core.h": "#pragma once\nint f();\n"}))
    assert setup.start().code == 3
    resolve_conflicts(setup)
    (setup.wt / "build-pyext").rmdir()
    outcome = setup.resume("--no-build")
    assert refusal(outcome) == (
        "merge_master: origin/master changed C++ (include/core.h); the suite needs a rebuilt"
        " _core and this run may not build C++. The merge is left in progress; hand back"
    )
    assert builds(outcome) == []
    assert_left_in_progress(setup)


def test_no_build_and_no_build_pyext_on_the_first_form_is_the_no_build_refusal(
    setup: Setup,
) -> None:
    (setup.wt / "build-pyext").rmdir()
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start("--no-build")
    assert refusal(outcome) == (
        "merge_master: origin/master changed C++ (include/core.h); the suite needs a rebuilt"
        " _core and this run may not build C++. No merge was started; hand back"
    )
    assert_no_merge_started(setup, outcome)
    assert builds(outcome) == []


# --- 8. MERGE_HEAD gone ------------------------------------------------------------------------


def test_continue_after_merge_head_is_gone_is_refused_and_commits_nothing(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    resolve_conflicts(setup)
    setup.merge_head.unlink()
    outcome = setup.resume()
    line = refusal(outcome)
    assert line.startswith("merge_master: the merge in ")
    assert " was left (MERGE_HEAD is gone: a checkout, reset or commit during it). " in line
    assert line.endswith("Nothing was committed by this tool; stop and hand back")
    assert head(setup.wt) == setup.old_head
    assert "commit" not in outcome.recorder.git_subs()


# --- 9. what the tool runs -----------------------------------------------------------------------

FORBIDDEN = {"stash", "reset", "rebase", "checkout", "restore", "push", "add", "rm"}


def test_the_tool_runs_nothing_guarded_and_commits_once_with_a_file(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    stopped = setup.start()
    assert stopped.code == 3
    resolve_conflicts(setup)
    done = setup.resume()
    assert done.code == 0, done.text
    argvs = stopped.recorder.argvs() + done.recorder.argvs()
    assert push_guard_findings(argvs) == []
    assert FORBIDDEN.isdisjoint(git_sub(a) for a in argvs), argvs
    commits = [a for a in argvs if git_sub(a) == "commit"]
    assert len(commits) == 1
    assert "-F" in commits[0] and "-m" not in commits[0]
    assert "--amend" not in commits[0] and "--no-verify" not in commits[0]


def test_a_clean_merge_runs_nothing_guarded_either(setup: Setup) -> None:
    put_core(setup, f"_core{SUFFIX}")
    bring_in(setup, (12, {"include/core.h": "#pragma once\nint f();\n"}))
    outcome = setup.start()
    assert outcome.code == 0, outcome.text
    argvs = outcome.recorder.argvs()
    assert push_guard_findings(argvs) == []
    assert FORBIDDEN.isdisjoint(git_sub(a) for a in argvs), argvs
    assert [git_sub(a) for a in argvs].count("commit") == 1


# --- 10. --abort --------------------------------------------------------------------------------


def test_abort_returns_the_tree_to_the_old_head_and_removes_the_state(setup: Setup) -> None:
    conflict_on_a_py_and_claude_md(setup)
    assert setup.start().code == 3
    outcome = setup.run([str(setup.wt), "--abort"])
    assert outcome.code == 0, outcome.text
    assert not setup.merge_head.exists()
    assert not setup.state.exists() and not setup.message_file.exists()
    assert head(setup.wt) == setup.old_head
    assert git(setup.wt, "status", "--porcelain") == ""
    assert (setup.wt / "a.py").read_text() == "x = 'branch'\ny = 2\n"
    assert push_guard_findings(outcome.recorder.argvs()) == []


def test_abort_with_no_merge_in_progress_is_refused(setup: Setup) -> None:
    outcome = setup.run([str(setup.wt), "--abort"])
    line = refusal(outcome)
    prefix, suffix = "merge_master: no merge in progress in ", "; nothing to abort"
    assert line.startswith(prefix) and line.endswith(suffix), line
    assert same(line.removeprefix(prefix).removesuffix(suffix), setup.wt), line
    assert head(setup.wt) == setup.old_head
    assert "merge" not in outcome.recorder.git_subs()


def test_the_default_run_answers_a_missing_program_with_127(setup: Setup, tmp_path: Path) -> None:
    tool = load(setup.repos.main, TOOL)
    done = tool.run(["rasputin-h18-no-such-program"], tmp_path)
    assert done.returncode == 127
    assert "not found" in done.stderr
