"""The after-commit gate checks the tree that was committed (h19).

Spec: `docs/increments/h19-commit-gate-tree.md` §5. The hook,
`.claude/hooks/gates_after_commit.py`, runs after every Bash call. A subagent's
shell returns to the session directory (the main checkout) between calls, so a
subagent works in a worktree only through `cd <wt> && ...` or `git -C <wt> ...`
inside the command. The hook must follow those to the committed tree, fire on
any committing git command (and on nothing else), and report a directory it
cannot work out as NOT CHECKED (exit 2) instead of gating some other tree.

The fixture is a temporary main repository with a linked worktree. The five
gates are stubs that append `<gate> <cwd>` to a log outside both trees, so
"which gates ran, where" is read from the log, never inferred from an exit
status. The worktree's detria stub is red and prints `RED IN <cwd>`. Every
event is subagent-shaped (`agent_id` set, `cwd` the main repository) unless a
test says otherwise. No network.

The pure functions of §3.5, `git_dirs` and `committed_trees`, are loaded from
the hook's path; a hook without them fails the test that needs them, naming
the function.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from harness_fixtures import (
    REAL,
    SUBAGENT,
    add_worktree,
    bash_event,
    clean_env,
    git,
    is_work_tree_top,
    make_repo,
    run_script,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

HOOK = ".claude/hooks/gates_after_commit.py"
SHELL_SCAN = "tools/shell_scan.py"
CHECKS = ("check_prohibited_deps", "check_detria_boundary", "check_citations")
CWD = Path("/c")


# --------------------------------------------------------------------------- fixture


def green_stub(name: str, log: Path) -> str:
    return (
        "import os\n"
        f"with open({str(log)!r}, 'a') as f:\n"
        f"    f.write({name!r} + ' ' + os.getcwd() + '\\n')\n"
    )


def red_stub(name: str, log: Path) -> str:
    return green_stub(name, log) + "print('RED IN ' + os.getcwd())\nraise SystemExit(1)\n"


def ruff_stub(log: Path) -> str:
    return f'#!/bin/sh\necho "ruff $(pwd -P)" >> {shlex.quote(str(log))}\nexit 0\n'


class Repos:
    """A main repository, a linked worktree of it, and the gates' log."""

    def __init__(self, root: Path) -> None:
        self.log = root / "gates.log"
        self.main = make_repo(root / "main")
        for relative in (HOOK, SHELL_SCAN):
            target = self.main / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text((REAL / relative).read_text(encoding="utf-8"))
        for name in CHECKS:
            (self.main / "tools" / f"{name}.py").write_text(green_stub(name, self.log))
        ruff = self.main / ".venv" / "bin" / "ruff"
        ruff.parent.mkdir(parents=True)
        ruff.write_text(ruff_stub(self.log))
        ruff.chmod(0o755)
        git(self.main, "add", "-A")
        git(self.main, "commit", "-q", "-m", "hook, stub gates, stub ruff")
        self.wt = add_worktree(self.main, root / "wt", "work")
        self.make_red(self.wt)
        git(self.wt, "commit", "-q", "-am", "red detria stub in the worktree")

    def make_red(self, tree: Path) -> None:
        name = "check_detria_boundary"
        (tree / "tools" / f"{name}.py").write_text(red_stub(name, self.log))

    def ran(self) -> list[tuple[str, Path]]:
        """(gate, directory) for every gate run so far, in order."""
        if not self.log.exists():
            return []
        lines = self.log.read_text().splitlines()
        return [(gate, Path(where)) for gate, where in (line.split(" ", 1) for line in lines)]

    def ran_in(self, tree: Path) -> list[str]:
        return [gate for gate, where in self.ran() if where == tree]

    def hook(self, command: str, **extra: Any) -> subprocess.CompletedProcess[str]:
        """Run the main repository's hook copy on a PostToolUse Bash event."""
        event = bash_event(self.main, command, hook_event_name="PostToolUse", **extra)
        return run_script(self.main, HOOK, event)

    def subagent(self, command: str) -> subprocess.CompletedProcess[str]:
        return self.hook(command, **SUBAGENT)


@pytest.fixture
def repos(tmp_path: Path) -> Repos:
    made = Repos(tmp_path.resolve())
    # Premises: two real work trees, both stubs in place, the worktree's detria
    # stub red and the main one green, and no gate has run yet.
    assert is_work_tree_top(made.main) and is_work_tree_top(made.wt)
    assert "RED IN" in (made.wt / "tools" / "check_detria_boundary.py").read_text()
    assert "RED IN" not in (made.main / "tools" / "check_detria_boundary.py").read_text()
    assert os.access(made.wt / ".venv" / "bin" / "ruff", os.X_OK)
    assert made.ran() == []
    return made


def load_hook(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"h19_hook_{abs(hash(path))}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def function(module: ModuleType, name: str) -> Any:
    if not hasattr(module, name):
        pytest.fail(f"{HOOK} has no {name}() (h19 §3.5)")
    return getattr(module, name)


@pytest.fixture(scope="module")
def git_dirs() -> Any:
    return function(load_hook(REAL / HOOK), "git_dirs")


# --------------------------------------------------------------------------- 1, 2


def run_in(tree: Path, *args: str) -> None:
    """Really run the git command the event describes, so the event is honest."""
    subprocess.run(
        ["git", "-C", str(tree), *args], check=True, capture_output=True, env=clean_env()
    )


def assert_red_in(result: subprocess.CompletedProcess[str], tree: Path) -> None:
    assert result.returncode == 2, f"exit {result.returncode}; stderr: {result.stderr}"
    assert f"RED IN {tree}" in result.stderr
    # Once in the stub's own output, and at least once more where the hook
    # heads the failure with the tree it concerns (§3.4).
    assert result.stderr.count(str(tree)) >= 2, result.stderr


def test_cd_into_the_worktree_gates_the_worktree(repos: Repos) -> None:
    run_in(repos.wt, "commit", "-q", "--allow-empty", "-m", "x")
    result = repos.subagent(f"cd {shlex.quote(str(repos.wt))} && git commit -q --allow-empty -m x")
    assert_red_in(result, repos.wt)
    assert repos.ran_in(repos.main) == [], "the main checkout was gated in the worktree's place"


def test_git_dash_c_into_the_worktree_gates_the_worktree(repos: Repos) -> None:
    run_in(repos.wt, "commit", "-q", "--allow-empty", "-m", "x")
    result = repos.subagent(f"git -C {shlex.quote(str(repos.wt))} commit -q --allow-empty -m x")
    assert_red_in(result, repos.wt)
    assert repos.ran_in(repos.main) == []


def test_merge_with_a_config_option_after_cd_gates_the_worktree(repos: Repos) -> None:
    run_in(repos.wt, "checkout", "-q", "-b", "y")
    run_in(repos.wt, "commit", "-q", "--allow-empty", "-m", "on y")
    run_in(repos.wt, "checkout", "-q", "work")
    run_in(repos.wt, "-c", "user.name=x", "merge", "--no-ff", "--no-edit", "y")
    assert "on y" in git(repos.wt, "log", "-3", "--format=%s"), "premise: y was merged"
    command = f"cd {shlex.quote(str(repos.wt))}; git -c user.name=x merge --no-ff --no-edit y"
    result = repos.subagent(command)
    assert_red_in(result, repos.wt)
    assert repos.ran_in(repos.main) == []


def test_main_session_commit_in_the_main_checkout_stays_gated(repos: Repos) -> None:
    run_in(repos.main, "commit", "-q", "--allow-empty", "-m", "x")
    green = repos.hook("git commit -m x")  # main-session shape: no agent_id
    assert green.returncode == 0, green.stderr
    assert sorted(repos.ran_in(repos.main)) == sorted([*CHECKS, "ruff", "ruff"])

    repos.make_red(repos.main)
    red = repos.hook("git commit -m x")
    assert_red_in(red, repos.main)
    assert repos.ran_in(repos.wt) == []


# --------------------------------------------------------------------------- 3, 4


TRIGGERS = [
    ("git -C /a commit", [Path("/a")]),
    ("git -c a=b commit", [CWD]),
    ("git --no-pager merge x", [CWD]),
    ("/usr/bin/git commit", [CWD]),
    ("env X=1 git commit", [CWD]),
    ("timeout 9 git commit", [CWD]),
    ('bash -c "cd /a && git commit"', [Path("/a")]),
    ("git pull", [CWD]),
    ("git cherry-pick abc", [CWD]),
    ("git revert abc", [CWD]),
    ("git am p.patch", [CWD]),
    ("git rebase master", [CWD]),
    ('echo "git commit"', []),
    ('grep -n "git merge" f', []),
    ('git log --grep="git commit"', []),
    ("git merge-base a b", []),
    ("git commit-tree T", []),
    ("git status", []),
]


@pytest.mark.parametrize(("command", "expected"), TRIGGERS, ids=[c for c, _ in TRIGGERS])
def test_which_commands_trigger(git_dirs: Any, command: str, expected: list[Path]) -> None:
    assert git_dirs(command, CWD) == expected


DIRECTORIES = [
    ("cd /a && git commit", [Path("/a")]),
    ("cd /a && cd b && git commit", [Path("/a/b")]),
    ("git -C /a -C b commit", [Path("/a/b")]),
    ("cd /a && git -C b commit", [Path("/a/b")]),
    ("cd b && git commit", [CWD / "b"]),
    ("cd ~/x && git commit", [Path(os.path.expanduser("~/x"))]),
    ("cd /a && git commit && cd /b && git commit", [Path("/a"), Path("/b")]),
    ('cd "$WT" && git commit', ["$WT"]),
    ("cd - && git commit", ["-"]),
    ("git --work-tree=/a commit", ["--work-tree=/a"]),
    # A relative step from an unknown directory stays unknown, holding the first unknown word.
    ("cd - && cd b && git commit", ["-"]),
    ("cd - && git -C b commit", ["-"]),
    ("cd - && cd /a && git commit", [Path("/a")]),
]


@pytest.mark.parametrize(("command", "expected"), DIRECTORIES, ids=[c for c, _ in DIRECTORIES])
def test_which_directory(git_dirs: Any, command: str, expected: list[Path | str]) -> None:
    found = git_dirs(command, CWD)
    assert found == expected
    # Path("/a") != "/a", but say which kind was wrong when it is.
    assert [type(x) for x in found] == [type(x) for x in expected]


@pytest.mark.parametrize("option", ["--work-tree /a", "--git-dir /a/.git", "--git-dir=/a/.git"])
def test_git_dir_and_work_tree_options_make_the_directory_unknown(
    git_dirs: Any, option: str
) -> None:
    found = git_dirs(f"cd /b && git {option} commit", CWD)
    assert len(found) == 1 and isinstance(found[0], str), found
    assert found[0].startswith(option.split()[0].split("=")[0])


# --------------------------------------------------------------------------- 5


def test_unresolvable_directory_is_not_checked_and_nothing_is_gated(repos: Repos) -> None:
    result = repos.subagent('cd "$WT" && git commit -m x')
    assert result.returncode == 2, result.stderr
    assert "NOT CHECKED" in result.stderr
    assert "$WT" in result.stderr
    assert repos.ran() == [], "a gate ran in some other tree in the unknown one's place"


def test_nonexistent_directory_is_not_checked_and_nothing_is_gated(repos: Repos) -> None:
    missing = repos.main.parent / "nonexistent"
    assert not missing.exists(), "premise"
    result = repos.subagent(f"cd {shlex.quote(str(missing))} && git commit -m x")
    assert result.returncode == 2, result.stderr
    assert "NOT CHECKED" in result.stderr
    assert "nonexistent" in result.stderr
    assert repos.ran() == []


def test_committed_trees_maps_directories_to_tops(repos: Repos) -> None:
    committed_trees = function(load_hook(repos.main / HOOK), "committed_trees")
    sub = repos.wt / "tools"
    assert sub.is_dir() and not is_work_tree_top(sub), "premise: a directory inside the worktree"
    assert not Path("/nonexistent").exists(), "premise"
    # The third commit runs in `sub` again (a `-C` does not move the shell), so
    # the worktree appears twice and must be listed once.
    command = (
        f"cd {sub} && git commit && git -C {repos.main} commit && git commit"
        " && cd /nonexistent && git commit"
    )
    event = bash_event(repos.main, command, hook_event_name="PostToolUse", **SUBAGENT)
    trees = committed_trees(event)
    assert trees[:2] == [repos.wt, repos.main]
    assert len(trees) == 3 and isinstance(trees[2], str), trees


# --------------------------------------------------------------------------- 6


def test_two_trees_are_both_gated_and_a_repeated_tree_once(repos: Repos) -> None:
    main, wt = shlex.quote(str(repos.main)), shlex.quote(str(repos.wt))
    command = (
        f"git -C {main} commit -q --allow-empty -m a && "
        f"git -C {wt} commit -q --allow-empty -m b && "
        "git commit -q --allow-empty -m c"  # cwd: the main repository, a second time
    )
    result = repos.subagent(command)
    assert_red_in(result, repos.wt)
    assert sorted(repos.ran_in(repos.main)) == sorted([*CHECKS, "ruff", "ruff"])
    assert sorted(repos.ran_in(repos.wt)) == sorted([*CHECKS, "ruff", "ruff"])


# --------------------------------------------------------------------------- 7


def test_unreadable_line_returns_none(git_dirs: Any) -> None:
    assert git_dirs('git commit -m "unclosed', CWD) is None


def test_unreadable_line_gates_cwd_and_says_it_could_not_be_read(repos: Repos) -> None:
    result = repos.subagent('git commit -m "unclosed')
    assert result.returncode == 2, result.stderr
    assert "NOT CHECKED" in result.stderr
    assert "could not be read" in result.stderr
    assert sorted(repos.ran_in(repos.main)) == sorted([*CHECKS, "ruff", "ruff"])


# --------------------------------------------------------------------------- 8


NO_COMMIT = ["git status", 'echo "git commit"', 'grep -n "git merge" notes.txt']


@pytest.mark.parametrize("command", [*NO_COMMIT, 'git log --grep="git commit"'])
def test_no_commit_no_gate(repos: Repos, command: str) -> None:
    result = repos.subagent(command)
    assert result.returncode == 0, result.stderr
    assert repos.ran() == []


def test_a_non_bash_event_runs_no_gate(repos: Repos) -> None:
    event = bash_event(repos.main, "git commit -m x", hook_event_name="PostToolUse", **SUBAGENT)
    event["tool_name"] = "Edit"
    result = run_script(repos.main, HOOK, event)
    assert result.returncode == 0, result.stderr
    assert repos.ran() == []


# --------------------------------------------------------------------------- 9


def test_hook_appends_tools_to_sys_path_and_never_inserts_first() -> None:
    source = (REAL / HOOK).read_text(encoding="utf-8")
    assert "sys.path.append(" in source
    assert "sys.path.insert(0," not in source


PROBE = """\
import importlib.util, json, sys
from pathlib import Path
before = importlib.util.find_spec("shell_scan") is not None
spec = importlib.util.spec_from_file_location("h19_hook", sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
dirs = "absent"
if hasattr(module, "git_dirs"):
    dirs = module.git_dirs("git commit -m x", Path(sys.argv[2]))
print(json.dumps({
    "findable_before": before,
    "in_modules": "shell_scan" in sys.modules,
    "attribute_is_none": getattr(module, "shell_scan", "absent") is None,
    "git_dirs": repr(dirs),
}))
"""


def test_without_shell_scan_git_dirs_is_none_and_the_text_test_gates_cwd(repos: Repos) -> None:
    (repos.main / SHELL_SCAN).unlink()
    assert not (repos.main / SHELL_SCAN).exists(), "premise"
    # A fresh, isolated interpreter (-I: no PYTHONPATH, no script dir on the
    # path), so no shell_scan already imported or findable can stand in for the
    # deleted copy.
    probe = subprocess.run(
        [sys.executable, "-I", "-c", PROBE, str(repos.main / HOOK), str(repos.main)],
        capture_output=True, text=True, cwd=repos.main, env=clean_env(), timeout=60, check=False,
    )  # fmt: skip
    assert probe.returncode == 0, probe.stderr
    seen = json.loads(probe.stdout)
    assert seen["findable_before"] is False, "premise: shell_scan is importable from elsewhere"
    assert seen["in_modules"] is False
    assert seen["attribute_is_none"] is True
    assert seen["git_dirs"] == "None"

    result = repos.subagent("git commit -m x")
    assert "NOT CHECKED" in result.stderr, result.stderr
    assert result.returncode == 2
    assert sorted(repos.ran_in(repos.main)) == sorted([*CHECKS, "ruff", "ruff"])
