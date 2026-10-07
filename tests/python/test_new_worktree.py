"""Tests for `tools/new_worktree.py` (h18, T4).

The spec is `docs/increments/h18-worktree-and-merge-tools.md` §3, and the
tests are its §3.6, numbered here as there. No network, no uv, no C++ build:
`git` runs for real against a temporary repository whose `origin` is a bare
repository on disk, and `uv`, `cmake` and the venv's Python are answered by
`NewTable` (`worktree_fixtures.Recorder`). Paths are compared after
`Path.resolve()`; the repository is made behind a symbolic link, so the paths
the tool is given are never the resolved ones.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import pytest

from worktree_fixtures import (
    Answer,
    Outcome,
    Recorder,
    Repos,
    failing_fetch,
    git,
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

TOOL = "tools/new_worktree.py"
#: What the table answers for `<py> -m pybind11 --cmakedir`.
CMAKEDIR = "/opt/h18/pybind11/share/cmake/pybind11"


@dataclass
class NewTable:
    """Answers uv, cmake and `<wt>/.venv/bin/python` for new_worktree.py.

    `uv venv` makes a stub `<target>/bin/python`, and the `build-pyext`
    configure writes `CMakeCache.txt`, so the tool's "unless it exists" checks
    see what a real run would leave. The import check (§3.2 step 7) answers as
    `before` until a `uv pip install` has succeeded, then as `after`:

    - `missing`: no `tin_engine` in the venv (exit 1);
    - `worktree`: both modules from `<wt>/src_python`, `_core` from `<wt>/.venv`;
    - `cli-main`: the package from the worktree, `tin_engine.cli` from `<main>`;
    - `main`: everything from the main checkout (a venv of another checkout).

    `<py> -m pybind11 --cmakedir` answers `CMAKEDIR` while `pybind11` is true,
    and exit 1 (no module) otherwise; a successful `uv pip install` whose
    arguments name `pybind11` (a version pin too) makes it true, whichever other packages it names.
    """

    main: Path
    before: str = "missing"
    after: str = "worktree"
    uv_version: Answer = field(default_factory=lambda: Answer(0, "uv 0.12.15\n"))
    cmake_version: Answer = field(default_factory=lambda: Answer(0, "cmake version 4.1.2\n"))
    install: Answer = field(default_factory=Answer)
    installed: bool = False
    pybind11: bool = True

    def __call__(self, argv: list[str], cwd: Path | None) -> Answer:
        name, rest = program(argv), argv[1:]
        if name == "uv":
            if rest == ["--version"]:
                return self.uv_version
            if rest[:1] == ["venv"]:
                stub_venv(Path(positionals(rest[1:])[-1]).parent, "3.14")
                return Answer()
            if rest[:2] == ["pip", "install"]:
                self.installed = self.install.returncode == 0
                if self.installed and any(w.startswith("pybind11") for w in positionals(rest[2:])):
                    self.pybind11 = True
                return self.install
        if name == "cmake":
            if rest == ["--version"]:
                return self.cmake_version
            if "-B" in rest and "--build" not in rest:
                build = Path(option(argv, "-B"))
                build.mkdir(parents=True, exist_ok=True)
                (build / "CMakeCache.txt").write_text("# h18 stub\n")
                return Answer(0, "-- Configuring done\n")
        if is_venv_python(argv):
            wt = Path(argv[0]).parent.parent.parent.resolve()
            if rest[:1] == ["-c"] and "tin_engine" in rest[1]:
                return self.check(wt)
            if rest == ["-m", "ruff", "--version"]:
                return Answer(0, "ruff 0.12.0\n")
            if rest == ["-m", "mypy", "--version"]:
                return Answer(0, "mypy 1.18.2 (compiled: yes)\n")
            if rest == ["-m", "pybind11", "--cmakedir"]:
                if not self.pybind11:
                    return Answer(1, "", f"{argv[0]}: No module named pybind11\n")
                return Answer(0, CMAKEDIR + "\n")
            if rest == ["--version"]:
                return Answer(0, "Python 3.14.7\n")
            if rest[:1] == ["-c"] and "version" in rest[1]:
                return Answer(0, "3.14.7\n")
        pytest.fail(f"new_worktree ran {argv} (cwd {cwd}), which the table does not answer")

    def check(self, wt: Path) -> Answer:
        mode = self.after if self.installed else self.before
        main = self.main.resolve()
        if mode == "missing":
            return Answer(1, "", "ModuleNotFoundError: No module named 'tin_engine'\n")
        package = (main if mode == "main" else wt) / "src_python" / "tin_engine"
        cli = (main if mode in ("main", "cli-main") else wt) / "src_python" / "tin_engine"
        core = (main if mode == "main" else wt) / ".venv" / "lib" / "python3.14"
        core = core / "site-packages" / "tin_engine" / "_core.cpython-314-darwin.so"
        return Answer(0, f"{package / '__init__.py'}\n{cli / 'cli.py'}\n{core}\n")


def positionals(args: list[str]) -> list[str]:
    """The arguments that are neither options nor the value of `--python`/`-p`."""
    found, skip = [], False
    for word in args:
        if skip:
            skip = False
        elif word in ("--python", "-p"):
            skip = True
        elif not word.startswith("-"):
            found.append(word)
    return found


def option(argv: list[str], flag: str) -> str:
    assert flag in argv, f"{flag} missing from {argv}"
    return argv[argv.index(flag) + 1]


def defines(argv: list[str]) -> dict[str, str]:
    """The `-DNAME[:TYPE]=VALUE` cache entries of a cmake argv, `-D NAME=VALUE` too."""
    found: dict[str, str] = {}
    for i, word in enumerate(argv):
        entry = argv[i + 1] if word == "-D" else word[2:] if word.startswith("-D") else None
        if entry and "=" in entry:
            key, value = entry.split("=", 1)
            found[key.split(":", 1)[0]] = value
    return found


def kind(argv: list[str]) -> str | None:
    """Which setup step an argv is: venv, install, configure, build or check."""
    name, rest = program(argv), argv[1:]
    if name == "uv" and rest[:1] == ["venv"]:
        return "venv"
    if name == "uv" and rest[:2] == ["pip", "install"]:
        return "install"
    if name == "cmake" and "--build" in rest:
        return "build"
    if name == "cmake" and rest != ["--version"]:
        return "configure"
    if is_venv_python(argv) and rest[:1] == ["-c"] and "tin_engine" in rest[1]:
        return "check"
    return None


def kinds(recorder: Recorder, *, checks: bool = False) -> list[str]:
    seen = [kind(argv) for argv in recorder.argvs()]
    return [k for k in seen if k is not None and (checks or k != "check")]


def same(a: str | Path, b: Path) -> bool:
    return Path(a).resolve() == b.resolve()


@pytest.fixture
def repos(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Repos:
    """The repositories, with origin's master one commit ahead of main's origin/master."""
    made = make_repos(tmp_path)
    made.master_commit({"notes.txt": "master moved on\n"}, "master moves on")
    made.publish()
    monkeypatch.chdir(made.main)
    return made


def run_tool(
    repos: Repos,
    args: list[str],
    capsys: pytest.CaptureFixture[str],
    table: NewTable | None = None,
    *,
    checkout: Path | None = None,
    fail_fetch: bool = False,
) -> Outcome:
    tool: ModuleType = load(checkout or repos.main, TOOL)
    recorder = Recorder(table or NewTable(repos.main), failing_fetch if fail_fetch else None)
    return invoke(tool, args, recorder, capsys)


def existing_worktree(repos: Repos, name: str, *, venv: bool, configured: bool) -> Path:
    """A worktree made by hand, with or without a stub `.venv` and a configured build-pyext."""
    wt = repos.worktree(name)
    git(repos.main, "worktree", "add", "-q", "-b", f"worktree-{name}", str(wt))
    if venv:
        stub_venv(wt, "3.14")
    if configured:
        (wt / "build-pyext").mkdir()
        (wt / "build-pyext" / "CMakeCache.txt").write_text("# h18 stub\n")
    return wt


def assert_nothing_made(repos: Repos, name: str, outcome: Outcome, *, path: bool = True) -> None:
    """No branch, no worktree add, no venv; and no directory unless the test made it."""
    if path:
        assert not repos.worktree(name).exists()
    assert git(repos.main, "branch", "--list", f"worktree-{name}").strip() == ""
    assert "worktree" not in outcome.recorder.git_subs()
    assert "venv" not in kinds(outcome.recorder)


def refusal(outcome: Outcome) -> str:
    """The refusal line: exit 2, and the last line on stderr."""
    assert outcome.code == 2, outcome.text
    assert outcome.err_lines(), "a refusal writes its reason on stderr"
    return outcome.err_lines()[-1]


# --- 1. where the new worktree lands -----------------------------------------------------


def test_new_worktree_is_under_main_on_its_branch_at_origin_master(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha"], capsys)
    assert outcome.code == 0, outcome.text
    wt = repos.worktree("alpha")
    assert same(git(wt, "rev-parse", "--show-toplevel").strip(), wt)
    assert git(wt, "rev-parse", "--abbrev-ref", "HEAD").strip() == "worktree-alpha"
    # origin's master, which only the tool's fetch brought into main's refs
    assert head(wt) == git(repos.origin, "rev-parse", "master").strip()
    fetches = [a for a in outcome.recorder.argvs() if git_sub(a) == "fetch"]
    assert fetches and all(a[-2:] == ["origin", "master"] for a in fetches), fetches


def test_run_from_inside_another_worktree_still_lands_under_the_main_checkout(
    repos: Repos, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    other = existing_worktree(repos, "other", venv=False, configured=False)
    monkeypatch.chdir(other)
    outcome = run_tool(repos, ["beta"], capsys, checkout=other)
    assert outcome.code == 0, outcome.text
    assert repos.worktree("beta").is_dir()
    assert not (other / ".claude" / "worktrees" / "beta").exists()
    assert git(repos.worktree("beta"), "rev-parse", "--abbrev-ref", "HEAD").strip() == (
        "worktree-beta"
    )


# --- 2. the setup commands -------------------------------------------------------------


def test_setup_commands_are_venv_install_configure_in_that_order_and_no_build(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha"], capsys)
    assert outcome.code == 0, outcome.text
    wt = repos.worktree("alpha")
    py = wt / ".venv" / "bin" / "python"
    assert kinds(outcome.recorder) == ["venv", "install", "configure"]
    by_kind = {kind(c.argv): c for c in outcome.recorder.calls}

    venv = by_kind["venv"].argv
    assert option(venv, "--python") == "3.14"
    assert same(positionals(venv[2:])[-1], wt / ".venv")

    install = by_kind["install"]
    assert same(option(install.argv, "--python"), py)
    assert option(install.argv, "-e") == ".[dev,codecs]"
    assert "pybind11" in install.argv[3:]
    assert install.cwd is not None and same(install.cwd, wt)

    configure = by_kind["configure"].argv
    assert same(option(configure, "-S"), wt)
    assert same(option(configure, "-B"), wt / "build-pyext")
    entries = defines(configure)
    assert same(entries["PYTHON_EXECUTABLE"], py)
    assert same(entries["Python_EXECUTABLE"], py)
    assert entries["pybind11_DIR"] == CMAKEDIR
    assert entries["CMAKE_BUILD_TYPE"] == "Release"
    assert entries["RASPUTIN_BUILD_PYTHON"] == "ON"
    assert entries["RASPUTIN_BUILD_TESTS"] == "OFF"
    assert entries["RASPUTIN_HARDENING"] == "ON"


def test_the_import_check_runs_in_the_worktree_with_its_venv(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha"], capsys)
    assert outcome.code == 0, outcome.text
    wt = repos.worktree("alpha")
    checks = [c for c in outcome.recorder.calls if kind(c.argv) == "check"]
    assert checks, "the step-7 check was not run"
    for call in checks:
        assert same(call.argv[0], wt / ".venv" / "bin" / "python")
        assert call.cwd is not None and same(call.cwd, wt)
        assert "tin_engine.cli" in call.argv[2] and "tin_engine._core" in call.argv[2]
    tail = [a[1:] for a in outcome.recorder.argvs() if is_venv_python(a)]
    assert ["-m", "ruff", "--version"] in tail
    assert ["-m", "mypy", "--version"] in tail


def test_a_pass_prints_the_block_naming_worktree_python_code_core_and_rebuild(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha"], capsys)
    assert outcome.code == 0, outcome.text
    lines = {line.split()[0]: line for line in outcome.out.splitlines() if line.split()}
    for label in ("worktree", "python", "code", "_core", "rebuild"):
        assert label in lines, f"no {label!r} line in:\n{outcome.out}"
    assert "worktree-alpha" in lines["worktree"]
    assert ".venv/bin/python" in lines["python"] and "3.14" in lines["python"]
    assert "src_python" in lines["code"]
    assert "_core.cpython-314-darwin.so" in lines["_core"]
    assert "cmake --build" in lines["rebuild"] and "--target _core" in lines["rebuild"]
    assert "build-pyext" in lines["rebuild"]


# --- 3. refusals ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["Alpha", "a/b", "_a", ".hidden", "a b", ""])
def test_a_bad_name_is_refused_and_nothing_is_made(
    repos: Repos, capsys: pytest.CaptureFixture[str], name: str
) -> None:
    outcome = run_tool(repos, [name], capsys)
    assert outcome.code == 2, outcome.text
    assert outcome.err_lines() == [
        f"new_worktree: name \"{name}\": use lower-case letters, digits, '.', '_' and '-' only"
    ]
    assert "worktree" not in outcome.recorder.git_subs()
    assert not (repos.main / ".claude" / "worktrees").exists()


def test_an_existing_path_is_refused_naming_the_existing_route(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    wt = repos.worktree("alpha")
    wt.mkdir(parents=True)
    outcome = run_tool(repos, ["alpha"], capsys)
    line = refusal(outcome)
    assert outcome.err_lines() == [line]
    hint = " already exists; to give it a venv run: python3 tools/new_worktree.py --existing "
    assert line.startswith("new_worktree: ") and hint in line
    assert same(line.rsplit("--existing ", 1)[1], wt)
    assert list(wt.iterdir()) == []
    assert_nothing_made(repos, "alpha", outcome, path=False)


def test_an_existing_branch_is_refused(repos: Repos, capsys: pytest.CaptureFixture[str]) -> None:
    git(repos.main, "branch", "worktree-alpha")
    outcome = run_tool(repos, ["alpha"], capsys)
    assert (
        refusal(outcome) == "new_worktree: branch worktree-alpha already exists; pick another name"
    )
    assert outcome.err_lines() == [refusal(outcome)]
    assert not repos.worktree("alpha").exists()
    assert "worktree" not in outcome.recorder.git_subs()


@pytest.mark.parametrize("which", ["uv", "cmake"])
def test_a_missing_or_broken_uv_or_cmake_is_refused_before_anything_is_made(
    repos: Repos, capsys: pytest.CaptureFixture[str], which: str
) -> None:
    table = NewTable(repos.main)
    broken = Answer(127, "", f"{which}: not found (h18 broken)\n")
    if which == "uv":
        table.uv_version = broken
    else:
        table.cmake_version = broken
    outcome = run_tool(repos, ["alpha"], capsys, table)
    line = refusal(outcome)
    assert line.startswith(f"new_worktree: {which} is not installed or does not run: ")
    assert "h18 broken" in line
    assert_nothing_made(repos, "alpha", outcome)


def test_a_failed_fetch_is_refused_as_blocked_on_network(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha"], capsys, fail_fetch=True)
    line = refusal(outcome)
    assert line.startswith("new_worktree: could not fetch origin/master (")
    assert "h18 test network down" in line
    assert line.endswith("); blocked on network: stop and hand back")
    assert_nothing_made(repos, "alpha", outcome)


def test_a_failed_worktree_add_makes_no_worktree_and_no_branch(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = run_tool(repos, ["alpha", "--from", "h18-no-such-revision"], capsys)
    line = refusal(outcome)
    assert line.startswith("new_worktree: git worktree add failed (exit ")
    assert line.endswith("its output is above. No worktree was made")
    assert not repos.worktree("alpha").exists()
    assert git(repos.main, "branch", "--list", "worktree-alpha").strip() == ""
    assert "venv" not in kinds(outcome.recorder)
    # --from is not origin/master, so there is nothing to fetch (§3.2 step 2)
    assert "fetch" not in outcome.recorder.git_subs()


def test_existing_refuses_a_directory_inside_a_checkout(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    wt = existing_worktree(repos, "alpha", venv=False, configured=False)
    outcome = run_tool(repos, ["--existing", str(wt / "tools")], capsys)
    line = refusal(outcome)
    assert line.startswith("new_worktree: ")
    assert line.endswith(" is not the top of a checkout of this repository")
    assert kinds(outcome.recorder, checks=True) == []
    assert not (wt / ".venv").exists() and not (wt / "tools" / ".venv").exists()


def test_existing_refuses_a_checkout_of_another_repository(
    repos: Repos, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    stranger = tmp_path / "link" / "stranger"
    stranger.mkdir()
    git(stranger, "init", "-q", "-b", "master")
    git(stranger, "commit", "-q", "--allow-empty", "-m", "elsewhere")
    outcome = run_tool(repos, ["--existing", str(stranger)], capsys)
    line = refusal(outcome)
    assert line.startswith("new_worktree: ")
    assert line.endswith(" is not the top of a checkout of this repository")
    assert kinds(outcome.recorder, checks=True) == []
    assert not (stranger / ".venv").exists()


# --- 4. the import check -----------------------------------------------------------------


def test_a_submodule_from_another_checkout_is_refused_naming_both_paths(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    table = NewTable(repos.main, after="cli-main")
    outcome = run_tool(repos, ["alpha"], capsys, table)
    line = refusal(outcome)
    wt = repos.worktree("alpha")
    elsewhere = repos.main.resolve() / "src_python" / "tin_engine" / "cli.py"
    assert line.startswith("new_worktree: ")
    assert ".venv imports tin_engine.cli from " in line
    assert str(elsewhere) in line
    assert line.endswith("; this venv runs another checkout's code")
    named = line.split(", not from ", 1)[1].split(";", 1)[0]
    assert same(named, wt)
    assert wt.is_dir(), "the worktree is kept"
    assert git(repos.main, "branch", "--list", "worktree-alpha").strip() != ""


# --- 5. --existing ---------------------------------------------------------------------------


def test_existing_with_a_good_venv_and_build_dir_runs_no_setup_step(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    wt = existing_worktree(repos, "alpha", venv=True, configured=True)
    table = NewTable(repos.main, before="worktree")
    outcome = run_tool(repos, ["--existing", str(wt)], capsys, table)
    assert outcome.code == 0, outcome.text
    assert kinds(outcome.recorder) == []
    assert [a for a in outcome.recorder.argvs() if program(a) == "cmake"] in (
        [],
        [["cmake", "--version"]],
    )


def test_existing_repairs_a_venv_that_runs_another_checkout(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    wt = existing_worktree(repos, "alpha", venv=True, configured=True)
    table = NewTable(repos.main, before="main", after="worktree")
    outcome = run_tool(repos, ["--existing", str(wt)], capsys, table)
    assert outcome.code == 0, outcome.text
    assert kinds(outcome.recorder, checks=True) == ["check", "install", "check"]


def test_existing_gets_pybind11_when_the_check_passes_but_build_pyext_is_missing(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    """Code review round 1, finding 1: the venv imports this worktree's code, so
    the check passes, but it has no pybind11 and there is no build-pyext. The
    run must end with pybind11 installed and build-pyext configured, not with a
    refusal that names `--existing` again (which would fail the same way).
    Either fix the reviewer offered passes: install because build-pyext is not
    configured, or install pybind11 alone before the configure."""
    wt = existing_worktree(repos, "alpha", venv=True, configured=False)
    table = NewTable(repos.main, before="worktree", after="worktree", pybind11=False)
    outcome = run_tool(repos, ["--existing", str(wt)], capsys, table)
    assert outcome.code == 0, outcome.text
    assert table.pybind11, "no successful install named pybind11"
    configures = [a for a in outcome.recorder.argvs() if kind(a) == "configure"]
    assert len(configures) == 1, configures
    assert defines(configures[0])["pybind11_DIR"] == CMAKEDIR
    assert (wt / "build-pyext" / "CMakeCache.txt").is_file()
    assert "venv" not in kinds(outcome.recorder)


def test_existing_sets_up_what_is_missing(repos: Repos, capsys: pytest.CaptureFixture[str]) -> None:
    wt = existing_worktree(repos, "alpha", venv=False, configured=False)
    outcome = run_tool(repos, ["--existing", str(wt)], capsys)
    assert outcome.code == 0, outcome.text
    assert kinds(outcome.recorder) == ["venv", "install", "configure"]
    assert "worktree" not in outcome.recorder.git_subs()


# --- 6. a failing step keeps the worktree ---------------------------------------------------


def test_a_failing_install_keeps_the_worktree_and_names_the_existing_route(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    table = NewTable(repos.main, install=Answer(1, "", "error: h18 build failed\n"))
    outcome = run_tool(repos, ["alpha"], capsys, table)
    line = refusal(outcome)
    wt = repos.worktree("alpha")
    assert line.startswith("new_worktree: ")
    tail = "failed (exit 1); its output is above. The worktree is kept; finish with: "
    assert tail + "python3 tools/new_worktree.py --existing " in line
    assert same(line.rsplit("--existing ", 1)[1], wt)
    assert "error: h18 build failed" in outcome.text, "the step's output is passed through"
    assert wt.is_dir()
    assert "configure" not in kinds(outcome.recorder)


# --- 7. nothing the push guard would ask about ----------------------------------------------


def test_no_command_the_tool_runs_is_one_the_push_guard_asks_about(
    repos: Repos, capsys: pytest.CaptureFixture[str]
) -> None:
    fresh = run_tool(repos, ["alpha"], capsys)
    wt = existing_worktree(repos, "beta", venv=True, configured=False)
    repair = run_tool(repos, ["--existing", str(wt)], capsys, NewTable(repos.main, before="main"))
    assert fresh.code == 0 and repair.code == 0, fresh.text + repair.text
    argvs = fresh.recorder.argvs() + repair.recorder.argvs()
    assert push_guard_findings(argvs) == []


def test_the_default_run_answers_a_missing_program_with_127(repos: Repos, tmp_path: Path) -> None:
    tool = load(repos.main, TOOL)
    done = tool.run(["rasputin-h18-no-such-program"], tmp_path)
    assert done.returncode == 127
    assert "not found" in done.stderr
