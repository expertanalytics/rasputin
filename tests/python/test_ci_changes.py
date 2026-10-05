"""h11: skip the code jobs on prose-only pull requests.

`docs/increments/h11-ci-path-filter.md` §6 specifies these tests:

- T1, T2: the classifier's pure functions in `tools/ci_changes.py`.
- T3: `main` end to end on a fixture git repository.
- T4: the shape of `.github/workflows/main.yaml`, read as text (PyYAML is not
  installed, and not added for this). T4 also runs the `CI result` step's shell
  against the rule of §2.3, since that step is the one required check.
- T5: the suite's own reads. The hook lives in `tests/python/conftest.py`; the
  pure part, `prose_reads`, is unit-tested here and the wiring is proved by a
  planted pytest run in a subprocess.

`docs/increments/h15-ci-speed.md` §7 PR 1 adds two kinds of T4 test: every
gating job is checked by the `CI result` step (each slot of its OTHERS list
fails the step when red), and each Python leg is split into a main-suite job
and a `python-extras` job (the `test_h15_*` tests).

The tool is loaded lazily through `harness_fixtures.Tool`, so while it is
missing every test that touches it fails naming `tools/ci_changes.py`.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from harness_fixtures import REAL, Tool, clean_env, git

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

ci = Tool("ci_changes")

WORKFLOW = REAL / ".github" / "workflows" / "main.yaml"
CONFTEST = REAL / "tests" / "python" / "conftest.py"

# ---------------------------------------------------------------------------
# T1: is_prose
# ---------------------------------------------------------------------------

PROSE = (
    "ROADMAP.md",
    "INSTALL.md",
    "testing.md",
    "docs/increments/23c-x.md",
    "docs/retrospectives/2026-10-04-x.md",
    "docs/GLOSSARY.md",
)

CODE = (
    "docs/increments/h4-probes/h4_probe.py",
    "docs/benchmarks/2026-09-26/quarter.geojson",
    ".claude/agents/tester.md",
    ".claude/briefs/common.md",
    ".github/pull_request_template.md",
    ".github/workflows/main.yaml",
    "lib/detria/README.md",
    "tests/python/notes.md",
    "src_python/tin_engine/cli.py",
    "LICENSE",
    "pyproject.toml",
    "CMakeLists.txt",
    ".gitignore",
    "tools/ci_changes.py",
    # Case-sensitive: neither the suffix nor the directory folds case.
    "docs/x.MD",
    "Docs/x.md",
)

#: The Markdown files the build or the suites read (§2.1's table). Every one
#: must be in NOT_PROSE, so the table and the constant cannot drift apart.
READ_BY_CODE = (
    "README.md",
    "CLAUDE.md",
    "docs/increments/README.md",
    "docs/PRINCIPLES.md",
)


@pytest.mark.parametrize("path", PROSE)
def test_t1_prose_paths_are_prose(path: str) -> None:
    assert ci.is_prose(path) is True


@pytest.mark.parametrize("path", CODE)
def test_t1_code_paths_are_not_prose(path: str) -> None:
    assert ci.is_prose(path) is False


@pytest.mark.parametrize("path", READ_BY_CODE)
def test_t1_markdown_read_by_code_is_in_not_prose(path: str) -> None:
    assert path in ci.NOT_PROSE


def test_t1_every_not_prose_entry_is_code() -> None:
    assert ci.NOT_PROSE, "NOT_PROSE is empty"
    for path in ci.NOT_PROSE:
        assert ci.is_prose(path) is False, path


# ---------------------------------------------------------------------------
# T2: needs_full_ci
# ---------------------------------------------------------------------------


def test_t2_empty_list_needs_full_ci() -> None:
    assert ci.needs_full_ci([]) is True


def test_t2_all_prose_does_not_need_full_ci() -> None:
    assert ci.needs_full_ci(list(PROSE)) is False


@pytest.mark.parametrize("code_path", ["src/a.cpp", "README.md", ".gitignore"])
def test_t2_one_code_path_among_prose_needs_full_ci(code_path: str) -> None:
    assert ci.needs_full_ci([*PROSE[:3], code_path, *PROSE[3:]]) is True


# ---------------------------------------------------------------------------
# T3: main end to end on a fixture repository
# ---------------------------------------------------------------------------


def _commit(repo: Path, message: str) -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "rev-parse", "HEAD").strip()


@pytest.fixture
def repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A repository with one code file and one prose file, committed.

    `main` takes only the two SHAs, so it diffs the repository it runs in: the
    test changes into the fixture, and strips GIT_* variables that would
    redirect git elsewhere.
    """
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    (root / "src").mkdir()
    (root / "src" / "a.cpp").write_text("int main() { return 0; }\n" * 20)
    (root / "docs").mkdir()
    (root / "docs" / "notes.md").write_text("notes\n")
    (root / "ROADMAP.md").write_text("# roadmap\n")
    _commit(root, "root")
    for name in [k for k in os.environ if k.startswith("GIT_")]:
        monkeypatch.delenv(name)
    monkeypatch.chdir(root)
    return root


def _head(repo: Path) -> str:
    return git(repo, "rev-parse", "HEAD").strip()


def _run_main(capsys: pytest.CaptureFixture[str], base: str, head: str) -> tuple[int, str, str]:
    code = ci.main([base, head])
    out, err = capsys.readouterr()
    return code, out, err


def test_t3_prose_only_commit_gives_code_false(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    (repo / "ROADMAP.md").write_text("# roadmap\n\nnew row\n")
    (repo / "docs" / "increments").mkdir()
    (repo / "docs" / "increments" / "h99.md").write_text("design\n")
    head = _commit(repo, "prose")
    assert _run_main(capsys, base, head)[:2] == (0, "code=false\n")


def test_t3_code_commit_gives_code_true(repo: Path, capsys: pytest.CaptureFixture[str]) -> None:
    base = _head(repo)
    (repo / "src" / "a.cpp").write_text("int main() { return 1; }\n")
    (repo / "ROADMAP.md").write_text("# roadmap\n\nrow\n")
    head = _commit(repo, "code and prose")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_markdown_read_by_code_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    (repo / "docs" / "increments").mkdir()
    (repo / "docs" / "increments" / "README.md").write_text("protocol\n")
    head = _commit(repo, "protocol text the session-state suite reads")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_rename_from_code_to_prose_path_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Red without `--no-renames`: rename detection lists only `docs/a.md`."""
    base = _head(repo)
    git(repo, "mv", "src/a.cpp", "docs/a.md")
    head = _commit(repo, "move a source file to a prose path")
    # Precondition: git itself sees this as a rename, so the probe can fail.
    status = git(repo, "diff", "--name-status", "-M", base, head)
    assert status.startswith("R"), status
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_changed_paths_lists_both_sides_of_a_rename(repo: Path) -> None:
    base = _head(repo)
    git(repo, "mv", "src/a.cpp", "docs/a.md")
    head = _commit(repo, "rename")
    assert set(ci.changed_paths(base, head, repo)) == {"src/a.cpp", "docs/a.md"}


def test_t3_deleting_a_source_file_alone_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    git(repo, "rm", "-q", "src/a.cpp")
    head = _commit(repo, "delete")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_unusual_prose_file_name_is_read_unquoted(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Red without `-z`: git quotes the name, and `"..."` does not end in `.md`."""
    base = _head(repo)
    (repo / "docs" / "naïve notes.md").write_text("prose\n")
    head = _commit(repo, "unusual name")
    assert _run_main(capsys, base, head)[:2] == (0, "code=false\n")


@pytest.mark.parametrize("base_kind", ["unknown SHA", "empty base", "base equals head"])
def test_t3_unusable_base_falls_back_to_a_full_run(
    repo: Path, capsys: pytest.CaptureFixture[str], base_kind: str
) -> None:
    (repo / "ROADMAP.md").write_text("# roadmap\n\nrow\n")
    head = _commit(repo, "prose")
    base = {
        "unknown SHA": "0123456789abcdef0123456789abcdef01234567",
        "empty base": "",
        "base equals head": head,
    }[base_kind]
    code, out, err = _run_main(capsys, base, head)
    assert (code, out) == (0, "code=true\n")
    assert err.strip(), f"{base_kind}: no sentence on stderr saying why"


# ---------------------------------------------------------------------------
# T4: the workflow's shape
# ---------------------------------------------------------------------------

#: Jobs that do not gate on `changes`: the classifier itself, governance (runs
#: on every event), and the aggregate.
UNGATED = ("changes", "governance", "result")
#: Jobs left out of `CI result`'s needs on purpose (§9, ruling 1).
EXEMPT_FROM_RESULT = ("tsan",)
#: Jobs on master before h11; the parser must find every one, or it is parsing
#: a layout it does not understand.
KNOWN_JOBS = {"cpp", "sanitizers", "tsan", "governance", "python"}

GATE = "if: needs.changes.outputs.code != 'false'"


def _top_level_block(text: str, key: str) -> list[str]:
    """The lines under a top-level `key:` up to the next top-level key."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.rstrip() == f"{key}:"]
    assert len(starts) == 1, f"main.yaml: expected one top-level '{key}:' line"
    block = []
    for line in lines[starts[0] + 1 :]:
        if line and not line[0].isspace() and not line.startswith("#"):
            break
        block.append(line)
    return block


def _jobs(text: str) -> dict[str, list[str]]:
    """Each job's lines, keyed by job id (two-space indent under `jobs:`)."""
    jobs: dict[str, list[str]] = {}
    current: str | None = None
    for line in _top_level_block(text, "jobs"):
        match = re.fullmatch(r"  ([A-Za-z_][\w-]*):\s*", line)
        if match:
            current = match.group(1)
            jobs[current] = []
        elif current is not None:
            jobs[current].append(line)
    missing = KNOWN_JOBS - jobs.keys()
    assert not missing, f"main.yaml: jobs not found by the parser: {sorted(missing)}"
    return jobs


def _job_keys(lines: list[str]) -> dict[str, str]:
    """A job's own keys (four-space indent) and their inline values."""
    keys = {}
    for line in lines:
        match = re.fullmatch(r"    ([\w-]+):\s*(.*?)\s*", line)
        if match:
            keys[match.group(1)] = match.group(2)
    return keys


def _flow_list(value: str) -> set[str]:
    match = re.fullmatch(r"\[(.*)\]", value)
    assert match, f"expected a one-line [a, b] list, got {value!r}"
    return {item.strip() for item in match.group(1).split(",") if item.strip()}


@pytest.fixture(scope="module")
def workflow() -> str:
    return WORKFLOW.read_text()


@pytest.fixture(scope="module")
def jobs(workflow: str) -> dict[str, list[str]]:
    return _jobs(workflow)


def test_t4_triggers_have_no_push(workflow: str) -> None:
    triggers = {
        m.group(1)
        for line in _top_level_block(workflow, "on")
        if (m := re.fullmatch(r"  ([\w-]+):.*", line))
    }
    assert {"pull_request", "merge_group"} <= triggers, triggers
    assert "workflow_dispatch" in triggers
    assert "push" not in triggers


def test_t4_changes_and_result_jobs_exist_under_their_names(
    jobs: dict[str, list[str]],
) -> None:
    assert _job_keys(jobs.get("changes", [])).get("name") == "Changed files"
    # `CI result` is the one required check (§2.3); its name is what branch
    # protection matches, so it is pinned exactly.
    assert _job_keys(jobs.get("result", [])).get("name") == "CI result"


def test_t4_changes_job_runs_the_classifier_with_full_history(
    jobs: dict[str, list[str]],
) -> None:
    body = "\n".join(jobs.get("changes", []))
    assert "fetch-depth: 0" in body
    assert "tools/ci_changes.py" in body
    assert '>> "$GITHUB_OUTPUT"' in body
    assert re.search(r"^    outputs:\s*$", body, re.M), "changes has no outputs:"
    assert re.search(r"^      code: \$\{\{ steps\.[\w-]+\.outputs\.code \}\}", body, re.M)


def test_t4_every_other_job_is_gated_on_changes(jobs: dict[str, list[str]]) -> None:
    gated = [job for job in jobs if job not in UNGATED]
    assert {"cpp", "sanitizers", "tsan", "python"} <= set(gated)
    for job in gated:
        keys = _job_keys(jobs[job])
        assert keys.get("needs") == "changes", f"{job}: needs is {keys.get('needs')!r}"
        assert f"if: {keys.get('if')}" == GATE, f"{job}: if is {keys.get('if')!r}"


def test_t4_governance_runs_on_every_event(jobs: dict[str, list[str]]) -> None:
    keys = _job_keys(jobs["governance"])
    assert "needs" not in keys and "if" not in keys, keys


def test_t4_result_always_runs_and_needs_every_gating_job(
    jobs: dict[str, list[str]],
) -> None:
    keys = _job_keys(jobs.get("result", []))
    assert keys.get("if") == "always()", keys
    expected = set(jobs) - {"result"} - set(EXEMPT_FROM_RESULT)
    assert _flow_list(keys.get("needs", "")) == expected, (
        "a job must be in CI result's needs or in EXEMPT_FROM_RESULT"
    )


def test_t4_unchecked_cpp_leg_does_not_block(jobs: dict[str, list[str]]) -> None:
    keys = _job_keys(jobs["cpp"])
    assert keys.get("continue-on-error") == "${{ matrix.hardening == 'OFF' }}", keys


# T4, continued: the `CI result` step's shell, run against §2.3's rule.


def _result_script(lines: list[str]) -> str:
    """The `run: |` block of the result job's one step, dedented."""
    starts = [i for i, line in enumerate(lines) if re.fullmatch(r"\s+run: \|.*", line)]
    assert len(starts) == 1, "result: expected one 'run: |' block"
    indent = len(lines[starts[0]]) - len(lines[starts[0]].lstrip()) + 2
    script = []
    for line in lines[starts[0] + 1 :]:
        if line.strip() and len(line) - len(line.lstrip()) < indent:
            break
        script.append(line[indent:])
    text = "\n".join(script).strip()
    assert text, "result: the run block is empty"
    assert "${{" not in text, "result: inputs come through env, not ${{ }} in run"
    return text + "\n"


S, K, F, C = "success", "skipped", "failure", "cancelled"


def _verdict(script: str, code: str, changes: str, governance: str, others: str) -> tuple[int, str]:
    bash = shutil.which("bash")
    assert bash is not None
    env = {
        **clean_env(),
        "CODE": code,
        "CHANGES": changes,
        "GOVERNANCE": governance,
        "OTHERS": others,
    }
    # GitHub runs `run:` with `bash --noprofile --norc -eo pipefail {0}`.
    done = subprocess.run(
        [bash, "--noprofile", "--norc", "-eo", "pipefail", "-c", script],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )
    return done.returncode, done.stdout + done.stderr


def _result_env(lines: list[str], key: str) -> str:
    """The value of one `env:` entry of the result job's step."""
    values = [m.group(1) for line in lines if (m := re.fullmatch(rf"\s+{key}: (.*?)\s*", line))]
    assert len(values) == 1, f"result: expected one env entry {key}, got {values}"
    return values[0]


def _needs_results(value: str) -> list[str]:
    """The job ids of every `needs.<id>.result` in `value`, in order."""
    pattern = r"needs(?:\.([\w-]+)|\['([\w-]+)'\])\.result"
    return [dotted or indexed for dotted, indexed in re.findall(pattern, value)]


def _shell_others(jobs: dict[str, list[str]]) -> list[str]:
    """The jobs the result step checks through OTHERS, in its order."""
    return _needs_results(_result_env(jobs.get("result", []), "OTHERS"))


#: An OTHERS value as (fill, position, value): every slot `fill`, except the
#: one at `position` (None: none), which is `value`. Written this way so the
#: cases do not assume how many jobs OTHERS lists; h15 PR 1 adds one.
Others = tuple[str, int | None, str | None]


def _others(n: int, shape: Others) -> str:
    fill, position, value = shape
    slots = [fill] * n
    if position is not None and value is not None:
        slots[position] = value
    return " ".join(slots)


@pytest.mark.parametrize(
    ("code", "changes", "governance", "shape"),
    [
        ("true", S, S, (S, None, None)),
        ("false", S, S, (K, None, None)),
        ("false", S, S, (S, 1, K)),
    ],
    ids=["code change, all green", "prose, all skipped", "prose, mixed"],
)
def test_t4_result_passes(
    jobs: dict[str, list[str]], code: str, changes: str, governance: str, shape: Others
) -> None:
    result_script = _result_script(jobs.get("result", []))
    others = _others(len(_shell_others(jobs)), shape)
    status, output = _verdict(result_script, code, changes, governance, others)
    assert status == 0, output


@pytest.mark.parametrize(
    ("code", "changes", "governance", "shape", "named"),
    [
        ("true", S, S, (S, 1, K), K),
        ("", S, S, (K, None, None), K),
        ("garbage", S, S, (K, None, None), K),
        ("false", S, K, (K, None, None), K),
        ("false", S, F, (K, None, None), F),
        ("", F, S, (K, None, None), F),
        ("true", S, S, (S, 1, F), F),
        ("false", S, S, (K, -1, F), F),
        ("true", S, S, (S, -1, C), C),
        ("false", C, S, (K, None, None), C),
    ],
    ids=[
        "code change, a job skipped",
        "empty code output, jobs skipped",
        "unexpected code output, jobs skipped",
        "governance skipped",
        "governance failed",
        "changes job crashed",
        "code change, a job failed",
        "prose, a job failed",
        "a job cancelled",
        "changes cancelled",
    ],
)
def test_t4_result_fails_and_says_how(
    jobs: dict[str, list[str]],
    code: str,
    changes: str,
    governance: str,
    shape: Others,
    named: str,
) -> None:
    result_script = _result_script(jobs.get("result", []))
    others = _others(len(_shell_others(jobs)), shape)
    status, output = _verdict(result_script, code, changes, governance, others)
    assert status == 1, output
    assert named in output, f"the failure does not say which job ended {named!r}"


# T4, continued (h15 §7 PR 1): no job gates nothing. A job left out of the
# result step's OTHERS, or listed there without a label in its loop, is never
# checked, and `CI result` passes with it red.


def test_t4_every_gating_job_is_in_the_result_shell(jobs: dict[str, list[str]]) -> None:
    result = jobs.get("result", [])
    assert _needs_results(_result_env(result, "CHANGES")) == ["changes"]
    assert _needs_results(_result_env(result, "GOVERNANCE")) == ["governance"]
    others = _shell_others(jobs)
    assert len(others) == len(set(others)), f"OTHERS lists a job twice: {others}"
    expected = set(jobs) - {"changes", "governance", "result"} - set(EXEMPT_FROM_RESULT)
    assert set(others) == expected, "a gating job must be in OTHERS or in EXEMPT_FROM_RESULT"


def test_t4_result_fails_on_each_job_in_others(jobs: dict[str, list[str]]) -> None:
    result_script = _result_script(jobs.get("result", []))
    count = len(_shell_others(jobs))
    for position in range(count):
        others = _others(count, (S, position, F))
        status, output = _verdict(result_script, "true", S, S, others)
        assert status == 1, f"a failure in OTHERS slot {position} passed:\n{output}"
        assert F in output, f"slot {position}: the failure is not named:\n{output}"


# ---------------------------------------------------------------------------
# h15 §5 P and §7 PR 1: each Python leg split into a main-suite job
# (`python`) and an extras job (`python-extras`). Steps are found by their
# `run:` text, which the design keeps unchanged. Only the positions §5 P names
# are tested (Install first, the trap right after it, codecs before viewer
# before the static gates); that every step's text and the full order are
# today's is `@reviewer`'s check by `git diff`, not a test (§7).
# ---------------------------------------------------------------------------

EXTRAS = "python-extras"
MAIN_SUITE = "pytest"
TRAP = "RASPUTIN_INSTALL_TRAP=1"
CODECS = 'pip install -e ".[dev,codecs]"'
VIEWER = 'pip install -e ".[dev,codecs,viewer]"'
STATIC_GATES = ("mypy", "ruff check .", "ruff format --check .")
FLOOR = "if: matrix.python-version == '3.12'"


def _steps(lines: list[str]) -> list[dict[str, str]]:
    """A job's steps, each as its `name`, `uses`, `if` and `run` (block text joined)."""
    steps: list[dict[str, str]] = []
    block: str | None = None
    for line in lines:
        if re.fullmatch(r"      - .*", line):
            steps.append({})
            line = "        " + line[8:]
            block = None
        if not steps:
            continue
        if block is not None and (not line.strip() or len(line) - len(line.lstrip()) > 8):
            steps[-1][block] += ("\n" if steps[-1][block] else "") + line.strip()
            continue
        block = None
        match = re.fullmatch(r"        ([\w-]+):\s*(.*?)\s*", line)
        if match:
            key, value = match.groups()
            if value in ("|", ">"):
                block, value = key, ""
            steps[-1][key] = value
    for step in steps:
        step["run"] = step.get("run", "").strip()
    return steps


def _index(steps: list[dict[str, str]], what: str, test: Callable[[str], bool]) -> int:
    found = [i for i, step in enumerate(steps) if test(step["run"])]
    assert len(found) == 1, f"expected one {what} step, found {len(found)}"
    return found[0]


def _install(steps: list[dict[str, str]]) -> int:
    return _index(steps, "Install", lambda run: 'pip install -e ".[dev]"' in run)


def _python_versions(lines: list[str]) -> list[str]:
    found = [
        m.group(1) for line in lines if (m := re.fullmatch(r"\s+python-version: (\[.*\])\s*", line))
    ]
    assert len(found) == 1, f"expected one matrix python-version line, got {found}"
    return [v.strip().strip('"') for v in found[0][1:-1].split(",") if v.strip()]


@pytest.fixture(scope="module")
def extras(jobs: dict[str, list[str]]) -> list[str]:
    """The extras job's lines; empty while it is missing, so each test fails on its own claim."""
    return jobs.get(EXTRAS, [])


def _assert_exists(extras: list[str]) -> None:
    assert extras, f"main.yaml has no job {EXTRAS!r} (h15 §7 PR 1)"


def test_h15_extras_job_runs_for_every_python_version(
    jobs: dict[str, list[str]], extras: list[str]
) -> None:
    _assert_exists(extras)
    assert _job_keys(extras).get("name") == "Python ${{ matrix.python-version }}, extras"
    assert _python_versions(extras) == _python_versions(jobs["python"])
    assert re.search(r"^      fail-fast: false\s*$", "\n".join(extras), re.M), (
        "one failing leg must not cancel the others, as in the python job"
    )


def test_h15_extras_job_installs_as_the_main_job_does(
    jobs: dict[str, list[str]], extras: list[str]
) -> None:
    _assert_exists(extras)
    main_steps, extras_steps = _steps(jobs["python"]), _steps(extras)
    main_install = main_steps[_install(main_steps)]
    extras_install = extras_steps[_install(extras_steps)]
    # §5 P: "Install (the same step, pip upgrade included)".
    assert "python -m pip install --upgrade pip" in extras_install["run"]
    assert extras_install["run"] == main_install["run"]
    uses = [step.get("uses", "") for step in extras_steps[: _install(extras_steps)]]
    assert uses == [step.get("uses", "") for step in main_steps[: _install(main_steps)]], (
        "checkout and setup-python before Install, as in the python job"
    )


def test_h15_extras_job_runs_the_install_trap_right_after_install(extras: list[str]) -> None:
    _assert_exists(extras)
    steps = _steps(extras)
    trap = _index(steps, "install trap", lambda run: TRAP in run)
    assert trap == _install(steps) + 1, "§5 P: the trap runs right after Install"
    assert f"if: {steps[trap].get('if')}" == FLOOR, "the trap runs on the floor leg only"


def test_h15_extras_job_holds_codecs_and_viewer_steps_after_install(
    extras: list[str],
) -> None:
    _assert_exists(extras)
    steps = _steps(extras)
    install = _install(steps)
    codecs = _index(steps, "codecs", lambda run: CODECS in run)
    viewer = _index(steps, "viewer", lambda run: VIEWER in run)
    assert install < codecs < viewer, (install, codecs, viewer)
    for gate in STATIC_GATES:
        index = _index(steps, gate, lambda run, gate=gate: run == gate)
        assert index > viewer, f"{gate} runs after the viewer step, as today"
        assert f"if: {steps[index].get('if')}" == FLOOR, f"{gate} runs on the floor leg only"


def test_h15_each_python_step_runs_in_one_job_only(
    jobs: dict[str, list[str]], extras: list[str]
) -> None:
    _assert_exists(extras)
    main_runs = [step["run"] for step in _steps(jobs["python"])]
    extras_runs = [step["run"] for step in _steps(extras)]
    # h17 §4b: the main suite's line may carry a marker expression (`-m "not harness"`).
    main_suite = re.compile(rf"{MAIN_SUITE}( -m .*)?")
    assert sum(bool(main_suite.fullmatch(run)) for run in main_runs) == 1, (
        "the main suite runs in the python job"
    )
    assert not any(main_suite.fullmatch(run) for run in extras_runs), (
        "the main suite runs once, not again in extras"
    )
    moved = (TRAP, CODECS, VIEWER, *STATIC_GATES)
    stayed = [what for what in moved if any(what in run for run in main_runs)]
    assert not stayed, f"steps still in the python job, not moved to extras: {stayed}"


def test_h15_extras_job_gates_ci_result(jobs: dict[str, list[str]], extras: list[str]) -> None:
    _assert_exists(extras)
    keys = _job_keys(extras)
    assert keys.get("needs") == "changes" and f"if: {keys.get('if')}" == GATE, keys
    assert EXTRAS in _flow_list(_job_keys(jobs["result"]).get("needs", "")), (
        "CI result does not need the extras job"
    )
    assert EXTRAS in _shell_others(jobs), "the result step does not check the extras job"


# ---------------------------------------------------------------------------
# h17 §4d (docs/increments/h17-ci-test-time.md): the harness tests in a job of
# their own (H1, H2), and the TSan job building exactly the suites it runs
# (H3). That the TSan list drops exactly the four thread-free suites, and that
# the parallel loop reports a failing suite, is `@reviewer`'s check by
# `git diff` and by reading the step, not a test (h15 §7's rule). The harness
# job's gate, its place in CI result's needs and in the result step's OTHERS
# are the T4 tests above, which cover every job.
# ---------------------------------------------------------------------------

HARNESS = "harness"
PYPROJECT = REAL / "pyproject.toml"
CMAKE_TESTS = REAL / "tests" / "cpp" / "CMakeLists.txt"
#: The dev extra's entries the harness job installs (§4b), the bounds read
#: from pyproject.toml rather than restated here.
HARNESS_REQUIREMENTS = ("pytest", "pytest-asyncio", "pytest-cov")


@pytest.fixture(scope="module")
def harness(jobs: dict[str, list[str]]) -> list[str]:
    """The harness job's lines; empty while it is missing, so each test fails on its own claim."""
    return jobs.get(HARNESS, [])


def _assert_harness_exists(harness: list[str]) -> None:
    assert harness, f"main.yaml has no job {HARNESS!r} (h17 §4b)"


def _pytest_step(steps: list[dict[str, str]], job: str) -> str:
    runs = [step["run"] for step in steps if re.match(r"pytest(\s|$)", step["run"])]
    assert len(runs) == 1, f"{job}: expected one step running pytest, found {runs}"
    return runs[0]


def _marker_expression(run: str) -> str | None:
    """The `-m` argument of a pytest command line, unquoted; None when absent."""
    match = re.search(r"\s-m\s+(\"[^\"]*\"|'[^']*'|\S+)", run)
    return match.group(1).strip("\"'") if match else None


def _dev_requirement(name: str) -> str:
    """`name`'s entry in pyproject.toml's dev extra, bound included."""
    block = re.search(r"^dev = \[(.*?)^\]", PYPROJECT.read_text(), re.M | re.S)
    assert block, "pyproject.toml: no dev extra"
    entries = re.findall(r"\"([^\"]+)\"", block.group(1))
    found = [e for e in entries if re.fullmatch(rf"{re.escape(name)}\s*[<>=~!].*", e)]
    assert len(found) == 1, f"dev extra: expected one entry for {name}, got {found}"
    return found[0]


def test_h1_harness_job_runs_the_harness_marker_without_coverage(harness: list[str]) -> None:
    _assert_harness_exists(harness)
    run = _pytest_step(_steps(harness), HARNESS)
    assert _marker_expression(run) == HARNESS, f"harness job's pytest: {run!r}"
    assert "--no-cov" in run.split(), (
        "the package is not installed, so coverage of tin_engine cannot be measured"
    )


def test_h1_harness_job_does_not_install_the_package(harness: list[str]) -> None:
    _assert_harness_exists(harness)
    for step in _steps(harness):
        for line in step["run"].splitlines():
            if "pip install" not in line:
                continue
            args = line.split("pip install", 1)[1].split()
            assert "-e" not in args and "--editable" not in args, f"installs editable: {line!r}"
            local = [a for a in args if a.strip("\"'") == "." or a.strip("\"'").startswith(".[")]
            assert not local, f"installs the package: {line!r}"


def test_h1_harness_job_installs_the_test_tools_at_the_dev_bounds(harness: list[str]) -> None:
    _assert_harness_exists(harness)
    installs = " ".join(
        line
        for step in _steps(harness)
        for line in step["run"].splitlines()
        if "pip install" in line
    )
    for name in HARNESS_REQUIREMENTS:
        requirement = _dev_requirement(name)
        assert requirement in installs, f"harness job does not install {requirement!r}"


def test_h1_harness_job_checks_out_full_history(harness: list[str]) -> None:
    """`test_count_loc` recounts recorded PRs, so the clone needs their commits."""
    _assert_harness_exists(harness)
    steps = _steps(harness)
    checkout = [
        i for i, step in enumerate(steps) if step.get("uses", "").startswith("actions/checkout")
    ]
    assert len(checkout) == 1, f"harness: expected one checkout step, got {len(checkout)}"
    body = "\n".join(harness)
    assert re.search(r"^          fetch-depth: 0\s*$", body, re.M), (
        "harness: checkout without fetch-depth: 0"
    )


def test_h1_harness_job_name_and_python(harness: list[str]) -> None:
    """§4b: named `Python harness tools`, on 3.12 (§8 question 1's default)."""
    _assert_harness_exists(harness)
    assert _job_keys(harness).get("name") == "Python harness tools"
    versions = [
        m.group(1)
        for line in harness
        if (m := re.fullmatch(r"\s+python-version:\s*\"?([\d.]+)\"?\s*", line))
    ]
    assert versions == ["3.12"], f"harness: python-version lines {versions}"


def test_h2_python_job_runs_the_complement_of_the_harness_marker(
    jobs: dict[str, list[str]], harness: list[str]
) -> None:
    """With H1's `-m harness`, every collected test runs in exactly one of the two jobs."""
    main = _marker_expression(_pytest_step(_steps(jobs["python"]), "python"))
    assert main == f"not {HARNESS}", f"python job's pytest -m is {main!r}"
    if harness:
        tools = _marker_expression(_pytest_step(_steps(harness), HARNESS))
        assert main == f"not {tools}", f"{main!r} is not the complement of {tools!r}"


def _tsan_suites(run: str) -> list[str]:
    """Every word of a step that names a C++ test suite (`test_*` or `prop_*`)."""
    return re.findall(r"(?<![\w/$.-])((?:test|prop)_\w+)(?![\w/.-])", run)


def _cmake_test_targets() -> set[str]:
    text = CMAKE_TESTS.read_text()
    return set(re.findall(r"^\s*add_terrain\w*_test\(\s*(\w+)", text, re.M))


def test_h3_tsan_builds_exactly_the_suites_it_runs(jobs: dict[str, list[str]]) -> None:
    steps = _steps(jobs["tsan"])
    build = [step["run"] for step in steps if "--target" in step["run"]]
    test = [step["run"] for step in steps if step.get("name") == "Test"]
    assert len(build) == 1 and len(test) == 1, (
        f"tsan: build steps {len(build)}, test steps {len(test)}"
    )
    built = _tsan_suites(build[0].split("--target", 1)[1])
    ran = _tsan_suites(test[0])
    assert built, "tsan: no suite after --target"
    assert len(built) == len(set(built)), f"tsan builds a suite twice: {built}"
    assert len(ran) == len(set(ran)), f"tsan runs a suite twice: {ran}"
    assert set(built) == set(ran), (
        f"built, not run: {sorted(set(built) - set(ran))}; "
        f"run, not built: {sorted(set(ran) - set(built))}"
    )
    targets = _cmake_test_targets()
    assert len(targets) > 20, f"the CMake parse found only {len(targets)} targets"
    unknown = sorted(set(built) - targets)
    assert not unknown, f"tsan names suites that are not add_terrain*_test targets: {unknown}"


# ---------------------------------------------------------------------------
# T5: the suite's own reads
# ---------------------------------------------------------------------------


def test_t5_prose_reads_names_prose_inside_the_root(tmp_path: Path) -> None:
    root = tmp_path / "root"
    paths = [
        str(root / "ROADMAP.md"),
        str(root / "docs" / "increments" / "h99.md"),
        str(root / "README.md"),  # in NOT_PROSE
        str(root / "docs" / "increments" / "README.md"),  # in NOT_PROSE
        str(root / "tools" / "ci_changes.py"),
        str(root / ".claude" / "agents" / "tester.md"),
        str(tmp_path / "elsewhere" / "ROADMAP.md"),  # outside the root
        str(tmp_path / "ROADMAP.md"),  # outside the root
    ]
    assert set(ci.prose_reads(paths, root)) == {"ROADMAP.md", "docs/increments/h99.md"}


def test_t5_prose_reads_of_nothing_is_empty(tmp_path: Path) -> None:
    assert list(ci.prose_reads([], tmp_path)) == []


PLANTED = """\
from pathlib import Path

REAL = Path({real!r})


def test_planted_read() -> None:
    assert (REAL / {target!r}).read_text()
"""


@pytest.fixture
def planted(tmp_path: Path) -> Callable[[str], subprocess.CompletedProcess[str]]:
    """Run one planted test in a subprocess, with the suite's hook loaded.

    The planted file sits outside `tests/python`, so the conftest is loaded as a
    plugin (`-p conftest`, `tests/python` on the path); `addopts` is cleared so
    coverage settings do not apply to the scratch run.
    """

    def run(target: str) -> subprocess.CompletedProcess[str]:
        assert CONFTEST.exists(), "tests/python/conftest.py (the T5 hook) does not exist"
        test = tmp_path / f"test_planted_{target.replace('/', '_').replace('.', '_')}.py"
        test.write_text(PLANTED.format(real=str(REAL), target=target))
        env = {**clean_env(), "PYTHONPATH": str(CONFTEST.parent)}
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "-p",
                "conftest",
                "-o",
                "addopts=",
                "--rootdir",
                str(tmp_path),
                str(test),
            ],
            cwd=tmp_path,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )

    return run


def test_t5_a_suite_that_reads_prose_fails_the_session(
    planted: Callable[[str], subprocess.CompletedProcess[str]],
) -> None:
    done = planted("ROADMAP.md")
    output = done.stdout + done.stderr
    assert done.returncode != 0, output
    assert "ROADMAP.md" in output, output
    assert "NOT_PROSE" in output, output


def test_t5_a_suite_that_reads_only_code_passes(
    planted: Callable[[str], subprocess.CompletedProcess[str]],
) -> None:
    """The control: the planted run itself works, so the failure above is the hook's."""
    done = planted("pyproject.toml")
    assert done.returncode == 0, done.stdout + done.stderr
